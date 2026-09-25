"""pytest 测试。

参考实现对 n<=5 的全部排列（n!）与全部方向（2^n）进行朴素穷举，
与子集动态规划结果逐项对比。
"""

from __future__ import annotations

import itertools
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

from dna_assembler.assembler import (
    NO_TWO_CONTIGS,
    AssemblyError,
    Read,
    assemble,
    assemble_two_contigs,
    max_overlap,
    parse_payload,
    reverse_complement,
    requested_contig_count,
)

REPO_ROOT = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------
# 独立的全枚举参考实现（刻意不依赖 assembler 内部逻辑，只复用 reverse_complement）
# --------------------------------------------------------------------------


def brute_force_assemble(reads: list[Read]) -> tuple[str, list[tuple[str, str]]]:
    """枚举所有排列与方向，返回 (装配串, [(id, 方向), ...])。"""
    best: tuple[int, str, tuple[tuple[str, str], ...]] | None = None
    n = len(reads)
    for order in itertools.permutations(range(n)):
        for bits in range(1 << n):
            seqs: list[str] = []
            path: list[tuple[str, str]] = []
            for pos, idx in enumerate(order):
                ori = "+" if (bits >> pos) & 1 == 0 else "-"
                seq = reads[idx].seq if ori == "+" else reverse_complement(
                    reads[idx].seq
                )
                seqs.append(seq)
                path.append((reads[idx].id, ori))
            assembly = seqs[0]
            for prev, nxt in zip(seqs, seqs[1:]):
                shift = len(nxt) - max_overlap(prev, nxt)
                assembly += nxt[len(nxt) - shift:]
            key = (len(assembly), assembly, tuple(path))
            if best is None or key < best:
                best = key
    assert best is not None
    return best[1], list(best[2])


def expected_layout(
    reads: list[Read], assembly: str, path: list[tuple[str, str]]
) -> list[dict]:
    """按路径逐步最大重叠，重放各读段 0 基半开区间的起止位置。"""
    by_id = {r.id: r for r in reads}
    placed: list[tuple[str, str, str]] = []
    for rid, sign in path:
        seq = by_id[rid].seq if sign == "+" else reverse_complement(by_id[rid].seq)
        placed.append((rid, sign, seq))
    layout: list[dict] = []
    cursor = 0
    prev_seq: str | None = None
    for rid, sign, seq in placed:
        if prev_seq is None:
            start = 0
        else:
            start = cursor - max_overlap(prev_seq, seq)
        layout.append(
            {"id": rid, "orientation": sign, "start": start, "end": start + len(seq)}
        )
        cursor = start + len(seq)
        prev_seq = seq
    assert cursor == len(assembly)
    return layout


def assert_matches_brute(reads: list[Read]) -> dict:
    """核心结果与全枚举参考逐项一致，并校验位置/切片自洽。"""
    result = assemble(reads)
    exp_assembly, exp_path = brute_force_assemble(reads)
    assert result["assembly"] == exp_assembly
    assert result["length"] == len(exp_assembly)
    assert [(item["id"], item["orientation"]) for item in result["path"]] == exp_path
    assert result["layout"] == expected_layout(reads, exp_assembly, exp_path)
    assert result["orientations"] == {rid: sign for rid, sign in exp_path}

    # 切片自洽：每条读段所选方向上的序列必须精确落在装配串对应区间。
    by_id = {r.id: r for r in reads}
    for item in result["layout"]:
        seq = by_id[item["id"]].seq
        if item["orientation"] == "-":
            seq = reverse_complement(seq)
        assert result["assembly"][item["start"]:item["end"]] == seq
        assert 0 <= item["start"] < item["end"] <= len(result["assembly"])
    return result


# --------------------------------------------------------------------------
# 恰好两条 contig 的全枚举参考实现
# --------------------------------------------------------------------------


def brute_force_two_contigs(
    reads: list[Read],
) -> tuple[
    tuple[
        tuple[str, tuple[tuple[str, str], ...]],
        tuple[str, tuple[tuple[str, str], ...]],
    ],
    int,
] | None:
    """枚举所有二分、排列和方向，返回规范 contig 对与总长度。"""
    n = len(reads)
    best: tuple[
        int,
        tuple[
            tuple[str, tuple[tuple[str, str], ...]],
            tuple[str, tuple[tuple[str, str], ...]],
        ],
    ] | None = None

    for assignment in range(1, 1 << n):
        groups = ([], [])
        for idx in range(n):
            groups[(assignment >> idx) & 1].append(idx)
        if not groups[0] or not groups[1]:
            continue

        descriptors = []
        valid = True
        for group in groups:
            group_best: tuple[
                int, str, tuple[tuple[str, str], ...]
            ] | None = None
            for order in itertools.permutations(group):
                for bits in range(1 << len(order)):
                    seqs = []
                    path = []
                    for pos, idx in enumerate(order):
                        sign = "+" if (bits >> pos) & 1 == 0 else "-"
                        seqs.append(
                            reads[idx].seq
                            if sign == "+"
                            else reverse_complement(reads[idx].seq)
                        )
                        path.append((reads[idx].id, sign))
                    if any(
                        max_overlap(prev, nxt) == 0
                        for prev, nxt in zip(seqs, seqs[1:])
                    ):
                        continue
                    assembly = seqs[0]
                    for prev, nxt in zip(seqs, seqs[1:]):
                        assembly += nxt[max_overlap(prev, nxt):]
                    candidate = (len(assembly), assembly, tuple(path))
                    if group_best is None or candidate < group_best:
                        group_best = candidate
            if group_best is None:
                valid = False
                break
            descriptors.append((group_best[1], group_best[2], group_best[0]))
        if not valid:
            continue

        descriptors.sort(key=lambda item: (item[0], item[1]))
        canonical = (
            (descriptors[0][0], descriptors[0][1]),
            (descriptors[1][0], descriptors[1][1]),
        )
        candidate = (
            descriptors[0][2] + descriptors[1][2],
            canonical,
        )
        if best is None or candidate < best:
            best = candidate

    return None if best is None else (best[1], best[0])


def assert_two_contigs_matches_brute(reads: list[Read]) -> dict | str:
    """双 contig 结果与独立穷举逐项一致，并校验局部坐标。"""
    result = assemble_two_contigs(reads)
    expected = brute_force_two_contigs(reads)
    if expected is None:
        assert result == NO_TWO_CONTIGS
        return result

    assert isinstance(result, dict)
    canonical, total_length = expected
    assert result["total_length"] == total_length
    assert set(result) == {"total_length", "contigs"}
    assert len(result["contigs"]) == 2
    assert all(contig["length"] >= 1 for contig in result["contigs"])

    by_id = {r.id: r for r in reads}
    observed = []
    all_ids = []
    for contig, (assembly, exp_path) in zip(result["contigs"], canonical):
        assert set(contig) == {"length", "path", "layout"}
        obs_path = tuple(
            (item["id"], item["orientation"]) for item in contig["path"]
        )
        assert obs_path == exp_path
        assert contig["length"] == len(assembly)
        assert contig["layout"] == expected_layout(
            reads, assembly, list(exp_path)
        )
        observed.append((assembly, obs_path))
        all_ids.extend(rid for rid, _sign in exp_path)

        for item in contig["layout"]:
            seq = by_id[item["id"]].seq
            if item["orientation"] == "-":
                seq = reverse_complement(seq)
            assert assembly[item["start"]:item["end"]] == seq
            assert 0 <= item["start"] < item["end"] <= len(assembly)

    assert tuple(observed) == canonical
    assert sorted(all_ids) == sorted(r.id for r in reads)
    assert len(all_ids) == len(set(all_ids)) == len(reads)
    return result


def payload(*pairs: tuple[str, str]) -> dict:
    return {"reads": [{"id": rid, "seq": seq} for rid, seq in pairs]}


# --------------------------------------------------------------------------
# 固定场景
# --------------------------------------------------------------------------


def test_simple_overlap():
    # 正向后缀/前缀直接重叠 "GT"
    reads = parse_payload(payload(("a", "AGCGT"), ("b", "GTTCA")))
    result = assert_matches_brute(reads)
    assert result["assembly"] == "AGCGTTCA"
    assert [(x["id"], x["orientation"]) for x in result["path"]] == [
        ("a", "+"),
        ("b", "+"),
    ]
    assert result["length"] == 8
    starts = {x["id"]: (x["start"], x["end"]) for x in result["layout"]}
    assert starts == {"a": (0, 5), "b": (3, 8)}


def test_needs_reverse_complement():
    # b 的反向互补是 GTTCA，与 a(AGCGT) 形成 "GT" 重叠；正向上无法装配得同样短。
    reads = parse_payload(payload(("a", "AGCGT"), ("b", "TGAAC")))
    result = assert_matches_brute(reads)
    assert result["assembly"] == "AGCGTTCA"
    chosen = {x["id"]: x["orientation"] for x in result["path"]}
    assert chosen == {"a": "+", "b": "-"}


def test_no_overlap_pair():
    # AAAA 与 CCCC 的任一对方向后缀/前缀均无公共碱基，直接拼接。
    reads = parse_payload(payload(("x", "AAAA"), ("y", "CCCC")))
    result = assert_matches_brute(reads)
    # 两种拼接长度相同，字典序 AAAA<CCCC，故 AAAA 在前
    assert result["assembly"] == "AAAACCCC"
    assert result["length"] == 8
    assert [x["id"] for x in result["path"]] == ["x", "y"]
    layout = {x["id"]: x for x in result["layout"]}
    assert (layout["x"]["start"], layout["x"]["end"]) == (0, 4)
    assert (layout["y"]["start"], layout["y"]["end"]) == (4, 8)


def test_reverse_complement_palindrome_read():
    # ACGT 与 GTAC 都是反向互补回文（自身 == 反向互补），
    # 两种排列分别以 "GT" / "AC" 重叠，长度同为 6，字典序 ACGTAC 更小。
    reads = parse_payload(payload(("p", "ACGT"), ("q", "GTAC")))
    result = assert_matches_brute(reads)
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("GTAC") == "GTAC"
    assert result["assembly"] == "ACGTAC"
    # 装配串相同的并列裁决到路径字典序：全部取 '+'（'+' < '-'）
    assert all(item["orientation"] == "+" for item in result["path"])
    assert result["orientations"] == {"p": "+", "q": "+"}
    assert {x["id"]: (x["start"], x["end"]) for x in result["layout"]} == {
        "p": (0, 4),
        "q": (2, 6),
    }


def test_plain_reverse_palindrome_read():
    # ATTA 是普通反转回文但不是反向互补回文；其反向互补 TAAT 与
    # AATCG 的前缀以 "AAT"（长度 3）重叠。反向排列（s 反向）同样长度，
    # 两个串 CGATTA / TAATCG 并列，按字典序取 CGATTA。
    reads = parse_payload(payload(("r", "ATTA"), ("s", "AATCG")))
    assert reads[0].seq == reads[0].seq[::-1]
    assert reverse_complement(reads[0].seq) == "TAAT"
    result = assert_matches_brute(reads)
    assert result["length"] == 6
    assert result["assembly"] == "CGATTA"
    # s 反向得 CGATT，与正向 ATTA 以 "ATT" 重叠
    assert [(i["id"], i["orientation"]) for i in result["path"]] == [
        ("s", "-"),
        ("r", "+"),
    ]
    layout = {x["id"]: (x["start"], x["end"]) for x in result["layout"]}
    assert layout == {"s": (0, 5), "r": (2, 6)}


def enumerate_all_minima(reads: list[Read]) -> list[tuple[str, tuple[tuple[str, str], ...]]]:
    """枚举所有取得全局最短长度的 (装配串, 路径)。"""
    minima: dict[str, list[tuple[tuple[str, str], ...]]] = {}
    best_len: int | None = None
    n = len(reads)
    for order in itertools.permutations(range(n)):
        for bits in range(1 << n):
            seqs: list[str] = []
            path: list[tuple[str, str]] = []
            for pos, idx in enumerate(order):
                sign = "+" if (bits >> pos) & 1 == 0 else "-"
                seqs.append(
                    reads[idx].seq
                    if sign == "+"
                    else reverse_complement(reads[idx].seq)
                )
                path.append((reads[idx].id, sign))
            assembly = seqs[0]
            for prev, nxt in zip(seqs, seqs[1:]):
                shift = len(nxt) - max_overlap(prev, nxt)
                assembly += nxt[len(nxt) - shift:]
            if best_len is None or len(assembly) < best_len:
                best_len = len(assembly)
                minima = {assembly: [tuple(path)]}
            elif len(assembly) == best_len:
                minima.setdefault(assembly, []).append(tuple(path))
    assert best_len is not None
    return [(assembly, min(paths)) for assembly, paths in minima.items()]


def test_no_overlap_then_lexicographic_order():
    # AAAA 与 CCCC 在所有 4 个方向组合下后缀/前缀重叠均为 0，
    # 8 种排列×方向装配长度同为 8，取字典序最小的拼接 AAAACCCC。
    reads = parse_payload(payload(("x", "AAAA"), ("y", "CCCC")))
    minima = enumerate_all_minima(reads)
    assert len(minima) == 8
    result = assert_matches_brute(reads)
    assert result["assembly"] == "AAAACCCC"
    assert [(i["id"], i["orientation"]) for i in result["path"]] == [
        ("x", "+"),
        ("y", "+"),
    ]
    layout = {i["id"]: (i["start"], i["end"]) for i in result["layout"]}
    assert layout == {"x": (0, 4), "y": (4, 8)}


def test_tied_optimum_picks_lexicographic_assembly():
    # a=AAAAG, rc(a)=CTTTT；b=GCCCC, rc(b)=GGGGC。
    # 四个方向/排列组合各以 1 个碱基重叠，产生 4 个长度同为 9 的不同装配串，
    # 必须按装配串字典序裁决出 AAAAGCCCC。
    reads = parse_payload(payload(("a", "AAAAG"), ("b", "GCCCC")))
    assert reverse_complement("AAAAG") == "CTTTT"
    assert reverse_complement("GCCCC") == "GGGGC"
    minima = enumerate_all_minima(reads)
    assert {assembly for assembly, _ in minima} == {
        "AAAAGCCCC",
        "AAAAGGGGC",
        "GCCCCTTTT",
        "GGGGCTTTT",
    }
    result = assert_matches_brute(reads)
    assert result["assembly"] == "AAAAGCCCC"
    assert [(i["id"], i["orientation"]) for i in result["path"]] == [
        ("a", "+"),
        ("b", "+"),
    ]
    layout = {i["id"]: (i["start"], i["end"]) for i in result["layout"]}
    assert layout == {"a": (0, 5), "b": (4, 9)}


def test_tied_optimum_same_string_path_tiebreak():
    # RC 回文读段 ACGT 与 GTAC：两个方向产生同一批装配串，
    # 串相同时必须进一步按 id/方向路径字典序裁决（'+' 恒小于 '-'）。
    reads = parse_payload(payload(("p", "ACGT"), ("q", "GTAC")))
    result = assert_matches_brute(reads)
    assert result["assembly"] == "ACGTAC"
    assert [(i["id"], i["orientation"]) for i in result["path"]] == [
        ("p", "+"),
        ("q", "+"),
    ]


def test_chained_three_reads():
    reads = parse_payload(
        payload(
            ("r1", "GATTACA"),
            ("r2", "ACATAGA"),
            ("r3", "TAGACAT"),
        )
    )
    assert_matches_brute(reads)


def test_five_reads_exhaustive():
    reads = parse_payload(
        payload(
            ("id1", "TGCATGCA"),
            ("id2", "GCATCGAT"),
            ("id3", "CGATTTAA"),
            ("id4", "TTAAGGCC"),
            ("id5", "GGCCATGC"),
        )
    )
    result = assert_matches_brute(reads)
    # 形成闭环状重叠，最短长度必短于总长 40
    assert result["length"] < 40
    # 同一组五条读段也穷举所有二分、排列和方向。
    assert_two_contigs_matches_brute(reads)


# --------------------------------------------------------------------------
# 恰好两条 contig：固定场景
# --------------------------------------------------------------------------


def test_two_contigs_zero_overlap_isolation():
    # a/b 在某方向下可以连成一条区域；c 与它们的任一方向均无后缀/前缀重叠。
    reads = parse_payload(
        payload(
            ("a", "AGCGT"),
            ("b", "TGAAC"),
            ("c", "GATC"),
        )
    )
    result = assert_two_contigs_matches_brute(reads)
    assert isinstance(result, dict)
    assert result["total_length"] == 12
    first, second = result["contigs"]

    # 规范排序：装配串 AGCGTTCA < GATC
    assert first["length"] == 8
    assert [(x["id"], x["orientation"]) for x in first["path"]] == [
        ("a", "+"),
        ("b", "-"),
    ]
    assert {(x["id"], x["start"], x["end"]) for x in first["layout"]} == {
        ("a", 0, 5),
        ("b", 3, 8),
    }

    assert second["length"] == 4
    assert [(x["id"], x["orientation"]) for x in second["path"]] == [
        ("c", "+")
    ]
    assert second["layout"] == [
        {"id": "c", "orientation": "+", "start": 0, "end": 4}
    ]
    assert "assembly" not in result
    assert all("assembly" not in contig for contig in result["contigs"])


def test_two_contigs_positive_overlap_and_local_coordinates():
    reads = parse_payload(
        payload(
            ("r1", "GATTACA"),
            ("r2", "ACATAGA"),
            ("z", "GGATTCGCC"),
        )
    )
    result = assert_two_contigs_matches_brute(reads)
    assert isinstance(result, dict)
    assert result["total_length"] == 20
    connected, isolated = result["contigs"]
    assert connected["length"] == 11
    assert [(x["id"], x["orientation"]) for x in connected["path"]] == [
        ("r1", "+"),
        ("r2", "+"),
    ]
    assert connected["layout"] == [
        {"id": "r1", "orientation": "+", "start": 0, "end": 7},
        {"id": "r2", "orientation": "+", "start": 4, "end": 11},
    ]
    assert isolated["length"] == 9
    assert isolated["layout"] == [
        {"id": "z", "orientation": "+", "start": 0, "end": 9}
    ]


def test_two_contigs_same_string_path_tiebreak():
    # 两条 RC 回文读段串成 ACGTAC；隔离读段为 GATC。
    # 相邻读段的 +/- 产生同一装配串，规范路径必须取 '+'。
    reads = parse_payload(
        payload(
            ("p", "ACGT"),
            ("q", "GTAC"),
            ("z", "GATC"),
        )
    )
    result = assert_two_contigs_matches_brute(reads)
    assert isinstance(result, dict)
    connected, isolated = result["contigs"]
    assert connected["length"] == 6
    assert [(x["id"], x["orientation"]) for x in connected["path"]] == [
        ("p", "+"),
        ("q", "+"),
    ]
    assert connected["layout"] == [
        {"id": "p", "orientation": "+", "start": 0, "end": 4},
        {"id": "q", "orientation": "+", "start": 2, "end": 6},
    ]
    assert isolated["length"] == 4
    assert isolated["layout"] == [
        {"id": "z", "orientation": "+", "start": 0, "end": 4}
    ]


def test_two_contigs_canonical_pair_tiebreak():
    # 两条 RC 回文读段可产生同长且同串的两条 singleton/connected 方案；
    # 规范结果应先按 contig 串、再按路径排序，不依赖上传次序。
    reads = parse_payload(payload(("p", "ACGT"), ("q", "GTAC")))
    reversed_upload = parse_payload(payload(("q", "GTAC"), ("p", "ACGT")))
    result = assert_two_contigs_matches_brute(reads)
    assert assert_two_contigs_matches_brute(reversed_upload) == result
    assert isinstance(result, dict)
    assert result["total_length"] == 8
    assert [
        [x["id"] for x in contig["path"]] for contig in result["contigs"]
    ] == [["p"], ["q"]]
    assert [contig["layout"] for contig in result["contigs"]] == [
        [{"id": "p", "orientation": "+", "start": 0, "end": 4}],
        [{"id": "q", "orientation": "+", "start": 0, "end": 4}],
    ]


def test_two_contigs_no_valid_split():
    # 正/反向边都只经过 a，无法选出两个互不共享 id 的非空内部重叠路径。
    reads = parse_payload(
        payload(
            ("a", "AAAA"),
            ("b", "AAAC"),
            ("c", "AAAT"),
            ("d", "AAGC"),
        )
    )
    assert assert_two_contigs_matches_brute(reads) == NO_TWO_CONTIGS
    assert assemble(reads, contigs=2) == NO_TWO_CONTIGS


def test_two_contigs_requires_two_reads():
    reads = [Read("only", "ACGT")]
    with pytest.raises(AssemblyError) as exc:
        assemble_two_contigs(reads)
    assert exc.value.code == "bad_read_count"


def test_requested_contig_count_modes():
    assert requested_contig_count({}) == 1
    assert requested_contig_count({"contigs": 2}) == 2
    assert requested_contig_count({"mode": "two_contigs"}) == 2
    assert requested_contig_count(
        {"contigs": 2, "mode": "two_contigs"}
    ) == 2
    with pytest.raises(AssemblyError) as exc:
        requested_contig_count({"contigs": 3})
    assert exc.value.code == "bad_contig_count"
    with pytest.raises(AssemblyError) as exc:
        requested_contig_count({"contigs": 2, "mode": "single"})
    assert exc.value.code == "bad_mode"


# --------------------------------------------------------------------------
# 随机合法输入：5 条以内，全部与全枚举参考比对
# --------------------------------------------------------------------------


def generate_valid_reads(
    rng: random.Random, n: int
) -> list[Read] | None:
    """生成 n 条通过全部校验的随机读段；失败返回 None。"""
    candidates: list[Read] = []
    attempts = 0
    while len(candidates) < n and attempts < 3000:
        attempts += 1
        length = rng.randint(4, 10)
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        rid = f"s{len(candidates)}"
        trial = candidates + [Read(rid, seq)]
        if len(trial) >= 2:
            try:
                parse_payload(
                    {"reads": [{"id": r.id, "seq": r.seq} for r in trial]}
                )
            except AssemblyError:
                continue
        candidates.append(Read(rid, seq))
    if len(candidates) < n:
        return None
    return candidates


@pytest.mark.parametrize("seed", range(30))
def test_random_against_brute_force(seed: int):
    rng = random.Random(1000 + seed)
    n = rng.randint(2, 5)
    reads = generate_valid_reads(rng, n)
    if reads is None:
        pytest.skip("随机生成器未能构造出合法读段集合")
    assert_matches_brute(reads)


@pytest.mark.parametrize("seed", range(30))
def test_random_two_contigs_against_brute_force(seed: int):
    rng = random.Random(2000 + seed)
    n = rng.randint(2, 5)
    reads = generate_valid_reads(rng, n)
    if reads is None:
        pytest.skip("随机生成器未能构造出合法读段集合")
    assert_two_contigs_matches_brute(reads)


# --------------------------------------------------------------------------
# 拒绝场景
# --------------------------------------------------------------------------


def test_reject_reverse_complement_pair():
    with pytest.raises(AssemblyError) as exc:
        parse_payload(payload(("a", "AACGT"), ("b", "ACGTT")))
    assert reverse_complement("AACGT") == "ACGTT"
    assert exc.value.code == "reverse_complement_pair"


def test_reject_containment_same_orientation():
    with pytest.raises(AssemblyError) as exc:
        parse_payload(payload(("a", "ACGTACGT"), ("b", "CGTA")))
    assert exc.value.code == "contained_read"


def test_reject_containment_after_reverse_complement():
    # rc(CCACGT) = ACGTGG，完整包含于 AACGTGG 中（此前校验漏过的方向）
    assert reverse_complement("CCACGT") == "ACGTGG"
    assert "ACGTGG" in "AACGTGG"
    reads_obj = payload(("a", "AACGTGG"), ("b", "CCACGT"))
    with pytest.raises(AssemblyError) as exc:
        parse_payload(reads_obj)
    assert exc.value.code == "contained_read"


@pytest.mark.parametrize(
    "bad",
    [
        {},
        {"reads": "not-a-list"},
        {"reads": []},
        {"reads": [{"id": "a", "seq": "ACGT"}]},  # 只有 1 条
        {"reads": [{"id": "a", "seq": "ACGT"} for _ in range(11)]},  # 11 条
        {"reads": [{"id": "a", "seq": "ACG"}]},  # 长度 3
        {"reads": [{"id": "a", "seq": "ACGTN"}]},  # 非法碱基
        {"reads": [{"id": "a", "seq": "acgt"}]},  # 小写
        {"reads": [{"id": "", "seq": "ACGT"}]},  # 空 id
        {"reads": [{"id": "a\x00", "seq": "ACGT"}]},  # 非 ASCII
        {"reads": [
            {"id": "a", "seq": "ACGT"},
            {"id": "a", "seq": "TTGG"},
        ]},  # 重复 id
        {"reads": [{"id": "a", "seq": 31337}]},
        {"reads": ["not-an-object"]},
    ],
)
def test_reject_malformed_payload(bad):
    with pytest.raises(AssemblyError):
        parse_payload(bad)


# --------------------------------------------------------------------------
# 输出完整性：不能只给长度
# --------------------------------------------------------------------------


def test_result_contains_full_layout_and_orientations():
    reads = parse_payload(payload(("a", "AAAAG"), ("b", "GCCCC"), ("c", "CCCTT")))
    result = assemble(reads)
    assert set(result) >= {"assembly", "length", "layout", "path", "orientations"}
    assert isinstance(result["assembly"], str)
    assert len(result["layout"]) == 3
    assert set(result["orientations"]) == {"a", "b", "c"}
    for item in result["layout"]:
        assert set(item) == {"id", "orientation", "start", "end"}
        assert item["orientation"] in "+-"


# --------------------------------------------------------------------------
# CLI 端到端
# --------------------------------------------------------------------------


def run_cli(data: dict | str) -> subprocess.CompletedProcess:
    if isinstance(data, dict):
        raw = json.dumps(data)
    else:
        raw = data
    return subprocess.run(
        [sys.executable, "-m", "dna_assembler"],
        input=raw,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )


def test_cli_success():
    proc = run_cli(payload(("a", "AGCGT"), ("b", "TGAAC")))
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    assert out["assembly"] == "AGCGTTCA"
    assert out["length"] == 8
    assert proc.stderr == ""


def test_cli_two_contigs_success():
    request = payload(
        ("a", "AGCGT"),
        ("b", "TGAAC"),
        ("c", "GATC"),
    )
    request["contigs"] = 2
    proc = run_cli(request)
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    assert out["total_length"] == 12
    assert "assembly" not in out
    assert [
        [item["id"] for item in contig["path"]]
        for contig in out["contigs"]
    ] == [["a", "b"], ["c"]]


def test_cli_two_contigs_mode_alias():
    request = payload(("a", "AAAA"), ("b", "CCGG"))
    request["mode"] = "two_contigs"
    proc = run_cli(request)
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    assert out["total_length"] == 8


def test_cli_two_contigs_no_split():
    request = payload(
        ("a", "AAAA"),
        ("b", "AAAC"),
        ("c", "AAAT"),
        ("d", "AAGC"),
    )
    request["contigs"] = 2
    proc = run_cli(request)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == NO_TWO_CONTIGS
    assert proc.stderr == ""


def test_cli_invalid_json():
    proc = run_cli("{not json")
    assert proc.returncode == 2
    err = json.loads(proc.stderr)
    assert err["code"] == "invalid_json"
    assert proc.stdout == ""


def test_cli_rejected_data():
    proc = run_cli(payload(("a", "AACGT"), ("b", "ACGTT")))
    assert proc.returncode == 2
    err = json.loads(proc.stderr)
    assert err["code"] == "reverse_complement_pair"


def test_cli_file_argument(tmp_path: Path):
    p = tmp_path / "in.json"
    p.write_text(json.dumps(payload(("x", "AAAA"), ("y", "CCCC"))), encoding="utf-8")
    proc = subprocess.run(
        [sys.executable, "-m", "dna_assembler", str(p)],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert proc.returncode == 0, proc.stderr
    assert json.loads(proc.stdout)["assembly"] == "AAAACCCC"
