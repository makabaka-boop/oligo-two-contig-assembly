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


def path_sequence_and_layout(
    reads: list[Read], path: list[tuple[str, str]]
) -> tuple[str, list[dict]]:
    """按定向路径生成串并计算局部 0 基半开坐标。"""

    by_id = {r.id: r for r in reads}
    placed: list[tuple[str, str, str]] = []
    for rid, sign in path:
        seq = by_id[rid].seq if sign == "+" else reverse_complement(by_id[rid].seq)
        placed.append((rid, sign, seq))

    assembly = placed[0][2]
    layout: list[dict] = []
    cursor = 0
    prev_seq: str | None = None
    for rid, sign, seq in placed:
        if prev_seq is None:
            start = 0
        else:
            overlap = max_overlap(prev_seq, seq)
            assert overlap >= 1
            start = cursor - overlap
        end = start + len(seq)
        layout.append(
            {"id": rid, "orientation": sign, "start": start, "end": end}
        )
        cursor = end
        if prev_seq is not None:
            assembly += seq[max_overlap(prev_seq, seq):]
        prev_seq = seq
    return assembly, layout


def brute_force_two_contigs(
    reads: list[Read],
) -> tuple[int, list[tuple[str, list[tuple[str, str]]]]]:
    """穷举所有非空二分、组内排列和方向，返回二 contig 规范结果。"""

    n = len(reads)
    full_mask = (1 << n) - 1
    best_by_mask: dict[
        int, tuple[int, str, tuple[tuple[str, str], ...]]
    ] = {}

    for mask in range(1, full_mask + 1):
        members = [i for i in range(n) if mask & (1 << i)]
        for order in itertools.permutations(members):
            for bits in range(1 << n):
                seqs: list[str] = []
                path: list[tuple[str, str]] = []
                for pos, idx in enumerate(order):
                    sign = "+" if (bits >> idx) & 1 == 0 else "-"
                    seqs.append(
                        reads[idx].seq
                        if sign == "+"
                        else reverse_complement(reads[idx].seq)
                    )
                    path.append((reads[idx].id, sign))

                overlaps = [
                    max_overlap(prev, nxt)
                    for prev, nxt in zip(seqs, seqs[1:])
                ]
                if overlaps and any(overlap == 0 for overlap in overlaps):
                    continue

                assembly = seqs[0]
                for nxt, overlap in zip(seqs[1:], overlaps):
                    assembly += nxt[overlap:]
                candidate = (len(assembly), assembly, tuple(path))
                current = best_by_mask.get(mask)
                if current is None or candidate < current:
                    best_by_mask[mask] = candidate

    best_total: int | None = None
    best_pair: list[tuple[str, list[tuple[str, str]]]] | None = None
    for mask_a in range(1, full_mask):
        if not mask_a & 1:
            continue
        mask_b = full_mask ^ mask_a
        state_a = best_by_mask.get(mask_a)
        state_b = best_by_mask.get(mask_b)
        if state_a is None or state_b is None:
            continue
        total = state_a[0] + state_b[0]
        canonical = sorted(
            (state_a, state_b), key=lambda state: (state[1], state[2])
        )
        pair = [
            (state[1], list(state[2]))
            for state in canonical
        ]
        pair_key = [(assembly, tuple(path)) for assembly, path in pair]
        if best_total is None or (total, pair_key) < (best_total, [
            (assembly, tuple(path)) for assembly, path in best_pair
        ]):
            best_total = total
            best_pair = pair

    assert best_total is not None and best_pair is not None
    return best_total, best_pair


def assert_matches_brute_two_contigs(reads: list[Read]) -> dict:
    """二 contig 结果与独立穷举参考逐项一致。"""

    result = assemble_two_contigs(reads)
    exp_total, exp_pair = brute_force_two_contigs(reads)
    assert isinstance(result, dict)
    assert result["length"] == exp_total
    assert len(result["contigs"]) == 2
    assert "assembly" not in result

    seen: set[str] = set()
    actual_pair = []
    for contig, (exp_assembly, exp_path) in zip(result["contigs"], exp_pair):
        assert set(contig) == {"length", "path", "layout"}
        assert "assembly" not in contig
        actual_path = [
            (item["id"], item["orientation"]) for item in contig["path"]
        ]
        assert actual_path == exp_path
        assert contig["length"] == len(exp_assembly)
        exp_assembly_again, exp_layout = path_sequence_and_layout(reads, exp_path)
        assert exp_assembly_again == exp_assembly
        assert contig["layout"] == exp_layout
        assert len(actual_path) >= 1
        assert not (seen & {rid for rid, _ in actual_path})
        seen.update(rid for rid, _ in actual_path)
        actual_pair.append((exp_assembly, tuple(actual_path)))

    assert seen == {r.id for r in reads}
    assert actual_pair == [
        (assembly, tuple(path)) for assembly, path in exp_pair
    ]
    assert sum(contig["length"] for contig in result["contigs"]) == exp_total
    return result


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


# --------------------------------------------------------------------------
# 恰有两条 contig 的固定场景
# --------------------------------------------------------------------------


def test_two_contigs_zero_overlap_pair_is_isolated():
    reads = parse_payload(payload(("x", "AAAA"), ("y", "CCCC")))
    result = assert_matches_brute_two_contigs(reads)
    assert result["length"] == 8
    assert [[p["id"] for p in c["path"]] for c in result["contigs"]] == [
        ["x"],
        ["y"],
    ]
    assert result["contigs"][0]["layout"] == [
        {"id": "x", "orientation": "+", "start": 0, "end": 4}
    ]
    assert result["contigs"][1]["layout"] == [
        {"id": "y", "orientation": "+", "start": 0, "end": 4}
    ]


def test_two_contigs_two_isolated_overlap_regions():
    # 两组各含内部重叠；穷举参考确认最优分配是这两个重叠区域，而非单读段混排。
    reads = parse_payload(
        payload(
            ("a", "AATT"),
            ("b", "TTAA"),
            ("c", "GCCG"),
            ("d", "CGCC"),
        )
    )
    result = assert_matches_brute_two_contigs(reads)
    assert result["length"] == 11
    assert [[p["id"] for p in c["path"]] for c in result["contigs"]] == [
        ["a", "b"],
        ["d", "c"],
    ]
    assert result["contigs"][0]["layout"] == [
        {"id": "a", "orientation": "+", "start": 0, "end": 4},
        {"id": "b", "orientation": "+", "start": 2, "end": 6},
    ]
    assert result["contigs"][1]["path"] == [
        {"id": "d", "orientation": "+"},
        {"id": "c", "orientation": "+"},
    ]
    assert result["contigs"][1]["layout"] == [
        {"id": "d", "orientation": "+", "start": 0, "end": 4},
        {"id": "c", "orientation": "+", "start": 1, "end": 5},
    ]


def test_two_contigs_reverse_complement_selection_and_path_tie():
    # p/q 是反向互补回文，产生相同 contig 串时必须选择 + 路径；
    # x 与它们没有跨组连接，整体规范为 [x]、[p,q]。
    reads = parse_payload(
        payload(("p", "ACGT"), ("q", "GTAC"), ("x", "AAAA"))
    )
    result = assert_matches_brute_two_contigs(reads)
    assert result["length"] == 10
    assert [[p["id"] for p in c["path"]] for c in result["contigs"]] == [
        ["x"],
        ["p", "q"],
    ]
    assert [
        (p["id"], p["orientation"])
        for p in result["contigs"][1]["path"]
    ] == [("p", "+"), ("q", "+")]
    assert result["contigs"][1]["layout"] == [
        {"id": "p", "orientation": "+", "start": 0, "end": 4},
        {"id": "q", "orientation": "+", "start": 2, "end": 6},
    ]


def test_two_contigs_reverse_complement_selection():
    # b 必须取反向互补 GTTCA，才能与 a 的 AGCGT 以 GT 重叠；
    # x 作为第二条 contig，其他二分总长度更差。
    reads = parse_payload(
        payload(("a", "AGCGT"), ("b", "TGAAC"), ("x", "AAAA"))
    )
    result = assert_matches_brute_two_contigs(reads)
    assert result["length"] == 12
    assert [
        [(p["id"], p["orientation"]) for p in c["path"]]
        for c in result["contigs"]
    ] == [[("x", "+")], [("a", "+"), ("b", "-")]]
    assert result["contigs"][1]["layout"] == [
        {"id": "a", "orientation": "+", "start": 0, "end": 5},
        {"id": "b", "orientation": "-", "start": 3, "end": 8},
    ]


def test_two_contigs_is_not_merely_cut_from_optimal_single_path():
    # 单串最优路径为 a+,b+,c+，但二 contig 最优会重新分配方向和次序：
    # [c+] 与 [b-,a+]；后者不是原单串路径上的一个连续切片。
    reads = parse_payload(
        payload(("a", "TTAA"), ("b", "AATC"), ("c", "CCAA"))
    )
    single = assert_matches_brute(reads)
    assert [
        (item["id"], item["orientation"]) for item in single["path"]
    ] == [("a", "+"), ("b", "+"), ("c", "+")]

    result = assert_matches_brute_two_contigs(reads)
    assert result["length"] == 10
    assert [
        [(p["id"], p["orientation"]) for p in c["path"]]
        for c in result["contigs"]
    ] == [[("c", "+")], [("b", "-"), ("a", "+")]]
    assert result["contigs"][1]["layout"] == [
        {"id": "b", "orientation": "-", "start": 0, "end": 4},
        {"id": "a", "orientation": "+", "start": 2, "end": 6},
    ]


def test_two_contigs_local_coordinates_independent():
    reads = parse_payload(
        payload(
            ("a", "AGCGT"),
            ("b", "GTTCA"),
            ("x", "AATT"),
            ("y", "TTAA"),
        )
    )
    result = assert_matches_brute_two_contigs(reads)
    paths = [
        [(p["id"], p["orientation"]) for p in c["path"]]
        for c in result["contigs"]
    ]
    assert paths == [[("x", "+"), ("y", "+")], [("a", "+"), ("b", "+")]]
    layouts = {
        p["id"]: (p["start"], p["end"])
        for c in result["contigs"]
        for p in c["layout"]
    }
    assert layouts == {
        "x": (0, 4),
        "y": (2, 6),
        "a": (0, 5),
        "b": (3, 8),
    }


def test_two_contigs_too_few_reads_sentinel():
    assert assemble_two_contigs([]) == NO_TWO_CONTIGS
    assert assemble_two_contigs([Read("only", "ACGT")]) == NO_TWO_CONTIGS


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
    assert_matches_brute_two_contigs(reads)


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


def test_cli_two_contigs_payload_mode():
    data = payload(("x", "AAAA"), ("y", "CCCC"))
    data["mode"] = "two_contigs"
    proc = run_cli(data)
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    assert "assembly" not in out
    assert out["length"] == 8
    assert [[p["id"] for p in c["path"]] for c in out["contigs"]] == [
        ["x"],
        ["y"],
    ]


def test_cli_two_contigs_flag():
    raw = json.dumps(payload(("x", "AAAA"), ("y", "CCCC")))
    proc = subprocess.run(
        [sys.executable, "-m", "dna_assembler", "--two-contigs"],
        input=raw,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    assert out["length"] == 8
    assert len(out["contigs"]) == 2


def test_cli_bad_mode():
    data = payload(("x", "AAAA"), ("y", "CCCC"))
    data["mode"] = "three"
    proc = run_cli(data)
    assert proc.returncode == 2
    err = json.loads(proc.stderr)
    assert err["code"] == "bad_mode"
    assert proc.stdout == ""


def test_cli_legacy_too_many_positional_args_error_is_unchanged():
    proc = subprocess.run(
        [sys.executable, "-m", "dna_assembler", "first.json", "second.json"],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert proc.returncode == 2
    err = json.loads(proc.stderr)
    assert err == {"error": "至多接受一个 JSON 文件参数"}
    assert proc.stdout == ""


def test_cli_explicit_single_mode_response_is_unchanged():
    data = payload(("a", "AGCGT"), ("b", "TGAAC"))
    data["mode"] = "single"
    proc = run_cli(data)
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    assert set(out) == {"assembly", "length", "layout", "path", "orientations"}
    assert out["assembly"] == "AGCGTTCA"


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
