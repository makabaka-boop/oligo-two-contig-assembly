"""核心装配逻辑。

规则摘要
--------
* 每条读段恰选择正向或反向互补一次，寻找使最终串最短的排列与方向组合。
* 相邻两条只按"前一条后缀 == 后一条前缀"的最大精确重叠连接。
* 单串平局裁决：先取装配串字典序最小；装配串仍相同（如回文读段）时，再按
  ``[(id, 方向), ...]`` 路径字典序最小。
* 双 contig 模式要求每个非空 contig 内部的相邻读段至少有 1 个碱基重叠；
  两条 contig 之间不拼接、不要求重叠。
* 双 contig 并列时，分别用 ``(装配串, 路径)`` 将两条 contig 规范排序后，再
  比较这一对结果。
* 方向标记：``"+"`` 正向，``"-"`` 反向互补。

算法：带子集的动态规划（Held-Karp 式）。10 条读段时
``10 * 2^10 * 2`` 个状态、每状态至多 20 条转移，完全枚举等价于穷举
所有排列与方向。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_BASES = frozenset("ACGT")

NO_TWO_CONTIGS = "NO_TWO_CONTIGS"


class AssemblyError(ValueError):
    """输入载荷不满足装配约定。"""

    def __init__(self, message: str, code: str = "invalid_input") -> None:
        super().__init__(message)
        self.message = message
        self.code = code


@dataclass(frozen=True)
class Read:
    """一条已校验的输入读段。"""

    id: str
    seq: str


def reverse_complement(seq: str) -> str:
    """返回 ``seq`` 的反向互补串。"""

    return seq.translate(_COMPLEMENT)[::-1]


def max_overlap(left: str, right: str) -> int:
    """``left`` 后缀与 ``right`` 前缀的最大精确重叠长度。

    重叠必须严格短于两条读段本身（完整包含的输入在校验阶段已拒绝）。
    """

    limit = min(len(left), len(right))
    for size in range(limit, 0, -1):
        if left[-size:] == right[:size]:
            return size
    return 0


def parse_payload(payload: object) -> list[Read]:
    """校验 JSON 载荷并返回读段列表（保持上传次序）。"""

    if not isinstance(payload, dict):
        raise AssemblyError("载荷必须是 JSON 对象", "bad_payload")

    if "reads" not in payload:
        raise AssemblyError("载荷必须是包含 reads 字段的 JSON 对象", "bad_payload")
    raw_reads = payload["reads"]
    if not isinstance(raw_reads, list):
        raise AssemblyError("reads 必须是数组", "bad_payload")
    if not 2 <= len(raw_reads) <= 10:
        raise AssemblyError("读段数量必须在 2 至 10 之间", "bad_read_count")

    reads: list[Read] = []
    seen_ids: set[str] = set()
    for index, item in enumerate(raw_reads):
        if not isinstance(item, dict):
            raise AssemblyError(f"第 {index} 条读段必须是对象", "bad_read")
        rid = item.get("id")
        seq = item.get("seq")
        if not isinstance(rid, str) or not rid:
            raise AssemblyError(f"第 {index} 条读段缺少非空字符串 id", "bad_id")
        if not rid.isascii():
            raise AssemblyError(f"id {rid!r} 必须全部为 ASCII 字符", "bad_id")
        if rid in seen_ids:
            raise AssemblyError(f"id {rid!r} 重复", "duplicate_id")
        if not isinstance(seq, str) or not 4 <= len(seq) <= 30:
            raise AssemblyError(
                f"读段 {rid!r} 的 seq 必须是长度 4 至 30 的字符串", "bad_seq"
            )
        if any(ch not in _BASES for ch in seq):
            raise AssemblyError(f"读段 {rid!r} 的 seq 只能含 A/C/G/T", "bad_seq")
        seen_ids.add(rid)
        reads.append(Read(id=rid, seq=seq))

    _reject_containment_or_complement(reads)
    return reads


def requested_contig_count(payload: object) -> Literal[1, 2]:
    """从请求载荷读取装配模式。

    支持两种等价写法：``{"contigs": 2}`` 或
    ``{"mode": "two_contigs"}``；省略时为原单串模式。
    """

    if not isinstance(payload, dict):
        raise AssemblyError("载荷必须是 JSON 对象", "bad_payload")

    mode = payload.get("mode", "single")
    if mode not in ("single", "two_contigs"):
        raise AssemblyError(
            "mode 只能是 'single' 或 'two_contigs'", "bad_mode"
        )
    mode_count: Literal[1, 2] = 2 if mode == "two_contigs" else 1

    if "contigs" not in payload:
        return mode_count

    count = payload["contigs"]
    if (
        isinstance(count, bool)
        or not isinstance(count, int)
        or count not in (1, 2)
    ):
        raise AssemblyError("contigs 只能为 1 或 2", "bad_contig_count")
    if "mode" in payload and count != mode_count:
        raise AssemblyError("mode 与 contigs 指定的装配模式不一致", "bad_mode")
    return 1 if count == 1 else 2


def _reject_containment_or_complement(reads: list[Read]) -> None:
    """拒绝反向互补对，以及任一方向上的完整包含。"""

    for i in range(len(reads)):
        rc_i = reverse_complement(reads[i].seq)
        for j in range(i + 1, len(reads)):
            sj = reads[j].seq
            if reads[i].seq == reverse_complement(sj):
                raise AssemblyError(
                    f"读段 {reads[i].id!r} 与 {reads[j].id!r} 互为反向互补",
                    "reverse_complement_pair",
                )
            # 取任一方为正向，检查两个方向上的完整包含。
            # 利用 u 是 v 的子串 <=> rc(u) 是 rc(v) 的子串，
            # 四项即覆盖 {si, rc_i} × {sj, rc_j} 的全部包含关系。
            si = reads[i].seq
            if si in sj or rc_i in sj or sj in si or sj in rc_i:
                raise AssemblyError(
                    f"读段 {reads[i].id!r} 与 {reads[j].id!r} 在某个方向上完整包含",
                    "contained_read",
                )


def _subset_dp(
    reads: list[Read], *, require_positive_overlap: bool
) -> tuple[
    list[tuple[str, str]],
    dict[tuple[int, int, int, int], int],
    dict[tuple[int, int, int], tuple[int, str, tuple[tuple[str, str], ...]]],
    dict[int, tuple[int, str, tuple[tuple[str, str], ...]]],
]:
    """运行子集 DP。

    ``require_positive_overlap=True`` 时，只允许至少一个碱基的相邻重叠；
    单读段状态不受此限制，因此可作为只含一条读段的 contig。
    """

    n = len(reads)
    variants = [(r.seq, reverse_complement(r.seq)) for r in reads]

    # ov[(i, oi, j, oj)] = 读段 i 取方向 oi 时的后缀与读段 j 取方向 oj
    # 时的前缀的最大精确重叠。
    ov: dict[tuple[int, int, int, int], int] = {}
    for i in range(n):
        for oi in range(2):
            left = variants[i][oi]
            for j in range(n):
                if i == j:
                    continue
                for oj in range(2):
                    ov[(i, oi, j, oj)] = max_overlap(left, variants[j][oj])

    # dp[(mask, last, ori)] -> 最优部分状态：
    #   长度更小者优；长度相同装配串字典序更小者优；再相同路径字典序更小者优。
    # 长度相同时未来转移完全一致，故按此三元组剪枝不影响最优性。
    State = tuple[int, str, tuple[tuple[str, str], ...]]
    dp: dict[tuple[int, int, int], State] = {}
    for i in range(n):
        for ori in range(2):
            dp[(1 << i, i, ori)] = (
                len(variants[i][ori]),
                variants[i][ori],
                ((reads[i].id, "+" if ori == 0 else "-"),),
            )

    for mask in range(1 << n):
        for last in range(n):
            for ori in range(2):
                state = dp.get((mask, last, ori))
                if state is None:
                    continue
                length, text, path = state
                for nxt in range(n):
                    if mask & (1 << nxt):
                        continue
                    for nori in range(2):
                        edge = (last, ori, nxt, nori)
                        overlap = ov[edge]
                        if require_positive_overlap and overlap == 0:
                            continue
                        nxt_seq = variants[nxt][nori]
                        candidate: State = (
                            length + len(nxt_seq) - overlap,
                            text + nxt_seq[overlap:],
                            path
                            + ((reads[nxt].id, "+" if nori == 0 else "-"),),
                        )
                        key = (mask | (1 << nxt), nxt, nori)
                        current = dp.get(key)
                        if current is None or candidate < current:
                            dp[key] = candidate

    best_by_mask: dict[int, State] = {}
    for (mask, _last, _ori), state in dp.items():
        current = best_by_mask.get(mask)
        if current is None or state < current:
            best_by_mask[mask] = state

    return variants, ov, dp, best_by_mask


def _build_layout(
    reads: list[Read],
    variants: list[tuple[str, str]],
    overlaps: dict[tuple[int, int, int, int], int],
    path: tuple[tuple[str, str], ...],
) -> list[dict]:
    """由路径重放每条读段在所属 contig 上的局部起止位置。"""

    index_by_id = {r.id: idx for idx, r in enumerate(reads)}
    layout: list[dict] = []
    cursor = 0
    prev_idx: int | None = None
    prev_ori: int | None = None
    for rid, sign in path:
        idx = index_by_id[rid]
        ori = 0 if sign == "+" else 1
        if prev_idx is None:
            start = 0
        else:
            start = cursor - overlaps[(prev_idx, prev_ori, idx, ori)]
        end = start + len(variants[idx][ori])
        cursor = end
        layout.append(
            {
                "id": rid,
                "orientation": sign,
                "start": start,
                "end": end,
            }
        )
        prev_idx, prev_ori = idx, ori
    return layout


def assemble(
    reads: list[Read],
    contigs: int = 1,
    *,
    mode: str | None = None,
) -> dict | Literal["NO_TWO_CONTIGS"]:
    """按请求的 contig 数量装配。

    可用 ``contigs=2`` 或 ``mode="two_contigs"`` 请求双 contig；两者都省略时，
    默认为保持向后兼容的单串模式。``contigs`` 也可作为第二个位置参数传入。
    """

    if isinstance(contigs, bool) or contigs not in (1, 2):
        raise AssemblyError("contigs 只能为 1 或 2", "bad_contig_count")
    if mode is not None and mode not in ("single", "two_contigs"):
        raise AssemblyError(
            "mode 只能是 'single' 或 'two_contigs'", "bad_mode"
        )

    mode_count = 1 if mode in (None, "single") else 2
    if mode is not None and contigs != mode_count:
        raise AssemblyError("mode 与 contigs 指定的装配模式不一致", "bad_mode")

    if contigs == 1:
        return _assemble_single(reads)
    return assemble_two_contigs(reads)


def _assemble_single(reads: list[Read]) -> dict:
    """求最短单串装配。"""

    n = len(reads)
    if n == 0:
        raise AssemblyError("没有读段可供装配")

    variants, overlaps, dp, _best_by_mask = _subset_dp(
        reads, require_positive_overlap=False
    )
    full_mask = (1 << n) - 1

    best = None
    for last in range(n):
        for ori in range(2):
            state = dp.get((full_mask, last, ori))
            if state is not None and (best is None or state < best):
                best = state
    assert best is not None
    length, assembly, path = best
    layout = _build_layout(reads, variants, overlaps, path)
    orientation_by_id = {rid: sign for rid, sign in path}

    return {
        "assembly": assembly,
        "length": length,
        "layout": layout,
        "path": [
            {"id": rid, "orientation": sign} for rid, sign in path
        ],
        "orientations": {
            r.id: orientation_by_id[r.id] for r in reads
        },
    }


def assemble_two_contigs(
    reads: list[Read],
) -> dict | Literal["NO_TWO_CONTIGS"]:
    """将读段恰好分入两条非空 contig，求长度之和最短的方案。

    成功时只返回各 contig 的长度、局部坐标和方向路径，不返回拼接片段；
    不存在满足正重叠约束的二分时返回 :data:`NO_TWO_CONTIGS`。
    """

    n = len(reads)
    if n < 2:
        raise AssemblyError("恰好两条 contig 至少需要两条读段", "bad_read_count")

    variants, overlaps, _dp, best_by_mask = _subset_dp(
        reads, require_positive_overlap=True
    )
    full_mask = (1 << n) - 1

    best_key: tuple[
        int, tuple[tuple[str, tuple[tuple[str, str], ...]], ...]
    ] | None = None
    best_descriptors: tuple[
        tuple[str, tuple[tuple[str, str], ...], int, list[dict]],
        tuple[str, tuple[tuple[str, str], ...], int, list[dict]],
    ] | None = None

    # 枚举所有非空真子集；A/B 互换会在下方按 (装配串, 路径) 归一化排序。
    for mask in range(1, full_mask):
        state_a = best_by_mask.get(mask)
        state_b = best_by_mask.get(full_mask ^ mask)
        if state_a is None or state_b is None:
            continue

        descriptors = []
        for length, assembly, path in (state_a, state_b):
            descriptors.append(
                (
                    assembly,
                    path,
                    length,
                    _build_layout(reads, variants, overlaps, path),
                )
            )
        descriptors.sort(key=lambda item: (item[0], item[1]))

        total_length = descriptors[0][2] + descriptors[1][2]
        canonical_pair = (
            (descriptors[0][0], descriptors[0][1]),
            (descriptors[1][0], descriptors[1][1]),
        )
        candidate_key = (total_length, canonical_pair)
        if best_key is None or candidate_key < best_key:
            best_key = candidate_key
            best_descriptors = (descriptors[0], descriptors[1])

    if best_descriptors is None or best_key is None:
        return NO_TWO_CONTIGS

    return {
        "total_length": best_key[0],
        "contigs": [
            {
                "length": length,
                "path": [
                    {"id": rid, "orientation": sign}
                    for rid, sign in path
                ],
                "layout": layout,
            }
            for _assembly, path, length, layout in best_descriptors
        ],
    }
