"""核心装配逻辑。

规则摘要
--------
* 每条读段恰选择正向或反向互补一次，寻找使最终串最短的排列与方向组合。
* 相邻两条只按"前一条后缀 == 后一条前缀"的最大精确重叠连接。
* 平局裁决：先取装配串字典序最小；装配串仍相同（如回文读段）时，再按
  ``[(id, 方向), ...]`` 路径字典序最小。
* 方向标记：``"+"`` 正向，``"-"`` 反向互补。

算法：带子集的动态规划（Held-Karp 式）。10 条读段时
``10 * 2^10 * 2`` 个状态、每状态至多 20 条转移，完全枚举等价于穷举
所有排列与方向。
"""

from __future__ import annotations

from dataclasses import dataclass

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_BASES = frozenset("ACGT")


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

    if not isinstance(payload, dict) or "reads" not in payload:
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


def assemble(reads: list[Read]) -> dict:
    """求最短装配，返回装配串与每条读段的方向、起止位置。

    位置采用 0 基、左闭右开区间（Python 切片语义），即
    ``assembly[start:end]`` 恰为该读段所选方向上的序列。
    """

    n = len(reads)
    if n == 0:
        raise AssemblyError("没有读段可供装配")

    variants = [(r.seq, reverse_complement(r.seq)) for r in reads]

    # ov[(i, oi, j, oj)] = 读段 i 取方向 oi 时的后缀与读段 j 取方向 oj
    # 时的前缀的最大精确重叠；tail 为随之需要追加到当前串末尾的片段。
    ov: dict[tuple[int, int, int, int], int] = {}
    tail: dict[tuple[int, int, int, int], str] = {}
    for i in range(n):
        for oi in range(2):
            left = variants[i][oi]
            for j in range(n):
                if i == j:
                    continue
                for oj in range(2):
                    nxt_seq = variants[j][oj]
                    size = max_overlap(left, nxt_seq)
                    edge = (i, oi, j, oj)
                    ov[edge] = size
                    tail[edge] = nxt_seq[size:]

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

    full_mask = (1 << n) - 1
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
                    new_mask = mask | (1 << nxt)
                    for nori in range(2):
                        edge = (last, ori, nxt, nori)
                        shift = len(variants[nxt][nori]) - ov[edge]
                        candidate: State = (
                            length + shift,
                            text + tail[edge],
                            path
                            + ((reads[nxt].id, "+" if nori == 0 else "-"),),
                        )
                        key = (new_mask, nxt, nori)
                        current = dp.get(key)
                        if current is None or candidate < current:
                            dp[key] = candidate

    best: State | None = None
    for last in range(n):
        for ori in range(2):
            state = dp.get((full_mask, last, ori))
            if state is not None and (best is None or state < best):
                best = state
    assert best is not None
    length, assembly, path = best

    # 由路径重放每条读段在装配串上的起止位置。
    index_by_id = {r.id: idx for idx, r in enumerate(reads)}
    orientation_by_id = {rid: sign for rid, sign in path}
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
            start = cursor - ov[(prev_idx, prev_ori, idx, ori)]
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
