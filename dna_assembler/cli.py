"""命令行入口：从 stdin（或文件参数）读取 JSON，输出装配结果 JSON。

用法::

    python -m dna_assembler < input.json
    python -m dna_assembler input.json

正常结果打印到 stdout；双 contig 无合法二分时 stdout 输出
``NO_TWO_CONTIGS``；输入问题以 JSON 形式打印到 stderr，退出码 2。
"""

from __future__ import annotations

import json
import sys

from .assembler import (
    NO_TWO_CONTIGS,
    AssemblyError,
    assemble,
    parse_payload,
    requested_contig_count,
)


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) > 1:
        print(
            json.dumps({"error": "至多接受一个 JSON 文件参数"}, ensure_ascii=False),
            file=sys.stderr,
        )
        return 2

    try:
        if argv:
            with open(argv[0], "r", encoding="utf-8") as fh:
                raw = fh.read()
        else:
            raw = sys.stdin.read()
        try:
            payload = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise AssemblyError(f"JSON 解析失败: {exc}", "invalid_json") from exc
        reads = parse_payload(payload)
        contig_count = requested_contig_count(payload)
        result = assemble(reads, contigs=contig_count)
    except AssemblyError as exc:
        print(
            json.dumps({"error": exc.message, "code": exc.code}, ensure_ascii=False),
            file=sys.stderr,
        )
        return 2
    except OSError as exc:
        print(
            json.dumps({"error": f"读取输入失败: {exc}", "code": "io_error"},
                       ensure_ascii=False),
            file=sys.stderr,
        )
        return 2

    if result == NO_TWO_CONTIGS:
        sys.stdout.write(NO_TWO_CONTIGS)
        sys.stdout.write("\n")
        return 0

    json.dump(result, sys.stdout, ensure_ascii=False, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
