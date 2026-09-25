"""命令行入口：从 stdin（或文件参数）读取 JSON，输出装配结果 JSON。

用法::

    python -m dna_assembler < input.json
    python -m dna_assembler input.json

正常结果打印到 stdout；输入问题以 JSON 形式打印到 stderr，退出码 2。
"""

from __future__ import annotations

import json
import sys

from .assembler import (
    AssemblyError,
    assemble,
    assemble_two_contigs,
    parse_payload,
)

TWO_CONTIGS_MODES = frozenset({"two_contigs", "two-contigs"})
SUPPORTED_MODES = frozenset({"single", "single_contig", *TWO_CONTIGS_MODES})


class _Options:
    def __init__(
        self,
        file: str | None = None,
        mode: str | None = None,
        two_contigs: bool = False,
    ) -> None:
        self.file = file
        self.mode = mode
        self.two_contigs = two_contigs


def _parse_argv(argv: list[str]) -> _Options | None:
    """解析参数；不带新模式参数时保持旧 CLI 的行为。"""

    if not argv:
        return _Options()
    if argv == ["--help"] or argv == ["-h"]:
        print("用法: python -m dna_assembler [--two-contigs | --mode MODE] [FILE]")
        raise SystemExit(0)

    options = _Options()
    positionals: list[str] = []
    index = 0
    while index < len(argv):
        arg = argv[index]
        if arg == "--two-contigs":
            options.two_contigs = True
        elif arg == "--mode":
            index += 1
            if index >= len(argv):
                raise AssemblyError("--mode 需要一个参数值", "bad_mode")
            options.mode = argv[index]
        elif arg.startswith("--mode="):
            options.mode = arg.split("=", 1)[1]
        elif arg != "-" and arg.startswith("-"):
            raise AssemblyError(f"未知参数 {arg}", "bad_arguments")
        else:
            positionals.append(arg)
        index += 1

    if len(positionals) > 1:
        print(
            json.dumps({"error": "至多接受一个 JSON 文件参数"}, ensure_ascii=False),
            file=sys.stderr,
        )
        return None

    options.file = positionals[0] if positionals else None
    if options.mode is not None and options.mode not in SUPPORTED_MODES:
        raise AssemblyError("mode 必须是 single 或 two_contigs", "bad_mode")
    return options


def _requested_mode(payload: object, cli_mode: str | None) -> str:
    mode = "single"
    if isinstance(payload, dict) and "mode" in payload:
        raw_mode = payload["mode"]
        if not isinstance(raw_mode, str) or raw_mode not in SUPPORTED_MODES:
            raise AssemblyError(
                "mode 必须是 single 或 two_contigs",
                "bad_mode",
            )
        mode = raw_mode
    if cli_mode is not None:
        mode = cli_mode
    return mode


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv

    try:
        args = _parse_argv(argv)
        if args is None:
            return 2
        if args.file:
            with open(args.file, "r", encoding="utf-8") as fh:
                raw = fh.read()
        else:
            raw = sys.stdin.read()
        try:
            payload = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise AssemblyError(f"JSON 解析失败: {exc}", "invalid_json") from exc
        mode = _requested_mode(
            payload,
            "two_contigs" if args.two_contigs else args.mode,
        )
        reads = parse_payload(payload)
        if mode in TWO_CONTIGS_MODES:
            result = assemble_two_contigs(reads)
        else:
            result = assemble(reads)
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

    json.dump(result, sys.stdout, ensure_ascii=False, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
