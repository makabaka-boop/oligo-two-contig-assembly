FROM python:3.11-slim

WORKDIR /app

# 仅用标准库，无需安装第三方依赖；pytest 仅在测试时需要。
COPY dna_assembler ./dna_assembler

# 默认从 stdin 读取 JSON；可带一个文件参数。
ENTRYPOINT ["python", "-m", "dna_assembler"]
