# DNA 短读段最短串装配器

每条 DNA 读段可能来自正链或反向互补链。本工具在**全部排列 × 每条读段两个方向**的
组合空间中，寻找最终装配串最短的方案；相邻两条只用「前一条后缀 == 后一条前缀」的
最大精确重叠连接。

## 输入

JSON 对象，2 至 10 条读段；`id` 为唯一 ASCII 非空字符串，`seq` 为 4 至 30 个
`A/C/G/T` 碱基。默认是单串装配；可选 `mode` 为 `"single"`（默认）或
`"two_contigs"`：

```json
{
  "reads": [
    {"id": "r1", "seq": "AGCGT"},
    {"id": "r2", "seq": "TGAAC"}
  ]
}
```

以下数据会**整份拒绝**（stderr 输出 JSON 错误，退出码 2）：

- 数量不在 2–10、id 为空/非 ASCII/重复、序列长度或碱基不合法、JSON 损坏；
- 两条读段**互为反向互补**；
- 任一方向上一条读段被另一条**完整包含**（正向或反向互补后是其子串）。

## 使用

```bash
# stdin / 文件参数
python -m dna_assembler < request.json
python -m dna_assembler request.json

# 恰有两条 contig：也可在请求 JSON 中设置 "mode": "two_contigs"
python -m dna_assembler --two-contigs < request.json

# Compose 的 assembler 服务
docker compose build
docker compose run --rm -T assembler < request.json
```

## 输出

```json
{
  "assembly": "AGCGTTCA",
  "length": 8,
  "layout": [
    {"id": "r1", "orientation": "+", "start": 0, "end": 5},
    {"id": "r2", "orientation": "-", "start": 3, "end": 8}
  ],
  "path": [
    {"id": "r1", "orientation": "+"},
    {"id": "r2", "orientation": "-"}
  ],
  "orientations": {"r1": "+", "r2": "-"}
}
```

- `orientation`：`"+"` 正向，`"-"` 反向互补；
- `start`/`end`：0 基、左闭右开区间，`assembly[start:end]` 恰为该读段所选
  方向上的序列；
- `layout`/`path` 按装配次序排列；`orientations` 以 id 为键，便于按 id 查询。

## 二 contig 模式

在请求对象中加入 `"mode": "two_contigs"`，或使用 CLI 参数
`--two-contigs`/`--mode two_contigs`。该模式把全部读段恰好分配到两条非空
contig；每条 contig 内部相邻读段必须有至少一个碱基的精确重叠，但两条 contig
之间不连接，因此不会产生零重叠“假接缝”。若无法二分，输出裸 JSON 字符串
`"NO_TWO_CONTIGS"`。

合法结果不输出任何 contig 拼接片段，只返回总长度、每条 contig 的局部路径和
局部起止位置；局部坐标各自从 0 开始：

```json
{
  "length": 8,
  "contigs": [
    {
      "length": 4,
      "path": [{"id": "r1", "orientation": "+"}],
      "layout": [{"id": "r1", "orientation": "+", "start": 0, "end": 4}]
    },
    {
      "length": 4,
      "path": [{"id": "r2", "orientation": "+"}],
      "layout": [{"id": "r2", "orientation": "+", "start": 0, "end": 4}]
    }
  ]
}
```

优化目标是两条 contig 长度之和最小。平局时，先在每条 contig 内按
`(装配串, [(id, orientation), ...])` 取规范结果，再把两条 contig 按该二元组
排序，最后比较这一对规范结果；结果因此与上传次序无关。

## 平局裁决

1. 装配串长度最短；
2. 长度相同，取**装配串字典序最小**；
3. 装配串也相同（典型情形：反向互补回文读段，方向不改变序列）时，再按
   `[(id, orientation), ...]` 路径字典序裁决（`"+" < "-"`）。

## 算法

子集动态规划（Held-Karp 式）：状态为「已用读段集合 × 最后一条读段 × 其方向」，
每条有向边的权值为两条定向读段间的最大精确重叠。状态数
`2^n · n · 2`，n=10 时约 2 万状态、18 万次转移，与穷举所有排列、方向等价，
运行时间在百毫秒内。平局键直接采用 `(长度, 装配串, 路径)` 三元组字典序——
相同长度的部分串追加同一未来序列后相对字典序不变，因此提前剪枝不影响最终裁决。
二 contig 模式沿用同一子集 DP，但只允许重叠长度大于 0 的转移；再枚举互补的两个
非空子集，选择规范后总长最短的一对。

## 测试

```bash
python -m pytest
```

测试中包含一个**独立的全枚举参考实现**，对 5 条以内的全部 `n! · 2^n` 排列方向
组合逐项比对，并覆盖：无重叠拼接、反向互补回文与普通反转回文、不同串并列与同串
路径并列、反向互补对/包含拒绝、随机合法输入、CLI 端到端。二 contig 模式还会
穷举所有非空二分及各组排列、方向，覆盖零重叠隔离、局部坐标、跨 contig 不输出
拼接片段以及二 contig 平局规范。
