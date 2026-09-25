# DNA 短读段最短串装配器

每条 DNA 读段可能来自正链或反向互补链。本工具在**全部排列 × 每条读段两个方向**的
组合空间中，寻找最终装配串最短的方案；相邻两条只用「前一条后缀 == 后一条前缀」的
最大精确重叠连接。

默认将全部读段装配为一条串；当文件中混入两段独立 DNA 区域时，也可以请求**恰好两条
contig**：每条 contig 至少包含一条读段，contig 内部相邻读段必须有至少一个碱基的
精确重叠，两条 contig 之间不拼接、不要求重叠。

## 输入

JSON 对象，2 至 10 条读段；`id` 为唯一 ASCII 非空字符串，`seq` 为 4 至 30 个
`A/C/G/T` 碱基：

```json
{
  "reads": [
    {"id": "r1", "seq": "AGCGT"},
    {"id": "r2", "seq": "TGAAC"}
  ]
}
```

请求恰好两条 contig 时，可使用下列任一等价字段；省略时保持原单串请求：

```json
{"contigs": 2, "reads": []}
{"mode": "two_contigs", "reads": []}
```

以下数据会**整份拒绝**（stderr 输出 JSON 错误，退出码 2）：

- 数量不在 2–10、id 为空/非 ASCII/重复、序列长度或碱基不合法、JSON 损坏；
- 两条读段**互为反向互补**；
- 任一方向上一条读段被另一条**完整包含**（正向或反向互补后是其子串）；
- `contigs` 不是 `1`/`2`，或 `mode` 与 `contigs` 冲突。

## 使用

```bash
# stdin / 文件参数，默认单串
python -m dna_assembler < request.json
python -m dna_assembler request.json

# Compose 的 assembler 服务
docker compose build
docker compose run --rm -T assembler < request.json
```

## 单串输出

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

## 恰好两条 contig 输出

成功时不输出 contig 拼接片段，只输出各 contig 的局部坐标、方向路径与长度：

```json
{
  "total_length": 12,
  "contigs": [
    {
      "length": 8,
      "path": [
        {"id": "r1", "orientation": "+"},
        {"id": "r2", "orientation": "-"}
      ],
      "layout": [
        {"id": "r1", "orientation": "+", "start": 0, "end": 5},
        {"id": "r2", "orientation": "-", "start": 3, "end": 8}
      ]
    },
    {
      "length": 4,
      "path": [{"id": "r3", "orientation": "+"}],
      "layout": [{"id": "r3", "orientation": "+", "start": 0, "end": 4}]
    }
  ]
}
```

每个 `start`/`end` 都是所属 contig 内的局部 0 基半开区间。若不存在满足条件的
二分，stdout 返回纯文本哨兵：

```text
NO_TWO_CONTIGS
```

## 平局裁决

1. 装配总长度最短；
2. 单串模式取**装配串字典序最小**；双 contig 模式先分别取各 contig 的最优结果，
   再将两条 contig 按 `(装配串, [(id, orientation), ...])` 规范排序，比较这一对
   规范结果；
3. 装配串也相同（典型情形：反向互补回文读段，方向不改变序列）时，再按
   `[(id, orientation), ...]` 路径字典序裁决（`"+" < "-"`）。

双 contig 的候选来自全部非空子集二分，而不是把原单串最优路径切开，因此结果不依赖
上传顺序。

## 算法

子集动态规划（Held-Karp 式）：状态为「已用读段集合 × 最后一条读段 × 其方向」，
每条有向边的权值为两条定向读段间的最大精确重叠。单串状态数
`2^n · n · 2`；双 contig 模式只保留正重叠边，再枚举所有互补子集对。平局键直接采用
`(长度, 装配串, 路径)` 三元组字典序——相同长度的部分串追加同一未来序列后相对字典序
不变，因此提前剪枝不影响最终裁决。

## 测试

```bash
python -m pytest
```

测试中包含独立的全枚举参考实现：单串模式对 5 条以内的全部 `n! · 2^n` 排列方向
组合逐项比对；双 contig 模式额外穷举所有二分，并覆盖：零重叠隔离、反向互补、
不同串与同串路径并列、contig 对规范排序、局部坐标和无合法二分时的
`NO_TWO_CONTIGS`。
