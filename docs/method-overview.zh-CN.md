# TagForge：核心原理与适用范围

[English](method-overview.md) | 简体中文

TagForge 将已知设计的条码测序数据转化为分子证据与计数。通用模型为
**身份条码 × 特征条码 × UMI**：身份条码区分细胞、样本等计数对象，特征条码识别检测标签；实际聚类按校正后的身份条码与注释特征名称分组。

[通用流程图](../figures/TagForge_core_principle.zh-CN.page-1.png) ·
[关键机制图](../figures/TagForge_core_principle.zh-CN.page-2.png) ·
[两页 PDF](../figures/TagForge_core_principle.zh-CN.pdf) ·
[自包含 HTML](../figures/TagForge_core_principle.zh-CN.html)

## 从 reads 到分子

1. **描述结构。** 配置片段角色、来源 read、长度、坐标和／或 linker。条码及 UMI 可以由多个片段组成，分布在 R1/R2。
2. **提取片段。** 使用固定位置或 linker 定位。配置 `methods: [linker, fixed]` 时，linker 提取失败后回到原始 read 的固定坐标；后续白名单校正失败不会再次触发提取回退。
3. **校正身份。** 各条码片段分别使用自己的白名单校正，再按顺序拼接。组合条码无需枚举所有完整组合的白名单。
4. **归并分子。** 在校正后的身份条码–特征组合内聚合原始 UMI 支持数，再调用 UMI-tools 聚类。一个校正后的 UMI 簇计为一个分子。
5. **保留证据与输出。** 分子记录保留代表 UMI、支持 reads 数和原始 UMI 类型数，形成稀疏计数、降采样诊断及样本／批次报告。

## 核心判定

### 条码校正

原始窗口精确命中白名单时直接接受。否则在配置允许的窗口偏移与 Hamming 错配范围内搜索，并按以下四元组依次比较候选：

`(I[偏移非零] + I[存在错配], |偏移| + 错配数, 错配数, |偏移|)`

即优先修改类型更少、总幅度更小、错配更少、偏移更小的候选。唯一最优条码身份被接受；多个不同身份并列最优时按歧义排除；无合法候选也被排除。多个位置指向相同的最优身份不属于身份歧义。该搜索处理位置偏移与碱基替换，不等同于任意 indel 比对。

### UMI 聚类

默认使用 UMI-tools `directional`，距离为 1。候选有向边 `u → v` 要求距离不超过 1，且满足 `n(u) >= 2*n(v) - 1`，其中 n 为 read 支持数；最终簇和代表 UMI 由 UMI-tools 确定。

在仅含两个相差 1 bp 的 UMI 的组内，支持数 8:1 可归并为一个分子，8:7 则保留为两个分子。图中两个例子属于独立组，所有结构、序列和计数均为示意。

矩阵统计 UMI 簇数，而非支持 reads 总数。例如支持数为 9 和 3 的两条分子记录，对计数矩阵贡献 **2 个分子**。

### 降采样诊断

对全深度分子中支持数为 n 的分子，在比例 p 下执行 `Binomial(n,p)` 抽样。当前每个比例一次带种子抽样，并提供 reads、分子和单例计数的解析期望；各深度不重跑片段提取、UMI 聚类或实验特异过滤。

- 测序饱和度：`1 − 分子数 / reads 数`。
- 非单例分子比例：`至少 2 条 reads 支持的分子数 / 分子数`。

所选诊断比例使观测非单例比例最大，平局时选择较低比例。主要分析使用全深度结果；该选择不直接构成测序深度建议。

## 工程能力

- **规模化处理与计算优化：** 分块读取、Parquet、DuckDB 磁盘聚合、稀疏输出，以及精确匹配快速路径、邻域查找和校正缓存。
- **按独立分组并行：** 完整身份条码–特征组合进入同一分桶；独立 UMI 连通分量可分别处理，聚类边界不被拆断。
- **断点恢复与结果校验：** 复用经校验的已完成阶段；提取和校正支持恢复已提交批次。输入／配置指纹、上游状态与输出校验共同决定是否复用；FASTQ 默认使用元数据指纹，可配置完整 SHA256。
- **批量运行与过程追踪：** Slurm 数组任务、实时进度、阶段日志、分段质控及样本／批次 HTML 和 Excel 报告。

这些是已实现的设计能力，不代表对其他工具的普遍速度或内存领先。实际资源需求仍受数据分布和 worker 处理的最大数据范围影响。

## 适用场景与扩展边界

| 场景 | 计数关系 | 适用范围 |
| --- | --- | --- |
| 抗体标签与组合条码 | 身份条码 × 抗体／特征标签 | 通用流程的主要应用。 |
| Sample Tags | 细胞 × 样本标签 | 适配方向；后续样本归属需使用对应分析方法。 |
| CRISPR Guide Capture | 细胞 × sgRNA | sgRNA 身份可直接解析时的适配方向；扰动归属分析另行完成。 |
| Antigen Capture | 细胞 × 抗原标签 | 适配方向；结合特异性解释另行完成。 |
| 启动子／报告基因条码 | 身份条码 × 设计标签 | 带 UMI 设计的适配方向；活性或表达解释取决于具体实验。 |

扩展应用要求输入为流程支持的双端 reads，片段和 UMI 可提取，特征身份可确定，且开启白名单校正的片段具备相应白名单。各文库均需结构配置和数据验证；通用计数模型不等于已经验证所有实验平台或完整下游流程。

## PI-seq 实验的两类定制文库

[第三页示意图](../figures/TagForge_core_principle.zh-CN.page-3.png)并列展示两类文库。它们复用通用核心，在 UMI 聚类后执行各自的归属规则。PB 是两类文库共享的身份条码；FB 是抗体特征条码；CB 是细胞条码。

| 文库 | 需要解决的问题 | 判定范围 | 主要结果 |
| --- | --- | --- | --- |
| PI-seq | 一条抗体标签分子归属于哪个 PB？ | 在同一 `(FB_name, corrected_umi)` 内比较各 PB | 过滤后的 PB × FB 分子明细与计数。 |
| PB–CB | 一个 PB 归属于哪个细胞 CB？ | 对每个 PB 汇总各 CB 的 read 支持 | PB → CB 映射、每个 CB 的 PB 数及分布。 |

### PI-seq：保留有明确 PB 归属的分子

开发示例中，R1 的五段 PB 各为 8 nt，经逐片段白名单校正后拼接；R2 为 `U1 (8 nt)–FB (12 nt)–U2 (8 nt)`，UMI 按配置顺序拼接。Linker 定位和固定位置回退复用通用流程。该布局是示例，具体实验可配置。

先在 `PB × FB_name` 内完成 UMI 聚类，再按 `(FB_name, corrected_umi)` 检查跨 PB 的归属。仅出现一个 PB 的组直接保留；出现多个 PB 时，保留支持 reads 最多的 PB，前提是其支持比例 **严格大于** `pi_seq.dominance_threshold`（默认 0.8），否则排除该整组。比例为该 PB 的支持 reads 除以组内所有 PB 的支持 reads。

例如 9:1 的 read 支持保留占 90% 的 PB；8:2 的支持恰好为 80%，默认排除。最高 read 支持并列时，程序按 PB 字典序选择候选，再检查阈值。过滤依据是 read 支持，最终计数仍以保留的 UMI 分子为单位。

主要分子结果为 `03_molecules/{sample}.molecule_detail.rmMP.parquet`，供 PI-seq 计数与下游逐 PB 分析使用；多 PB、优势归属和移除比例记录于 `pi_seq_qc.tsv`。

### PB–CB：建立 PB 到细胞的映射

该文库读取 `PB × CB × UMI`。CB 以实测序列保留，要求 `barcode2.sequence_only: true`、关闭其白名单校正且不配置 CB 白名单；PB 片段仍可按配置使用白名单校正。

默认启用可关闭的观测 CB 保守纠错：全样本中仅使用不含 N、至少 5 条 reads 的候选母序列；低支持 CB 与其相差 1 bp，或仅含 1 个可匹配的 N，且母序列支持严格更多、至少为子序列的 10 倍时，只有唯一符合条件的母序列才建立别名。多候选或不符合条件的序列不强制归并。参数由 `pb_cb` 配置。

CB 别名处理后，相同 `(PB, CB, corrected_umi)` 的分子记录重新合并，避免 CB 归并后重复计数。对每个 PB 按 CB 汇总分子数与支持 reads，用 **read 支持比例** 选择优势 CB，比例必须 **严格大于** `pb_cb.dominance_threshold`（默认 0.8）。最高 read 支持并列时先按 CB 字典序选择候选，再检查阈值。默认配置下，9:1 保留映射，8:2 不分配 CB。

主要输出在 `04_matrix/`：`pb_cb_mapping.tsv.gz` 保留原始 PB–CB 配对支持；`pb_cb_map.tsv.gz` 保存通过阈值的映射及分子数、read 支持和优势比例；`cb_pb_counts.tsv.gz` 与 `cb_pb_count_distribution.tsv` 描述细胞的 PB 数；`cb_observed_correction.tsv.gz` 记录 CB 别名。

### 通过共享 PB 连接两类结果

PI-seq 提供 `PB × 抗体特征` 的分子证据，PB–CB 提供 `PB → 细胞` 的映射。两类输出可按共享且定义一致的 PB 关联，作为细胞层面抗体特征整合的输入；关联和下游解释另行完成。多个 PB 可以映射到同一 CB，每个通过筛选的 PB 仅保留一个 CB 归属。

## 实现依据

- [结构配置](../src/tagforge/config.py)
- [片段提取](../src/tagforge/extract.py)与[条码校正](../src/tagforge/barcode_correct.py)
- [UMI 聚类与聚合](../src/tagforge/umi_correct.py)
- [PI-seq 归属过滤](../src/tagforge/pi_seq.py)与[PB–CB 映射](../src/tagforge/matrix.py)
- [稀疏计数](../src/tagforge/sparse_matrix.py)与[降采样](../src/tagforge/downsample.py)
- [检查点与依赖校验](../src/tagforge/checkpoint.py)
- [Slurm 集成](../src/tagforge/slurm.py)

Linker 匹配使用 Cutadapt 与精确匹配路径，UMI 聚类使用 UMI-tools；共同工作流、数据组织和执行逻辑由 TagForge 实现。
