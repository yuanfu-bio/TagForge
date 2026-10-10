# TagForge 2.0.0 实施记录

更新日期：2026-10-09。当前源码版本：2.0.0，本地修复、原计划验收与最终安装包验证已完成，尚未正式发布。

下文按批次保留实施时的状态；当前结论以文末“2.0.0 最终本地验收与下一步”为准。历史候选与首轮构建不能替代最终产物。

## 第一批：M1 已完成

| 任务 | 实现结果 |
| --- | --- |
| F01 阶段指纹与依赖失效 | JSON 完成 manifest 记录规范化配置、直接输入指纹、上游代次、输出清单与行数。白名单及 annotation 校验内容；FASTQ 默认校验元数据，可配置 SHA256。参数或输入变化后拒绝旧结果，下游独立命令检查前置阶段，汇总排除失效样本。 |
| F02 阶段发布与恢复 | 同一样本互斥运行。阶段在临时目录计算，使用发布日志恢复多产物替换，最后原子发布完成 manifest。提取和纠错保留批次恢复状态直至发布成功，核对已提交 fragment 和行数，恢复时丢弃未提交 fragment 并重建提取 preview。其他阶段中断后安全整步重算。 |
| F03 selected 引用 | `selection.json` 明确记录 ratio、repeat、seed、criterion、指标及引用产物。全深度直接引用本次 raw 产物，报告统一解析 selected 引用；重算移走旧 optimal 与可能遮蔽新计数的旧 TSV。选择结果及引用都参与完成检查。 |

新增逻辑集中在 `src/tagforge/checkpoint.py` 和 `src/tagforge/selected.py`。README 和 CHANGELOG 已同步当前恢复、迁移与产物行为。原有工作区修改保留，尚未提交或发布。

## 验证结果

完整回归通过 42 项测试，其中新增 18 项检查点与选择结果测试。覆盖白名单/annotation/UMI/seed 变化、旧完成标记、最终统计写入失败、未提交 fragment、发布中断、同一样本互斥、缺失 selected fragment、零分子，以及 partial/full-depth 切换。测试和真实流程均在 `/tmp` 的独立副本运行。

四个 PI-seq fixture 共 200,000 对 reads 完整运行成功。以下计数以及全部分子明细、PI 过滤后明细、矩阵 pairs、矩阵汇总、纠错统计、PI QC 和降采样指标均与修复前基线一致。

| 样本 | 输入 reads | 有效 reads | 去重分子 | PI 过滤后分子 | PB 数 | feature 数 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 100-1-N1 | 50,000 | 47,281 | 46,423 | 46,413 | 46,037 | 102 |
| 100-2-SLE1 | 50,000 | 47,231 | 46,198 | 46,184 | 45,577 | 102 |
| 100-3-SLE2 | 50,000 | 47,604 | 46,430 | 46,422 | 45,731 | 103 |
| 100-5-N2 | 50,000 | 47,425 | 46,454 | 46,450 | 45,866 | 102 |

相同配置再次运行成功，24 个阶段的 generation ID 全部保持不变。最终完整运行耗时 189.87 秒，峰值 RSS 850,456 KiB；这是当前环境的运行记录，尚未建立可比的性能提升结论。

机器可读证据保存在 [tests/validation/m1.json](../validation/m1.json)，包括基线与验证源码 SHA256、fixture 配置摘要、运行依赖、计数和重复运行结果。基线取自实施前包含已有未提交修改的工作区，不能只用 Git HEAD 复原，需结合源码摘要识别。

完整基线和最终验证副本已持久归档到 `docs/archive/testing/m1-20261009/{baseline,validation}`，包含源码、配置、输入 fixture、全部输出、报告和运行日志；`archive-manifest.json` 保存逐文件 SHA256，归档时已与原始副本核对。另将最终输出和批次报告复制到原 fixture 的 `tests/fixtures/pi_seq_dev/{02_output,00_report}`，便于直接查看。归档索引见 [tests/validation/README.md](../../validation.md)。这些生成目录由现有 `.gitignore` 排除，不会加入源码版本控制；不再依赖 `/tmp` 保留结果。

前两个 fixture 为固定 seed 的随机子集，后两个为前缀子集。本次验证证明这些输入的计数保持一致及故障恢复行为，不证明纠错算法准确率、完整样本分布或百万级性能。Parquet 检查覆盖 footer、schema 和 fragment 元数据，不扫描所有数据页校验摘要；恢复保证针对进程中断，未测试断电持久性。

## 升级与下一批

旧版只写程序版本的 `.done` 不会复用。升级后应先运行完整流程重建完成状态，再执行独立下游阶段；现有白名单/annotation/输入改变也会按依赖使结果失效。

M0 已建立四个 fixture 的可追踪基线；合成真值、完整过滤链指标定义及开发变更的发布审阅仍需继续。M2 已实现并完成一轮验证：

1. F04：quick-test 与正式纠错共用 read 判定；输出逐 read QC 与互斥失败原因，覆盖部分片段失败、无效 UMI、缺失注释和零有效 reads。只有显式配置 `quick_test.min_valid_rate` / `max_n_base_rate` 时才阻止流程。
2. F05：集中表来源识别、fragment 校验和稳定 schema；空表可读，计数采用整数类型。完整多段 feature 的注释缺口由纠错阶段直接产出，trace 开/关一致，旧目录 trace 可读取。
3. F06：实现禁用降采样时的全深度报告、禁用恢复时全部重算；拒绝未知配置项、负距离、错误白名单长度和输出名称冲突。PB–CB 生成独立映射报告与汇总，alias 后按最终 PB–CB–UMI 合并分子并保留 reads 总和。

M2 独立快照通过 **63 项测试**；四个真实 fixture 全部运行成功，分子明细、矩阵 pairs、矩阵汇总、PI QC、降采样指标及原有纠错统计与 M1 一致。新增失败计数与有效 reads 之和等于输入 reads。再次运行复用全部 24 个阶段。完整运行耗时 148.16 秒，GNU time 记录峰值 RSS 868,732 KiB；记录用于验收，性能归因将在 M3 使用相同基准核对。

证据见 [tests/validation/m2.json](../validation/m2.json)，源码快照、完整输出、日志和比较结果持久保存在 `docs/archive/testing/m2-20261009/`。M1 日期归档保持不变。

M3 已完成：exact 快路、分批矩阵汇总/导出、资源设置、一次批次汇总发布、UMI 预算、降采样分析语义和 linker 策略评估均已实现并验证。

### M3 第二轮实现与证据

最新冻结快照 `docs/archive/testing/m3-20261009/validation-r2` 已通过 **81 项完整回归**（16.929 秒），源码摘要见 `source-r2-sha256.json`，日志见 `tests-r2.log`。首轮出现旧测试 mock 未接受新 linker 参数，已修正并保留首轮失败日志。

- UMI：两个后端均记录 `peak_group_umis` / `peak_component_umis` / `peak_batch_umis`，支持可选 `performance.umi_max_group_umis` / `umi_max_component_umis`；超限说明 barcode/feature，不拆连通纠错关系。components 逐个生成，减少整组复制。SQLite 的 `total_raw_umis` 改为聚合后的唯一 UMI 数，进度区分 groups 与 components。
- 降采样：ratio 按 repeat 均值选择；代表 repeat 取最接近均值者，再按 repeat 编号打破平局。保存所有 repeats、均值/样本 SD、独立解析期望计数和 criterion。HTML 误差条及 Excel 说明明确这是固定全深度分子组上的 reads 抽样诊断；不推导推荐测序深度。
- Linker：8 个合成真值场景通过。四个真实 fixture 共 200,000 对 reads 的新策略新增 841 有效 reads、丢失 3、改变 21 条有效标签，所有 865 条差异均保存逐条 trace 与坐标参考证据。两条丢失 read 的原固定坐标精确命中白名单。按 A01 保留 `legacy_shortest` 默认，`length_position` 显式可选，不宣称真实样本准确率已证明；决策见 `strategy/decision.json`。
- 已生成 1,000,000 对真实配对 reads 的前缀基准，输入 SHA256、原始来源及前缀采样局限记录在 `million/fixture/fixture.json`。正在顺序测量 M1 基线与新快照在合成、四个 fixture 和百万对 reads 上的耗时、进程树 RSS、临时磁盘峰值；采样峰值不保证捕获瞬间尖峰，进程树 RSS 会重复计入共享页。
- 版本来源已统一到 `tagforge.__version__`，setup shim 与测试不再硬编码版本；仍为 0.1.14。新增线程/chunk/trace/固定 seed 一致性测试尚待下一轮执行。M3 尚未完成性能与计数验收。

机器可读进度见 [tests/validation/m3.json](../validation/m3.json)。所有实验输出仍保存在项目内的忽略目录，M1/M2 日期归档不变。

### M3 最终验收与 M4 候选验证

六次顺序基准均已完成。单次完整 CLI 运行（包含 quick-test）的耗时为：合成样本 8.02 → 7.52 秒，四个真实 fixture 336.94 → 260.59 秒，百万对 reads 556.55 → 528.85 秒。百万级 GNU time 峰值 RSS 为 1,477,692 → 1,027,056 KiB；采样进程树 RSS 为 5.09 → 4.69 GB，临时逻辑磁盘峰值为 167.29 → 165.26 MB。共享机器、单次测量及 RSS 共享页重复计数等限制记录在 `performance/summary.json`，不将这些结果解释为稳定倍数提升。

四个 fixture 及百万对 reads 的有效 reads、分子明细、PI 过滤后明细、矩阵 pairs、矩阵汇总、PI QC 和降采样指标逐项一致。百万对 reads 得到 947,722 有效 reads、705,912 分子、PI 过滤后 704,259 分子；逐表对称差异均为 0，见 `performance/million/count-comparison.json`。

M4 开发包预检已完成 wheel/sdist 构建及真实数据排除检查。独立 Python 3.9.23 与 3.11.15 环境的已安装 wheel 均通过三类文库、多样本、配置开关、三个导出格式及两种 Slurm 脚本检查；Python 3.9 完整回归 **82 项通过**（18.191 秒），因此保留最低 Python 3.9 声明。UMI-tools 旧构建脚本的 pip 安装问题已复现并通过兼容构建依赖解决，README 记录可重复安装步骤。

真实 Chrome 渲染 7 个 HTML、openpyxl 打开 13 个 workbook；仅保存 Arial PDF 预览。打开检查发现并修复默认 Excel 样式、冻结窗格 OOXML 和大小写工作表名冲突；报告算法版本已单独提升，避免为样式重算 FASTQ。现已进入 `2.0.0rc1`，候选源码、构建及完整回归保存在 `docs/archive/testing/m4-20261009/`。

新增 `init-example` 让已安装 wheel 独立生成公开合成项目；MIGRATION.md 记录配置、完成状态、表格式、selected 指标和资源约定。候选包安装/报告复验、最终文档与发布审查仍未完成，当前不能按 2.0.0 正式发布；本轮持续目标仍包含全部已发现问题的修复与优化。

### rc2 验收记录（历史批次）

本轮补齐坐标、长度、线程、chunk、seed 和计数字段的整数校验，拒绝小数、布尔值及非有限数，防止 `int()` 静默截断；PB–CB alias 比率也拒绝 NaN/Infinity。真实 PI-seq FASTQ 加入 Git 忽略规则，冻结源码与 wheel/sdist 均包含仓库已有完整 MIT LICENSE。

新增测试中一条布尔断言位置错误已修正。rc2 冻结源码在 Python 3.11 和 3.9 上均通过 **85 项完整回归**（18.871 / 18.268 秒）。首轮两套回归共用输出目录触发样本互斥，保留原日志；Python 3.9 使用逐文件摘要一致的独立副本重新通过。

**Python 3.9.23 / 3.11.15 / 3.13.12 已安装 rc2 wheel 均通过 31 条 CLI 检查**，覆盖三类文库、多样本、配置开关、三个矩阵导出格式及两种 Slurm 脚本。已安装 `init-config` 与仓库示例一致，覆盖拒绝覆盖现有配置。13 个 Excel 和 7 个 HTML 实际打开检查通过；仅生成 Arial PDF 预览。构建检查确认源码与 wheel 一致、完整许可证存在、真实 FASTQ 不进入发布包。证据见 [m4-rc2.json](../validation/m4-rc2.json)，完整副本、运行日志、构建包和预览保存在 `docs/archive/testing/m4-20261009/`。

rc2 验收时需要收尾的发布问题（已在最终验收关闭）：

1. 对照 M0–M4 原验收条件逐项建立最终证据清单，统一审阅实施前已有修改、合成真值覆盖、指标分母与迁移拒绝路径。
2. 整理最终 CHANGELOG、迁移说明和 Release 说明，修正历史阶段记录与当前状态混用。
3. 最终审查通过后将版本改为 `2.0.0`，准备本地最终源码、构建产物及 SHA256；不自动推送、打 tag 或上传。

已知边界继续明确保留：`length_position` 在真实差异中仍有丢失 read，默认保留 `legacy_shortest`；降采样是固定全深度分子组的诊断，不用于宣称最佳测序深度；超大 UMI 组目前提供预算和超限诊断，完整算法重构作为后续优化。UMI-tools 的旧 pip 构建问题使用 README 已验证步骤处理。上述边界应进入最终发布说明。


### 2.0.0 最终本地验收与下一步

原 M0–M4、F01–F07、A01/A02 与 11 项发布检查已逐项审阅，30 个条款的映射见 [RELEASE_READINESS.md](../../RELEASE_READINESS.md)，最终机器可读证据见 [release.json](../../../tests/validation/release.json)。实施前已有 5 个修改模块的原始基线和审阅依据保留在独立记录中。

最终审查补齐 quick-test duplication 百分数单位、UMI 不适用配置/冲突 segment role 拒绝、PB–CB 等支持数 reciprocal alias 防护、迁移命令 threads 校验、集成测试隔离，以及开发安装 metadata 与源码版本同步。

最后复现并修复矩阵导出未检查配置依赖的问题：将 `correction_umi.max_distance` 从 1 改为 0 后，旧包仍导出 6 个分子，而重算为 7。现在 export 检查 matrix 及上游完成状态，拒绝过期计数并保留已有导出文件；重算后正确导出 7。实际 CLI 拒绝/重算/新计数导出及完整回归通过，复现记录位于 M4 归档的 `stale-export-repro/summary.json`。

- Python 3.9.23、3.11.15、3.13.12 各通过 **91 项完整回归**，无跳过。3.9/3.11 使用独立冻结源码，3.13 使用已安装最终 wheel。
- 三套环境各完成 31 条安装包 CLI 检查，加上最终 60 条独立阶段、指标、迁移、阈值、并发和导出检查，共 **153 条 CLI 检查**通过。
- 13 个 Excel 与 13 个 HTML 实际打开通过；仅保存 13 个 PDF 预览，PDF 字体条目核对使用用户指定 Arial。
- 快照的 62 个文件与工作区逐文件摘要一致，回归后输入未变；wheel 全部运行模块与快照一致。wheel/sdist 含完整 MIT LICENSE，sdist 仅含两个公开合成 FASTQ，真实 FASTQ 与实验结果未分发。
- M1 历史归档 552 个文件重新核验 SHA256；M3 四个 fixture 及百万对 reads 的逐表计数比较保留完整证据。性能数字属于 M3 快照的单次共享机器测量，不作为最终快照重新测量的结果。

最终安装包和摘要位于 `docs/archive/testing/m4-20261009/dist-release-r3/`，源码为 `release-source-r3/`、摘要为 `release-source-r3-sha256.json`。`dist-release/` 与 `dist-release-r2/` 缺少最后的导出修复，保留历史但不能用于发布。最终回归、三套安装流程、补充 CLI 与报告使用 `*-r3` 归档名称，所有完整结果保存在项目内。

当前没有遗留的已发现必修项；本地 2.0.0 准备完成。正式发布下一步是审阅与已验证快照一致的源码变更，形成对应提交，再按实际发布渠道准备 `v2.0.0` 和 Release/安装包；本轮未提交、推送、打 tag 或上传。

后续按原计划保留研究范围：优先改善超大 UMI 组资源使用，再验证 linker 定位策略、大白名单邻居索引及质量值纠错。全深度模拟与测序深度推荐需要独立真值及真实样本证据，不阻塞本次 2.0.0。
