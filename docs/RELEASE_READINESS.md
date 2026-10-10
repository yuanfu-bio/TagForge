# TagForge 2.0.0 本地发布验收

本记录按原 [2.0.0 计划](archive/2.0.0/TagForge_2.0.0_plan.md) 核对范围。最终状态、源码摘要、
安装包 SHA256 和日志索引以 [release.json](../tests/validation/release.json) 为准。
验收产物保存于 `docs/archive/testing/m4-20261009/`；正式上传、推送和 tag
不在本次本地准备范围内。发布时必须使用与已验证快照一致的源码提交和安装包。

## 当前验收环境

日常修改和发布验收均使用 `/data/home/yfu/miniconda3/envs/tagforge/bin/python`（Python 3.11）。下表的三版本结果是历史证据，不是当前发布门槛；历史测试环境已移除。实验快照和日志归档到 `docs/archive/testing/`。当前包仍需在该环境中重建并完成与发布相关的验证。

## 最新流程目录与结果入口

样本输出已按处理流程重新编号，checkpoints 合并到 `00_logs/checkpoints/` 并使用编号清单，有效 reads 移入 `02_corrected/`，分子结果集中在 `03_molecules/`。三种文库的主结果在 HTML 和 Excel 中直接标识。现有 Python 3.11 的 12 个相关场景通过；真实四样本完整重跑后，24 张逻辑表、QC 和降采样值与备份一致，全部 24 个阶段有效，5 份工作簿的原有统计值保持一致（耗时、速度随重跑更新）。证据见 [本次验收](archive/testing/pipeline-layout-20261010-183901/summary.json)。两个样本的 quick-test 分段 linker 诊断略有变化，总体 QC 和正式数据一致，但仍需在发布前排查可重复性。当前发布安装包仍需重建，以下历史验收不代表新目录已打包。

## 输出文件精简

当前源码采用精简产物契约，四样本迁移核对通过；仅在当前 Python 3.11 环境检查 17 个直接相关场景。证据见 [精简输出验收](archive/testing/compact-output-20261010-005722/summary.json)。此前完整回归和跨版本记录属于对应历史快照；发布前需重新构建当前安装包。

## 批次工作簿合并

当前源码只生成一份批次 Excel，区分全深度和选择计数，并保留 QC、降采样及注释详情。四样本报告核对通过，证据见 [合并验收](archive/testing/merged-batch-report-20261010-003646/summary.json)。本次已完成三版本各 100 项回归；日常修改不再默认执行跨版本回归，后续发布验收也仅使用已有 Python 3.11 环境，暂不创建跨版本测试环境。当前安装包仍需重建。

## 批次报告目录统一

最新源码将全部批次报告集中在 `{output_dir}/00_report/`，通过 Python 3.11 的 97 项完整回归及四样本实际报告刷新。报告内容、数据和 24 个阶段状态核对通过，见 [目录变更验收](archive/testing/report-layout-20261010-002219/summary.json)。此变更未构建新的安装包；之前的跨版本、已安装 wheel 验收保留为对应快照的历史证据，发布前需重建并验证最新包。

## Snappy 默认压缩变更

Snappy 变更快照已通过 Python 3.9/3.11/3.13 各 96 项回归；四个真实样本全部重跑，92 个 Parquet 文件均为 Snappy，24 张逻辑表与上一轮一致。证据见 [本次运行摘要](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838/summary.json)。原 `dist-release-r3` 保留为历史验收包，不包含此次变更；正式发布前需重新构建和验证安装包。

## 原计划逐项证据

| 条款                      | 实现与验收证据                                                                                                                                                                                                              |
| ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| M0.1 已有修改与基线       | M1 baseline 的 22 个模块逐文件 SHA256 已重新核对；`initial-dirty-inventory.json` 与 `initial-dirty.patch` 区分 Git HEAD 和实施前已有的 5 个源码修改。最终处理见下表。                                                   |
| M0.2 合成及四个真实基线   | M1/M2 独立基线与结果归档；M3 合成、四样本及百万级的逻辑计数比较、耗时/RSS/临时空间记录。                                                                                                                                    |
| M0.3 合成真值             | `test_qc`：多段、部分失败、歧义、N/非 DNA；`test_linker_strategy` 和 8 场景 JSON：重复 linker/位置/歧义/回退；`test_core`：移位、UMI 邻居、PI 多 PB；`test_workflow`：alias 后共享 UMI；支持数比率为 1 的独立真值。 |
| M0.4 过滤链与分母         | `evaluate_extracted_read` 的互斥失败原因；quick/formal 逐 read 对照，M2 四样本失败计数之和核对，最终 CLI 零有效结果与阈值检查。                                                                                           |
| M0.5 指标约定             | `calculate_metrics`、独立支持数公式、穷举抽样期望测试；最终 CLI 审查核对 TSV/JSON/workbook/HTML 数值。quick duplication 与降采样均为百分数，其他 QC rate 为小数。                                                         |
| F01 配置/输入与下游失效   | 白名单、同元数据 annotation 内容、坐标、UMI 参数、seed、FASTQ SHA256、线程和报告版本用例；独立命令及 export 拒绝失效上游，拒绝时保留现有导出，重算后导出新计数；汇总排除失效样本。                                          |
| F02 发布/恢复/互斥        | 提取及纠错最后统计失败、未提交 fragment/preview、发布第二产物中断、损坏 fragment、同一样本锁用例；实际两套 CLI 共用样本曾被锁拒绝；最终 CLI 两个不同样本并发后保留两行汇总。                                                |
| F03 selected 一致性       | subset → full → subset、seed 变化、删除 fragment、零记录、旧 optimal/TSV 清理；JSON 引用与 point、分子行数、矩阵总和核对；最终 CLI 支持 reads 总和核对。                                                                  |
| F04 QC                    | 共用完整 read 判定，部分 PB/FB/UMI 失败、annotation 缺失、barcode 歧义逐条预期断言；默认诊断与显式阈值经实际 CLI 检查。                                                                                                     |
| F05 数据集                | 单文件、多 fragment、空 schema、schema 不符、trace 开/关、多段注释缺口、旧目录 trace、旧 TSV 迁移计数；不完整当前 Parquet 拒绝回退。                                                                                        |
| F06 流程与配置            | 三类文库完整及独立阶段、重复执行、禁用 downsample/resume、线程/chunk/trace、多样本与并发；PB–CB 去 alias 后重新按最终 PB–CB–UMI 合并并保留 reads，strict dominance 真值；未知/无效/不生效设置明确拒绝。                  |
| F07 保守性能              | 原位 exact 命中不调用 distance search；矩阵`fetchmany` 与 DuckDB 资源；每次 invocation 汇总一次；两 UMI 后端的完整/分量/预算对照；六次顺序基准及全部逻辑行对称差异为 0。                                                  |
| A01 条件策略              | 八个独立 tag 真值场景、200,000 对真实 reads、865 条差异完整审阅与摘要；按原条件保留`legacy_shortest`，`length_position` 显式可选。                                                                                      |
| A02 分析语义              | 每个 ratio 单次采样、按 observed 指标选 ratio、tie 规则、观测值/期望计数分别存储；穷举 16 个 read 子集证明期望计数；不把期望计数比值作为期望比例。                                                                |
| M4.1 候选到正式源码       | rc1/rc2 安装包验证、rc3 全量与补充 CLI 审查均留档；最终源码版本 2.0.0，经最终回归/包检查记录。                                                                                                                              |
| M4.2 单一版本来源         | `tagforge.__version__`、setuptools dynamic version、setup shim；源码、CLI、日志、已安装 wheel metadata 和构建版本核对。                                                                                                   |
| M4.3 版本约定             | README 改为 semantic versioning；CHANGELOG、MIGRATION、RELEASE_NOTES 解释完成状态、配置、表和指标契约变化。                                                                                                                 |
| M4.4 支持范围             | Linux/Python 3.9、3.11、3.13；完整回归和独立 wheel 环境；Cutadapt/UMI-tools/Arrow/DuckDB/NumPy/Pandas/PyYAML/GNU sort 版本由最终记录保存。                                                                                  |
| M4.5 迁移                 | legacy flat 配置实际跑通；旧`.done` 被下游 CLI 拒绝并重建；旧 TSV 显式读取、当前缺失/损坏无静默回退；文档标出必须重算与能力限制。                                                                                         |
| C01–C06 正确性与回归门禁 | 上述 F01–F07、三类文库、四真实 fixture、故障/配置/空数据/并发、独立真值、线程/chunk/固定 seed 的证据组合。                                                                                                                 |
| C07 已安装包门禁          | 三个隔离环境使用构建 wheel，各完成 31 条 CLI 检查，包括三文库、多样本、开关、三个导出格式、两种 Slurm 的`bash -n`；补充独立命令/QC/迁移/并发审查。                                                                        |
| C08 文档与打开检查        | 安装包`init-config` 与仓库 YAML 规范化相等，覆盖覆盖保护；真实 Chrome/openpyxl 打开 sample/batch/downsample/quick-test HTML 与 workbook；预览仅 Arial PDF。                                                               |
| C09 历史迁移与拒绝        | 同 M4.5/F03/F05；CLI 错误信息与旧状态重建、当前数据损坏检查。                                                                                                                                                               |
| C10 分发边界              | wheel 无 FASTQ，sdist 仅九记录的两个公开合成 FASTQ；真实 PI-seq 与实验结果排除且本地 raw 被 Git 忽略；`init-example` 可从已安装包生成。                                                                                   |
| C11 发布资料与产物        | 最终 CHANGELOG/MIGRATION/RELEASE_NOTES、本记录、源码 SHA256、测试/CLI/报告日志、wheel/sdist 和 SHA256SUMS 由最终 release.json 索引。                                                                                        |

## 实施前已有修改的处理

| 原模块     | 最终处理与审阅依据                                                                                                                                                                |
| ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| cli        | 保留 fixture 路径迁移、report 自动更新诊断及按全部已完成样本汇总；增加 disabled 条件、前置依赖、显式 QC 阈值与 migration threads 校验。已安装独立命令、旧状态拒绝及并发汇总验证。 |
| downsample | 保留 read-level saturation 与 non-singleton 区分；按单次观测指标选择 ratio，并用明确 selected 引用。数学真值、全深度切换与跨格式数值验证。                              |
| io_utils   | 原 schema-version helper 保留为 legacy 辅助，当前正式完成判定统一进入 checkpoint；version-only marker 不再复用。旧状态用例及 CLI 迁移验证。                                       |
| pipeline   | 保留 PB–CB 诊断和批次重扫描方向，替换为按配置的阶段计划、互斥及 journal 发布；三文库/disabled/多样本完整和独立流程验证。                                                         |
| reports    | 保留原指标区分与布局，改为 current manifests/selected 来源、单次采样观测值、PB–CB mapping 报告与正确 Excel OOXML。实际打开及跨格式数值验证。                                            |

旧 examples 数据迁移到 `tests/fixtures/small_fastq`；公开生成器与 fixture 的配对 records
保留，配置和教程引用更新。真实数据及历史日期归档保留在项目内的忽略目录。

## 公开配置生效审阅

| 配置组                            | 使用方与约定                                                                                                                                                                |
| --------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| project/samples                   | 路径、输入和 output 根目录由 config 解析；name 是描述性标签，不改变计数/阶段身份。样本名和模板受校验。                                                                      |
| barcode1/barcode2/segments/linker | 提取方法、坐标、linker 策略、白名单与纠错进入 extract/correct 身份和共用判定；旧分组及 flat alias 有兼容测试；segment target 冲突拒绝。                                     |
| annotation/sequence_only          | 注释列、重复名称与 full-feature 映射作用于共用有效性和 QC；sequence-only 直接保留观测序列，PB–CB 禁用 CB 白名单纠错。                                                      |
| umi/correction_umi                | method/distance 作用于 UMI-tools，UMI 片段只提取和完整性校验；不生效的 whitelist/barcode correction 激活拒绝。                                                              |
| downsample/quick_test             | 开关、ratios/seed、前缀 read 数及显式 QC 阈值都有对应运行/判定；threads/chunk 的输出一致性有真值比较。                                                              |
| performance                       | workers/backend/sort/sqlite/matrix 预算进入实际执行资源；batch/chunk/preview/compression 控制缓冲、诊断或编码；group/component 预算拒绝连通范围超限。资源变化无需改变计数。 |
| resume/output                     | reuse/batch resume/overwrite/fingerprint 模式进入阶段路径；trace 关闭保留紧凑 annotation QC。                                                                               |
| library/pi_seq/pb_cb              | 文库类型决定 stage plan，dominance/alias 参数作用于各自过滤与映射；reads/molecules 分别记录。                                                                               |

## 证据边界和后续研究

- 两个真实 fixture 是随机子集，两个是前缀；百万基准是一个真实样本的前缀。逻辑计数一致与性能记录不能证明总体生物学准确率。
- 新 linker 策略有 841 rescued、3 lost、21 assignment changes；保留默认是 A01 明确允许的决策。
- 断点恢复测试针对进程中断；FASTQ 默认元数据校验、Parquet footer/schema/fragment 校验的范围在迁移说明中明确，不扩张为断电或全部 payload 摘要保证。
- RSS 为 GNU time 与进程树采样两种口径，临时空间为逻辑字节采样；共享页及短暂峰值限制见性能 summary。
- 原计划明确排除的深度模拟、新深度推荐、quality-aware/大白名单索引、超大 UMI 完整解法和 dominan
- ce 新默认仍作为后续研究，未缩减原发布修复范围。
