# 验证结果与归档

完整流程可在隔离副本中运行，但验收完成后应将结果和日志保存到项目目录，不只保留 `/tmp` 副本。每批次使用独立归档名称，保留对应源码、配置、输入和基线；不要用下一次运行覆盖历史验收快照。

## 当前验证环境与结果位置（2026-10-10）

日常修改和 2.0.0 发布验收均使用 `/data/home/yfu/miniconda3/envs/tagforge/bin/python`（Python 3.11），只运行本次变更相关的检查。暂不执行跨版本验收或创建新的测试环境。`tests/validation/` 保留验证工具和当前索引；12 批历史实验结果移至 [docs/archive/testing](archive/testing/)，三套历史 Python 环境及专用 Conda 缓存已删除。当前真实数据及结果仍在 `tests/fixtures/pi_seq_dev/`。快照内部的原执行路径按 [迁移表](archive/testing/relocation.json) 查找，不改写历史文件内容。详情见 [清理记录](archive/validation/testing-cleanup-20261010.json)。

## 下采样图片宽度减半（2026-10-10）

HTML 图与 PDF 画布从 1000 改为 500 像素宽，高度仍为 540。逐项比对四个样本页面，仅宽度设置变动，曲线、数值标注、Arial 字体及数据均不变。详见 [验收记录](archive/testing/downsample-half-width-20261010-202545/summary.json)。

## 下采样图内添加数值（2026-10-10）

上方面板内标记样本名、全量（ratio = 1）Saturation 和 Non-singleton molecule fraction 百分比，保留两位小数。100-1-N1 分别为 1.82% 和 1.80%。原来的两个面板、四条曲线和全部数据不变；5 项相关测试通过，浏览器确认数字完整显示且字体为 Arial。详见 [验收记录](archive/testing/downsample-numeric-labels-20261010-202408/summary.json) 与 [PDF 预览](archive/testing/downsample-numeric-labels-20261010-202408/100-1-N1.downsample.pdf)。

## 下采样页面保留原来的双图（2026-10-10）

根据最终要求，恢复原图上方 Saturation / Non-singleton molecule fraction 和下方 UMI Types / UMI detected once 两个面板。四条曲线、配色、坐标轴范围、面板位置及高度均与原图一致，只删除图外标题和小字说明，图内加样本名，统一 Arial。5 项相关测试与浏览器检查通过，四样本曲线数据与原图逐项一致，批次工作簿全部单元格不变。详见 [验收记录](archive/testing/two-panel-downsample-figure-20261010-201817/summary.json) 和 [PDF 预览](archive/testing/two-panel-downsample-figure-20261010-201817/100-1-N1.downsample.pdf)。

## 先前单图方案（已由双图方案替代，2026-10-10）

下采样 HTML 只展示一张图，图内标记样本名，保留 Saturation 和 Non-singleton molecule fraction 两条曲线。去掉页面标题与小字说明，统一 Arial，纵轴按数据自动缩放。5 项相关测试通过；浏览器验证实际字体为 Arial、两条曲线与图例完整，并生成 [PDF 预览](archive/testing/single-downsample-figure-20261010-200716/100-1-N1.downsample.pdf)。四样本仅更新下采样和报告，16 个上游清单、全部数值及批次工作簿单元格不变，详见 [验收记录](archive/testing/single-downsample-figure-20261010-200716/summary.json)。

## 下采样指标去重（2026-10-10）

下采样只保留 `sequencing_saturation` 和 `non_singleton_molecule_fraction` 两个比例指标，移除原始指标表、曲线表、选中结果及报告中的重复 `duplication_ratio`。批次 HTML 表头改为 `Sat.` 和 `Non-singleton (%)`，曲线使用 Saturation 名称；`Main results` 和 Downsample 小字说明继续隐藏。

105 项测试通过，四样本仅重算下采样和报告 8 个阶段，16 个上游清单原字节不变。除删除的重复字段外，所有原始指标和曲线数值完全一致；批次 Excel 全部单元格不变，始终含 `Counts_raw` 和 `Counts_final`。9 份 HTML 脚本语法检查通过，旧输出、源码、测试日志和核对摘要见 [验收记录](archive/testing/two-downsample-metrics-20261010-195556/summary.json)。

## 下采样改为单次采样（2026-10-10）

每个比例仅采样一次，移除重复次数配置、重复编号、重复选择规则、均值/标准差和曲线误差线。选择记录使用 schema 2；每个比例的观测指标和解析期望计数分别保留。当前 Python 3.11 环境的 105 项测试通过，包含 CLI 集成运行；9 份更新的 HTML 内嵌脚本语法检查通过，5 份工作簿 XML、Arial 样式和重复字段移除检查通过。

四个真实样本仅重算下采样与报告共 8 个阶段，16 个上游阶段清单保持原字节。各比例的观测计数保持不变，最终比例均为 1。批次工作簿除 `Downsample` 外的工作表全部单元格不变；`Downsample` 从 13 列精简为 10 列。原报告和下采样输出已备份，详见 [验收摘要](archive/testing/single-sampling-20261010-193032/summary.json)。`latest_run.json` 已更新，之前的完整运行记录保留为基线。

## 最新流程目录与结果入口（2026-10-10）

样本目录按流程统一为 `00_logs/`、`01_extracted/`、`02_corrected/`、`03_molecules/`、`04_matrix/`、`05_downsample/`、`06_report/`。有效 reads 移至纠错目录；阶段清单放在 `00_logs/checkpoints/` 并编号为 `01_extract.manifest.json` 至 `06_report.manifest.json`。HTML 和 Excel 直接标识标准、PI-seq、PB–CB 各自的全深度主结果，不增加重复数据文件。

现有 Python 3.11 环境的 12 个直接相关测试通过，覆盖三种文库入口、阶段恢复、纠错及分子读写、导出和 CLI。真实四样本 200,000 对 reads、24 个阶段全部重新计算，含 quick-test 耗时 256.50 秒。24 张逻辑 Parquet 表的 schema 和逐行内容与旧结果完全一致，QC、降采样及选择值保持一致；24 个阶段状态有效，5 份工作簿结构、样式和原有统计值核对通过（耗时、速度随重跑更新），主结果链接可达。旧结果已完整备份。源码、日志、输入摘要、新旧结果和核对记录见 [本次验收](archive/testing/pipeline-layout-20261010-183901/summary.json)。另有两个样本的 quick-test 分段 linker 诊断略有变化，总体有效 reads、分子数及纠错 QC 一致，原因尚未确认，已记为发布前排查项。当前安装包仍待发布前重建；以下条目属于各自历史快照。

## 阶段清单改名（2026-10-10）

该次快照使用 `01_checkpoint/{step}.manifest.json`。四样本的 24 份旧 JSON `.done` 清单已原样改名，SHA256、修改时间和 generation 全部不变，24 个阶段仍有效，未重跑真实数据。旧名称的有效 JSON 可读取并在阶段复用时迁移；损坏的新清单不会回退到旧文件。指定 Python 3.11 环境的 7 个相关场景通过，见 [改名验收](archive/testing/manifest-naming-20261010-173845/summary.json)。安装包仍待发布前重建。

## 移除独立临时目录（2026-10-10）

流程不再创建 `08_tmp`。阶段暂存迁入 `01_checkpoint/.work`，成功后清理，失败时保留续跑数据；旧暂存中的纠错批次和发布日志迁移后可继续恢复。UMI 排序、SQLite 和 PI 过滤的临时文件用完即清理。已有 TagForge Python 3.11 环境的 6 个相关场景通过；当前四样本 24 个阶段仍有效，未重跑真实数据，未更改当前结果运行记录。验证与日志见 [本次摘要](archive/testing/remove-tmp-directory-20261010-012400/summary.json)。安装包仍待发布前重建。

## 前次 Snappy 四样本完整运行（2.0.0，2026-10-10）

该次快照将 `pi_seq_dev` 输出更新为默认 Snappy 的 2.0.0：200,000 对 reads、24 个阶段全部重算，完整 CLI 耗时 265.03 秒（包含 prefix quick-test）。有效 reads、分子明细、PI 过滤、矩阵、PI QC 和降采样指标与上一轮 2.0.0 逐项一致；完成状态、计数和 selected 引用核对通过。7 个 Excel、13 个 HTML 实际打开通过，预览只保存 Arial PDF。

| 样本 | 有效 reads | 去重分子 | PI 过滤后分子 |
| --- | ---: | ---: | ---: |
| 100-1-N1 | 47,281 | 46,423 | 46,413 |
| 100-2-SLE1 | 47,231 | 46,198 | 46,184 |
| 100-3-SLE2 | 47,604 | 46,430 | 46,422 |
| 100-5-N2 | 47,425 | 46,454 | 46,450 |

- 当前结果：[02_output](../tests/fixtures/pi_seq_dev/02_output)，[latest_run.json](../tests/fixtures/pi_seq_dev/02_output/latest_run.json)。
- 当前批次报告：[HTML](../tests/fixtures/pi_seq_dev/02_output/00_report/TagForge_batch_report.html)、[Excel](../tests/fixtures/pi_seq_dev/02_output/00_report/TagForge_batch_report.xlsx)（已合并汇总表）。
- 完整独立归档：[pi-seq-dev-2.0.0-20261009-220831](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838)，[运行摘要](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838/summary.json)。
- [计数对比](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838/count-comparison.json)、[产物验收](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838/validation.json)、[运行日志](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838/run.log)、[资源记录](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838/resources.json)。
- [逐文件 SHA256](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838/archive-manifest.json)、[报告打开与 PDF 预览](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838/report-open)。

下文保留各历史阶段的证据，M1 链接指向当时独立快照。新归档保存本次源码、完整输入、旧/新输出及日志；复制的完成 manifest 保留原运行路径，重放归档输入时需生成新的完成状态。

92 个 Parquet 文件全部使用 Snappy；24 张逻辑表逐行、schema 对比一致。Python 3.9/3.11/3.13 各通过 96 项完整回归。详见 [压缩验收](archive/testing/pi-seq-dev-2.0.0-snappy-20261009-234838/compression-verification.json)。旧 `dist-release-r3` 不包含此次压缩变更，发布前需重新构建并验证安装包。

## 输出文件精简（2026-10-10）

当前四样本已迁移到精简输出，每个样本 TSV（含 gzip）从 17 份减少到 6 份。进度、纠错 QC 和选择记录使用 JSON；下采样曲线合并 observed 均值/SD 与 expected 计数；逐 read 调试导出默认关闭。24 张逻辑表、全部 QC 和曲线值与原输出一致，24 个阶段有效，批次工作簿所有单元格不变，5 份改动的 Excel 实际打开通过。

仅在当前 Python 3.11 环境核对 17 个直接相关场景，修复后仅重跑失败场景。迁移复用四个提取阶段，重算纠错及其下游共 20 个阶段；未重跑无变化的 prefix quick-test，既有聚合 QC 保留。详见 [迁移摘要](archive/testing/compact-output-20261010-005722/summary.json)、[原输出](archive/testing/compact-output-20261010-005722/previous)、[当前输出快照](archive/testing/compact-output-20261010-005722/results/02_output)。旧完整运行及跨版本证据保持原样，当前安装包仍需发布前重建。

## 批次工作簿合并（2026-10-10）

当前批次只保留 `TagForge_batch_report.xlsx` 和 HTML。工作簿包含 `Summary`、`Counts_raw`、`Downsample` 和 `Warnings`；仅在存在选择比例小于 1 的样本时增加 `Counts_final`。四个真实样本计数、QC、200 条注释缺失记录与原工作簿逐项核对通过；220 个样本文件保持原身份，HTML 字节不变，样本阶段未重算。旧工作簿已备份，见 [合并验收](archive/testing/merged-batch-report-20261010-003646/summary.json)。

本次已完成 Python 3.9/3.11/3.13 各 100 项回归以及并发样本更新核对；这些是本次保留的证据，不是日常修改的默认门槛。根据用户反馈，日常变更只在当前开发环境执行与变更相关的测试，后续发布验收也使用已有 Python 3.11 环境，不再默认创建多版本环境。安装包仍需在发布时重建验证。

## 批次报告目录统一（2026-10-10）

所有当前批次报告集中到 `02_output/00_report/`。通过实际 `tagforge report` 刷新四样本报告，未重算样本阶段；新旧 HTML 字节及三份 XLSX 内部内容一致，除四份追加运行日志外的 220 个样本文件保持原 SHA256、大小和修改时间，24 个阶段仍有效。旧报告位置和手工 `03_export/` 副本已清理，原文件备份在本次归档中。

最新目录变更通过 Python 3.11 的 97 项完整回归；前述 Python 3.9/3.13 的 96 项记录对应 Snappy 变更快照。目录变更的安装包仍需重建与验证。详见 [验收摘要](archive/testing/report-layout-20261010-002219/summary.json)、[回归日志](archive/testing/report-layout-20261010-002219/tests.log)、[备份](archive/testing/report-layout-20261010-002219/previous)。最终结果索引仍处于方案讨论阶段，尚未实现。

## M1（2026-10-09）

- 验证摘要：[m1.json](archive/validation/m1.json)。
- 当时 fixture 输出：[pi_seq_dev/02_output](archive/testing/m1-20261009/validation/tests/fixtures/pi_seq_dev/02_output)，包含各样本的明细、矩阵、降采样和报告。
- 当时汇总表：[00_summary.xlsx](archive/testing/m1-20261009/validation/tests/fixtures/pi_seq_dev/02_output/00_summary.xlsx)。
- 当时批次报告：[HTML](archive/testing/m1-20261009/validation/tests/fixtures/pi_seq_dev/00_report/TagForge_batch_report.html)、[Excel](archive/testing/m1-20261009/validation/tests/fixtures/pi_seq_dev/00_report/TagForge_batch_report.xlsx)。
- 完整基线快照：[baseline](archive/testing/m1-20261009/baseline)。
- 完整最终验证快照：[validation](archive/testing/m1-20261009/validation)。
- 最终运行日志：[final-real-run.log](archive/testing/m1-20261009/validation/final-real-run.log)。
- 耗时与内存：[final-real-run.time.txt](archive/testing/m1-20261009/validation/final-real-run.time.txt)。
- 42 项测试日志：[final-tests.txt](archive/testing/m1-20261009/validation/final-tests.txt)。
- 基线比较：[comparison.json](archive/testing/m1-20261009/validation/comparison.json)。
- 重复运行验证：[rerun-check.json](archive/testing/m1-20261009/validation/rerun-check.json)。
- 归档文件 SHA256：[archive-manifest.json](archive/testing/m1-20261009/archive-manifest.json)。

原 fixture 的输出副本便于查看，后续运行可以更新；日期归档保留本次验收结果。归档时逐文件核对 SHA256，保留原始日志、manifest 和文件时间，未修改历史输入身份。完成 manifest 中记录的 `/tmp` 路径属于当时运行的审计信息；从新位置重新运行时路径指纹会改变，需要重新计算，不能直接当作新位置的可复用检查点。

`02_output/` 和 `00_report/` 已由项目 `.gitignore` 排除。这些结果保存在本地项目目录，尚未上传远端或纳入发布包。

## M2（2026-10-09）

- 验证摘要：[m2.json](archive/validation/m2.json)，63 项测试及四个真实 fixture 首轮验收通过。
- 源码及输入快照：[validation](archive/testing/m2-20261009/validation)。
- 完整输出：[PI-seq 02_output](archive/testing/m2-20261009/validation/tests/fixtures/pi_seq_dev/02_output)。
- 批次报告：[HTML](archive/testing/m2-20261009/validation/tests/fixtures/pi_seq_dev/00_report/TagForge_batch_report.html)。
- [回归日志](archive/testing/m2-20261009/tests.log)、[真实运行日志](archive/testing/m2-20261009/real-run.log)、[耗时和内存](archive/testing/m2-20261009/real-run.time.txt)。
- [源码摘要](archive/testing/m2-20261009/source-sha256.json)、[M1 对比](archive/testing/m2-20261009/comparison.json)、[重复运行](archive/testing/m2-20261009/rerun-check.json)。
# M3（已完成）

[M3 历史摘要](archive/validation/m3.json) 记录 81 项冻结快照完整回归、后续 Python 3.9 的 82 项回归及 linker 差异审查。
完整源码、测试日志、865 条差异 read 审阅、默认策略决策和顺序性能基准
保存于 `docs/archive/testing/m3-20261009/`。百万对 reads 的来源、前缀采样方式和
输入 SHA256 见该目录的 `million/fixture/fixture.json`。六次基准均已完成，
`performance/summary.json` 保存测量与限制，`performance/{synthetic,pi4,million}/count-comparison.json`
证明逐表计数一致。

## M4（本地准备已完成）

[M4 预检摘要](archive/validation/m4-preflight.json) 记录开发 wheel 在独立 Python 3.9 和 3.11 中的三类文库、
多样本、导出与 Slurm 检查。候选源码及构建保存在 `docs/archive/testing/m4-20261009/`。
报告打开检查的 PDF 预览使用用户指定 Arial 字体；历史候选记录保留对应批次事实。

rc2 历史验收：[m4-rc2.json](archive/validation/m4-rc2.json)。85 项回归在 Python 3.9/3.11 通过，三套 Python 的安装包各通过 31 条 CLI 检查；13 个 Excel、7 个 HTML 打开检查通过。完整候选包、SHA256、源码及 Arial PDF 预览保存于上述 M4 日期归档。最终状态见下文。


### 最终 2.0.0

最终状态以 [release.json](../tests/validation/release.json) 为准；原计划逐项证据见 [RELEASE_READINESS.md](RELEASE_READINESS.md)，既有修改审阅见 [source-review.json](../tests/validation/source-review.json)。本地准备完成，尚未提交、推送、打 tag 或上传。

- 最终安装包：[dist-release-r3](archive/testing/m4-20261009/dist-release-r3)，[SHA256SUMS](archive/testing/m4-20261009/dist-release-r3/SHA256SUMS)。
- 源码：[release-source-r3](archive/testing/m4-20261009/release-source-r3)，[逐文件 SHA256](archive/testing/m4-20261009/release-source-r3-sha256.json)。易变的验证索引 JSON/Markdown 不进入冻结快照。
- 91 项完整回归：[Python 3.9](archive/testing/m4-20261009/tests-release-39-r3.log)、[Python 3.11](archive/testing/m4-20261009/tests-release-311-r3.log)、[已安装 wheel/Python 3.13](archive/testing/m4-20261009/tests-release-installed-313-r3.log)，均无跳过。
- 三套安装流程、配置、导出与 Slurm 检查：[smoke39](archive/testing/m4-20261009/smoke39-release-r3)、[smoke311](archive/testing/m4-20261009/smoke311-release-r3)、[smoke313](archive/testing/m4-20261009/smoke313-release-r3)，各 31 条命令。
- 60 条独立阶段、迁移、指标、阈值、并发与过期导出检查：[release-audit-final-r3](archive/testing/m4-20261009/release-audit-final-r3)。
- 13 个 Excel/13 个 HTML 的打开结果及仅 Arial PDF 预览：[report-open-release-r3](archive/testing/m4-20261009/report-open-release-r3)。
- 配置模板/覆盖保护：[init-config-release-r3.json](archive/testing/m4-20261009/init-config-release-r3.json)。
- M1 552 个文件重新摘要核验：[historical-archive-verification-r3.json](archive/testing/m4-20261009/historical-archive-verification-r3.json)。

`dist-release/` 与 `dist-release-r2/` 是历史构建，尚未包含最后的过期矩阵导出修复，不能用于发布。实际复现（旧计数 6，重算 7）见 [stale-export-repro](archive/testing/m4-20261009/stale-export-repro/summary.json)；最终单元及已安装 CLI 验证已通过。
