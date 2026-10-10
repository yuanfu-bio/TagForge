# 文档与历史归档

整理日期：2026-10-09。根目录仅保留 README；使用说明、发布资料、教程和验证索引集中在 `docs/`，历史版本、设计与实施过程记录保留在 `docs/archive/`。

## 当前文档

- [README](../README.md)：安装、配置、命令与输出。
- [CHANGELOG](CHANGELOG.md)：当前 2.0.0 变更。
- [MIGRATION](MIGRATION.md)：升级与结果重算约定。
- [RELEASE_NOTES](RELEASE_NOTES.md)：2.0.0 发布说明。
- [RELEASE_READINESS](RELEASE_READINESS.md)：发布验收清单。
- [示例教程](tutorial.md)与[PI-seq fixture 说明](pi_seq_dev.md)。
- [最终验收记录](../tests/validation/release.json)与[完整结果索引](validation.md)。

## 历史版本与设计

- [历史 CHANGELOG](archive/CHANGELOG.md)：2.0.0rc1–rc3 与 0.1.1–0.1.14 记录；原日期及 Unreleased 状态保留。
- [旧设计说明](archive/TagForge_design.md)：历史架构与接口草案。
- [2.0.0 优化与发布计划](archive/2.0.0/TagForge_2.0.0_plan.md)。
- [2.0.0 分批实施记录](archive/2.0.0/TagForge_2.0.0_progress.md)。

## 阶段验收记录

| 阶段 | 记录 |
| --- | --- |
| M1 | [m1.json](archive/validation/m1.json) |
| M2 | [m2.json](archive/validation/m2.json) |
| M3 | [m3.json](archive/validation/m3.json) |
| M4 开发包预检 | [m4-preflight.json](archive/validation/m4-preflight.json) |
| 2.0.0rc1 | [m4-rc1.json](archive/validation/m4-rc1.json) |
| 2.0.0rc2 | [m4-rc2.json](archive/validation/m4-rc2.json) |

JSON 内的证据路径仍相对于项目根目录。完整实验输入、输出、冻结源码、报告和安装包已归档到 `docs/archive/testing/`；上述 JSON 是这些结果的索引。历史快照内部仍保留当时的执行路径，通过 [路径迁移表](archive/testing/relocation.json) 定位现存文件。三套历史 Python 环境和专用 Conda 缓存已删除。

## 本地旧笔记

[旧命令](archive/legacy/commands.md)、[Git 操作笔记](archive/legacy/git.md)、[工作日志](archive/legacy/journal.md)及 `archive/legacy/worklog/` 保留原内容，供本地查阅。这些笔记沿用原 Git 忽略策略，不进入发布包；当前操作以项目 README 为准。

[归档清单](archive/manifest.json)记录原位置、新位置以及迁移前后文件摘要。历史候选及冻结构建保留原样；当前工作区的文档路径以本索引为准。[测试归档清理记录](archive/validation/testing-cleanup-20261010.json) 记录本次环境删除和输出迁移。
