# 项目级 Agent 说明

## 1. 范围与基本原则

本文件适用于整个仓库，任何 Agent 开始工作前先读本文件。

本项目研究频域相位预编码矩阵 `V`、CDD、DMRS 与接收机信道估计之间的关系。结论必须由可复现配置、原始数据和量化判据支持；事实、推断、假设和待验证解释必须明确区分。数学符号和物理量定义以 `DESIGN.md` 为准。

只读取完成当前任务所需的最小文档集合。除本文件外，不得因为“可能有用”而默认加载全部全局文档、历史 result 或原始输出。

## 2. 按任务读取

| 任务 | 必读内容 |
|---|---|
| 理论问答、设计讨论 | `DESIGN.md` 的相关章节；仅在需要实验依据时读取 `KNOWLEDGE.md` 或相关 `result-NNN-text.md` |
| 小型文档修改 | 目标文件及其直接引用的规范；不自动读取 `GOALS.md`、全部 `DESIGN.md` 或历史 result |
| 代码诊断或局部实现 | 目标代码、相关配置、测试和日志；涉及物理定义时再读 `DESIGN.md` 相关章节 |
| 制定正式实验计划 | `GOALS.md` 相关问题、`DESIGN.md` 相关章节、`KNOWLEDGE.md` 相关结论、最新相关 `result-NNN-text.md` 和 `docs/agent/EXPERIMENT_PLAN.md` |
| 按 plan 执行 | 指定 plan、plan 引用的代码与配置，以及 `docs/agent/EXPERIMENT_PLAN.md` |
| 分析正式实验 | 对应 plan、展开配置、精确输出与日志，以及 `docs/agent/RESULT_SPEC.md` |
| 更新全局结论 | 已确认的 result、`GOALS.md`/`KNOWLEDGE.md` 的相关章节，以及 `docs/agent/RESULT_SPEC.md` |

先用目录、标题、搜索或窄范围片段定位内容；只有跨章节关系确有必要时才阅读全文。`docs/design/DESIGN_ANNOTATED.md` 和 `docs/archive/` 只用于按明确指针追溯。

## 3. 文档与目录职责

- `GOALS.md`：研究问题、关键进展、阶段验收和跨实验全局变更。
- `DESIGN.md`：当前理论框架、公式、物理量定义和符号规范。
- `KNOWLEDGE.md`：有结果支持的结论、已排除方向和开放问题。
- `research/plan-NNN.md`：一轮代码修改或正式实验的自包含执行规格。
- `research/result-NNN.md` 与 `research/result-NNN-text.md`：同一轮正式实验的图文版与可独立核验的无图版。
- `cdd_lls/`：可复用算法；`run.py`：通用 YAML 入口；`tools/`：专题扫描或分析。
- `tools/run_bler_curves.py`：static TDL-A 多 candidate BLER 曲线的固定配置入口；用法和配置字段见 `README.md` 的“固定 BLER 曲线入口”。
- `configs/`：稳定配置；`outputs/`：展开配置、原始数据、日志和未引用图。
- `docs/reports/`：面向研究者的报告；`docs/figures/`：被 Markdown 引用的图。

完整结构见 `README.md`。不要在本文件硬编码“最新/下一”实验编号；从 `research/README.md` 和实际文件名核对，以未占用编号为准。

## 4. 工作流边界

新的正式仿真实验必须先有 plan，并执行“确认问题与判据 → 实现与测试 → smoke → 正式运行 → 成对 result → 研究者确认 → 更新全局文档”的流程。详细要求见 `docs/agent/EXPERIMENT_PLAN.md` 和 `docs/agent/RESULT_SPEC.md`。

理论问答、小型文档修改、只读诊断和范围明确的局部修复不要求机械创建 plan；按任务风险进行必要验证。result 未经研究者确认，不得写入 `KNOWLEDGE.md` 作为已验证结论，也不得据此改变 `GOALS.md` 的阶段验收状态。

Git checkpoint 以完成的任务为单位，不以单次文件编辑为单位。任务完成并验证通过后，先向研究者汇报；只有研究者明确确认本轮交付并允许创建 checkpoint，才创建仅包含本任务相关修改的语义完整 Git commit。在此之前以当前工作区文件为内容事实源，不自动 stash、暂存或提交既有用户改动。

## 5. 实现与证据原则

- 可复用算法放 `cdd_lls/`，不要复制到多个实验脚本。
- 每次增加或追加符合固定范围的 BLER 曲线时，优先为 `tools/run_bler_curves.py` 新建 YAML 配置并依次执行 `--stage validate`、`--stage run`；不得为只增加场景、candidate、SNR 或 trial 预算而复制链路主循环。已有点追加必须通过 `base_csv` 和绝对 trial 区间继续，禁止从 trial 1 重跑后与旧数据相加。超出 README 声明的固定物理范围时，才扩展入口 schema、实现和测试。
- 新行为补充测试；无法自动测试时说明原因并执行可复现的 smoke。
- 实际运行保存展开配置，并核对单位、维度、索引约定、归一化、随机种子和输出字段。
- 仿真图必须在本地工作区内由可复现脚本读取已保存的本地原始数据生成；不得把数据交给云端绘图服务或模型侧图像生成能力代画。绘图脚本、配置或完整复现命令必须随结果保留。
- `research/` 不保存大规模原始数组；正式产物写入 `outputs/<experiment_name>/<run_id>/`。
- 物理量定义、单位、相位分母、循环移位/真实时延等约定必须显式记录，不得凭名称推断。
- 保留工作区中与当前任务无关的用户改动；不得用破坏性 Git 命令清理工作区。

## 6. 上下文与工具开销

- 当前工作区文件是内容事实源；先读取目标文件的相关部分。用 `git status --short -- <path>` 判断目标是否已有未提交改动，只有需要区分既有改动边界时才查看相关 diff。
- 同一内容不重复读取；已有可靠摘要时，只补读完成任务所需的变化或原文片段。
- 目标文件干净时，不为回溯历史而读取 diff；目标文件非干净时，只查看会与当前任务重叠的 diff 片段。避免同时输出完整文件与完整 diff。
- 工具输出默认限制在可审阅范围内；结果过长时改用筛选、分页或统计摘要。
- 分析仿真结果时只向模型提供 `result-NNN-text.md`、CSV/TSV/JSON、日志、绘图输入数据或由本地脚本生成的数值摘要，不得向模型传输或加载结果图片。即使需要新增或修改图，也应由本地脚本完成，并通过数据、脚本日志、文件存在性和尺寸等本地可核验信息检查；图片文件仅保存到仓库或输出目录供研究者查看。
- 非仿真任务需要处理图片时仍遵循最小传输原则：先使用对应数据或文字版；只有任务本身明确要求分析图片内容且不存在可替代的结构化信息时，才传输图片。
- 单个范围明确的文件修改通常不启用多个 Agent。并行任务必须独立、边界不重叠，并要求返回简短结论和精确路径。

## 7. 编辑与写作

编码、补丁、Windows 执行、Markdown 与数学公式规则见 `docs/agent/EDITING_RULES.md`。修改文本后必须以 UTF-8 回读并执行定向 `git diff --check`。

## 8. 完成条件

交付前检查任务范围内的 diff、编码、路径、引用和编号；运行与风险相称的测试或 smoke。最终回复说明改了什么、验证结果、未执行项及其原因。正式实验还必须按 `docs/agent/RESULT_SPEC.md` 检查证据和文档配对。
