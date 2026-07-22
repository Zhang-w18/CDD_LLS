# 项目级 Agent 说明

## 1. 适用范围

本文件适用于整个仓库。任何 Agent 开始工作前先读本文件。

本项目研究频域相位预编码矩阵 `V`、CDD、DMRS 与接收机信道估计之间的关系。结论必须由可复现配置、原始数据和量化判据支持。

## 2. 文档职责

| 文件或目录 | 职责 | 读者 |
|---|---|---|
| `AGENTS.md` | 结构、读取顺序、工作流、输出规则 | 所有 Agent |
| `GOALS.md` | 目标、背景、进展、阶段验收 | 人与 Agent |
| `DESIGN.md` | 无注释理论框架；制定 plan 时必读 | 人与 Agent |
| `docs/design/DESIGN_ANNOTATED.md` | 推导、解释和历史修订 | 按需追溯 |
| `KNOWLEDGE.md` | 已验证结论、已排除方向、开放问题 | Agent |
| `research/plan-NNN.md` | 一轮实验的自包含执行规格 | 实现 Agent |
| `research/result-NNN.md` | 给人审阅的结果，包含图和关键表 | 人 |
| `research/result-NNN-text.md` | 与同编号结果等价的无图版 | Agent |
| `docs/reports/` | 阶段报告和专题报告 | 人 |
| `docs/figures/` | 被 Markdown 引用的图 | 人与 Agent |
| `docs/archive/` | 被当前文档替代的历史材料 | 按需追溯 |
| `outputs/` | 原始数据、日志、未引用图 | 实验执行者 |

`KNOWLEDGE.md` 只记录有结果支持的内容。未验证判断必须标记为假设或推断。

## 3. 项目结构与目录用法

```text
.
├── AGENTS.md / README.md
├── GOALS.md / DESIGN.md / KNOWLEDGE.md
├── cdd_lls/                 # 可复用仿真库
│   ├── core/                # 配置、MCS 等基础定义
│   ├── phy/                 # 信道、估计、LDPC、预编码、调制、资源网格
│   ├── sim/                 # 仿真编排与统计
│   └── utils/               # 通用工具
├── configs/                 # run.py 使用的 YAML
├── tools/                   # 专题实验和分析入口
├── tests/                   # 单元测试与回归测试
├── research/                # plan/result 与模板
├── docs/                    # design、reports、figures、archive
└── outputs/                 # 可再生运行产物，默认不提交
```

- 可复用算法放 `cdd_lls/`，不要复制到多个实验脚本。
- `run.py` 是通用 YAML 入口；专题扫描或分析放 `tools/`。
- 稳定配置放 `configs/`；实际运行必须保存展开配置。
- 新行为补充测试，或说明无法自动测试的原因。
- `research/` 不保存大规模原始数组。
- `docs/archive/` 不是当前规范来源。

## 4. 按任务读取文档

### 制定实验计划

按顺序读取 `AGENTS.md`、`GOALS.md`、无注释版 `DESIGN.md`、`KNOWLEDGE.md`、最新 `result-NNN-text.md`。只在这些文档给出指针时读取更早材料。

制定 plan 前，与研究者确认问题、基线、公平性、判据、预算和停止条件。未确认的 plan 标记为“草案”。

### 按 plan 实现实验

读取 `AGENTS.md`、指定 plan 及 plan 引用的代码和配置。plan 必须自包含；缺少关键参数、判据、路径或复现入口时先补全，不得猜测。

### 分析与更新

读取对应 plan、展开配置、原始输出和日志，逐条回答验收标准。result 经研究者确认后，再按要求更新报告、`KNOWLEDGE.md` 和 `GOALS.md`。

## 5. 标准工作流

1. 与研究者确认问题和判据。
2. 新建 `research/plan-NNN.md`。
3. 根据 plan 修改代码、配置和测试。
4. 运行测试或 smoke，验证单位、维度、种子和输出。
5. 正式产物写入 `outputs/<experiment_name>/<run_id>/`。
6. 成对生成 `result-NNN.md` 和 `result-NNN-text.md`。
7. 研究者确认后，更新报告或 `KNOWLEDGE.md`。
8. 更新 `GOALS.md` 和 `research/README.md`。

最新完成编号是 `024`；下一默认编号是 `025`。

## 6. plan 必需内容

1. 研究问题和可证伪假设；
2. 与 GOALS、DESIGN 的对应关系；
3. 系统模型、符号、全部参数；
4. 候选、基线和公平性约束；
5. 按执行顺序编写的步骤；
6. 种子、trial、SNR 网格和预算；
7. smoke 与正式运行条件；
8. 输出目录、文件、字段和单位；
9. 统计方法、置信区间和样本不足判定；
10. 量化验收、停止条件和后续动作；
11. 复现命令及预计修改文件；
12. 两个 result 版本的必需内容。

不得把关键定义只写成“同上一轮”或“见历史文档”。

## 7. result、图和数据规则

每轮正式实验必须成对输出：

- `result-NNN.md`：给人阅读，可以嵌入图；
- `result-NNN-text.md`：禁止嵌入图片，仅靠文字、表格、原始数据或采样数据核验结论。

两版必须有相同结论、配置、数字、异常和证据路径，并包含配置回执、逐项判定、统计不确定性、异常、数据路径、复现命令和适用范围。

- result 文档必须完整说明仿真条件；一般使用参数表列出全部仿真参数，并注明参数含义、取值和单位（如适用）。
- 需要向模型上传仿真结果图片时，优先上传该图片对应的文字或数据描述版本；仅当不存在文字或数据描述时，才上传原图。
- 被 Markdown 引用的图只放 `docs/figures/<topic-or-result>/`。
- 其他图留在 `outputs/`。
- 每张图必须有对应文字数据；禁止用“见图可知”替代分析。
- 采样数据说明规则并覆盖有利、持平和不利条件。
- 大规模数据留在 `outputs/`，result 给出精确路径。

## 8. Windows execution rules

- 所有源代码、Markdown、YAML、JSON、CSV 和文本文件统一使用 UTF-8。
- 仅使用 Python `pathlib` 修改文本文件，并显式指定 `encoding="utf-8"`。
- 不要通过 PowerShell 管道或 Base64 传递补丁。
- 不要通过 Windows PowerShell 5.1 管道传递非 ASCII 补丁；如确实无法避免，必须先将 `$OutputEncoding` 显式设置为 UTF-8。
- 不要使用 `cmd.exe echo`、shell 重定向或 ANSI 代码页写入中文文本。
- 不要通过 `powershell.exe -Command` 传递大型 Base64 补丁。
- 大型修改应先写入 UTF-8 patch 文件，再对该文件运行 `git apply`。
- 将修改拆分为小型、可独立验证的补丁。
- 每完成一组文件修改后立即运行 `git diff --check`。
- 单条简单文件系统命令运行超过 30 秒时，停止命令并诊断命令运行器；不要切换到 PowerShell、`cmd`、Node 或其他备用通道重复执行相同操作。
- 修改中文文件后，必须按 UTF-8 回读，并检查 ASCII 问号 `?` 和 Unicode replacement character `U+FFFD`。

## 9. 文档写作规则

- 禁止使用比喻、拟人和含义不确定的修辞。
- 禁止跳过推导、实验或判定步骤。
- 首次出现的缩写、符号和指标必须定义。
- 明确区分事实、假设、推断和待验证解释。
- 比较说明基线、方向、单位、样本数和不确定性。
- 记录负结果。
- 使用 UTF-8 和仓库相对路径。
- 数学符号以 `DESIGN.md` 为准；新增符号先加入符号表。

## 10. 完成检查

检查根目录入口、最新 result 无图版、GOALS/KNOWLEDGE 状态、引用图位置、路径和编号。实验任务还必须报告测试或 smoke 结果。
