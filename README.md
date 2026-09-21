# CDD LLS 仿真项目

本项目研究频域相位预编码矩阵 `V`、CDD、DMRS 和信道估计对链路级 BLER 的联合影响。

任何 Agent 先读 `AGENTS.md`。制定实验计划时，再按顺序读 `GOALS.md`、`DESIGN.md`、`KNOWLEDGE.md` 和最新 `research/result-NNN-text.md`。

当前最新完成实验是 `024`。result-024 已确认 Sidon delay 集在限定的 48 PRB 平坦信道、V-aware matched LMMSE 条件下优于等差 QC CDD；下一阶段验证物理 PDP、失配、定时误差和带宽变化。

## 代码入口

```powershell
python run.py --config configs/smoke.yaml
python -m unittest discover -s tests
```

专题实验入口位于 `tools/`，命令以对应 plan 为准。

## 固定 BLER 曲线入口

static TDL-A、8 Tx、16QAM MCS 8（码率 553/1024，取自 NR 256QAM MCS table）、V-aware matched LMMSE 或真实 data-RE CSI 接收机的多 candidate BLER 曲线，统一调用：

```powershell
python tools/run_bler_curves.py --config <bler-config.yaml> --stage validate
python tools/run_bler_curves.py --config <bler-config.yaml> --stage run
python tools/run_bler_curves.py --config <bler-config.yaml> --stage merge
```

`run` 会在运行结束后自动执行 merge 和出图；单独的 `merge` 用于只重建合并 CSV 与图片。可直接复制 `configs/bler_curves_result028_augmentation.yaml` 作为配置起点，无需修改入口代码。配置中指定：

- manifest、SHA-256 与审批回执路径；
- `scenario_id` 和需要运行的 `candidate_ids`，`all` 表示 manifest 中全部 candidate；
- `estimated`、`ideal` 或两种 receiver；
- 新曲线使用 `snr_db + fixed_total_trials`；已有曲线追加使用 `base_csv + target_errors + max_total_trials`；
- `stop_below_bler`：每条曲线保留首个低于该门限的点，之后不再追加或绘制，适合避免 ideal-CSI 在远低于 1% 的高 SNR 浪费计算量。
- `candidate_overrides`：按 candidate ID 覆盖上述预算或截断字段；适合只提高个别抖动尾点的 `target_errors`/`max_total_trials`，而不扩大整组曲线的高 SNR 计算量。
- `mode: transparent_prg_dft`：plan-028 第 10 节的 A100 comb-6 透明基线模式；由 `transparent_prg_baselines` 声明 4-RB/6-RB 的 8×8 空域 DFT cycling，4-RB 尾部映射按 absolute trial 无放回重抽，接收端在 PRG 内使用底层 PDP 形状的 matched LMMSE。
- `plan_section: 15`：仅用于 plan-028 第 15 节的 1Rx transparent PRG 6-RB estimated-CSI SNR 延伸；候选、receiver 和 SNR 网格均由入口严格校验，不生成 ideal-CSI 数据。
- `paired_policy`：透明 PRG 模式中两条基线及 ideal/estimated 接收机共同使用的 `target_errors`、`min_total_trials` 和 `max_total_trials`；`max_interval_trials` 只限制单次可恢复落盘区间，不改变总预算。`max_snr_tasks_per_run` 可限制一次命令处理的 SNR 区间数，便于外层编排安全续跑。
- `mode: transparent_cdd_physical_covariance`：plan-028 第 11 节的 A100 comb-6 透明 CDD 接收机模式；由 `transparent_cdd_baselines` 把新曲线 ID 映射到既有 CDD 发射 candidate，真实等效信道仍使用 source delay set，但五条 estimated 接收机共用 $8R_{phy}$，不生成 ideal-CSI 数据。

追加运行按绝对 trial 编号派生随机数。若 `base_csv` 某点已有 3000 trials，新数据从 trial 3001 开始；错误块数、Wilson 95% 区间和 estimated-CSI 的线性域 CE NMSE 会与原数据合并。每个完成的 SNR/trial 区间立即保存，重复调用同一命令会校验区间连续性并从未完成部分继续。不得删除或改写 `supplemental_points.csv` 后继续沿用同一输出目录。

透明 PRG 模式还保存逐 trial 4-RB 尾部映射及 SHA-256、CE 误差能量、真实信道能量和 trial 级线性 NMSE；ideal/estimated 与两条基线必须追加相同 absolute trial 区间。可重复调用同一配置，或使用 `tools/run_plan028_transparent_prg_until_complete.ps1` 逐区间运行至预算条件满足。

透明 CDD 模式保存五条曲线共同的 absolute-trial 区间、estimated error flags、CE 三类逐 trial 数组，以及按真实 CDD 协方差计算的解析失配 trace-NMSE。正式配置为 `configs/bler_curves_result028_a100_transparent_cdd.yaml`，可用 `tools/run_plan028_transparent_cdd_until_complete.ps1` 安全续跑。

入口输出包括 `supplemental_points.csv`、逐 trial error flags、`final/<receiver>_csi_bler_points.csv`、共享 `curve_styles.json` 及 candidate/缩写两套图。当前入口的固定物理范围写在本节首句；超出该范围时先扩展配置 schema 和测试，不得只凭字段名称假定已经支持。

## plan-031 严格 Sidon 候选预扫

运行 2026-09-08 新增规划中的 6 个 `(Tx, AL)` 场景和 48 个冻结候选时，在仓库根目录执行一次：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py
```

入口先输出当前完成量和总预算，再校验配置与 `sidon_shortlist.json` 完全一致、模 $K$ 严格 Sidon、fold residue、pilot rank/condition、零噪声 CE floor 和无噪声 PDCCH 译码。校验通过后最多并行运行三个 candidate，每 100 trials 输出一次进度，每完成一个 candidate/SNR 点立即保存 BLER CSV、逐 trial error flags 和 CE NMSE 数组。重复执行相同命令会跳过已经完成的点并从剩余 candidate 继续。

所有结果、校验回执、环境信息、编排状态和完整日志保存在 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/`。只检查进度而不运行时增加 `--summary-only`。预扫完成后先按 plan 冻结各 candidate 的正式 SNR 网格，再生成正式配置；该入口不会以后验规则自动选择正式网格。

若原始预扫完成后，冻结清单中的 5 个 C300 候选尚未对 1% BLER（其中一条还包括 10% BLER）形成双侧 bracket，执行一次固定的高 SNR 补扫：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py --phase extension
```

该阶段只增加 54 个 candidate/SNR 点（16,200 candidate-trials），继续写入原 candidate 目录，并将独立状态、验证回执和日志写入 `prescan_extension_orchestration/`。重复执行同一命令可续跑。

排除有明确 CE/BLER floor 的 5 个候选后，对其余 43 个候选执行 0.5 dB、3,000 trials/点的粗确认：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py --phase coarse-confirmation
```

该阶段共 254 个 candidate/SNR 点（762,000 candidate-trials），每个候选分别覆盖原预扫 10% 和 1% BLER 首个下降 bracket 的两端及中点。状态和日志写入 `coarse_confirmation_orchestration/`，BLER、error flags、CE NMSE 和展开配置写入 `coarse_confirmation/`；同一命令可安全续跑。

粗确认筛选后，对 18 个仍可能最优的候选执行最终 0.25 dB 细扫描：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py --phase fine
```

该阶段共 149 个 candidate/SNR 点。每点至少运行 10,000 trials，达到 200 errors 后停止，最多运行 50,000 trials；最低总预算 1,490,000 candidate-trials，最大总预算 7,450,000。入口最多并行三个候选，自动输出逐点进度并保存 BLER、error flags、CE NMSE、展开配置、验证回执、日志和环境信息。重复执行同一命令会跳过已完成点并安全续跑。

若正式统计审计提示目标交点位于网格边界，执行冻结的 23 点边界补扫：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py --phase fine-supplement
```

该命令只运行 10 个受影响候选的新 SNR 点，输出独立保存到 `fine_supplement/`，不会重跑或覆盖原 149 点；同一命令可安全续跑。

## result-032 移动性续跑

将 result-032 的 14–20 dB 曲线续跑至 result-028 的样本量级时，只需在仓库根目录执行一次：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_result032_result028_scale.py
```

入口首先输出已完成/目标/剩余的共同 trial 摘要，随后最多同时运行三个互不重叠的 SNR 分片，并将每个分片的实时进度输出到同一窗口。每完成一个 SNR 点，内层入口立即保存该点的 `intervals.csv`、逐 trial error flags、CE NMSE 数组和分片合并 CSV；中断时只丢失尚未完成的当前 SNR 计算。再次执行完全相同的命令会核对 absolute-trial 连续性、跳过已经完成的 SNR，并从第一个未完成点继续。

编排状态保存在 `outputs/experiment032_tdl_mobility/20260906_main/result028_scale_orchestration/status.json`，每次运行的完整控制台日志保存在同目录的 `logs/`。若只想检查当前进度而不启动仿真，可增加 `--summary-only`。不要同时在两个窗口运行该命令；入口具有单实例锁，会拒绝第二个进程以避免重复写入。

## 项目结构

```text
.
├── AGENTS.md               # Agent 工作流、文档规则和执行约束
├── README.md               # 项目入口和结构说明
├── GOALS.md                # 研究目标、关键进展和跨实验全局变更
├── DESIGN.md               # 当前理论框架和符号规范
├── KNOWLEDGE.md            # 已验证结论、已排除方向和开放问题
├── run.py                  # 通用 YAML 仿真入口
├── cdd_lls/                # 可复用仿真库
│   ├── core/               # 配置、MCS 等基础定义
│   ├── phy/                # 信道、估计、LDPC、预编码、调制和资源网格
│   ├── sim/                # 仿真编排与统计
│   └── utils/              # 通用工具
├── configs/                # run.py 使用的稳定 YAML 配置
├── tools/                  # 专题实验、参数扫描和分析入口
├── tests/                  # 单元测试与回归测试
├── research/               # plan、result、无图 result 和研究索引
├── docs/
│   ├── design/             # 带解释的设计材料
│   ├── reports/            # 阶段报告和专题报告
│   ├── figures/            # 被 Markdown 文档引用的图
│   └── archive/            # 已被当前文档替代的历史材料
└── outputs/                # 可再生的原始数据、日志和未引用图，默认不提交
```

各轮实验的执行规格位于 `research/plan-NNN.md`，结果成对记录在 `research/result-NNN.md` 和 `research/result-NNN-text.md`。正式运行的完整产物保存在 `outputs/<experiment_name>/<run_id>/`。
