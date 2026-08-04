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

static TDL-A、8 Tx / 1 Rx、256QAM MCS 8、V-aware matched LMMSE 或真实 data-RE CSI 接收机的多 candidate BLER 曲线，统一调用：

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

追加运行按绝对 trial 编号派生随机数。若 `base_csv` 某点已有 3000 trials，新数据从 trial 3001 开始；错误块数、Wilson 95% 区间和 estimated-CSI 的线性域 CE NMSE 会与原数据合并。每个完成的 SNR/trial 区间立即保存，重复调用同一命令会校验区间连续性并从未完成部分继续。不得删除或改写 `supplemental_points.csv` 后继续沿用同一输出目录。

入口输出包括 `supplemental_points.csv`、逐 trial error flags、`final/<receiver>_csi_bler_points.csv`、共享 `curve_styles.json` 及 candidate/缩写两套图。当前入口的固定物理范围写在本节首句；超出该范围时先扩展配置 schema 和测试，不得只凭字段名称假定已经支持。

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
