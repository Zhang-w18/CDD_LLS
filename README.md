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

目录：`cdd_lls/` 为仿真库，`configs/` 为配置，`tools/` 为专题入口，`tests/` 为测试，`research/` 为研究记录，`docs/` 为文档和引用图，`outputs/` 为可再生输出。
