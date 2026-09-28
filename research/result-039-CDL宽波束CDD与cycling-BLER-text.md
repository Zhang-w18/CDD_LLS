# result-039：CDL 宽波束 CDD 与 cycling BLER（当前已完成批次）

> 状态：阶段性记录，待研究者确认。阶段 1B 的 `cycling_initial` 已完成；CDD 已完成全网格
> 初始 1,000 trials，并只在 4 个 SNR 点完成第一批追加。补充的 60 km/h、40 ms outdated-CSI
> MRT 闭环曲线已完成冻结的 26,000 trials。原三曲线的自适应追加、paired bootstrap 和最终
> 验收尚未完成，因此本文件不把原三曲线的当前 crossing 当作最终结论。

## 1. 仿真条件与复现

- 场景：`C_UE2R_DUAL_ASD25`，Sionna CDL-C，RMS delay spread 100 ns，4 GHz，3 km/h，
  双极化 32 TXRU 发射、2 Rx、rank 1。
- 资源：48 PRB、576 个 active subcarriers、30 kHz SCS、10 个 PDSCH symbols；DMRS symbols
  `[2,7]`、comb-6；16QAM，NR 256QAM table MCS 8，目标码率 553/1024。
- 接收机：`B0_QC` 与 `SIDON_SELECTED` 使用 `PLAN039_COMMON_REFERENCE_PDP`；
  `PRG6_CYCLING` 使用每个 6-RB PRG 内独立的 `PLAN039_PRG_COMMON_REFERENCE_PDP`。
- 基础种子 20260939；evaluation channel seed 2152137127。候选与 seed 冻结见
  `outputs/experiment039_cdl_beam_bler/stage0_20260925/estimated_confirm/C_UE2R_DUAL_ASD25/stage1b_method1_only_freeze.json`。
- 本次记录对应 Git `3d98527f46cbd4b305e6d44cbd6e1d15e40f56e0` 加未提交的 plan-039 工作区实现；
  工作区不是干净 checkpoint。
- 已执行命令：

  ```powershell
  python tools/run_plan039_cdl_beam_bler.py --config configs/plan039_cdl_wide_beam_bler.yaml --stage cycling_initial
  python tools/analyze_plan039_stage1b_initial.py --input outputs/experiment039_cdl_beam_bler/stage0_20260925/estimated_confirm --output outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis --figure-dir docs/figures/result-039-CDL宽波束CDD与cycling-BLER
  python tools/run_plan039_cdl_beam_bler.py --config configs/plan039_cdl_wide_beam_bler.yaml --stage aged_mrt_formal
  python tools/analyze_plan039_aged_mrt.py --input outputs/experiment039_cdl_beam_bler/stage0_20260925/aged_mrt_60kmh/C_UE2R_DUAL_ASD25 --output outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis/aged_mrt --figure-dir docs/figures/result-039-CDL宽波束CDD与cycling-BLER --stage1b-points outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis/stage1b_current_points.csv
  ```

分析脚本只纳入存在 `batch_receipt.json` 和 `summary.csv` 的完成批次；无 receipt 的残留目录不计入。

## 2. 当前完成范围

| 曲线 | SNR 点数 | SNR 范围 | 1,000 trials 点数 | 2,000 trials 点数 | 累计 trials |
|---|---:|---:|---:|---:|---:|
| `B0_QC` | 36 | 7.5--20.0 dB | 32 | 4 | 40,000 |
| `SIDON_SELECTED` | 36 | 7.5--20.0 dB | 32 | 4 | 40,000 |
| `PRG6_CYCLING` | 36 | 7.5--20.0 dB | 36 | 0 | 36,000 |
| **合计** | **108 个曲线-SNR 点** | — | **100** | **8** | **116,000** |

CDD 的 2,000-trial 点是 13.75、14.0、16.0 和 16.5 dB，B0 与 Sidon 在同一绝对 trial
区间 `[0,1000)`、`[1000,2000)` 配对。其余 CDD 点和全部 cycling 点当前均为 `[0,1000)`。
逐点 errors、trials、BLER、Wilson 95% 区间和线性域聚合后的 CE NMSE 位于
`outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis/stage1b_current_points.csv`；范围汇总位于
`outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis/stage1b_trial_coverage.csv`。

## 3. 当前数值

以下 crossing 仅按现有原始 BLER 点，在相邻真实双侧 bracket 内做 log-BLER 线性插值；尚未做
预声明的 1,000 次 paired bootstrap。

| 曲线 | 10% crossing / bracket | 1% crossing / bracket | 当前判定 |
|---|---|---|---|
| `B0_QC` | 14.051 dB / [14.0, 14.25] | 16.189 dB / [16.0, 16.5] | 两目标均有当前 bracket |
| `SIDON_SELECTED` | 14.039 dB / [14.0, 14.25] | 16.357 dB / [16.0, 16.5] | 两目标均有当前 bracket |
| `PRG6_CYCLING` | 18.168 dB / [18.0, 18.25] | 无；20 dB 为 22/1000 = 2.2% | 1% 未达到 |

当前 B0 与 Sidon 的 10% crossing 差只有约 0.012 dB；1% crossing 上 B0 暂时领先约
0.167 dB。相对于这一级别的差异，当前样本尚不足以支持“二者有确定 BLER 优劣”的结论。
cycling 的当前 10% crossing 比两条 CDD 曲线高约 4.1 dB，但最终 gain 仍须等自适应预算和
paired bootstrap 完成后报告。

CE NMSE 与 BLER 呈随 SNR 改善的趋势。一个需要明确记录的现象是 cycling 虽然 CE NMSE 更低，
BLER 却显著更差；例如 18.0 dB 时 cycling 为 111/1000、CE NMSE -23.37 dB，而 B0 为
1/1000、-20.27 dB。这是数据事实，说明单独的等效信道 NMSE 不能替代不同预编码方案的 BLER
比较；其物理分解尚未在本轮完成，不能把该现象直接解释为某一单一机制。

## 4. 前置 gate 与候选冻结

- validate 状态为 PASS。ASD25 的功率加权圆均值约 $0^\circ$，实测 ASD 24.936°；AoD
  `[-60°,60°]` 覆盖率 0.9484，AoD--ZoD 目标矩形覆盖率 0.9351，均通过预声明门限。
- ASD25 的 8 波束域非对角相关均值为 0.4530，低于 ASD10 的 0.6656；有效秩从 2.306
  提高到 3.428，相关性 gate 通过。
- stage-0 smoke 的 8 个 variant、两个 SNR、每点 20 trials 均通过数值有限性、配对重放和
  单位预编码功率检查；smoke 不作为性能证据。
- 阶段 1A 每条已执行曲线在 11 个 SNR 点各运行 400 trials。唯一实际执行的 Sidon 候选
  `[0,11,28,148,170,233,277,351]` 被冻结为 `SIDON_SELECTED`；其余 top-8 候选均未运行并标记
  `DEFERRED_NOT_RUN`。因此当前结果不能支持“在全部 top-8 中择优”的解释。

## 5. 结论边界与未完成项

事实：cycling 全部 36 个初始点已完成；两条 CDD 曲线各完成 36 个初始点及 4 个追加点。
当前三条曲线均随 SNR 总体下降，未发现足以否定本批数据的明显趋势错误。

推断：当前数据强烈提示 cycling 在该公共参考 PDP、PRG-local LMMSE 和冻结波束映射下的 BLER
明显劣于两条 CDD 曲线，而 B0 与 Sidon 目前接近。

尚不能支持：最终 10%/1% crossing、方案 gain、置信区间显著性、最终 capped 状态，以及
“Sidon 优于 B0”的结论。正式完成前仍需由研究者决定是否恢复自适应追加；若恢复，必须从现有
`absolute_trial_stop` 继续，随后完成 paired bootstrap 和最终验收。未经研究者确认，不更新
`KNOWLEDGE.md` 或 `GOALS.md`。

## 6. 补充闭环曲线：60 km/h、40 ms outdated-CSI MRT

### 6.1 冻结条件与完成性

该补充曲线的 variant 为 `aged_mrt_prg6__60kmh__csi_age40ms`。UE 速度为 60 km/h；发射端在
每个 6-RB PRG 内，使用同一 realization 在当前 slot 起点前 40 ms 的物理 CSI，由两根 Rx 和
该 PRG 全部子载波的 Gram 矩阵主特征向量形成未量化 rank-1 MRT。权值在一个 PRG 和整个
10-symbol slot 内固定。接收端仍为 `PLAN039_PRG_COMMON_REFERENCE_PDP` estimated-CSI LMMSE
与 2Rx MRC。SNR 网格为 7.5--20.0 dB、间隔 0.5 dB，共 26 点；每点固定 1,000 trials，未做
prescan 或自适应追加。

完成回执状态为 `COMPLETE`。汇总表包含 26 点，trial 表包含 26,000 行；每点绝对 trial 区间均为
`[0,1000)`。汇总表和 trial 表 SHA-256 分别为
`1fecac46775ad081f2b695521e3c3bb4404f77b64597600e7f9a94bf147ff7b7` 和
`62cbb7dbc20bdd5b06390b2cc5e9019f270ed942c6a68d5a9ed07d85f8c9e4ac`，与完成回执一致。
全部汇总行的 CSI age 均为 40 ms；同 realization aged/current replay 最大逐元素误差为 0；
MRT 原始权值功率范围为 `[0.9999999999999971, 1.0000000000000027]`，通过单位功率审计。

### 6.2 原始数值、crossing 与统计不确定性

下表抽样规则为：覆盖最低 SNR、10% bracket 两端、1% bracket 两端、首个零错误点和最高 SNR。
CE NMSE 先在线性域跨每点 1,000 trials 聚合，再转为 dB。

| SNR (dB) | TB errors / trials | BLER | Wilson 95% 区间 | CE NMSE (dB) |
|---:|---:|---:|---:|---:|
| 7.5 | 239 / 1,000 | 0.239 | [0.2136, 0.2664] | -17.25 |
| 10.0 | 105 / 1,000 | 0.105 | [0.0875, 0.1255] | -19.11 |
| 10.5 | 95 / 1,000 | 0.095 | [0.0783, 0.1148] | -19.51 |
| 14.5 | 14 / 1,000 | 0.014 | [0.00836, 0.02336] | -22.65 |
| 15.0 | 9 / 1,000 | 0.009 | [0.00474, 0.01702] | -23.10 |
| 18.0 | 0 / 1,000 | 0 | [0, 0.00383] | -25.53 |
| 20.0 | 0 / 1,000 | 0 | [0, 0.00383] | -27.21 |

基于相邻真实采样点的 log-BLER 线性插值，10% crossing 为 10.244 dB，真实 bracket 为
[10.0, 10.5] dB；1% crossing 为 14.881 dB，真实 bracket 为 [14.5, 15.0] dB。二者均为
固定预算下的描述性点估计，没有 bootstrap crossing 区间。18.0--20.0 dB 共 5 点为零错误；
原始 CSV 保留 BLER 0，绘图时按预声明规则置于 `0.5/1000` 并用空心圆标记。全部 26 点累计
1,437 个 TB errors。

事实：BLER 随 SNR 总体下降，CE NMSE 从 7.5 dB 的 -17.25 dB 改善到 20 dB 的 -27.21 dB；
10% 和 1% 均存在唯一真实下降 bracket，未发现多 crossing 或反常的大幅回升。

另生成四条曲线的 BLER/CE NMSE 合图，不绘制 Wilson 区间。该合图把 aged-MRT 与原三条曲线
并列显示，但不改变上述比较边界，不据此计算未配对 gain。

解释边界：该曲线与原三条曲线在速度和发射端知识上均不同，而且没有共享 realization，不能把
crossing 差值归因为 MRT 本身，也不能进行 paired gain 或显著性判定。它只支持“在冻结的
60 km/h、40 ms outdated-CSI、PRG6 MRT 和指定接收机配置下取得上述单曲线性能”的结论。

## 7. 证据路径

- 完成批次：`outputs/experiment039_cdl_beam_bler/stage0_20260925/estimated_confirm/batches/`
- cycling 合并原始表：`outputs/experiment039_cdl_beam_bler/stage0_20260925/estimated_confirm/C_UE2R_DUAL_ASD25/cycling_initial_summary.csv`
- 当前三曲线统一表：`outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis/stage1b_current_points.csv`
- 当前 crossing 表：`outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis/stage1b_current_crossings.csv`
- 绘图样式：`outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis/plot_style.json`
- validate：`outputs/experiment039_cdl_beam_bler/stage0_20260925/validate/C_UE2R_DUAL_ASD25/validate_report.json`
- selection freeze：`outputs/experiment039_cdl_beam_bler/stage0_20260925/selection/C_UE2R_DUAL_ASD25/selection_freeze.json`
- aged-MRT 冻结与完成回执：`outputs/experiment039_cdl_beam_bler/stage0_20260925/aged_mrt_60kmh/C_UE2R_DUAL_ASD25/aged_mrt_formal_freeze.json`、`aged_mrt_formal_report.json`
- aged-MRT 原始汇总与 trial 数据：同目录下 `aged_mrt_formal_summary.csv`、`aged_mrt_formal_trial_metrics.csv`
- aged-MRT 分析表与审计：`outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis/aged_mrt/aged_mrt_points.csv`、`aged_mrt_crossings.csv`、`aged_mrt_analysis_audit.json`
- aged-MRT 绘图脚本与样式：`tools/analyze_plan039_aged_mrt.py`、`outputs/experiment039_cdl_beam_bler/stage0_20260925/analysis/aged_mrt/aged_mrt_plot_style.json`
