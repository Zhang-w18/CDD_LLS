# result-026-text：CDD 协方差适配、厚 Sidon 潜力与近似平坦 TDL 复核（无图版）

> 对应 `research/plan-026.md` 和 `research/result-026.md`。状态：正式运行已完成并经研究者确认。本文不嵌入图片；所有图均有对应 CSV/JSON 文字数据。

## 1. 结论与验收判定

本轮按 plan-026 完成 E1、E2、E3。统一以“基线目标 SNR 减去候选目标 SNR”为增益，正值表示候选达到相同 outage/BLER 所需的 SNR 更低。

| 假设 | 判定 | 关键证据 |
|---|---|---|
| H1：完整协方差下存在同时通过 outage、16 dB CE、80 dB CE 地板和 rank 的非等差设计 | **不成立** | 10 个冻结候选中无一通过全部门槛。S0 的 10% ideal-CSI outage 增益为 `+0.246 dB`，但 16/80 dB CE 分别劣化 `3.030/40.571 dB`；新联合候选 E1_0028 的增益为 `+0.161 dB`，16 dB CE 劣化 `0.688 dB`，但 80 dB 地板劣化 `3.608 dB`。 |
| H2：至少两个预定场景存在满足厚 Sidon 硬约束、outage 增益不少于 0.10 dB、CE 劣化不超过 1 dB 的候选 | **成立** | TDL-A 10 ns、TDL-C 5 ns、TDL-C 10 ns 三个场景通过；增益为 `+0.160/+0.202/+0.174 dB`。 |
| H3a：TDL-A 0.1 ns 下原 Sidon 的 10% estimated-CSI BLER 增益不少于 0.15 dB，且区间不支持反向 | **成立** | 增益 `+0.304 dB`，保守 95% 区间 `[+0.169,+0.439] dB`。 |
| H3b：TDL-D/E LoS 诊断 | **D0.1、E0.1、E1 区间支持 Sidon；D1 未闭合** | 三个可拟合场景的 10% 增益为 `+6.576/+6.661/+6.704 dB`，区间全正；D1 的 QC 在 19.5 dB 仍为 41/400 错误，未在预定扫描上限内形成 10% 跨越。 |

事实层面的总判断是：

1. 在 TDL-A 5 ns 完整协方差下，放松 residue 并联合搜索 lift 确实能找到不同的 `J_CDD`–CE 折中，但没有候选通过 H1 的全部工程门槛；
2. 只使用 `T_epsilon` 的厚 Sidon 生成规则在三个预定 NLoS 场景中找到了值得进入后续 estimated-CSI 链路的候选，因此该方向不应被一般性否定；
3. 原严格 Sidon 在 TDL-A 近似平坦极限恢复了 10% BLER 增益；TDL-D/E 的大增益还同时包含 LoS 统计改变，不能单独归因于时延扩展趋零。

## 2. 仿真条件与复现标识

| 类别 | 实际取值 |
|---|---|
| 有效带宽 | 48 PRB，576 active subcarriers，30 kHz SCS |
| FFT / CP / symbols | 4096 / 288 samples / 10 |
| Tx / Rx / layer | 8 / 1 / 1 |
| DMRS | symbols `[2,7]`，comb 24，offset 0；24 个唯一频率导频，5712 data RE |
| CDD | `V[k,n]=exp(-j2π k j_n/576)`，恒模 1；一个 delay index 为 `57.870370 ns` |
| 接收机 | known-V、known-PDP、delay-matched 全带频率 LMMSE；两个静态 DMRS 等效平均 |
| 链路 | 16QAM，MCS 8，码率 553/1024，LDPC 最多 8 次迭代 |
| 信道 | Sionna 1.0.2 TDL，0 km/h，3.5 GHz；不做 per-realization normalization |
| 共同随机数 | 同场景、SNR、trial 内候选共享 TDL、payload、平均 LS noise 和 data noise |
| E1/E2 outage | 16QAM BICM ideal-CSI outage；每个候选/场景 200,000 个共同样本 |
| E3 BLER | 粗扫 400 paired trials/点；自动精扫 3000 paired trials/点 |
| seed | `20260726` |
| run id | `20260724_main` |
| 执行代码版本 | Git `ecfb6a95ef3582499f03571b41bd99a8654f1e63` 加本轮未提交实现 |
| Python 环境 | Python 3.11.9，NumPy 1.26.4，SciPy 1.15.3，TensorFlow 2.15.1，Sionna 1.0.2 |

固定历史设计为：

| 设计 | delay index `j` | 物理时延，ns |
|---|---|---|
| B0/QC | `[0,9,18,27,36,45,54,63]` | `[0,520.833,1041.667,1562.500,2083.333,2604.167,3125.000,3645.833]` |
| S0/Sidon | `[0,1,3,7,12,20,30,65]` | `[0,57.870,173.611,405.093,694.444,1157.407,1736.111,3761.574]` |

## 3. E1：完整协方差下的 residue/lift 联合设计

E1 在 TDL-A 5 ns 完整协方差下评估 203 个唯一候选，底层搜索缓存共评估 221,430 个候选状态，得到 9 个全候选 Pareto 点，并在看到 outage 前冻结 10 个候选。优化等差基线 B1 为：

`j=[0,1,2,3,4,5,6,7]`

`tau=[0,57.870,115.741,173.611,231.481,289.352,347.222,405.093] ns`

B1 的 10%/1% ideal-CSI outage 目标分别为 `10.2078/12.4117 dB`，16/80 dB matched CE NMSE 为 `-22.691/-74.655 dB`。关键候选如下。

| 候选 | delay index `j` | 10% outage 增益 vs B1 | 1% outage 增益 vs B1 | 16 dB CE 劣化 | 80 dB CE 地板劣化 | rank | H1 |
|---|---|---:|---:|---:|---:|---:|---|
| E1_0028，达到 outage 门槛的新联合候选 | `[0,1,2,3,4,5,252,569]` | +0.161 dB | +0.281 dB | +0.688 dB | +3.608 dB | 8 | 失败：80 dB 地板 |
| E1_0006，最接近全部门槛 | `[0,1,2,3,4,5,6,574]` | +0.032 dB | +0.058 dB | +0.274 dB | +1.639 dB | 8 | 失败：outage 与地板 |
| S0 | `[0,1,3,7,12,20,30,65]` | +0.246 dB | +0.430 dB | +3.030 dB | +40.571 dB | 8 | 失败：CE |
| B0/QC | `[0,9,18,27,36,45,54,63]` | +0.060 dB | +0.145 dB | +2.081 dB | +27.632 dB | 8 | 失败：CE |

E1_0028 的物理时延为

`[0,57.870,115.741,173.611,231.481,289.352,14583.333,32928.241] ns`。

因此 H1 的负结论不是“搜索不到 outage 更好的非等差集合”，而是“outage 改善与高 SNR CE 可辨识性之间的联合门槛未被满足”。完整候选及指标见 `outputs/experiment026_cdd_design/20260724_main/e1_search/e1_metrics.csv`，冻结候选的 outage 结果见 `e1_outage/outage_targets.csv`，最终判据表见 `final/e1_h1_candidates.csv`。

## 4. E2：只知道 `T_epsilon` 时的厚 Sidon 潜力

使用 `epsilon=0.01`，`T_margin=T_sync=0`。候选生成只接收 `T_epsilon` 和系统几何；实际 PDP/协方差只在候选冻结后的评价阶段进入。五个预定场景的生成统计为：

| 场景 | `T_epsilon` | 所需 pair-sum gap，index | 所需 fold gap，index | 200,000 个几何样本中的硬候选数 |
|---|---:|---:|---:|---:|
| TDL-A 10 ns | 47.966 ns | 2 | 1 | 17,065 |
| TDL-A 20 ns | 95.932 ns | 4 | 2 | 164 |
| TDL-A 30 ns | 143.898 ns | 5 | 3 | 3 |
| TDL-C 5 ns | 31.533 ns | 2 | 1 | 17,391 |
| TDL-C 10 ns | 63.065 ns | 3 | 2 | 1,515 |

各场景同知识等级的 T1 等差基线均为：

`j=[0,33,87,174,261,348,435,522]`

`tau=[0,1909.722,5034.722,10069.444,15104.167,20138.889,25173.611,30208.333] ns`。

通过 H2 的硬厚 Sidon 候选为：

| 场景 | delay index `j` | 物理时延，ns | 10% outage 增益 vs T1 | 1% outage 增益 vs T1 | 16 dB CE 劣化 |
|---|---|---|---:|---:|---:|
| TDL-A 10 ns | `[0,4,69,198,330,489,496,540]` | `[0,231.481,3993.056,11458.333,19097.222,28298.611,28703.704,31250.000]` | +0.160 dB | +0.234 dB | +0.152 dB |
| TDL-C 5 ns | `[0,3,36,101,320,402,526,567]` | `[0,173.611,2083.333,5844.907,18518.519,23263.889,30439.815,32812.500]` | +0.202 dB | +0.291 dB | +0.026 dB |
| TDL-C 10 ns | `[0,3,42,87,105,165,246,300]` | `[0,173.611,2430.556,5034.722,6076.389,9548.611,14236.111,17361.111]` | +0.174 dB | +0.248 dB | -0.065 dB |

1% ideal-CSI outage 使用同一轮 200,000 个共同信道样本和 0.5 dB SNR 网格的对数域插值；它不是 H2 的预定验收条件，本轮也没有为目标 SNR 计算置信区间。

TDL-A 20/30 ns 的几何池仍分别存在 164/3 个硬候选，但预冻结的 outage 子集没有包含硬候选；其被测试的最佳软候选的 10% 增益为 `+0.116/+0.091 dB`，1% 增益为 `+0.168/+0.124 dB`，不能标记为满足厚 Sidon。故本轮可以确认 H2 已由 A10、C5、C10 三个场景满足，但不能据此判定 A20/A30 的硬厚 Sidon 性能。该覆盖限制是候选冻结策略的限制，不是硬约束不可行证明。

附注：

1. “系统几何”指不读取实际 TDL profile/PDP 时已经确定的离散结构，包括 `K=576` 的循环 DFT delay 栅格、30 kHz 子载波间隔、8 Tx、每 index `57.870370 ns`、DMRS comb-24 的 residue/fold 结构、模 576 的 pair-sum 圆周距离以及候选等价化规则；不包含实际抽头功率或完整协方差。
2. “200,000 个几何样本中的硬候选数”是随机生成候选经公共移位和天线排列等价规范化去重后，同时满足所需 pair-sum gap 与 fold gap 的候选数，在最终最多保留 500 个候选之前统计；它不是 576 点全空间的穷举可行候选总数。
3. E2 表中的增益均相对 T1 几何等差基线 `[0,33,87,174,261,348,435,522]`。E1 的 `[0,1,2,3,4,5,6,7]` 不能直接替代：它在 TDL-A 5 ns 完整协方差和 E1 CE 门槛下选择，而 E2 包含不同 profile/delay spread，且要求发射端选择阶段不知道实际 PDP。若目标是不考虑知识公平性、比较各场景实际 outage 最优的等差集合，应另行计算每场景的 `B1_oracle_arithmetic`；本轮没有完成该 oracle 对比。

生成统计见 `e2_thick_sidon/e2_scenarios.csv`，全部场景指标见 `e2_thick_sidon/e2_metrics.csv`，冻结候选见 `e2_thick_sidon/frozen_outage_candidates.csv`，最终判据见 `final/e2_h2_candidates.csv`。

## 5. E3：近似平坦 TDL 的 estimated-CSI BLER

### 5.1 10% 与 1% BLER 合并结果

| 场景 | 10% QC/Sidon SNR | 10% Sidon 增益，保守 95% 区间 | 1% QC/Sidon SNR | 1% Sidon 增益，保守 95% 区间 |
|---|---:|---:|---:|---:|
| result-024 纯平坦 `h~CN(0,I_8)` | 14.65/14.32 dB | +0.33 dB，[+0.20,+0.45] | 17.23/16.38 dB | +0.85 dB，[+0.63,+1.07] |
| TDL-A 0.1 ns | 14.632/14.329 dB | +0.304 dB，[+0.169,+0.439] | 17.544/16.452 dB | +1.091 dB，[+0.560,+1.623]，QC 明显外推 |
| TDL-A 1 ns | 14.840/14.438 dB | +0.402 dB，[+0.256,+0.547] | 17.809/16.566 dB | +1.242 dB，[+0.846,+1.638]，QC 轻微外推 |
| TDL-A 5 ns | 15.271/15.097 dB | +0.174 dB，[+0.055,+0.292] | 17.606/17.610 dB | -0.003 dB，[-0.231,+0.224] |
| TDL-D 0.1 ns | 19.280/12.704 dB | +6.576 dB，[+6.521,+6.631] | — | 未闭合 |
| TDL-D 1 ns | — | 未闭合 | — | 未闭合 |
| TDL-E 0.1 ns | 19.314/12.652 dB | +6.661 dB，[+6.606,+6.716] | — | 未闭合 |
| TDL-E 1 ns | 19.452/12.748 dB | +6.704 dB，[+6.647,+6.761] | — | 未闭合 |

TDL-A 三点的 10% 增益不是随 RMS delay spread 单调变化，依次为 `0.304/0.402/0.174 dB`；plan 明确不以事后单调性为验收条件。A0.1 的 QC 1% 目标超过精扫上限且目标带累计错误为 0，只能作为先导外推；A1 的 QC 1% 目标轻微超过精扫上限。D1 粗扫在 19.0/19.5 dB 的 QC 错误分别为 63/400 和 41/400，而 Sidon 在 13.0 dB 已为 30/400，因此方向上强烈支持 Sidon 至少约 6.5 dB 的优势，但没有在预定网格内闭合 QC 的 10% 目标，不能给正式点估计和区间。

D/E 没有在预定粗扫范围内同时形成 QC/Sidon 的 1% 跨越，因此不报告 1% 目标。精扫网格在运行前由 QC/Sidon 各自相邻跨越区间的并集自动冻结，见 `e3_refinement_grid.json`。

## 6. CE 与机制诊断

E3 在共同 SNR 下的 QC/Sidon CE NMSE 很接近：TDL-A 0.1 ns 在 14.25 dB 为 `-21.483/-21.429 dB`；TDL-D 0.1 ns 在 15 dB 为 `-22.467/-22.555 dB`；TDL-E 0.1 ns 在 15 dB 为 `-22.729/-22.846 dB`。因此 D/E 的 6.6 dB 量级 BLER分离不是由 Sidon 的 CE NMSE 明显更好造成的。

在 16 dB，TDL-A 0.1/1/5 ns 下 QC 与 Sidon 的 `J_CDD` 分别约为 `1315.6/47.9`、`1233.5/47.9`、`530.3/47.9`；对应有效四阶相关和 `M4` 约为 `27286.9/9143.8`、`27180.4/9121.1`、`24996.3/8660.8`。事实是 Sidon 的高阶相关代理持续较低。

可作出的推断是：在 NLoS 近似平坦极限，人工时延结构重新成为主要频率相关来源，所以严格 Sidon 的优势恢复；在 D/E 中，LoS 确定性分量还会与 CDD 阵列因子相干叠加，可能放大 QC 的频率深衰落和 Sidon 的频率铺展差异。该解释与 CE 和相关代理一致，但本轮没有做 LoS K-factor、理想 CSI 或确定性分量消融，不能把 6.6 dB 唯一归因于某一机制。

## 7. 验证、异常与适用范围

1. Phase-0 preflight 通过：576 点相位/时延单位、24 个唯一导频、5712 data RE、pilot rank、`J_CDD` 与 matched CE 均为有限值；相关 13 个定向测试全部通过。
2. 为避免 Sionna 旧协方差辅助函数在 profile/spread 扫描时构造约 5.75 GiB 的 `23×4096×4096` 临时数组，新增了直接计算 576 个 active subcarrier 协方差的等价实现。A0.1 smoke 重放的 BLER 完全一致，CE NMSE 差异约 `1e-14 dB`。
3. E3 最初以 batch 200 运行 TDL-A 时内存需求约 22.46 GiB，进程在写出结果前失败；正式 A 场景改用 batch 100，trial 数、seed、共同随机数和统计口径不变。
4. 自动精扫网格的早期实现曾错误取两个候选跨越区间的凸包；在正式 refine 前修正为 plan 要求的离散并集并重新固化，正式结果没有使用错误网格。
5. E1/E2 输出是 ideal-CSI 16QAM BICM outage，不是 BLER。H2 只表示候选具有进入后续 estimated-CSI 链路的潜力。
6. CDD 仍作为理想频域相位矩阵施加；最大人工时延超过 CP 不代表本轮模拟了真实时域延迟引起的 ISI/ICI。
7. 结论限于 8×1、48 PRB、当前 DMRS/MCS/LDPC、零速度、matched/oracle 协方差与预定 TDL 场景。

## 8. 证据路径与复现命令

固定证据根目录：

`outputs/experiment026_cdd_design/20260724_main/`

主要证据：

- `validation/environment.json`、`validation/preflight.json`：环境和确定性回执；
- `e1_search/e1_metrics.csv`、`e1_search/e1_pareto.csv`、`e1_search/frozen_outage_candidates.csv`：E1 搜索、二维指标和冻结；
- `e1_outage/outage_curves.csv`、`e1_outage/outage_targets.csv`：E1 200,000 样本 outage；
- `e2_thick_sidon/e2_scenarios.csv`、`e2_metrics.csv`、`frozen_outage_candidates.csv`：E2 `T_epsilon`、几何池和冻结；
- `e2_outage/outage_curves.csv`、`e2_outage/outage_targets.csv`：E2 outage；
- `e3_prescan/`、`e3_refine_10pct/`、`e3_refine_1pct/`：各场景完整 BLER、CE、配对错误计数和展开配置；
- `e3_refinement_grid.json`：正式 refine 前固化的网格；
- `final/final_summary.json`、`e1_h1_candidates.csv`、`e2_h2_candidates.csv`、`e3_targets.csv`：验收汇总；
- `final/candidate_catalog.csv`：所有最终报告候选的时延目录。

核心复现命令：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_cdd_design.py --stage validate --run-id 20260724_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_cdd_design.py --stage e1-search --run-id 20260724_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_cdd_design.py --stage e1-outage --run-id 20260724_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_cdd_design.py --stage e2 --run-id 20260724_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_cdd_design.py --stage e2-outage --run-id 20260724_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage smoke --run-id 20260724_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage prescan --run-id 20260724_main --trials 400 --batch-size 100
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage refine-grid --run-id 20260724_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage refine-10pct --run-id 20260724_main --scenarios A_0p1ns,A_1ns,A_5ns --batch-size 100
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage refine-10pct --run-id 20260724_main --scenarios D_0p1ns,D_1ns,E_0p1ns,E_1ns --batch-size 200
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage refine-1pct --run-id 20260724_main --scenarios A_0p1ns,A_1ns,A_5ns --batch-size 100
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_plan026.py --run-id 20260724_main
```

## 9. 状态

026 的正式执行、分析、成对 result 和索引更新已完成；结论仍待研究者确认。按项目规范，本轮尚未把 026 结论写入 `KNOWLEDGE.md` 或据此改变 `GOALS.md` 的全局阶段判断，也未创建 Git checkpoint。
