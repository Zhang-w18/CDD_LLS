# result-028-text：comb-6 三类 CSI 曲线统一增补与 TDL-A 300 ns 对比（无图版）

> 对应 `research/plan-028.md`。A30/A100 的 estimated-CSI 数据由 result-027 正式点只读合并；A300 的搜索、outage、estimated-CSI 链路和三个场景的真实 ideal-CSI 链路均已完成。2026-08-04 按研究者要求直接在本轮结果上追加 trial 以降低曲线波动，未新增实验 plan。本结果尚待研究者确认。

## 1. 结论与状态

事实：A30、A100、A300 的物理候选数分别为 8、10、10，完整数据仍分别包含 72、85、91 个 `candidate_id + SNR` 点。增量运行保留原 SNR 网格，estimated-CSI 与 CE 使用同一批追加 trial；ideal-CSI 独立按 1%附近的停止规则追加，因此追加后 estimated 与 ideal 的逐点 trial 数不再强制相同。

事实：ideal-CSI 曲线来自真实 TDL-A 信道、payload、LDPC 编码、256QAM、AWGN、真实 data-RE 等效 CSI 均衡、软解调与译码，不是 outage 或互信息曲线。追加合并后三场景全部共同采样点的 ideal BLER 均不高于 estimated BLER；原 A30 `AP_RMS_T1 @ 16.5 dB` 的一错误块有限样本反向不再出现。

事实：两轮增量相对原数据累计追加 A30 `67,400/91,900`、A100 `129,400/302,300`、A300 `22,800/108,100` 个 estimated/ideal candidate-trials，合计 721,900；无 trial 区间重叠或缺口。A100 本次二次追加为 `91,800/217,500`。1%附近（合并 BLER 为 0.5%–2%）的 pointwise Wilson 95%区间中位宽度：A30 estimated `0.00762→0.00578`、ideal `0.01155→0.00532`；A100 estimated `0.00704→0.00374`、ideal `0.00910→0.00390`；A300 estimated `0.00755→0.00582`、ideal `0.01926→0.00505`。

事实：A300 的 99% PDP 支撑宽度为 `1438.98 ns`，硬厚要求为 pair/fold `50/25` 个 q 单位。A300 estimated-CSI 正式闭合者为：10% 目标的 `AP_RMS_T1`、`AP_TALIAS_NT`、`B0_QC`、`S0_SIDON`，以及 1% 目标的 `AP_RMS_T1`、`B0_QC`、`S0_SIDON`。其余候选在冻结的正式扩展范围内未形成双侧 bracket，不外推。

事实：A300 六基线中 `AP_RMS_T1` 在 10%/1% 都最优，目标 SNR 为 `14.152/15.365 dB`。`B0_QC` 分别差 `0.250 dB`（成对 95%区间 `[0.189,0.311]`）和 `0.210 dB`（`[0.080,0.340]`）。固定参考 `S0_SIDON` 分别差 `0.108 dB`（`[0.050,0.167]`）和 `0.156 dB`（`[0.008,0.305]`）。没有搜索 family 在 A300 的 estimated-CSI 正式目标上闭合并优于六基线最优者。

推断：A300 中大跨距候选同时呈现较高 CE NMSE 地板和 BLER 不闭合，支持“当前 comb-6 matched CE 是其链路瓶颈”的解释；这是本轮条件内的诊断关联，不是对所有信道、DMRS 或接收机的普遍因果结论。

## 2. 仿真条件、数据口径与复现

- 系统：48 PRB、$K=576$、30 kHz SCS、4096 FFT、288 CP、10 OFDM symbols、8 Tx / 1 Rx / 1 layer；static Sionna 1.0.2 TDL-A，0 km/h，3.5 GHz。
- DMRS 与负载：symbols `[2,7]`，comb-6，每 symbol 96 pilot、总计 192 pilot RE、5568 data RE；256QAM MCS 8，LDPC 最多 8 次迭代。
- CDD：$V_{k,n}=\exp(-j2\pi k j_n/576)$；相位参考首个 active subcarrier；`1q=57.870370370... ns`，连续时延不量化。
- estimated 接收机：two-DMRS averaged、known-$\mathbf V$/known-PDP、V-aware matched 全带频域 LMMSE；CE NMSE 在 data RE 上先逐 trial 在线性域求比并平均，再转 dB。
- A30/A100：原 estimated 点按“trial 数最大；相同时 `refine_1pct > refine_10pct > prescan`”合并。源 manifest SHA-256 分别为 `b20cca1e3082638ff53a029f3daf034f2aec34664fc140393816acd44097748a`、`33ee2b87b7da18b3d2bc5d1faeccd92b070f105d5c4c9f5c65b084ccd9a07eb1`。
- A300：T2/几何搜索各 200,000 状态；outage 200,000 共同信道样本、1,000 次 bootstrap；smoke 20 trials，prescan 每点 400，正式细扫每点 3,000；seed `20260727`。链路 manifest SHA-256 为 `6ebda32026eaf78138126f4aaf31f13978d648061dee9146969c398942a3a665`。
- 曲线增量：固定入口 `tools/run_bler_curves.py`。首轮配置 `configs/bler_curves_result028_augmentation.yaml` 以累计约 50 个错误块为目标、最多 5,000 trials。A100 二次配置 `configs/bler_curves_result028_a100_more_trials.yaml` 将常规 estimated/ideal 提高到 100 个错误块、最多 10,000 trials，截断分别为 0.1%/0.3%；candidate 02 ideal 通过 `candidate_overrides` 提高到目标 150 个错误块、最多 20,000 trials且不截断。每次最多三轮自适应补齐；seed 仍为 `20260727`，通过绝对 trial 编号确保新旧区间不重叠。
- 展开配置：`outputs/experiment028_csi_curves/20260803_main/a300_comb6/comb6/search/resolved_experiment.json`、`.../outage/resolved_experiment.json`、`.../link/A300/<stage>/resolved_experiment.json` 和三个场景各自的 `ideal_link/resolved_experiment.json`。

核心复现命令：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage import-existing --scenario A30 --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage import-existing --scenario A100 --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage search --scenario A300 --run-id 20260803_main --search-states 200000 --geometry-states 200000
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage outage --scenario A300 --run-id 20260803_main --outage-samples 200000 --bootstrap-repeats 1000
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage manifest --scenario A300 --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage estimated-<smoke|prescan|refine-10pct|refine-1pct> --scenario A300 --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage ideal --scenario <A30|A100|A300> --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage analyze --scenario <A30|A100|A300> --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_augmentation.yaml --stage validate
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_augmentation.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_more_trials.yaml --stage validate
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_more_trials.yaml --stage run
```

## 3. comb-6 三场景统一增补

本章集中交付用户要求的全部新增结果。下表中的 q 为 576 点 DFT 网格坐标，ns 为实际物理时延；显示值为便于阅读的舍入值，精确浮点数、family、represented families、角色、rank 和 condition number 均在各场景 `final/candidate_number_delay_table.csv`。

### 3.1 A30 comb-6

| 编号 | 原缩写 | delay q | delay ns | rank / cond |
|---|---|---|---|---|
| candidate 01 | `A30_B0_QC` | `[0,9,18,27,36,45,54,63]` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` | 8 / 1.000 |
| candidate 02 | `A30_AP_RMS_T1` | `[0,.5184,1.0368,1.5552,2.0736,2.592,3.1104,3.6288]` | `[0,30,60,90,120,150,180,210]` | 8 / 185.588 |
| candidate 03 | `A30_AP_TEPS_T1` | `[0,2.4866,4.9731,7.4597,9.9462,12.4328,14.9193,17.4059]` | `[0,143.898,287.796,431.694,575.592,719.49,863.388,1007.286]` | 8 / 1.225 |
| candidate 04 | `A30_AP_TU_NT` | `[0,72,144,216,288,360,432,504]` | `[0,4166.667,8333.333,12500,16666.667,20833.333,25000,29166.667]` | 4 / `2.254e13` |
| candidate 05 | `A30_AP_TALIAS_NT` | `[0,12,24,36,48,60,72,84]` | `[0,694.444,1388.889,2083.333,2777.778,3472.222,4166.667,4861.111]` | 8 / 1.000 |
| candidate 06 | `A30_S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0,57.87,173.611,405.093,694.444,1157.407,1736.111,3761.574]` | 8 / 1.000 |
| candidate 07 | `A30_AP_T2_00` | `[0,71,142,213,284,355,426,497]` | `[0,4108.796,8217.593,12326.389,16435.185,20543.981,24652.778,28761.574]` | 8 / 1.000 |
| candidate 08 | `A30_MEFF_T2_05` | `[0,49,130,197,303,390,462,524]` | `[0,2835.648,7523.148,11400.463,17534.722,22569.444,26736.111,30324.074]` | 8 / 1.000 |

A30 的 72 个点全部三网格对齐。estimated-CSI 的正式闭合目标沿用 result-027：`AP_RMS` 15.404/18.625、`AP_T2` 14.353/15.917、`AP_TALIAS` 14.190/15.732、`AP_TEPS` 14.202/16.047、`B0` 14.222/15.821、`S0` 14.066/15.619 dB（10%/1%）；`AP_TU`、`MEFF` 未闭合。完整区间、点数与错误计数见原 `target_summary.csv`。

增量合并曲线数据：`outputs/experiment028_csi_curves/20260803_main/curve_augmentation_v2/a30/final/estimated_csi_bler_points.csv`（同时含 CE）和同目录 `ideal_csi_bler_points.csv`；原始点与每个追加区间分别保留在原路径及 `a30/<receiver>/supplemental_points.csv`/`error_flags/`。

### 3.2 A100 comb-6

| 编号 | 原缩写 | delay q | delay ns | rank / cond |
|---|---|---|---|---|
| candidate 01 | `A100_B0_QC` | `[0,9,18,27,36,45,54,63]` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` | 8 / 1.000 |
| candidate 02 | `A100_AP_RMS_T1` | `[0,1.728,3.456,5.184,6.912,8.64,10.368,12.096]` | `[0,100,200,300,400,500,600,700]` | 8 / 1.407 |
| candidate 03 | `A100_AP_TEPS_T1` | `[0,8.2885,16.577,24.8656,33.1541,41.4426,49.7311,58.0197]` | `[0,479.66,959.32,1438.98,1918.64,2398.3,2877.96,3357.62]` | 8 / 1.061 |
| candidate 04 | `A100_AP_TU_NT` | `[0,72,144,216,288,360,432,504]` | `[0,4166.667,8333.333,12500,16666.667,20833.333,25000,29166.667]` | 4 / `2.254e13` |
| candidate 05 | `A100_AP_TU_NTM1` | `[0,82.2857,164.5714,246.8571,329.1429,411.4286,493.7143,576]` | `[0,4761.905,9523.81,14285.714,19047.619,23809.524,28571.429,33333.333]` | 7 / `1.117e13` |
| candidate 06 | `A100_AP_TALIAS_NT` | `[0,12,24,36,48,60,72,84]` | `[0,694.444,1388.889,2083.333,2777.778,3472.222,4166.667,4861.111]` | 8 / 1.000 |
| candidate 07 | `A100_S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0,57.87,173.611,405.093,694.444,1157.407,1736.111,3761.574]` | 8 / 1.000 |
| candidate 08 | `A100_AP_T2_01` | `[0,65,138,211,284,357,430,503]` | `[0,3761.574,7986.111,12210.648,16435.185,20659.722,24884.259,29108.796]` | 8 / 1.000 |
| candidate 09 | `A100_MEFF_T2_04` | `[0,40,115,176,225,324,395,484]` | `[0,2314.815,6655.093,10185.185,13020.833,18750,22858.796,28009.259]` | 8 / 1.000 |
| candidate 10 | `A100_MEFF_T2_06` | `[0,32,157,257,298,334,406,518]` | `[0,1851.852,9085.648,14872.685,17245.37,19328.704,23495.37,29976.852]` | 8 / 1.000 |

A100 的 85 个点全部三网格对齐。estimated-CSI 正式闭合者及 10%/1% 目标为：`AP_RMS` 14.232/16.039、`AP_TALIAS` 14.662/16.100、`AP_TEPS` 14.477/15.928、`AP_TU_NTM1` 14.843/16.505、`B0` 14.504/15.930、`S0` 14.187/15.543 dB；`AP_TU`、`AP_T2` 和两条 `MEFF` 未闭合。

增量合并曲线数据：`outputs/experiment028_csi_curves/20260803_main/curve_augmentation_v2/a100/final/estimated_csi_bler_points.csv` 和同目录 `ideal_csi_bler_points.csv`；追加区间见 `a100/<receiver>/supplemental_points.csv`/`error_flags/`。

A100 二次追加后，candidate 02 ideal 的重点尾点为：

| SNR (dB) | trials | 错误块 | BLER | Wilson 95%区间 |
|---:|---:|---:|---:|---:|
| 16.00 | 20,000 | 100 | 0.00500 | [0.00411, 0.00608] |
| 16.25 | 20,000 | 62 | 0.00310 | [0.00242, 0.00397] |
| 16.50 | 20,000 | 55 | 0.00275 | [0.00211, 0.00358] |
| 17.00 | 20,000 | 26 | 0.00130 | [0.00089, 0.00190] |

相对二次追加前，以上四点的 Wilson 95%区间宽度分别由 `0.00406/0.00313/0.00344/0.00951` 收窄到 `0.00196/0.00155/0.00146/0.00102`；17 dB 从 400 trials、0 错误扩充后得到 26 个真实误块，消除了零错误点造成的宽上界和视觉跳变。

### 3.3 A300 comb-6

| 编号 | 原缩写 | delay q | delay ns | rank / cond |
|---|---|---|---|---|
| candidate 01 | `A300_B0_QC` | `[0,9,18,27,36,45,54,63]` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` | 8 / 1.000 |
| candidate 02 | `A300_AP_RMS_T1` | `[0,5.184,10.368,15.552,20.736,25.92,31.104,36.288]` | `[0,300,600,900,1200,1500,1800,2100]` | 8 / 1.090 |
| candidate 03 | `A300_AP_TEPS_T1` | `[0,24.8656,49.7311,74.5967,99.4623,124.3279,149.1934,174.059]` | `[0,1438.98,2877.96,4316.94,5755.92,7194.9,8633.88,10072.86]` | 8 / 1.119 |
| candidate 04 | `A300_AP_TU_NT` | `[0,72,144,216,288,360,432,504]` | `[0,4166.667,8333.333,12500,16666.667,20833.333,25000,29166.667]` | 4 / `2.254e13` |
| candidate 05 | `A300_AP_TU_NTM1` | `[0,82.2857,164.5714,246.8571,329.1429,411.4286,493.7143,576]` | `[0,4761.905,9523.81,14285.714,19047.619,23809.524,28571.429,33333.333]` | 7 / `1.117e13` |
| candidate 06 | `A300_AP_TALIAS_NT` | `[0,12,24,36,48,60,72,84]` | `[0,694.444,1388.889,2083.333,2777.778,3472.222,4166.667,4861.111]` | 8 / 1.000 |
| candidate 07 | `A300_S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0,57.87,173.611,405.093,694.444,1157.407,1736.111,3761.574]` | 8 / 1.000 |
| candidate 08 | `A300_AP_T2_01` | `[0,65,138,211,284,357,430,503]` | `[0,3761.574,7986.111,12210.648,16435.185,20659.722,24884.259,29108.796]` | 8 / 1.000 |
| candidate 09 | `A300_MEFF_T2_05` | `[0,32,146,215,271,403,448,513]` | `[0,1851.852,8449.074,12442.13,15682.87,23321.759,25925.926,29687.5]` | 8 / 1.000 |
| candidate 10 | `A300_MEFF_T2_06` | `[0,34,149,194,234,303,366,485]` | `[0,1967.593,8622.685,11226.852,13541.667,17534.722,21180.556,28067.13]` | 8 / 1.000 |

正式 estimated-CSI 目标及不确定性：

| 目标 | candidate | SNR (dB) | 95%区间 (dB) | own CE (dB) | 目标带错误数 / trials |
|---|---|---:|---:|---:|---:|
| 10% | `AP_RMS_T1` | 14.152 | [14.109, 14.195] | -18.951 | 808 / 9000 |
| 10% | `AP_TALIAS_NT` | 15.845 | [15.782, 15.908] | -16.102 | 1023 / 9000 |
| 10% | `B0_QC` | 14.402 | [14.358, 14.445] | -18.102 | 797 / 12000 |
| 10% | `S0_SIDON` | 14.260 | [14.221, 14.299] | -18.674 | 960 / 9000 |
| 1% | `AP_RMS_T1` | 15.365 | [15.271, 15.459] | -20.043 | 63 / 12000 |
| 1% | `B0_QC` | 15.575 | [15.486, 15.664] | -19.130 | 55 / 12000 |
| 1% | `S0_SIDON` | 15.521 | [15.407, 15.636] | -19.767 | 112 / 12000 |

10% 未闭合：`AP_T2_01`、`AP_TEPS_T1`、`AP_TU_NT`、`AP_TU_NTM1`、`MEFF_T2_05`、`MEFF_T2_06`；1% 未闭合还包括 `AP_TALIAS_NT`。所有已报告拟合的目标带错误数均不少于 30，统计样本判据通过。

完整证据：原正式目标仍由 `a300_comb6/comb6/link/A300/final/target_summary.csv` 支撑；本次没有后验重定义正式 target fit。增量合并曲线为 `curve_augmentation_v2/a300/final/estimated_csi_bler_points.csv` 和 `ideal_csi_bler_points.csv`，追加区间及 flags 为 `a300/<receiver>/supplemental_points.csv`/`error_flags/`。原 stage、outage、CE 与散点数据保持不变。

### 3.4 曲线降波动增量

增量入口只运行未覆盖的绝对 trial 区间。例如原点为 3,000 trials 时，新样本从 3,001 开始；错误数直接相加，CE NMSE 在线性域按 trial 加权后再转 dB。每个完成区间立即写入 `supplemental_points.csv` 和独立 `.npy` flags；重复调用配置后六类 receiver/scenario 的剩余任务数均为 0。

首轮 ideal 图只使用至首个低于 0.5%的点；A100 二次追加后，常规 candidate 改为使用至首个低于 0.3%的点，candidate 02 则保留全部高 SNR 尾点。A30/A100/A300 最终分别绘制完整 CSV 中的 `48/72`、`45/85`、`33/91` 个 ideal 点。estimated 累计分别追加 67,400、129,400、22,800 candidate-trials；ideal 累计分别追加 91,900、302,300、108,100。

### 3.5 图形一致性与验收

三个场景分别只使用一个 `curve_styles.json`；同一物理方案在 estimated、ideal、CE、缩写图例、编号图例及 A300 散点中的 color/linestyle/marker 一致。图的标题/轴标签/图例不小于 16 pt，刻度不小于 14 pt，线宽不小于 2.5 pt，marker 不小于 8 pt；已按约 13 cm PPT 半页宽预览，坐标、标题、图例与曲线可辨认。零错误 BLER 点只在绘图时用 `0.5/trials` 显示下界，CSV 保留真实 0 和错误计数。

验收状态：plan-028 要求的 A30/A100 三张 candidate 图、A300 三张缩写图、三张 candidate 图、prescan/正式 estimated 图及两个 outage–CE–BLER 散点均已生成；A100 三张 candidate 图已再次替换为二次增量合并版本。真实 ideal 接收机、绝对 trial 续跑、区间连续性、candidate 级预算覆盖、动态 A300 时延、manifest hash、编号双射和样式下限均通过测试或运行时校验。

### 3.6 A100 五 candidate、SNR 不高于 16 dB 的无标题图

本节基于合并后的现有 CSV，只保留 candidate 01、02、04、07、10，不运行新 trial。三张交付图分别为 estimated-CSI BLER、ideal-CSI BLER 和 CE NMSE；均无标题，横轴显示范围为 13.75–16.25 dB，现有数据为 14–16 dB，左右各保留 0.25 dB 空白。竖向主网格为 0.5 dB、次网格为 0.1 dB，图例自动置于坐标区内不与曲线重叠的空白位置。同一 candidate 的 color、linestyle 和 marker 与 result-028 其他图一致。图路径为：

- `docs/figures/result-028/a100_comb6_subset_estimated_csi_bler_snr_le_16.png`；
- `docs/figures/result-028/a100_comb6_subset_ideal_csi_bler_snr_le_16.png`；
- `docs/figures/result-028/a100_comb6_subset_ce_nmse_snr_le_16.png`。

`delay q` 是当前 `/K` 相位定义中的无量纲 DFT 栅格坐标 $j_n$，满足 $V_{k,n}=\exp(-j2\pi k j_n/K)$。整数 q 对应 576 点有效带宽 DFT 栅格上的整数循环移位，`1q=57.87037 ns`；小数 q 是分数循环移位。$K=576$、$N_{\rm FFT}=4096$ 时，同一物理时延对应的 FFT sample 数为

$$
d_n^{\rm FFT}=\tau_nN_{\rm FFT}\Delta f=j_n\frac{N_{\rm FFT}}K=j_n\frac{64}{9}.
$$

该列是物理时延单位换算，不表示把仿真相位分母从 576 改成 4096。当前数字实现逐子载波施加复相位，因而支持任意实数 q，并按 q 模 576 周期等价；有限字长只造成相位量化。若实现只能使用 4096 点时域整数循环 buffer，则必须满足 $d_n^{\rm FFT}\in\mathbb Z$，等价于 q 为 $9/64$ 的整数倍；其他取值需要频域相位旋转或分数时延滤波。数字域循环移位不作为真实传播时延计入 CP。

| 编号 | candidate ID | delay q | delay ns | 等效 delay（FFT samples） |
|---|---|---|---|---|
| candidate 01 | `A100_B0_QC` | `[0,9,18,27,36,45,54,63]` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` | `[0,64,128,192,256,320,384,448]` |
| candidate 02 | `A100_AP_RMS_T1` | `[0,1.728,3.456,5.184,6.912,8.64,10.368,12.096]` | `[0,100,200,300,400,500,600,700]` | `[0,12.288,24.576,36.864,49.152,61.44,73.728,86.016]` |
| candidate 04 | `A100_AP_TU_NT` | `[0,72,144,216,288,360,432,504]` | `[0,4166.667,8333.333,12500,16666.667,20833.333,25000,29166.667]` | `[0,512,1024,1536,2048,2560,3072,3584]` |
| candidate 07 | `A100_S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0,57.87,173.611,405.093,694.444,1157.407,1736.111,3761.574]` | `[0,7.111,21.333,49.778,85.333,142.222,213.333,462.222]` |
| candidate 10 | `A100_MEFF_T2_06` | `[0,32,157,257,298,334,406,518]` | `[0,1851.852,9085.648,14872.685,17245.37,19328.704,23495.37,29976.852]` | `[0,227.556,1116.444,1827.556,2119.111,2375.111,2887.111,3683.556]` |

ideal 图使用与 3.2 相同的 `plot_included` 统计截断。原始 CSV 中仍有以下弱统计尾点，但修订图不再绘制：

- candidate 10 的 15.5/16 dB 均为 400 trials、0 错误；此前对数图以 `0.5/trials=1.25e-3` 显示两个零值，产生了人为平台；
- candidate 07 的 15.75/16 dB 均为 3000 trials、2 错误，经验 BLER 都是 `6.67e-4`，相同有限样本计数产生水平线。

因此不是“没有数据”，而是这些点的错误数不足，且已被正式作图规则排除；candidate 10 和 07 在修订 ideal 图中分别止于 15 和 15.5 dB。

candidate 04 的 CE NMSE 高不是待验证猜测，而是由三组一致证据支持的缺秩结果：

1. comb-6 下 $N_p=K/6=96$，导频相位只由 $j_n\bmod96$ 决定。candidate 04 的折叠位置为 `[0,72,48,24,0,72,48,24]`；第 1/5、2/6、3/7、4/8 列分别相同，所以 $\mathbf V_P$ 从 8 列降为 rank 4，condition number 为 `2.254e13`。
2. 在 matched TDL-A 等效协方差上，以 $10^{-10}$ 相对奇异值为秩门限，完整 $\mathbf R_g$ 为 rank 136，导频协方差 $\mathbf R_{PP}$ 为 rank 68，rank gap 为 68。数据 RE 无噪声预测残差给出 `-2.863 dB` NMSE 地板。
3. 16 dB 的理论 matched-LMMSE NMSE 为 `-2.839 dB`，实际 Monte Carlo 为 `-2.863 dB`；二者都贴近无噪声地板。estimated BLER 在 14–16 dB 均为 1，而 ideal BLER 已下降，故瓶颈是缺秩导频可观测性下的信道估计，而不是理想 CSI 链路。

零错误 ideal 点在对数图中仍显示为 `0.5/trials`，CSV 保留真实零值。复现命令为 `python tools/plot_result028_a100_subset.py`；精确诊断和时延换算分别保存为 `outputs/experiment028_csi_curves/20260803_main/curve_augmentation_v2/a100/final/a100_subset_matrix_diagnostics.csv` 与 `a100_subset_delay_table.csv`。

## 4. 异常、边界与待确认项

- A30/A100 estimated-CSI 的 base 是 result-027 原始正式数据，本次只追加未覆盖 trial；ideal-CSI 同样在 result-028 原数据上追加，均未从 trial 1 重跑后重复相加。
- 正式 target summary 仍采用追加前的预定 fit 口径；新增数据用于降低逐点曲线波动和 pointwise Wilson 区间，不据此后验改变闭合/未闭合规则。
- `AP_TU_NT` rank 4、`AP_TU_NTM1` rank 7 且 condition number 极大；rank 仍只作诊断，不单独替代 CE/BLER 判据。
- A300 未闭合曲线不做 logistic 外推；ideal 曲线也不额外声称目标门限，只报告实际采样点。
- 适用范围仅限当前 48 PRB、static TDL-A、30/100/300 ns、comb-6、零速度、当前 MCS 与 matched receiver；不推广至其他 TDL/CDL、移动性、相关性或失配。
- 本 result 尚待研究者确认，因此未更新 `KNOWLEDGE.md` 或 `GOALS.md` 的已验证结论/阶段状态，也未创建 Git checkpoint。
