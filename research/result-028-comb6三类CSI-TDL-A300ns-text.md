# result-028-text：comb-6 三类 CSI 曲线统一增补与 TDL-A 300 ns 对比（无图版）

> 对应 `research/plan-028-comb6三类CSI-TDL-A300ns.md`。A30/A100 的 estimated-CSI 数据由 result-027 正式点只读合并；A300 的搜索、outage、estimated-CSI 链路和三个场景的真实 ideal-CSI 链路均已完成。2026-08-04 按研究者要求直接追加 trial；2026-08-05 按第 10 节新增两条透明 PRG；2026-08-10 按第 11 节新增五组透明 CDD；2026-08-11 按第 12 节新增一条小时延透明 CDD；2026-08-12 按第 15 节扩展小时延透明 CDD 的 14–15 dB 段，并把 transparent PRG 6 RB 延伸到 1% BLER 以下；随后按第 16 节补充同一小时延 CDD 的 ideal CSI 和 matched-covariance estimated CSI/CE。本结果尚待研究者确认。

## 1. 结论与状态

事实：A30、A100、A300 的物理候选数分别为 8、10、10，完整数据仍分别包含 72、85、91 个 `candidate_id + SNR` 点。增量运行保留原 SNR 网格，estimated-CSI 与 CE 使用同一批追加 trial；ideal-CSI 独立按 1%附近的停止规则追加，因此追加后 estimated 与 ideal 的逐点 trial 数不再强制相同。

事实：A100 主候选仍为上述 10 条，另在第 3.6 节加入 `A100_PRG_DFT8_4RB`、`A100_PRG_DFT8_6RB` 两条透明基线，各含相同的 8 个 SNR 点。新增正式运行完成 84,740 个共同 paired trials，等价于 estimated/ideal 各 169,480 candidate-trials；全部 32 个 BLER 点至少 200 个误块，四条 BLER 曲线均无相邻 SNR 反向。

事实：第 3.6 节再加入与五组 delay set 一一对应的透明 CDD 接收机曲线；发射端不变，UE 对五组均使用 $8\mathbf R_{\rm phy}$，不知道 CDD delay。8 个 SNR 各 10,000 个共同 trials，共 400,000 个 estimated candidate-trials；40 个 BLER 点至少 9,987 个误块，CE 逐 trial 数组和解析失配 NMSE 审计通过，不生成 ideal-CSI 数据。

事实：第 3.6 节小时延透明 CDD `j=[0,0.25,...,1.75]` 已扩为 20 个 SNR，共 298,260 trials；每点至少 200 errors且无相邻反向，estimated BLER 同时闭合 10%与 1%，CE NMSE 为 `-19.368` 至 `-22.632 dB`，全点优于 `-15 dB`。逐 trial CE 数组和 307 个 absolute-trial 区间审计通过，不生成 ideal-CSI 数据。

事实：第 3.6 节 transparent PRG 6 RB 新增 `16.25–17.75 dB` 的 7 个 estimated-CSI/CE 点，共 143,980 trials。BLER 无相邻反向，1%由 17/17.25 dB 的 `0.010984/0.009038` 双侧夹住；CE NMSE 为 `-24.468` 至 `-25.803 dB`。148 个 absolute-trial 区间和全部 CE 数组审计通过，不扩展 ideal-CSI。

事实：同一小时延 CDD 的 matched-covariance estimated CSI 新增 13 个 SNR 点、247,660 trials，BLER 无相邻反向，10%与 1%分别由 `15/15.5 dB` 和 `18/18.5 dB` 双侧夹住，CE NMSE 为 `-24.213` 至 `-29.689 dB`。19.5/20 dB 达到 50,000-trial 上限时分别为 196/102 errors，保留并明确标注样本边界。新增 ideal CSI 共 11 点、170,960 trials，每点至少 200 errors，无相邻反向，1%同样由 `18/18.5 dB` 双侧夹住。

事实：第 3.7 节从第 3.6 节正式数据只读提取 large-delay CDD、small-delay CDD 的透明/非透明接收机和 transparent 6-RB precoder cycling，生成 14–20 dB 专用 estimated-CSI BLER/CE NMSE 子集图和精确 CSV；不新增仿真或插值。

事实：ideal-CSI 曲线来自真实 TDL-A 信道、payload、LDPC 编码、16QAM 调制、AWGN、真实 data-RE 等效 CSI 均衡、软解调与译码，不是 outage 或互信息曲线。追加合并后三场景全部共同采样点的 ideal BLER 均不高于 estimated BLER；原 A30 `AP_RMS_T1 @ 16.5 dB` 的一错误块有限样本反向不再出现。

事实：两轮增量相对原数据累计追加 A30 `67,400/91,900`、A100 `129,400/302,300`、A300 `22,800/108,100` 个 estimated/ideal candidate-trials，合计 721,900；无 trial 区间重叠或缺口。A100 本次二次追加为 `91,800/217,500`。1%附近（合并 BLER 为 0.5%–2%）的 pointwise Wilson 95%区间中位宽度：A30 estimated `0.00762→0.00578`、ideal `0.01155→0.00532`；A100 estimated `0.00704→0.00374`、ideal `0.00910→0.00390`；A300 estimated `0.00755→0.00582`、ideal `0.01926→0.00505`。

事实：A300 的 99% PDP 支撑宽度为 `1438.98 ns`，硬厚要求为 pair/fold `50/25` 个 q 单位。A300 estimated-CSI 正式闭合者为：10% 目标的 `AP_RMS_T1`、`AP_TALIAS_NT`、`B0_QC`、`S0_SIDON`，以及 1% 目标的 `AP_RMS_T1`、`B0_QC`、`S0_SIDON`。其余候选在冻结的正式扩展范围内未形成双侧 bracket，不外推。

事实：A300 六基线中 `AP_RMS_T1` 在 10%/1% 都最优，目标 SNR 为 `14.152/15.365 dB`。`B0_QC` 分别差 `0.250 dB`（成对 95%区间 `[0.189,0.311]`）和 `0.210 dB`（`[0.080,0.340]`）。固定参考 `S0_SIDON` 分别差 `0.108 dB`（`[0.050,0.167]`）和 `0.156 dB`（`[0.008,0.305]`）。没有搜索 family 在 A300 的 estimated-CSI 正式目标上闭合并优于六基线最优者。

推断：A300 中大跨距候选同时呈现较高 CE NMSE 地板和 BLER 不闭合，支持“当前 comb-6 matched CE 是其链路瓶颈”的解释；这是本轮条件内的诊断关联，不是对所有信道、DMRS 或接收机的普遍因果结论。

## 2. 仿真条件、数据口径与复现

- 系统：48 PRB、$K=576$、30 kHz SCS、4096 FFT、288 CP、10 OFDM symbols、8 Tx / 1 Rx / 1 layer；static Sionna 1.0.2 TDL-A，0 km/h，3.5 GHz。
- DMRS 与负载：symbols `[2,7]`，comb-6，每 symbol 96 pilot、总计 192 pilot RE、5568 data RE；16QAM、MCS 8、码率 553/1024（取自 NR 256QAM MCS table），LDPC 最多 8 次迭代。
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

### 3.6 A100 五组 delay set、透明 PRG 与小时延 CDD 接收机对比

原五组曲线基于合并后的现有 CSV，只保留 candidate 01、02、04、07、10，并按本节顺序重新编号为 `delay set 1` 至 `delay set 5`，不重跑旧 trial；随后在相同 A100 链路上正式运行两条透明 PRG precoder cycling、五条“发射端保持原 CDD、UE 不知道 CDD delay”的透明 CDD，以及显式小时延 CDD 的透明、matched-covariance 与 ideal 三种接收机结果。五条 source-mapped 透明 CDD 仍只加入 estimated-CSI BLER 和 CE NMSE；小时延 CDD 则按第 16 节补充 ideal CSI。三张交付图均无标题；ideal 图横轴扩展到 13.75–19.25 dB，estimated 与 CE 图扩展到 13.75–21.25 dB。竖向主网格为 0.5 dB、次网格为 0.1 dB；原五组样式不变，小时延透明、matched 和 ideal 曲线均使用青色，通过不同线型与 marker 区分接收机口径。图路径为：

- `docs/figures/result-028/a100_comb6_subset_estimated_csi_bler_snr_le_21.png`；
- `docs/figures/result-028/a100_comb6_subset_ideal_csi_bler_snr_le_19.png`；
- `docs/figures/result-028/a100_comb6_subset_ce_nmse_snr_le_21.png`。

`delay q` 是当前 `/K` 相位定义中的无量纲 DFT 栅格坐标 $j_n$，满足 $V_{k,n}=\exp(-j2\pi k j_n/K)$。整数 q 对应 576 点有效带宽 DFT 栅格上的整数循环移位，`1q=57.87037 ns`；小数 q 是分数循环移位。$K=576$、$N_{\rm FFT}=4096$ 时，同一物理时延对应的 FFT sample 数为

$$
d_n^{\rm FFT}=\tau_nN_{\rm FFT}\Delta f=j_n\frac{N_{\rm FFT}}K=j_n\frac{64}{9}.
$$

该列是物理时延单位换算，不表示把仿真相位分母从 576 改成 4096。当前数字实现逐子载波施加复相位，因而支持任意实数 q，并按 q 模 576 周期等价；有限字长只造成相位量化。若实现只能使用 4096 点时域整数循环 buffer，则必须满足 $d_n^{\rm FFT}\in\mathbb Z$，等价于 q 为 $9/64$ 的整数倍；其他取值需要频域相位旋转或分数时延滤波。数字域循环移位不作为真实传播时延计入 CP。

| 编号 | candidate ID | delay q | delay ns | 等效 delay（FFT samples） |
|---|---|---|---|---|
| delay set 1 | `A100_B0_QC` | `[0,9,18,27,36,45,54,63]` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` | `[0,64,128,192,256,320,384,448]` |
| delay set 2 | `A100_AP_RMS_T1` | `[0,1.728,3.456,5.184,6.912,8.64,10.368,12.096]` | `[0,100,200,300,400,500,600,700]` | `[0,12.288,24.576,36.864,49.152,61.44,73.728,86.016]` |
| delay set 3 | `A100_AP_TU_NT` | `[0,72,144,216,288,360,432,504]` | `[0,4166.667,8333.333,12500,16666.667,20833.333,25000,29166.667]` | `[0,512,1024,1536,2048,2560,3072,3584]` |
| delay set 4 | `A100_S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0,57.87,173.611,405.093,694.444,1157.407,1736.111,3761.574]` | `[0,7.111,21.333,49.778,85.333,142.222,213.333,462.222]` |
| delay set 5 | `A100_MEFF_T2_06` | `[0,32,157,257,298,334,406,518]` | `[0,1851.852,9085.648,14872.685,17245.37,19328.704,23495.37,29976.852]` | `[0,227.556,1116.444,1827.556,2119.111,2375.111,2887.111,3683.556]` |

| 图例 | candidate ID | PRG 划分 | DFT 向量映射 |
|---|---|---|---|
| transparent PRG 4 RB | `A100_PRG_DFT8_4RB` | 12 个 PRG，每个 48 子载波 | 前 8 个固定 `[0,1,2,3,4,5,6,7]`；尾部 4 个在每个 absolute trial 从 8 个向量中有序无放回重抽，trial 内跨全部 symbol、DMRS/data、ideal/estimated 固定 |
| transparent PRG 6 RB | `A100_PRG_DFT8_6RB` | 8 个 PRG，每个 72 子载波 | 固定 `[0,1,2,3,4,5,6,7]` |

透明基线实际发射向量为 $V_{n,m}=\exp(-j2\pi nm/8)$，每列功率为 8；DMRS 与 data 使用相同的 PRG 向量，并继续采用 `noise_variance=8/SNR_linear`。对 PRG $b$，等效标量信道为 $g_k=\sum_n h_{k,n}V_{n,m_b}$。在 8 根天线独立且 PDP 相同的当前 A100 模型下，PRG 内协方差为 $\mathbf R_g=8\mathbf R_{\rm phy}$，所以 PDP 形状与底层物理信道相同，只差功率因子 8。两个 static DMRS 的 LS 观测先平均，噪声方差为 $4/SNR_{\rm linear}$；接收机随后逐 PRG 独立计算

$$
\widehat{\mathbf g}_{D_b}
=8\mathbf R_{{\rm phy},D_bP_b}
\left(8\mathbf R_{{\rm phy},P_bP_b}+\frac{4}{SNR_{\rm linear}}\mathbf I\right)^{-1}
\bar{\mathbf z}_b,
$$

不跨 PRG 插值，也不联合非相邻、重复使用同一向量的 PRG。CE NMSE 与原五组的统计方法相同：每个 trial 先在全部 data RE 上计算 $\mathrm{NMSE}_t=\sum|\widehat g-g|^2/\sum|g|^2$，再在线性域跨 trial 求平均，最后转换为 dB；不是 PRG 等权平均、总能量比或 dB 域平均。

五条透明 CDD 曲线与原 delay set 一一对应：

| 图例 | transparent candidate ID | source CDD candidate ID |
|---|---|---|
| delay set 1, transparent CDD | `A100_B0_QC_TRANSPARENT_CDD` | `A100_B0_QC` |
| delay set 2, transparent CDD | `A100_AP_RMS_T1_TRANSPARENT_CDD` | `A100_AP_RMS_T1` |
| delay set 3, transparent CDD | `A100_AP_TU_NT_TRANSPARENT_CDD` | `A100_AP_TU_NT` |
| delay set 4, transparent CDD | `A100_S0_SIDON_TRANSPARENT_CDD` | `A100_S0_SIDON` |
| delay set 5, transparent CDD | `A100_MEFF_T2_06_TRANSPARENT_CDD` | `A100_MEFF_T2_06` |

透明 CDD 的发射端仍严格使用上表 source candidate 的 $V_{k,n}=\exp(-j2\pi k j_n/576)$，真实等效信道、DMRS 和 data 均为 $g_k=\sum_nV_{k,n}h_{k,n}$，噪声方差仍为 $8/SNR_{\rm linear}$。区别只在接收端：UE 不知道 $j_n$，五条曲线在同一 SNR 下共用同一个假设协方差 $\mathbf R_a=8\mathbf R_{\rm phy}$ 和同一个全带 LMMSE 滤波器

$$
\mathbf W_a=\mathbf R_{a,DP}
\left(\mathbf R_{a,PP}+\frac{4}{SNR_{\rm linear}}\mathbf I\right)^{-1},
\qquad
\widehat{\mathbf g}_D=\mathbf W_a\bar{\mathbf z}_P.
$$

这里保留功率因子 8，是为了只测未知 CDD delay 造成的协方差形状失配，不额外混入总接收功率标定误差。误差真值仍是每条 source CDD 的真实等效信道；CE 仍按每个 trial 的 $\sum|\widehat g-g|^2/\sum|g|^2$ 在线性域平均，不能用底层物理信道代替 $g$。

正式运行使用相同的 8 点 SNR 网格 `[14,14.25,14.5,14.75,15,15.5,15.75,16]` dB，共完成 84,740 个共同 paired trials，即每个接收机 169,480 candidate-trials。14–15.5 dB 各 10,000 trials；15.75/16 dB 为使四个 candidate×receiver 点均达到至少 200 个误块，自适应追加到 11,180/13,560 trials。精确结果如下；“errors / BLER”保留原始误块数，CE 区间是 trial 级线性 NMSE 均值的 95% Monte Carlo 区间转 dB：

| PRG | SNR dB | paired trials | estimated errors / BLER | ideal errors / BLER | CE NMSE dB | CE 95% CI dB |
|---|---:|---:|---:|---:|---:|---:|
| 4 RB | 14 | 10000 | 1858 / 0.1858 | 1157 / 0.1157 | -21.669 | [-21.702, -21.637] |
| 4 RB | 14.25 | 10000 | 1420 / 0.142 | 861 / 0.0861 | -21.899 | [-21.932, -21.867] |
| 4 RB | 14.5 | 10000 | 1263 / 0.1263 | 725 / 0.0725 | -22.101 | [-22.134, -22.069] |
| 4 RB | 14.75 | 10000 | 978 / 0.0978 | 575 / 0.0575 | -22.304 | [-22.337, -22.271] |
| 4 RB | 15 | 10000 | 792 / 0.0792 | 450 / 0.045 | -22.521 | [-22.553, -22.489] |
| 4 RB | 15.5 | 10000 | 473 / 0.0473 | 248 / 0.0248 | -22.970 | [-23.002, -22.938] |
| 4 RB | 15.75 | 11180 | 421 / 0.037657 | 201 / 0.017979 | -23.178 | [-23.209, -23.148] |
| 4 RB | 16 | 13560 | 412 / 0.030383 | 201 / 0.014823 | -23.375 | [-23.403, -23.348] |
| 6 RB | 14 | 10000 | 1807 / 0.1807 | 1218 / 0.1218 | -22.469 | [-22.504, -22.435] |
| 6 RB | 14.25 | 10000 | 1436 / 0.1436 | 990 / 0.099 | -22.709 | [-22.743, -22.674] |
| 6 RB | 14.5 | 10000 | 1232 / 0.1232 | 796 / 0.0796 | -22.961 | [-22.996, -22.926] |
| 6 RB | 14.75 | 10000 | 1003 / 0.1003 | 620 / 0.062 | -23.133 | [-23.168, -23.099] |
| 6 RB | 15 | 10000 | 845 / 0.0845 | 520 / 0.052 | -23.344 | [-23.379, -23.309] |
| 6 RB | 15.5 | 10000 | 513 / 0.0513 | 315 / 0.0315 | -23.824 | [-23.859, -23.789] |
| 6 RB | 15.75 | 11180 | 467 / 0.041771 | 290 / 0.025939 | -24.035 | [-24.067, -24.003] |
| 6 RB | 16 | 13560 | 470 / 0.034661 | 258 / 0.019027 | -24.265 | [-24.294, -24.236] |

稳定性与可恢复性审计通过：四条 BLER 的相邻 SNR 反向数为 0，所有 32 个 BLER 点均不少于 200 个误块；87 个 absolute-trial 区间在两基线、ideal/estimated 间完全配对且连续无缺口。4-RB 的 84,740 个尾部映射全部按 seed 精确重放、每个四元组内部无重复，共出现 1,680 种有序组合。169,480 个 estimated candidate-trials 的 CE 逐 trial 比值、和与平方和全部一致，因此后续可从各点当前 `trial_end+1` 直接追加，无需重跑旧 trial。

第 15 节只延伸 6-RB transparent PRG 的 estimated-CSI/CE，原 16 dB 及更低 SNR 点和全部 ideal-CSI 点保持不变。200-trial prescan 在 `16.5/17/17.5/18 dB` 得到 BLER `0.03/0.015/0.005/0`，随后在正式 trial 1 前冻结 `16.25–17.75 dB`、间隔 0.25 dB 的 7 点网格。正式结果为：

| SNR dB | trials | estimated errors / BLER | BLER Wilson 95% CI | CE NMSE dB | CE 95% CI dB |
|---:|---:|---:|---:|---:|---:|
| 16.25 | 10000 | 237 / 0.023700 | [0.020896, 0.026869] | -24.468 | [-24.502, -24.433] |
| 16.50 | 10000 | 209 / 0.020900 | [0.018275, 0.023893] | -24.718 | [-24.754, -24.684] |
| 16.75 | 14220 | 200 / 0.014065 | [0.012256, 0.016136] | -24.928 | [-24.957, -24.899] |
| 17.00 | 18300 | 201 / 0.010984 | [0.009573, 0.012600] | -25.132 | [-25.158, -25.107] |
| 17.25 | 22460 | 203 / 0.009038 | [0.007882, 0.010363] | -25.348 | [-25.371, -25.325] |
| 17.50 | 31000 | 203 / 0.006548 | [0.005710, 0.007509] | -25.573 | [-25.593, -25.554] |
| 17.75 | 38000 | 200 / 0.005263 | [0.004584, 0.006042] | -25.803 | [-25.821, -25.785] |

7 点累计 143,980 trials，每点至少 10,000 trials 和 200 errors，相邻 BLER 反向数为 0；17/17.25 dB 分别为 `0.010984/0.009038`，因此在实测点上双侧闭合 1%。CE NMSE 从 `-24.468 dB` 降至 `-25.803 dB`。148 个 absolute-trial 区间逐 SNR 从 trial 1 连续、无重叠，全部 143,980 组 CE trial arrays 已验证；后续仍可从各点 `trial_end+1` 继续追加。该延伸未构造或输出 ideal receiver 数据。

透明 CDD 正式运行在每个 SNR 使用 10,000 个共同 trials，共 80,000 个 SNR-level trials、400,000 个 estimated candidate-trials；不生成 ideal-CSI 数据。40 个 BLER 点的误块数均不少于 9,987，故全部在最小 trial 预算处满足停止规则。精确结果如下；CE 区间仍是 trial 级线性 NMSE 均值的 95% Monte Carlo 区间转 dB：

| delay set | SNR dB | estimated errors / BLER | CE NMSE dB | CE 95% CI dB |
|---:|---:|---:|---:|---:|
| 1 | 14 | 10000 / 1 | -0.942 | [-0.952, -0.932] |
| 1 | 14.25 | 10000 / 1 | -0.949 | [-0.959, -0.939] |
| 1 | 14.5 | 10000 / 1 | -0.940 | [-0.949, -0.930] |
| 1 | 14.75 | 10000 / 1 | -0.957 | [-0.967, -0.948] |
| 1 | 15 | 10000 / 1 | -0.950 | [-0.959, -0.940] |
| 1 | 15.5 | 10000 / 1 | -0.964 | [-0.974, -0.954] |
| 1 | 15.75 | 10000 / 1 | -0.951 | [-0.961, -0.941] |
| 1 | 16 | 10000 / 1 | -0.962 | [-0.972, -0.953] |
| 2 | 14 | 9998 / 0.9998 | -5.193 | [-5.225, -5.161] |
| 2 | 14.25 | 9995 / 0.9995 | -5.214 | [-5.246, -5.182] |
| 2 | 14.5 | 9996 / 0.9996 | -5.192 | [-5.224, -5.160] |
| 2 | 14.75 | 9991 / 0.9991 | -5.234 | [-5.266, -5.201] |
| 2 | 15 | 9997 / 0.9997 | -5.236 | [-5.268, -5.204] |
| 2 | 15.5 | 9996 / 0.9996 | -5.277 | [-5.310, -5.245] |
| 2 | 15.75 | 9996 / 0.9996 | -5.256 | [-5.289, -5.224] |
| 2 | 16 | 9987 / 0.9987 | -5.312 | [-5.344, -5.279] |
| 3 | 14 | 10000 / 1 | 0.060 | [0.050, 0.071] |
| 3 | 14.25 | 10000 / 1 | 0.045 | [0.034, 0.055] |
| 3 | 14.5 | 10000 / 1 | 0.059 | [0.049, 0.070] |
| 3 | 14.75 | 10000 / 1 | 0.044 | [0.033, 0.054] |
| 3 | 15 | 10000 / 1 | 0.058 | [0.047, 0.068] |
| 3 | 15.5 | 10000 / 1 | 0.057 | [0.046, 0.067] |
| 3 | 15.75 | 10000 / 1 | 0.063 | [0.052, 0.073] |
| 3 | 16 | 10000 / 1 | 0.062 | [0.052, 0.072] |
| 4 | 14 | 10000 / 1 | -2.803 | [-2.823, -2.782] |
| 4 | 14.25 | 10000 / 1 | -2.826 | [-2.847, -2.806] |
| 4 | 14.5 | 10000 / 1 | -2.798 | [-2.819, -2.778] |
| 4 | 14.75 | 10000 / 1 | -2.837 | [-2.858, -2.816] |
| 4 | 15 | 10000 / 1 | -2.825 | [-2.846, -2.805] |
| 4 | 15.5 | 10000 / 1 | -2.841 | [-2.862, -2.821] |
| 4 | 15.75 | 10000 / 1 | -2.813 | [-2.833, -2.793] |
| 4 | 16 | 10000 / 1 | -2.827 | [-2.848, -2.807] |
| 5 | 14 | 10000 / 1 | -0.503 | [-0.511, -0.496] |
| 5 | 14.25 | 10000 / 1 | -0.512 | [-0.520, -0.504] |
| 5 | 14.5 | 10000 / 1 | -0.501 | [-0.509, -0.494] |
| 5 | 14.75 | 10000 / 1 | -0.510 | [-0.518, -0.502] |
| 5 | 15 | 10000 / 1 | -0.499 | [-0.507, -0.492] |
| 5 | 15.5 | 10000 / 1 | -0.504 | [-0.511, -0.496] |
| 5 | 15.75 | 10000 / 1 | -0.493 | [-0.501, -0.486] |
| 5 | 16 | 10000 / 1 | -0.501 | [-0.508, -0.493] |

透明 CDD 的 CE 曲线在 14–16 dB 基本保持平台，说明当前范围由协方差形状失配主导：五组全网格范围分别为 `[-0.964,-0.940]`、`[-5.312,-5.192]`、`[0.044,0.063]`、`[-2.841,-2.798]`、`[-0.512,-0.493]` dB。逐 SNR 的解析失配 trace-NMSE 与 Monte Carlo trial-ratio 均值最大相差 `0.063 dB`。BLER 只有 delay set 2 出现成功块，且仍为 `0.9987–0.9998`；其 14.5 和 15 dB 相对前一点分别有 `0.0001/0.0006` 的反向，Wilson 区间重叠且图上不可见，因此未追加或平滑。其余 32 个透明 CDD BLER 点均精确为 1。

80 个 1000-trial absolute 区间在五条曲线间完全一致、逐 SNR 从 1 连续到 10,000 且无重叠。400,000 个 CE trial arrays 已逐项验证 `r_t=E_t/S_t` 以及区间和、平方和；后续可把最小 trial 或目标误块数提高后从 10,001 继续追加。两次外层 1 小时执行窗口结束时，未落盘的 15.5/16 dB 区间被丢弃并按相同 absolute seed 重算，不进入累计数据。

新增小时延基线 `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` 使用显式冻结的 $\mathbf j=[0,0.25,0.5,0.75,1,1.25,1.5,1.75]$，对应人工时延 `[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273] ns`。发射仍采用 $V_{k,n}=\exp(-j2\pi k j_n/576)$、总功率 8；UE 仍不知道 CDD delay，只使用 $8\mathbf R_{\rm phy}$ 的全带失配 LMMSE。101.273 ns 的最大人工时延明显小于五个原 delay set，用来降低人工频率选择性。第 12/15 节当时未运行 ideal CSI；第 16 节现已在相同发射候选上补充正式 ideal-BLER，结果见下文。

解析扫描首先选择的 $\Delta q=0.5$ 候选虽然 CE NMSE 有余量，但 13–19 dB 的 200-trial 链路 prescan 只能从 BLER `0.49` 降至 `0.075`，无法闭合 1%，因此在正式 trial 1 前淘汰。$\Delta q=0.25$ 在 prescan 中同时闭合 10%与 1%，且比同样合格的 $\Delta q=0.1$ 保留更大非零 delay span，故冻结为正式候选。第 12 节原网格为 `[15.5,15.75,16,16.25,16.5,17,17.5,18,18.5,19,19.25,19.5,19.75,20,20.25,20.5,21]` dB；第 15 节在低 SNR 侧补充 `[14,14.5,15]` dB。完整正式网格共 20 点，只运行 estimated CSI 与 CE，不生成 ideal 数据。

| SNR dB | trials | estimated errors / BLER | CE NMSE dB | CE 95% CI dB |
|---:|---:|---:|---:|---:|
| 14.0 | 10000 | 3023 / 0.302300 | -19.368 | [-19.425, -19.311] |
| 14.5 | 10000 | 2367 / 0.236700 | -19.732 | [-19.788, -19.676] |
| 15.0 | 10000 | 1993 / 0.199300 | -19.933 | [-19.995, -19.873] |
| 15.5 | 10000 | 1435 / 0.143500 | -20.321 | [-20.382, -20.262] |
| 15.75 | 10000 | 1281 / 0.128100 | -20.369 | [-20.434, -20.305] |
| 16.0 | 10000 | 1159 / 0.115900 | -20.493 | [-20.556, -20.431] |
| 16.25 | 10000 | 1011 / 0.101100 | -20.639 | [-20.702, -20.577] |
| 16.5 | 10000 | 870 / 0.087000 | -20.776 | [-20.842, -20.712] |
| 17.0 | 10000 | 653 / 0.065300 | -20.998 | [-21.065, -20.932] |
| 17.5 | 10000 | 477 / 0.047700 | -21.232 | [-21.300, -21.165] |
| 18.0 | 10000 | 367 / 0.036700 | -21.467 | [-21.537, -21.399] |
| 18.5 | 10000 | 242 / 0.024200 | -21.709 | [-21.778, -21.641] |
| 19.0 | 12040 | 201 / 0.016694 | -21.903 | [-21.969, -21.839] |
| 19.25 | 12400 | 200 / 0.016129 | -21.895 | [-21.960, -21.830] |
| 19.5 | 16340 | 202 / 0.012362 | -22.073 | [-22.129, -22.017] |
| 19.75 | 19100 | 201 / 0.010524 | -22.162 | [-22.215, -22.110] |
| 20.0 | 21780 | 200 / 0.009183 | -22.257 | [-22.307, -22.207] |
| 20.25 | 25260 | 202 / 0.007997 | -22.357 | [-22.404, -22.311] |
| 20.5 | 27100 | 200 / 0.007380 | -22.430 | [-22.477, -22.383] |
| 21.0 | 44240 | 200 / 0.004521 | -22.632 | [-22.668, -22.596] |

第 15 节新增三点的 BLER Wilson 95%区间分别为：14 dB `[0.293376,0.311376]`、14.5 dB `[0.228471,0.245131]`、15 dB `[0.191586,0.207245]`。

正式累计共 298,260 trials。每点至少 10,000 trials，之后逐点追加到至少 200 errors；21 dB 最终为 44,240 trials。20 个 BLER 点相邻反向数为 0；10%由 16.25/16.5 dB 的 `0.1011/0.0870` 双侧夹住，1%由 19.75/20 dB 的 `0.010524/0.009183` 双侧夹住，因此曲线同时闭合两个目标且没有平滑或修改原始点。CE NMSE 全网格为 `-19.368` 至 `-22.632 dB`，全部优于 `-15 dB`；解析 trace-NMSE 与 Monte Carlo trial-ratio 均值最大差 `1.046 dB`。

307 个 absolute-trial 区间逐 SNR 从 trial 1 连续、无重叠，区间长度均不超过 1,000；298,260 组 CE trial arrays 已逐项验证 $r_t=E_t/S_t$ 及区间和、平方和。后续可从各点当前 `trial_end+1` 继续追加，不需重跑现有 trials。

第 16 节保持上述同一发射候选不变，只增加两个接收机口径：ideal CSI 在 data RE 直接使用真实等效信道；matched-covariance estimated CSI 假设 UE 知道小时延 CDD，并用真实等效协方差

$$
\mathbf R_g=\mathbf R_{\rm phy}\odot(\mathbf V\mathbf V^H)
$$

构造全带 LMMSE。它与透明接收机的唯一区别是协方差知识由 $8\mathbf R_{\rm phy}$ 改为 $\mathbf R_g$；噪声、DMRS、data、信道 realizations 的生成规则、功率归一化和每 trial CE NMSE 定义均不改变。matched estimated 正式结果如下：

| SNR dB | trials | estimated errors / BLER | BLER Wilson 95% CI | CE NMSE dB | CE 95% CI dB |
|---:|---:|---:|---:|---:|---:|
| 14.0 | 10000 | 2353 / 0.235300 | [0.227089, 0.243715] | -24.213 | [-24.266, -24.162] |
| 14.5 | 10000 | 1778 / 0.177800 | [0.170430, 0.185417] | -24.709 | [-24.761, -24.657] |
| 15.0 | 10000 | 1432 / 0.143200 | [0.136472, 0.150202] | -25.048 | [-25.103, -24.994] |
| 15.5 | 10000 | 992 / 0.099200 | [0.093494, 0.105214] | -25.641 | [-25.693, -25.589] |
| 16.0 | 10000 | 786 / 0.078600 | [0.073486, 0.084038] | -25.996 | [-26.050, -25.943] |
| 16.5 | 10000 | 504 / 0.050400 | [0.046282, 0.054863] | -26.506 | [-26.559, -26.454] |
| 17.0 | 10000 | 355 / 0.035500 | [0.032048, 0.039309] | -26.938 | [-26.991, -26.886] |
| 17.5 | 10000 | 240 / 0.024000 | [0.021178, 0.027188] | -27.395 | [-27.448, -27.341] |
| 18.0 | 12540 | 201 / 0.016029 | [0.013974, 0.018380] | -27.875 | [-27.923, -27.828] |
| 18.5 | 21780 | 206 / 0.009458 | [0.008256, 0.010833] | -28.295 | [-28.331, -28.260] |
| 19.0 | 33340 | 202 / 0.006059 | [0.005281, 0.006951] | -28.784 | [-28.813, -28.755] |
| 19.5 | 50000 | 196 / 0.003920 | [0.003409, 0.004507] | -29.225 | [-29.249, -29.202] |
| 20.0 | 50000 | 102 / 0.002040 | [0.001681, 0.002476] | -29.689 | [-29.713, -29.666] |

累计 247,660 trials、250 个 absolute-trial 区间，13 个 BLER 点相邻反向数为 0。10%由 15/15.5 dB 的 `0.1432/0.0992` 双侧夹住，1%由 18/18.5 dB 的 `0.016029/0.009458` 双侧夹住。CE NMSE 从 `-24.213 dB` 单调改善至 `-29.689 dB`；247,660 组逐 trial CE 数组和线性域充分统计量完全一致，解析 matched trace-NMSE 与 Monte Carlo 最大相差 `1.166 dB`。19.5/20 dB 已达到每点 50,000-trial 上限，但只有 196/102 errors；这两点的 Wilson 区间已给出，不能声称达到预定 200-error 精度，不过不影响 18/18.5 dB 对 1% 的双侧闭合。

同一发射候选的 ideal-CSI 正式结果为：

| SNR dB | trials | ideal errors / BLER | BLER Wilson 95% CI |
|---:|---:|---:|---:|
| 14.0 | 10000 | 2028 / 0.202800 | [0.195034, 0.210794] |
| 14.5 | 10000 | 1527 / 0.152700 | [0.145783, 0.159883] |
| 15.0 | 10000 | 1220 / 0.122000 | [0.115730, 0.128560] |
| 15.5 | 10000 | 838 / 0.083800 | [0.078528, 0.089392] |
| 16.0 | 10000 | 630 / 0.063000 | [0.058404, 0.067932] |
| 16.5 | 10000 | 407 / 0.040700 | [0.037000, 0.044752] |
| 17.0 | 10000 | 261 / 0.026100 | [0.023152, 0.029412] |
| 17.5 | 11000 | 211 / 0.019182 | [0.016781, 0.021918] |
| 18.0 | 17000 | 201 / 0.011824 | [0.010305, 0.013562] |
| 18.5 | 26960 | 201 / 0.007455 | [0.006496, 0.008555] |
| 19.0 | 46000 | 200 / 0.004348 | [0.003787, 0.004992] |

ideal 累计 170,960 trials、172 个连续区间，每点至少 200 errors且无相邻反向；10%由 15/15.5 dB 的 `0.1220/0.0838` 夹住，1%由 18/18.5 dB 的 `0.011824/0.007455` 夹住。全部共同 SNR 点均满足 ideal BLER 不高于 matched estimated BLER。与透明小时延接收机相比，matched estimated 的 1%夹点从 19.75/20 dB 前移到 18/18.5 dB；这是当前固定发射、固定链路下消除协方差形状失配后的接收机收益，不推广为其他 CDD 或移动信道的普遍增益。

第 16 节源 manifest SHA-256 仍为 `33ee2b87b7da18b3d2bc5d1faeccd92b070f105d5c4c9f5c65b084ccd9a07eb1`。matched 策略 manifest / 最终配置 SHA-256 分别为 `7cfc862c7b268dc8006ef96fae7b100814890b4eadf0056e8d0f4688cabffed9` / `591e3a718795ef10a73a457b96641d7b5351b1a544ded5d854539d0e68a7d3d0`；ideal 分别为 `564c3713f67c06052be001a75d22f2c507e6e7eb8f74b2524b881627f1da4091` / `c77351e0a0dfa2446b58cdd8c148e313f6adba27d58ce4518b6e179404dd676e`。精确证据位于：

- `outputs/experiment028_csi_curves/20260803_main/small_delay_cdd_matched/a100/final/`：`estimated_csi_bler_points.csv`、`ce_nmse_uncertainty.csv`、`stability_audit.csv`、`bler_bracket_audit.csv`、`analytic_vs_monte_carlo_ce.csv`、`analysis_summary.json`；
- `outputs/experiment028_csi_curves/20260803_main/small_delay_cdd_ideal/a100/final/`：`ideal_csi_bler_points.csv`、`stability_audit.csv`、`bler_bracket_audit.csv`、`analysis_summary.json`；
- estimated/CE 同图精确 CSV 为 `outputs/experiment028_csi_curves/20260803_main/transparent_cdd_small_delay/a100/final/a100_subset_with_small_delay_transparent_cdd_<estimated|ce>_points.csv`，ideal 同图精确 CSV 为 `outputs/experiment028_csi_curves/20260803_main/transparent_prg_baselines/a100/final/a100_subset_with_transparent_ideal_points.csv`。

ideal 图使用与 3.2 相同的 `plot_included` 统计截断。原始 CSV 中仍有以下弱统计尾点，但修订图不再绘制：

- delay set 5（原 candidate 10）的 15.5/16 dB 均为 400 trials、0 错误；此前对数图以 `0.5/trials=1.25e-3` 显示两个零值，产生了人为平台；
- delay set 4（原 candidate 07）的 15.75/16 dB 均为 3000 trials、2 错误，经验 BLER 都是 `6.67e-4`，相同有限样本计数产生水平线。

因此不是“没有数据”，而是这些点的错误数不足，且已被正式作图规则排除；delay set 5 和 4 在修订 ideal 图中分别止于 15 和 15.5 dB。

delay set 3（原 candidate 04）的 CE NMSE 高不是待验证猜测，而是由三组一致证据支持的缺秩结果：

1. comb-6 下 $N_p=K/6=96$，导频相位只由 $j_n\bmod96$ 决定。delay set 3 的折叠位置为 `[0,72,48,24,0,72,48,24]`；第 1/5、2/6、3/7、4/8 列分别相同，所以 $\mathbf V_P$ 从 8 列降为 rank 4，condition number 为 `2.254e13`。
2. 在 matched TDL-A 等效协方差上，以 $10^{-10}$ 相对奇异值为秩门限，完整 $\mathbf R_g$ 为 rank 136，导频协方差 $\mathbf R_{PP}$ 为 rank 68，rank gap 为 68。数据 RE 无噪声预测残差给出 `-2.863 dB` NMSE 地板。
3. 16 dB 的理论 matched-LMMSE NMSE 为 `-2.839 dB`，实际 Monte Carlo 为 `-2.863 dB`；二者都贴近无噪声地板。estimated BLER 在 14–16 dB 均为 1，而 ideal BLER 已下降，故瓶颈是缺秩导频可观测性下的信道估计，而不是理想 CSI 链路。

零错误 ideal 点在对数图中仍显示为 `0.5/trials`，CSV 保留真实零值。透明基线源 A100 manifest SHA-256 为 `33ee2b87b7da18b3d2bc5d1faeccd92b070f105d5c4c9f5c65b084ccd9a07eb1`，透明策略 manifest SHA-256 为 `224ff43bf05aba42c5c9726cea8de12301f1df73479c54261b397f5b3883c7aa`，最终配置 SHA-256 为 `281069610f6b5b497fa404eb5e6668bc744bad4c15517c32dad4abbf95014ece`。正式证据路径为：

- `outputs/experiment028_csi_curves/20260803_main/transparent_prg_baselines/a100/final/estimated_csi_bler_points.csv`；
- `.../final/ideal_csi_bler_points.csv`、`stability_audit.csv`、`ce_nmse_uncertainty.csv`、`prg_mapping_audit.json`、`analysis_summary.json`；
- `.../final/a100_subset_with_transparent_<estimated|ideal|ce>_points.csv`；
- 逐 trial error flags、CE 三类数组、4-RB 尾部映射、filter diagnostics、resolved run 与进度日志均保存在同一 `a100/` 目录。

原五组精确诊断和时延换算仍保存为 `outputs/experiment028_csi_curves/20260803_main/curve_augmentation_v2/a100/final/a100_subset_matrix_diagnostics.csv` 与 `a100_subset_delay_table.csv`。透明基线复现命令为：

透明 CDD 使用同一个源 A100 manifest；接收机策略 manifest SHA-256 为 `89c89236e62f9983279436f2017b9a00de78c08f34ab72ab846c3f41cb7b2aa7`，最终配置 SHA-256 为 `b15bde6ff36f4d375f7b53fb9d535e0f75de6f0a5a44ee9a56dcc337f253b245`。正式证据位于 `outputs/experiment028_csi_curves/20260803_main/transparent_cdd_physical_covariance/a100/`：

- `final/estimated_csi_bler_points.csv`、`stability_audit.csv`、`ce_nmse_uncertainty.csv`；
- `final/analytic_vs_monte_carlo_ce.csv`、`source_mapping_audit.json`、`analysis_summary.json`；
- `final/a100_subset_with_transparent_cdd_<estimated|ce>_points.csv`；
- 逐 trial error flags、CE 三类数组、filter diagnostics、resolved run、80 个区间的进度和后台状态/日志均保存在同一正式目录；该目录不存在 ideal receiver 或 ideal-CSI CSV。

小时延透明 CDD 使用同一个源 manifest；扩展后显式 delay 策略 manifest SHA-256 为 `b202488173c9658164334f1b00044c0eb5f5d5b5827b342a0b9d96fbae886533`，最终配置 SHA-256 为 `7de44a586bbdbd87c507bf8d8f88fd86fd76acf191fd8c848fd1941a7df34bf1`。正式证据位于 `outputs/experiment028_csi_curves/20260803_main/transparent_cdd_small_delay/a100/`：

- `final/estimated_csi_bler_points.csv`、`stability_audit.csv`、`ce_nmse_uncertainty.csv`；
- `final/analytic_vs_monte_carlo_ce.csv`、`bler_bracket_audit.csv`、`analysis_summary.json`；
- `final/a100_subset_with_small_delay_transparent_cdd_<estimated|ce>_points.csv`；
- 逐 trial error flags、CE 三类数组、filter diagnostics、resolved run、307 个区间和后台状态/日志均保存在同一目录；该目录不存在 ideal receiver 或 ideal-CSI CSV。

6-RB transparent PRG 延伸继续使用源 manifest `33ee2b87b7da18b3d2bc5d1faeccd92b070f105d5c4c9f5c65b084ccd9a07eb1`；延伸策略 manifest SHA-256 为 `ccc3f402652948f57ee7e5fed31eb89cd897d16382aded43df6d7cf808a5587e`，最终配置 SHA-256 为 `1648ba6e349b9d765a042867d5e032383629d85a4bb11d7d36bbfe8227b43246`。正式证据位于 `outputs/experiment028_csi_curves/20260803_main/transparent_prg_6rb_1pct_extension/a100/`：

- `final/estimated_csi_bler_points.csv`、`stability_audit.csv`、`ce_nmse_uncertainty.csv`、`bler_bracket_audit.csv` 和 `analysis_summary.json`；
- 最终同图精确 estimated/CE CSV 由绘图脚本写入小时延目录的 `final/a100_subset_with_small_delay_transparent_cdd_<estimated|ce>_points.csv`；
- 逐 trial error flags、CE 三类数组、148 个区间、resolved run 和运行日志均保留；为并行调度建立的分片目录只改变 absolute-trial 区间分配，不改变物理配置或 seed。该目录不存在 ideal receiver 或 ideal-CSI CSV。

透明 PRG、透明 CDD 与小时延透明 CDD 的复现命令为：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_transparent_prg.yaml --stage validate
& powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools\run_plan028_transparent_prg_until_complete.ps1
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_result028_transparent_prg.py
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_transparent_cdd.yaml --stage validate
& powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools\run_plan028_transparent_cdd_until_complete.ps1
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_result028_transparent_cdd.py
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_transparent_cdd_small_delay.yaml --stage validate
& powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools\run_plan028_transparent_cdd_small_delay_until_complete.ps1
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_result028_transparent_cdd_small_delay.py
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_transparent_prg6_1pct_extension.yaml --stage validate
& powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools\run_plan028_transparent_prg6_extension_until_complete.ps1
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_result028_transparent_prg6_extension.py
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_small_delay_cdd_matched.yaml --stage validate
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_small_delay_cdd_matched.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_small_delay_cdd_ideal.yaml --stage validate
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_bler_curves.py --config configs\bler_curves_result028_a100_small_delay_cdd_ideal.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_result028_small_delay_cdd_matched_ideal.py
& D:\venvs\cdd-s102\Scripts\python.exe tools\plot_result028_a100_subset.py
```

### 3.7 A100 透明/非透明 CDD 与 6-RB precoder cycling 对比

本节只读复用第 3.6 节已经完成并审计的正式点，不新增、不重跑 trial，也不插值或外推。图例中的原 `delay set 2` 改称 `large-delay CDD`；“transparent”表示 UE 不知道 CDD delay 或 PRG 预编码索引，只使用底层物理信道协方差形状，“non-transparent”表示 CDD 接收机知道实际 delay，并使用匹配的真实等效协方差。横轴固定为 `13.75–20.25 dB`，即要求的 14–20 dB 范围两侧各保留 0.25 dB；为控制图面范围，两张图中的 6-RB cycling 截止 17 dB，small-delay CDD non-transparent 截止 18.5 dB，其余曲线只画自身已有的正式点。

| 图例方案 | candidate ID | 发射端参数 | 接收端信道估计 |
|---|---|---|---|
| large-delay CDD, transparent | `A100_AP_RMS_T1_TRANSPARENT_CDD` | $j=[0,1.728,3.456,5.184,6.912,8.64,10.368,12.096]$；人工时延 $[0,100,200,300,400,500,600,700]$ ns；等效 FFT samples $[0,12.288,24.576,36.864,49.152,61.44,73.728,86.016]$，扩展为 700 ns / 86.016 samples | UE 不知道 $j$；全带 LMMSE 使用 $8\mathbf R_{\rm phy}$ |
| large-delay CDD, non-transparent | `A100_AP_RMS_T1` | 与上一行完全相同 | UE 知道 $j$；全带 LMMSE 使用匹配的真实等效协方差 $\mathbf R_g$ |
| small-delay CDD, transparent | `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | $j=[0,0.25,0.5,0.75,1,1.25,1.5,1.75]$；人工时延 $[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273]$ ns；等效 FFT samples $[0,1.778,3.556,5.333,7.111,8.889,10.667,12.444]$，扩展为 101.273 ns / 12.444 samples | UE 不知道 $j$；全带 LMMSE 使用 $8\mathbf R_{\rm phy}$ |
| small-delay CDD, non-transparent | `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | 与上一行完全相同 | UE 知道 $j$；全带 LMMSE 使用匹配的真实等效协方差 $\mathbf R_g$ |
| precoder cycling 6 RB, transparent | `A100_PRG_DFT8_6RB` | 48 PRB 分为 8 个 PRG，每 PRG 6 RB = 72 子载波；PRG $m=0,\ldots,7$ 使用第 $m$ 个 8 维 DFT 向量，$w_m[n]=\exp(-j2\pi nm/8)$，$n=0,\ldots,7$；8 个向量固定轮询，无随机尾部，每向量总功率为 8 | UE 不需要知道 DFT 向量索引；每个 PRG 内独立 LMMSE，使用 $8\mathbf R_{\rm phy}$，不跨 PRG 插值 |

CDD 的相位均为 $V_{k,n}=\exp(-j2\pi k j_n/576)$；表中“等效 FFT samples”按 $d_n^{\rm FFT}=j_n(4096/576)=j_n64/9$ 换算，“扩展”是最大与最小人工时延之差，不把数字循环移位作为真实传播时延计入 CP。6-RB cycling 的 DFT 权值在 DMRS 与 data 上保持一致，因此 PRG 内等效 PDP 与底层 PDP 同形，仅相差总功率因子 8。

五条曲线共同具备的 14/16 dB 正式点如下，CE NMSE 均为每 trial data-RE 线性比值先平均再转 dB：

| 方案 | 14 dB BLER | 14 dB CE NMSE | 16 dB BLER | 16 dB CE NMSE |
|---|---:|---:|---:|---:|
| large-delay CDD, transparent | 0.9998 | -5.193 dB | 0.9987 | -5.312 dB |
| large-delay CDD, non-transparent | 0.1270 | -22.677 dB | 0.010667 | -24.450 dB |
| small-delay CDD, transparent | 0.3023 | -19.368 dB | 0.1159 | -20.493 dB |
| small-delay CDD, non-transparent | 0.2353 | -24.213 dB | 0.0786 | -25.996 dB |
| precoder cycling 6 RB, transparent | 0.1807 | -22.469 dB | 0.034661 | -24.265 dB |

事实：在两个共同 SNR 点，large-delay CDD 的透明接收机存在明显协方差失配地板；换为 matched non-transparent 接收机后，16 dB BLER 从 `0.9987` 降为 `0.010667`。small-delay CDD 的透明失配较小，但 matched receiver 仍把 16 dB CE NMSE 从 `-20.493 dB` 改善到 `-25.996 dB`，BLER 从 `0.1159` 降到 `0.0786`。6-RB cycling 的透明估计在 16 dB 为 `-24.265 dB / 0.034661`；这些是当前 A100、1Rx、static TDL-A、comb-6 配置内的 pointwise 事实，不外推缺少正式点的高 SNR 区间。

两张图的精确数据分别位于：

- `outputs/experiment028_csi_curves/20260803_main/section37_receiver_comparison/a100/final/a100_section37_estimated_csi_bler_points.csv`；
- `outputs/experiment028_csi_curves/20260803_main/section37_receiver_comparison/a100/final/a100_section37_ce_nmse_points.csv`。


图文件为 `docs/figures/result-028/a100_comb6_section37_estimated_csi_bler_snr_14_20.png` 与 `docs/figures/result-028/a100_comb6_section37_ce_nmse_snr_14_20.png`；本无图版不嵌入图片。

## 4. 异常、边界与待确认项

- A30/A100 estimated-CSI 的 base 是 result-027 原始正式数据，本次只追加未覆盖 trial；ideal-CSI 同样在 result-028 原数据上追加，均未从 trial 1 重跑后重复相加。
- 透明基线在首个正式区间落盘前把 batch size 从 100 调为 20、可恢复区间冻结为最多 1000 trials；首批 8 个区间后把每进程 SNR task 上限从 1 调为 8。三者只影响计算分块/调度，不改变物理参数、seed 或 absolute-trial 集合。
- 研究者要求暂停时，正在计算的 16 dB `5001–6000` 尚未落盘，因此被丢弃；恢复后从 5001 确定性重算。最终 87 个区间已验证连续、无重叠、两接收机配对一致。
- 透明 CDD 的隔离 batch-size benchmark 证明 batch 100 会在 Sionna 信道生成阶段申请约 11.23 GiB 单张量并 OOM，因此正式运行保持 batch 20。两次外层 1 小时窗口结束只丢弃未落盘的 15.5/16 dB 当前区间；恢复后按 absolute seed 重算，最终 80 个区间连续、无重叠、五条曲线配对一致。
- 正式 target summary 仍采用追加前的预定 fit 口径；新增数据用于降低逐点曲线波动和 pointwise Wilson 区间，不据此后验改变闭合/未闭合规则。
- `AP_TU_NT` rank 4、`AP_TU_NTM1` rank 7 且 condition number 极大；rank 仍只作诊断，不单独替代 CE/BLER 判据。
- A300 未闭合曲线不做 logistic 外推；ideal 曲线也不额外声称目标门限，只报告实际采样点。
- 小时延透明 CDD 的 307 个区间连续无重叠，全部 20 点至少 200 errors；未对曲线做平滑。第 12/15 节当时未运行 ideal CSI，第 16 节已用独立 absolute-trial 数据补充，不与透明 CE trial 混合。
- 第 16 节并行正式运行的若干一小时前台窗口只中断尚未落盘的当前区间；已完成区间通过 base CSV 和相同 absolute seed 继续，matched/ideal 最终 250/172 个区间均连续无重叠。matched 19.5/20 dB 达到 50,000-trial 上限但只有 196/102 errors，已保留原始点与 Wilson 区间，不作平滑或外推。
- 6-RB 延伸的第一次四进程调度在完整新区间落盘前因主机内存不足退出；随后使用至多三个独立 absolute-trial 分片完成。一次一小时前台窗口只中断未落盘区间，恢复时从各分片已完成区间的 `trial_end+1` 继续。最终合并的 148 个区间连续无重叠，无完成 trial 被丢失或重复累计。
- 适用范围仅限当前 48 PRB、static TDL-A、30/100/300 ns、comb-6、零速度与当前 MCS；主结果使用 matched receiver，第 11–12 节透明 CDD 明确使用 $8\mathbf R_{\rm phy}$ 失配接收机。不推广至其他 TDL/CDL、移动性、相关性或其他失配模型。
- 本 result 尚待研究者确认，因此未更新 `KNOWLEDGE.md` 或 `GOALS.md` 的已验证结论/阶段状态，也未创建 Git checkpoint。
