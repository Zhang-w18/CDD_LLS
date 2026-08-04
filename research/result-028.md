# result-028：comb-6 三类 CSI 曲线统一增补与 TDL-A 300 ns 对比

> 对应 `research/plan-028.md`。A30/A100 的 estimated-CSI 数据由 result-027 正式点只读合并；A300 的搜索、outage、estimated-CSI 链路和三个场景的真实 ideal-CSI 链路均已完成。2026-08-04 直接在本轮结果上追加 trial 以降低曲线波动，未新增实验 plan。本结果尚待研究者确认。精确表格与原始证据口径和 `research/result-028-text.md` 完全一致。

## 1. 结论与状态

事实：A30、A100、A300 分别有 8、10、10 条物理曲线和 72、85、91 个完整采样点。增量合计追加 412,600 candidate-trials；estimated 与 CE 共用追加 trial，ideal 按 1%附近独立停止。合并后三场景的共同采样点均满足 ideal BLER 不高于 estimated，原 A30 的一错误块反向不再出现。

事实：A300 的 10% estimated 目标只闭合 `AP_RMS_T1`、`AP_TALIAS_NT`、`B0_QC`、`S0_SIDON`，1% 只闭合 `AP_RMS_T1`、`B0_QC`、`S0_SIDON`。六基线最优均为 `AP_RMS_T1=14.152/15.365 dB`；没有搜索 family 闭合并优于该基线。大跨距候选的 CE 地板与不闭合相关，但只作为当前配置内的诊断推断。

## 2. 仿真条件与证据口径

系统固定为 48 PRB、$K=576$、30 kHz、8 Tx / 1 Rx、static TDL-A、comb-6、192 pilot/5568 data RE、256QAM MCS 8、LDPC 8 iterations。CDD 为 $V_{k,n}=\exp(-j2\pi k j_n/576)$，`1q=57.870370370... ns`。estimated 接收机为 V-aware matched LMMSE；ideal 接收机仍运行完整 TDL-A、编码、256QAM、AWGN、真实 data-RE CSI 均衡、软解调与译码。

A30/A100 源 manifest SHA-256 分别为 `b20cca1e...97748a` 和 `33ee2b87...07eb1`；A300 manifest 为 `6ebda320...3a665`。增量配置为 `configs/bler_curves_result028_augmentation.yaml`：目标约 50 个错误块、每点最多 5,000 trials，estimated/ideal 分别在首个低于 0.3%/0.5% 的点停止。完整参数、复现命令和证据路径见无图版第 2、3.4 节。

## 3. comb-6 三场景统一增补

### 3.1 A30 comb-6

| 编号 | 原缩写 | delay ns |
|---|---|---|
| candidate 01 | `A30_B0_QC` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` |
| candidate 02 | `A30_AP_RMS_T1` | `[0,30,60,90,120,150,180,210]` |
| candidate 03 | `A30_AP_TEPS_T1` | `[0,143.898,287.796,431.694,575.592,719.49,863.388,1007.286]` |
| candidate 04 | `A30_AP_TU_NT` | `[0,4166.667,8333.333,12500,16666.667,20833.333,25000,29166.667]` |
| candidate 05 | `A30_AP_TALIAS_NT` | `[0,694.444,1388.889,2083.333,2777.778,3472.222,4166.667,4861.111]` |
| candidate 06 | `A30_S0_SIDON` | `[0,57.87,173.611,405.093,694.444,1157.407,1736.111,3761.574]` |
| candidate 07 | `A30_AP_T2_00` | `[0,4108.796,8217.593,12326.389,16435.185,20543.981,24652.778,28761.574]` |
| candidate 08 | `A30_MEFF_T2_05` | `[0,2835.648,7523.148,11400.463,17534.722,22569.444,26736.111,30324.074]` |

精确 q/ns、rank 和 condition number 见 `outputs/experiment028_csi_curves/20260803_main/a30_comb6/final/candidate_number_delay_table.csv`。

![A30 comb-6 candidate estimated CSI BLER](../docs/figures/result-028/a30_comb6_candidate_estimated_csi_bler.png)

![A30 comb-6 candidate ideal CSI BLER](../docs/figures/result-028/a30_comb6_candidate_ideal_csi_bler.png)

![A30 comb-6 candidate CE NMSE](../docs/figures/result-028/a30_comb6_candidate_ce_nmse.png)

A30 estimated/ideal 分别追加 67,400/91,900 candidate-trials；ideal 图使用 48/72 个点。estimated 正式闭合目标仍沿用原预定 fit。增量 CSV 为 `curve_augmentation_v2/a30/final/<receiver>_csi_bler_points.csv`。

### 3.2 A100 comb-6

| 编号 | 原缩写 | delay ns |
|---|---|---|
| candidate 01 | `A100_B0_QC` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` |
| candidate 02 | `A100_AP_RMS_T1` | `[0,100,200,300,400,500,600,700]` |
| candidate 03 | `A100_AP_TEPS_T1` | `[0,479.66,959.32,1438.98,1918.64,2398.3,2877.96,3357.62]` |
| candidate 04 | `A100_AP_TU_NT` | `[0,4166.667,8333.333,12500,16666.667,20833.333,25000,29166.667]` |
| candidate 05 | `A100_AP_TU_NTM1` | `[0,4761.905,9523.81,14285.714,19047.619,23809.524,28571.429,33333.333]` |
| candidate 06 | `A100_AP_TALIAS_NT` | `[0,694.444,1388.889,2083.333,2777.778,3472.222,4166.667,4861.111]` |
| candidate 07 | `A100_S0_SIDON` | `[0,57.87,173.611,405.093,694.444,1157.407,1736.111,3761.574]` |
| candidate 08 | `A100_AP_T2_01` | `[0,3761.574,7986.111,12210.648,16435.185,20659.722,24884.259,29108.796]` |
| candidate 09 | `A100_MEFF_T2_04` | `[0,2314.815,6655.093,10185.185,13020.833,18750,22858.796,28009.259]` |
| candidate 10 | `A100_MEFF_T2_06` | `[0,1851.852,9085.648,14872.685,17245.37,19328.704,23495.37,29976.852]` |

精确 q/ns、rank 和 condition number 见 `a100_comb6/final/candidate_number_delay_table.csv`。

![A100 comb-6 candidate estimated CSI BLER](../docs/figures/result-028/a100_comb6_candidate_estimated_csi_bler.png)

![A100 comb-6 candidate ideal CSI BLER](../docs/figures/result-028/a100_comb6_candidate_ideal_csi_bler.png)

![A100 comb-6 candidate CE NMSE](../docs/figures/result-028/a100_comb6_candidate_ce_nmse.png)

A100 estimated/ideal 分别追加 37,600/84,800 candidate-trials；ideal 图使用 40/85 个点。正式闭合与目标表仍沿用原预定 fit。

### 3.3 A300 comb-6

A300 的 99% PDP 支撑为 1438.98 ns。动态重算后的候选如下：

| 编号 | 原缩写 | delay ns |
|---|---|---|
| candidate 01 | `A300_B0_QC` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` |
| candidate 02 | `A300_AP_RMS_T1` | `[0,300,600,900,1200,1500,1800,2100]` |
| candidate 03 | `A300_AP_TEPS_T1` | `[0,1438.98,2877.96,4316.94,5755.92,7194.9,8633.88,10072.86]` |
| candidate 04 | `A300_AP_TU_NT` | `[0,4166.667,8333.333,12500,16666.667,20833.333,25000,29166.667]` |
| candidate 05 | `A300_AP_TU_NTM1` | `[0,4761.905,9523.81,14285.714,19047.619,23809.524,28571.429,33333.333]` |
| candidate 06 | `A300_AP_TALIAS_NT` | `[0,694.444,1388.889,2083.333,2777.778,3472.222,4166.667,4861.111]` |
| candidate 07 | `A300_S0_SIDON` | `[0,57.87,173.611,405.093,694.444,1157.407,1736.111,3761.574]` |
| candidate 08 | `A300_AP_T2_01` | `[0,3761.574,7986.111,12210.648,16435.185,20659.722,24884.259,29108.796]` |
| candidate 09 | `A300_MEFF_T2_05` | `[0,1851.852,8449.074,12442.13,15682.87,23321.759,25925.926,29687.5]` |
| candidate 10 | `A300_MEFF_T2_06` | `[0,1967.593,8622.685,11226.852,13541.667,17534.722,21180.556,28067.13]` |

精确 q/ns、rank 和 condition number 见 `a300_comb6/final/candidate_number_delay_table.csv`。

缩写图例版：

![A300 comb-6 abbreviated estimated CSI BLER](../docs/figures/result-028/a300_comb6_abbrev_estimated_csi_bler.png)

![A300 comb-6 abbreviated ideal CSI BLER](../docs/figures/result-028/a300_comb6_abbrev_ideal_csi_bler.png)

![A300 comb-6 abbreviated CE NMSE](../docs/figures/result-028/a300_comb6_abbrev_ce_nmse.png)

candidate 编号版（样式与缩写版逐物理方案完全一致）：

![A300 comb-6 candidate estimated CSI BLER](../docs/figures/result-028/a300_comb6_candidate_estimated_csi_bler.png)

![A300 comb-6 candidate ideal CSI BLER](../docs/figures/result-028/a300_comb6_candidate_ideal_csi_bler.png)

![A300 comb-6 candidate CE NMSE](../docs/figures/result-028/a300_comb6_candidate_ce_nmse.png)

A100 同口径的链路增补图：

![A300 comb-6 all-candidate prescan](../docs/figures/result-028/a300_comb6_prescan_all_curves.png)

![A300 comb-6 formal estimated CSI BLER](../docs/figures/result-028/a300_comb6_bler_curves.png)

![A300 comb-6 10pct outage CE BLER link](../docs/figures/result-028/a300_comb6_10pct_outage_ce_link.png)

![A300 comb-6 1pct outage CE BLER link](../docs/figures/result-028/a300_comb6_1pct_outage_ce_link.png)

正式目标：`AP_RMS_T1` 为 14.152 `[14.109,14.195]` / 15.365 `[15.271,15.459]` dB；`B0_QC` 为 14.402 `[14.358,14.445]` / 15.575 `[15.486,15.664]`；`S0_SIDON` 为 14.260 `[14.221,14.299]` / 15.521 `[15.407,15.636]`；`AP_TALIAS_NT` 仅 10% 闭合，为 15.845 `[15.782,15.908]`。完整 own CE、错误数与 trials 见无图版表及 `a300_comb6/comb6/link/A300/final/target_summary.csv`。

A300 estimated/ideal 分别追加 22,800/108,100 candidate-trials；ideal 图使用 33/91 个完整数据点。上述正式目标保持原预定 fit，不因本次视觉降波动增量而后验改口径。

### 3.4 曲线降波动增量

增量入口按绝对 trial 编号续跑，完成区间立即保存；六类 receiver/scenario 的剩余任务数均为 0。1%附近的 Wilson 95%区间中位宽度均下降，其中 A300 ideal 从 `0.01926` 降至 `0.00505`。ideal 图不再绘制远低于 1% 的高 SNR 点；A30/A100/A300 分别保留 48、40、33 个点。正式 target summary 未后验重拟合。

### 3.5 作图与验收

每个场景使用一个 `curve_styles.json`，同一物理方案跨 CSI 类型、缩写/编号图例和散点保持相同 color/linestyle/marker。标题、轴标签、图例至少 16 pt，刻度至少 14 pt，线宽至少 2.5 pt，marker 至少 8 pt；16 张图均按约 13 cm PPT 半页宽预览通过。零错误点的图示下界为 `0.5/trials`，CSV 保留真实值。

## 4. 异常、边界与待确认项

- A30/A100 estimated 以 result-027 为 base；本次 estimated/ideal 都只追加未覆盖的绝对 trial 区间。
- 未闭合目标不外推；ideal 曲线不另做门限外推。
- 结论只适用于当前 48 PRB、static TDL-A、comb-6、零速度、当前 MCS 与 matched receiver。
- 本 result 尚待研究者确认，因此未更新 `KNOWLEDGE.md`/`GOALS.md`，也未创建 Git checkpoint。
