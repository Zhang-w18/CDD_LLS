# result-035-text：PDSCH 2Rx、60 km/h 下 4/8Tx 的 estimated/ideal CSI 与 CE NMSE

> 状态：正式实验已执行，结果待研究者确认。本文件为无图版；结论、配置和数字与 `research/result-035-PDSCH 2Rx仿真.md` 一致。

## 1. 配置、执行和审计

本轮按 `research/plan-035-PDSCH 2Rx仿真.md` 运行 `A100_NT4_NR2_V60` 与
`A100_NT8_NR2_V60`。系统采用 TDL-A、100 ns RMS delay spread、3.5 GHz、60 km/h、
48 PRB、30 kHz SCS、10 个 PDSCH symbols、DMRS symbols `[2,7]` comb-6、16QAM、
NR 256QAM table MCS 8、单层、2Rx MRC。estimated CSI 使用真实时频协方差的二维 LMMSE；
ideal CSI 仅替换接收端信道知识，aged MRT 的发射权值仍来自 4.994791667 ms 旧 CSI。
横轴为每根 Rx 分支 SNR，发射端保持单位总功率。

正式网格为 8--14 dB、间隔 0.25 dB。每点先运行 1,000 个共同 paired trials，再按冻结的
10%/1% bracket 端点规则以 1,000 trials 为单位追加，目标 200 errors、单点上限 50,000。
4Tx 共运行 321,000 paired trials，8Tx 共运行 333,000 paired trials；两场景均为
`complete`，无 capped endpoint。base seed 为 `20260727`，CPU-only。

- 4Tx 配置 SHA-256：`95e30cee2eb1079478a7ad9526505dc686a72d0305b852097ba96ca1c0580b42`。
- 8Tx 配置 SHA-256：`c6e4d5fd8a9342deac504eb30e081b832c6de6c12bad030a572adb19db60c0d4`。
- 4Tx `intervals.csv` SHA-256：`6b4f909b0648b60b5507b827b6b11e44f9097155a08a69e796563dc53411a07a`。
- 8Tx `intervals.csv` SHA-256：`dcbb4dd9df8158a64cb414264baa463aad6f260e33d44f31388d1be384d6c3c8`。
- 4Tx 审计：3,852 interval rows、5,778 个数组、300 个曲线点全部一致，峰值 RSS 3,229,732,864 bytes。
- 8Tx 审计：3,996 interval rows、5,994 个数组、300 个曲线点全部一致，峰值 RSS 3,423,227,904 bytes。

审计逐 interval 核对 absolute-trial 连续性、error flag 求和、NMSE 数组有限性及求和、
estimated/ideal 配对数量、NMSE 仅属于 estimated、CPU placement 和自适应停止状态。完整审计分别见
`outputs/experiment035_pdsch_2rx/20260917_fix1/formal/<scenario>/analysis/data_audit.json`。

首次 `20260917_main` prescan 因 estimated error-flag 保存键错误而作废；本 result 仅使用修复后
`20260917_fix1` 的 smoke、prescan 和 formal 数据。

## 2. 10%/1% BLER crossing

表中数值为 `crossing SNR dB [500 次 paired bootstrap 95% interval]`。crossing 仅在相邻真实点
双侧 bracket 内按 log-BLER 线性插值，不做外推或单调修正。

| 场景 | 方案 | estimated 10% | estimated 1% | ideal 10% | ideal 1% |
|---|---|---:|---:|---:|---:|
| 4Tx | B0_QC | 9.683 [9.650,9.711] | 10.985 [10.940,11.037] | 8.644 [8.585,8.693] | 9.917 [9.850,9.990] |
| 4Tx | S0_SIDON | 9.457 [9.391,9.535] | 11.068 [11.028,11.101] | 8.892 [8.837,8.958] | 10.433 [10.366,10.513] |
| 4Tx | transparent PRG6 | 9.665 [9.622,9.701] | 11.313 [11.239,11.368] | 9.003 [8.915,9.039] | 10.475 [10.415,10.531] |
| 4Tx | aged MRT PRG6 | 10.330 [10.275,10.368] | 12.597 [12.521,12.651] | 9.621 [9.540,9.663] | 11.683 [11.620,11.773] |
| 4Tx | small-delay transparent | 10.599 [10.568,10.634] | 12.830 [12.750,12.894] | 9.986 [9.949,10.021] | 12.189 [12.119,12.267] |
| 4Tx | small-delay matched | 10.450 [10.406,10.504] | 12.633 [12.523,12.745] | 9.986 [9.949,10.021] | 12.189 [12.119,12.267] |
| 8Tx | B0_QC | 9.701 [9.646,9.760] | 10.681 [10.646,10.721] | 8.215 [8.171,8.262] | 9.108 [9.078,9.137] |
| 8Tx | S0_SIDON | 9.364 [9.340,9.392] | 10.429 [10.394,10.471] | 8.241 [8.191,8.296] | 9.204 [9.177,9.239] |
| 8Tx | transparent PRG6 | 9.435 [9.380,9.516] | 10.868 [10.817,10.915] | 8.773 [8.704,8.821] | 10.194 [10.082,10.274] |
| 8Tx | aged MRT PRG6 | 9.885 [9.798,9.925] | 12.200 [12.104,12.285] | 9.165 [9.134,9.196] | 11.201 [11.147,11.280] |
| 8Tx | small-delay transparent | 10.451 [10.425,10.483] | 12.626 [12.563,12.680] | 9.454 [9.399,9.537] | 11.426 [11.371,11.500] |
| 8Tx | small-delay matched | 9.980 [9.918,10.024] | 11.921 [11.870,11.977] | 9.454 [9.399,9.537] | 11.426 [11.371,11.500] |

每个 crossing 的真实 bracket、两端 BLER、errors 和 trials 见
`analysis/target_crossings_bootstrap.csv`。例如 8Tx estimated S0_SIDON 的 1% bracket 为
10.25--10.5 dB，两端分别为 357/22,000 与 206/25,000 errors/trials；8Tx estimated B0_QC
为 10.5--10.75 dB，两端 403/25,000 与 200/24,000；4Tx estimated B0_QC 为
10.75--11 dB，两端 207/12,000 与 203/21,000。所有触发追加规则的端点均达到 200 errors。

## 3. 方案差值与接收机间隔

差值定义为“前者 crossing SNR 减后者”，负值表示前者所需 SNR 更低。

| 场景/接收机 | 比较 | 10%差值 dB [95%] | 1%差值 dB [95%] |
|---|---|---:|---:|
| 4Tx estimated | S0_SIDON − B0_QC | -0.226 [-0.291,-0.151] | +0.083 [+0.024,+0.127] |
| 8Tx estimated | S0_SIDON − B0_QC | -0.337 [-0.403,-0.282] | -0.252 [-0.302,-0.203] |
| 4Tx estimated | small matched − transparent | -0.149 [-0.190,-0.105] | -0.197 [-0.307,-0.112] |
| 8Tx estimated | small matched − transparent | -0.471 [-0.546,-0.423] | -0.705 [-0.773,-0.624] |
| 4Tx estimated | aged MRT − transparent PRG6 | +0.665 [+0.604,+0.722] | +1.284 [+1.200,+1.380] |
| 8Tx estimated | aged MRT − transparent PRG6 | +0.451 [+0.336,+0.517] | +1.333 [+1.219,+1.436] |

estimated 相对 ideal 的 penalty 在 4Tx 各方案约为 0.44--1.07 dB，在 8Tx 约为
0.50--1.57 dB。8Tx 的 B0_QC penalty 最大：10% 为 1.486 dB [1.419,1.571]，1% 为
1.573 dB [1.527,1.626]。完整逐方案表见 `analysis/target_snr_differences_bootstrap.csv`。

## 4. CE NMSE

NMSE 为每 trial 在两根 Rx 和全部 data RE 上先求线性比值，再跨 trial 线性平均并转 dB。
在 12 dB 时：

| 方案 | 4Tx NMSE dB | 8Tx NMSE dB |
|---|---:|---:|
| B0_QC | -17.574 | -15.835 |
| S0_SIDON | -20.345 | -17.024 |
| transparent PRG6 | -19.445 | -19.554 |
| aged MRT PRG6 | -19.334 | -19.693 |
| small-delay transparent | -19.780 | -17.294 |
| small-delay matched | -21.452 | -21.325 |

事实：matched small-delay 在两个场景均给出最低 NMSE；8Tx B0_QC、SIDON 与 transparent
small-delay 的 NMSE 明显差于各自 4Tx 值。推断：8Tx estimated 排序比 ideal 排序更受
频域协方差匹配和可估计性影响。NMSE 本身不作为方案胜出判据。

## 5. 结论与边界

1. 8Tx estimated 下 S0_SIDON 在 10%和1%均优于 B0_QC；4Tx 只在10%优于 B0_QC，
   到1%反而多需约0.083 dB。因此“SIDON 总是优于 B0”不能跨 4/8Tx 泛化。
2. small-delay matched 与 transparent 在 ideal CSI 下完全相同，符合其发射波形相同的定义；
   estimated 下 matched 明确更优，且8Tx增益大于4Tx，说明差异来自接收端协方差知识。
3. 60 km/h、约5 ms aged MRT 在 ideal 和 estimated 下都差于 transparent PRG6；本结果不支持
   旧 CSI MRT 带来性能增益。
4. ideal/estimated 间隔对所有方案均为正，说明 CE 损失不可忽略；8Tx B0_QC 尤其明显。
5. 4Tx 与8Tx不是 paired 场景；跨 Tx 数的差值只能结合各自独立区间解释。
6. plan-033 的对应4Rx正式结果尚未形成可直接引用的完整同口径汇总，因此本 result 不给出
   2Rx 对4Rx的定量差值。该问题仍待 plan-033 完成后分析。
7. 自适应采样在固定1,000-trial检查点停止；Wilson单点区间和500次 bootstrap 是最终数据的
   不确定性描述。BLER 图纵轴截取为 0.008--1.1，横轴主/次网格间隔分别为 0.5/0.25 dB；
   CSV 保留全部真实数值，零误块点不被改写。

## 6. 产物与复现

每个场景的 `analysis/` 包含三张 PNG、对应绘图 CSV、`style_map.json`、`figure_audit.json`、
`data_audit.json`、`target_crossings_bootstrap.csv` 和 `target_snr_differences_bootstrap.csv`。
PNG 均为 2160×1560，约 18.3×13.2 cm@300 dpi，满足约13 cm宽预览要求；视觉布局仍需研究者人工核对。

复现分析：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_result035_pdsch_2rx.py --scene outputs\experiment035_pdsch_2rx\20260917_fix1\formal\A100_NT4_NR2_V60
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_result035_pdsch_2rx.py --scene outputs\experiment035_pdsch_2rx\20260917_fix1\formal\A100_NT8_NR2_V60
```

结果尚未获研究者确认，因此未更新 `KNOWLEDGE.md`/`GOALS.md`，也未创建 Git checkpoint。
