# result-033-text：4Rx 下 4/8Tx、3/60 km/h 对比（部分正式结果，无图版）

对应计划：`research/plan-033-PDSCH-4Rx-MIMO移动性.md`。状态：**正式实验进行中**；本版包含已完成初始 10,000 trials/点的 `A100_NT8_NR4_V3`、`A100_NT8_NR4_V60` 与 `A100_NT4_NR4_V60`，以及第 4 节 8Tx/4Rx ideal-CSI Sidon/PRG6 轻量诊断。正式场景尚未完成 1%端点追加，诊断图也不进入正式门限验收，因此本版不构成 plan-033 最终验收。`A100_NT4_NR4_V3`、完整跨 Tx 比较、paired bootstrap 和最终结论均待补。

本版 estimated-CSI 表格和图文版中的图冻结于三个场景各点初始 10,000 trials 的 `analysis_partial/` 快照。快照生成后，`A100_NT8_NR4_V60` 的 5 个待追加 SNR 已从 absolute trial 10,001 连续推进到 11,000；8 个冻结端点仍少于 200 errors，因此这些新增 trial 尚未并入本版 crossing，待该场景端点闭合后统一重算。

## 1. 当前可报告的事实

- 8Tx/4Rx、3 km/h 下，5 ms aged MRT 的 10%/1%初始目标 SNR 为 `-1.336/-0.544 dB`，相对 transparent PRG6 左移 `6.256/6.360 dB`。同一 8Tx/4Rx 链路在 60 km/h 下，aged MRT 的初始目标变为 `5.354/6.947 dB`，相对 transparent PRG6 高 `0.297/0.966 dB`。这是两个场景的初始点估计；速度差异方向与旧 CSI 去相关一致，但端点追加和跨场景独立不确定性分析前不作最终效应量结论。
- 同一 8Tx/3 km/h 场景的非 aged 方案中，transparent PRG6 的 10%目标最好；Sidon CDD 在 1%目标仅比它高 `0.021 dB`。Sidon 相对 B0-QC 在 10%/1%分别节省 `0.292/0.284 dB`。
- 8Tx/60 km/h 下，transparent PRG6 的 10%/1%初始目标为 `5.056/5.981 dB`，在六条曲线中最低；Sidon CDD 为 `5.296/6.016 dB`，分别高 `0.240/0.035 dB`。1%差异很小，不能在追加和 paired bootstrap 前判定稳定优劣。
- 4Tx/4Rx、60 km/h 下，Sidon CDD 的 10%/1%初始目标 SNR 为 `5.152/6.153 dB`，在六条曲线中最低；相对 transparent PRG6 分别节省 `0.075/0.144 dB`。效应量较小，最终是否稳定需等待端点追加和 paired bootstrap。
- 8Tx/4Rx、3 km/h ideal-CSI 诊断中，Sidon 在 `3.5--5.75 dB` 的 10 个共同采样点上均有较低的 BLER，`6 dB` 两者均为 `1/1000`。按原始采样点，Sidon 在 `4.25 dB` 已低于 10% ，PRG6 到 `4.5 dB` 才低于 10%；低于 1% 分别首次出现在 `5.0 dB` 和 `5.25 dB`。这只是每点 1000 paired trials 的趋势诊断，未做 crossing 拟合、追加或 paired 不确定性分析，不视为已确认的 SNR 增益。
- small-delay matched Rx 在三个已完成场景均优于 transparent Rx：8Tx/3 km/h 的 10%/1%差为 `-0.376/-0.368 dB`，8Tx/60 km/h 为 `-0.424/-0.360 dB`，4Tx/60 km/h 为 `-0.126/-0.124 dB`。负值表示 matched Rx 所需 SNR 更低。
- 三个场景的目标绘图区均无相邻 SNR 的 BLER 反向，但 1% bracket 中共有 30 个 candidate-endpoint record 少于 200 errors，因此以下 crossing 仍是初始点估计，不能作为最终置信结论。

## 2. 冻结条件与当前样本

| 项目 | 取值 |
|---|---|
| 信道 | Sionna TDL-A，RMS delay spread 100 ns，3.5 GHz，20 sinusoids；天线分支独立、无空间相关 |
| 天线/层数 | `8Tx/4Rx/1 layer` 或 `4Tx/4Rx/1 layer` |
| 波形 | 48 PRB，576 active subcarriers，30 kHz SCS，FFT 4096，CP 288，10 个 PDSCH symbols |
| DMRS/接收机 | symbols `[2,7]`、comb-6；二维时频 LMMSE；各 Rx 独立估计，data RE 上 MRC |
| 调制编码 | 16QAM；NR 256QAM table MCS 8；目标码率 553/1024；Sionna LDPC 最多 8 次迭代 |
| 发射归一化 | 每频点预编码向量满足 $\lVert\mathbf W_k\rVert_2^2=1$；噪声方差 $1/\mathrm{SNR}_{linear}$ |
| CSI age | 10 slots，实际 4.994791667 ms；old/current 来自同一连续 realization |
| 随机性 | seed `20260727`；同一场景/SNR/absolute trial 的六曲线共享信道、payload 和 data/DMRS AWGN |
| 当前样本 | 8Tx/V3：33 点；8Tx/V60：17 点；4Tx/V60：17 点；每点六曲线共同 10,000 trials，每 interval 1,000 trials |

## 3. 初始目标 SNR 与原始 bracket

目标只在相邻真实采样点双侧 bracket 内按 log10(BLER) 线性插值，不做平滑、单调修正或外推。

表中 1%端点均报告原始错误数；其 Wilson 95%区间完整保存在 `target_summary_initial.csv`。这些端点区间总体覆盖约 `0.0020--0.0243`，例如 8Tx/V60 aged MRT 的上下端分别为 `[0.0115,0.0161]` 与 `[0.0075,0.0113]`，说明追加前的 1% crossing 仍有不可忽略的抽样不确定性。

| 场景 | 曲线 | 10% SNR (dB) | 1% SNR (dB) | 1% bracket：低 SNR 点 → 高 SNR 点 |
|---|---|---:|---:|---|
| 8Tx/V3 | aged MRT PRG6 | -1.336 | -0.544 | -0.75 dB: 0.0191 (191/10000) → -0.50 dB: 0.0087 (87/10000) |
| 8Tx/V3 | B0-QC | 5.430 | 6.120 | 6.00: 0.0145 (145) → 6.25: 0.0067 (67) |
| 8Tx/V3 | Sidon | 5.139 | 5.836 | 5.75: 0.0156 (156) → 6.00: 0.0043 (43) |
| 8Tx/V3 | small CDD matched | 5.278 | 6.565 | 6.50: 0.0115 (115) → 6.75: 0.0067 (67) |
| 8Tx/V3 | small CDD transparent | 5.655 | 6.932 | 6.75: 0.0141 (141) → 7.00: 0.0088 (88) |
| 8Tx/V3 | transparent PRG6 | 4.921 | 5.815 | 5.75: 0.0126 (126) → 6.00: 0.0052 (52) |
| 8Tx/V60 | aged MRT PRG6 | 5.354 | 6.947 | 6.75: 0.0136 (136) → 7.00: 0.0092 (92) |
| 8Tx/V60 | B0-QC | 5.539 | 6.274 | 6.25: 0.0111 (111) → 6.50: 0.0037 (37) |
| 8Tx/V60 | Sidon | 5.296 | 6.016 | 6.00: 0.0109 (109) → 6.25: 0.0029 (29) |
| 8Tx/V60 | small CDD matched | 5.405 | 6.792 | 6.75: 0.0116 (116) → 7.00: 0.0048 (48) |
| 8Tx/V60 | small CDD transparent | 5.829 | 7.152 | 7.00: 0.0153 (153) → 7.25: 0.0076 (76) |
| 8Tx/V60 | transparent PRG6 | 5.056 | 5.981 | 5.75: 0.0213 (213) → 6.00: 0.0094 (94) |
| 4Tx/V60 | aged MRT PRG6 | 5.679 | 7.250 | 7.00: 0.0174 (174) → 7.25: 0.0100 (100) |
| 4Tx/V60 | B0-QC | 5.404 | 6.338 | 6.25: 0.0138 (138) → 6.50: 0.0055 (55) |
| 4Tx/V60 | Sidon | 5.152 | 6.153 | 6.00: 0.0154 (154) → 6.25: 0.0076 (76) |
| 4Tx/V60 | small CDD matched | 5.751 | 7.422 | 7.25: 0.0120 (120) → 7.50: 0.0092 (92) |
| 4Tx/V60 | small CDD transparent | 5.878 | 7.546 | 7.50: 0.0119 (119) → 7.75: 0.0046 (46) |
| 4Tx/V60 | transparent PRG6 | 5.228 | 6.297 | 6.25: 0.0109 (109) → 6.50: 0.0069 (69) |

## 4. 8Tx/4Rx ideal-CSI Sidon 与 PRG6 轻量诊断

该诊断使用 `A100_NT8_NR4_V3`、TDL-A 100 ns、8Tx/4Rx/1 layer、单位总发射功率和每根 Rx 分支 SNR 口径。Sidon delay grid 为 `[0,1,3,7,12,20,30,65]`；PRG6 按 8 个 6-RB PRG 依次使用 DFT8 向量 `[0,1,2,3,4,5,6,7]`。接收机在每个 data RE 使用真实等效信道执行 4Rx MRC，不包含信道估计误差。两条曲线在每个 `SNR + absolute trial` 共享底层信道、payload 和 data noise；seed 为 `20260727`，batch size 为 `25`，每点 1000 trials。

| SNR (dB) | Sidon errors/trials | Sidon BLER [Wilson 95%] | PRG6 errors/trials | PRG6 BLER [Wilson 95%] |
|---:|---:|---:|---:|---:|
| 3.50 | 394/1000 | 0.394 [0.3642, 0.4246] | 428/1000 | 0.428 [0.3977, 0.4589] |
| 3.75 | 236/1000 | 0.236 [0.2107, 0.2633] | 303/1000 | 0.303 [0.2753, 0.3322] |
| 4.00 | 153/1000 | 0.153 [0.1320, 0.1766] | 212/1000 | 0.212 [0.1878, 0.2384] |
| 4.25 | 90/1000 | 0.090 [0.0738, 0.1093] | 127/1000 | 0.127 [0.1078, 0.1491] |
| 4.50 | 50/1000 | 0.050 [0.0381, 0.0653] | 88/1000 | 0.088 [0.0720, 0.1072] |
| 4.75 | 13/1000 | 0.013 [0.0076, 0.0221] | 48/1000 | 0.048 [0.0364, 0.0631] |
| 5.00 | 9/1000 | 0.009 [0.0047, 0.0170] | 20/1000 | 0.020 [0.0130, 0.0307] |
| 5.25 | 1/1000 | 0.001 [0.0002, 0.0056] | 8/1000 | 0.008 [0.0041, 0.0157] |
| 5.50 | 0/1000 | 0.000 [0.0000, 0.0038] | 5/1000 | 0.005 [0.0021, 0.0117] |
| 5.75 | 0/1000 | 0.000 [0.0000, 0.0038] | 3/1000 | 0.003 [0.0010, 0.0088] |
| 6.00 | 1/1000 | 0.001 [0.0002, 0.0056] | 1/1000 | 0.001 [0.0002, 0.0056] |

事实：Sidon 的点估计在 `3.5--5.75 dB` 全部低于 PRG6，其中 `3.75/4.0/4.5/4.75 dB` 两条曲线的单点 Wilson 95% 区间不重叠。但是 Wilson 区间是各曲线单点区间，不是 paired 差值区间，不可据此报告精确配对增益。`5.5/5.75/6.0 dB` 的 Sidon 误块数为 `0/0/1`，属于低误块数抽样反向，图上的零值仅用 `0.5/1000` 下界显示，CSV 保留真实零值。

推断：在本轮独立天线、TDL-A 100 ns 和 ideal-CSI 设定下，Sidon 的频率/空间分集趋势优于固定 PRG6 cycling，而 estimated-CSI 主实验中两者接近的排序还包含信道估计和接收机模型影响。该机理解释仍待更大 ideal-CSI 样本或 paired 差值分析验证；本节不拟合 10%/1% crossing，不将其写成全局结论。

证据与复现：

- 原始汇总：`outputs/experiment033_tdl_mobility_mimo/20260916_ideal_sidon_prg6/A100_NT8_NR4_V3/ideal_csi_bler_points.csv`，SHA-256 `656E97C9DE2E713FC654490082B71DCB84DD502839A5919DF1FE0F3586BCE85C`；
- 展开配置：`.../resolved_config.json`，SHA-256 `BBEA08F572A07CF06CD2B228F88FDB1719BA3AE6A4CE00EAEEA26BAC804D220B`；
- 同目录的 `intervals.csv` 共 22 条记录，对应 2 条曲线 × 11 个 SNR，合计 22,000 curve-trials；22 个 error-flag 数组引用全部存在；
- 绘图数据为同目录 `plot_data.csv`，样式为 `curve_styles.json`，图文版引用的常规图和 13 cm 预览分别为 `docs/figures/result-033/a100_nt8_nr4_v3_ideal_sidon_vs_prg6.png` 与 `..._preview_13cm.png`；
- 脚本：`tools/plot_plan033_ideal_sidon_prg6.py`，SHA-256 `47C2ACB600E87FF23C53B9CCB589190EE442B8E6E614533463C198D1FE40D83B`。

```powershell
python tools/plot_plan033_ideal_sidon_prg6.py
```

## 5. 图数据、逻辑检查与证据

图文版当前包含三个场景各一张 BLER 图和一张 data-RE CE-NMSE 图。BLER 图只显示 `BLER >= 0.005` 的目标区域；范围外点仍保留在源 CSV。所有 marker/线段均来自原始 Monte Carlo 点，同一方案在全部图中复用 `curve_styles.json`；另生成约 13 cm 宽预览供研究者人工检查，Agent 未加载图片。

机器可核验证据：

- `outputs/experiment033_tdl_mobility_mimo/20260915_main/formal/{A100_NT8_NR4_V3,A100_NT8_NR4_V60,A100_NT4_NR4_V60}/final/estimated_csi_bler_points.csv`：BLER、Wilson 95%区间、错误数、trials、CE NMSE；
- `outputs/experiment033_tdl_mobility_mimo/20260915_main/analysis_partial/available_formal_points.csv`：三场景当前全部 402 个曲线点；
- `.../target_summary_initial.csv`：36 个目标记录及实测 bracket；
- `.../append_requirements_initial.csv`：30 个需继续检查/追加的 candidate-endpoint record；
- `.../partial_analysis.json`：完整/缺失场景与图文件审计；
- `.../curve_styles.json`：共享颜色、线型和 marker；
- 绘图脚本：`tools/analyze_result033_tdl_mobility_mimo.py`，SHA-256 `F8AE5AC3FBC10870DE83D37F3AAE341BB4ABDB2AE699E5152DA4A21695029A88`。

复现命令：

```powershell
python tools/analyze_result033_tdl_mobility_mimo.py --formal-root outputs/experiment033_tdl_mobility_mimo/20260915_main/formal --analysis-dir outputs/experiment033_tdl_mobility_mimo/20260915_main/analysis_partial --figure-dir docs/figures/result-033
```

自动检查确认 8Tx/V3 有 `1980` 个 1000-trial interval，8Tx/V60 与 4Tx/V60 各有 `1020` 个，所有数组引用存在且 absolute-trial 区间连续、无重叠；三个场景在 `BLER >= 0.005` 的相邻正式点中均无反向。当前六张主图均成功生成，另有六张 13 cm 预览；`partial_analysis.json` 记录文件路径和大小，PNG 格式与像素尺寸另经本地元数据检查通过。

## 6. 未完成项与解释边界

- 必须完成 `A100_NT4_NR4_V3`，之后才能回答 plan 的完整跨 Tx 问题；8Tx 跨速度目前只有初始点估计，最终比较仍需端点追加和跨场景独立不确定性分析。
- 必须按 plan 对实际 1% bracket 中 0.5%--2% BLER 且少于 200 errors 的端点追加，之后重新计算 crossing、Wilson 区间和图。
- `A100_NT8_NR4_V60` 第一轮追加已落盘：5 个 SNR 各新增 1,000 paired trials，形成 30 条 curve interval record，数组引用均存在；调度状态和初始证据收据分别位于 `formal/A100_NT8_NR4_V60/append_scheduler/state.json` 与 `append_scheduler/baseline/manifest.json`。
- 当前没有 paired-bootstrap 区间；小于约 0.2 dB 的方案差异只作为点估计描述，不作稳定优劣结论。
- ideal-CSI Sidon/PRG6 诊断每点只有 1000 trials，尚无 paired 差值区间；它只支持趋势观察，不支持精确 crossing 或已确认 SNR 增益。
- 8Tx/V3 aged MRT 的巨大优势不能解释成“移动场景 MRT 总是更好”；它同时依赖低速度、约 5 ms CSI age、未量化精确 MRT和本轮无空间相关模型。
- 不模拟 ICI/CFO、反馈量化、反馈误码或多层检测。结果未经研究者最终确认，不更新 `KNOWLEDGE.md`/`GOALS.md`，不创建 Git checkpoint。
