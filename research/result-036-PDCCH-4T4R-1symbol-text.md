# result-036-text：PDCCH 4Tx、4Rx/2Rx、1-symbol、5 Hz 的 CDD 与 DFT cycling 对比

> 状态：**正式数据已完成，结果待研究者确认；本次仅执行 5 Hz，属于 plan-036 的缩减范围**。研究者明确要求暂不运行 1100 Hz，因此本文不回答跨 Doppler 不变性问题，也不把结果写入 `KNOWLEDGE.md` 或 `GOALS.md`。

## 1. 配置与复现

本次运行 AL1、AL2，并分别覆盖 4Tx/4Rx 与追加的 4Tx/2Rx。共同条件为单层、1-symbol CORESET、48 RB、non-interleaved CCE-to-REG、6-REG bundle、30 kHz SCS、4096 FFT、288-sample CP、TDL-C 300 ns、4 GHz、20 sinusoids、PDCCH Polar A=40/CRC24C/list size 8、QPSK 和 comb-4 DMRS。每个 RE 的 4Tx 总功率为 1；横轴 SNR 是总发射功率相对单个 Rx 分支噪声功率，独立 Rx 分支使用 coherent MRC。

请求 Doppler 为 5 Hz，通过 $v=f_Dc/f_c$ 换算为 `1.349066061 km/h`，继续调用平台已有的 Sionna TDL 实现；没有新增独立 Doppler、symbol 内时变或 ICI 模型。单 symbol 内只取一个信道时间采样。

候选为 Sidon non-transparent、QC-delay non-transparent、相同 QC 波形 transparent，以及 DFT4 cycling transparent。AL1 Sidon 坐标为 `[0,1,3,7]`，QC 为 `[0,0.24192,0.48384,0.72576]`，DFT order 为 `[0]`；AL2 Sidon 为 `[0,11,19,64]`，QC 为 `[0,0.41472,0.82944,1.24416]`，DFT order 为 `[0,1]`。ideal CSI 只保留 Sidon、QC 和 DFT 三个不重复波形。

正式配置及 SHA-256：

- `configs/pdcch_plan036_al1_fd5_formal.yaml`：`A99CFD454455EA3C7DC6E2AF99A97CF467592AE979B42945E16F19E8678D3B7D`
- `configs/pdcch_plan036_al2_fd5_formal.yaml`：`FEF5648BED69EE4DC1685ADE983BECDECED9F9495130B1A06131EDF97965EFBD`
- `configs/pdcch_plan036_al1_fd5_4t2r_formal.yaml`：`631C77F40E0F1BCEF932A0DDAB823F8BE4F80D08F60ED0C47C834B9EE212CA90`
- `configs/pdcch_plan036_al1_fd5_4t2r_snr8_extension.yaml`：`479CA421F64CC56C8E0F285D26EFB638AF2F9CB83992DFFAECFF061E2EA30222`
- `configs/pdcch_plan036_al2_fd5_4t2r_formal.yaml`：`4749FB2ACA987A137197666738382DFB55B66FFD510A0CB7F7C6A398D2989586`

代码基准提交为 `83b25bc2ea77324ae98f7dc67136fc1e93b1f32c`，仿真使用包含 plan-036 未提交实现的工作区，不能只靠该提交复现。执行命令为：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan036_pdcch.py --config configs\pdcch_plan036_al1_fd5_formal.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan036_pdcch.py --config configs\pdcch_plan036_al2_fd5_formal.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan036_pdcch.py --config configs\pdcch_plan036_al1_fd5_formal.yaml --stage analyze
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan036_pdcch.py --config configs\pdcch_plan036_al2_fd5_formal.yaml --stage analyze
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan036_pdcch.py --config configs\pdcch_plan036_al1_fd5_4t2r_snr8_extension.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan036_pdcch.py --config configs\pdcch_plan036_al1_fd5_4t2r_snr8_extension.yaml --stage analyze
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan036_pdcch.py --config configs\pdcch_plan036_al2_fd5_4t2r_formal.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan036_pdcch.py --config configs\pdcch_plan036_al2_fd5_4t2r_formal.yaml --stage analyze
```

每个 SNR 最少 2,000、最多 50,000 个共同 trial，每 1,000 trial 检查一次。七条曲线各自达到 200 errors、Wilson 95% 上限低于 0.01 或达到最大预算后，该 SNR 才停止。4Rx 单点样本范围：AL1 为 2,000–21,000，AL2 为 2,000–23,000；2Rx 单点样本范围：AL1 为 2,000–22,000，AL2 为 2,000–23,000。没有点达到 50,000 上限。

## 2. 原始数据与统计口径

完整逐点数据分别为：

- `outputs/experiment036_pdcch_4t4r_1symbol/20260918_formal/al1_fd5/bler_points.csv`
- `outputs/experiment036_pdcch_4t4r_1symbol/20260918_formal/al2_fd5/bler_points.csv`
- `outputs/experiment036_pdcch_4t2r_1symbol/20260919_formal/al1_fd5/bler_points.csv`（AL1 已追加至 8 dB）
- `outputs/experiment036_pdcch_4t2r_1symbol/20260919_formal/al2_fd5/bler_points.csv`

每份 CSV 有 112 行，即 16 个 SNR 点乘以 4 条 estimated 和 3 条 ideal 曲线；字段包含 trials、errors、原始 BLER、Wilson 95% 区间、estimated CE NMSE、pilot rank/condition 和停止原因。逐 trial error flags、estimated CE NMSE 和断点文件位于各场景的 `trial_data/`，逐轮 Wilson 与停止审计位于 `adaptive_status.json`，配置和配对审计位于 `run_metadata.json`。

目标 SNR 仅在相邻 1 dB 原始点形成双侧 bracket 时，对 `log10(BLER)` 线性插值。下表中的增益定义为 DFT baseline 目标 SNR 减候选目标 SNR，正值表示候选更好；CE 代价为 estimated 目标 SNR 减 ideal 目标 SNR，正值表示信道估计造成损失。

| AL | CSI | 方案 | 10% SNR / dB | 1% SNR / dB | 相对 DFT 10% / dB | 相对 DFT 1% / dB |
|---:|---|---|---:|---:|---:|---:|
| 1 | estimated | Sidon | -0.143 | 1.367 | -0.159 | +1.137 |
| 1 | estimated | QC non-transparent | -0.400 | 1.972 | +0.098 | +0.532 |
| 1 | estimated | QC transparent | -0.047 | 2.225 | -0.254 | +0.279 |
| 1 | estimated | DFT baseline | -0.302 | 2.504 | 0 | 0 |
| 1 | ideal | Sidon | -2.904 | -1.346 | +1.057 | +1.851 |
| 1 | ideal | QC | -2.114 | 0.006 | +0.267 | +0.499 |
| 1 | ideal | DFT baseline | -1.847 | 0.506 | 0 | 0 |
| 2 | estimated | Sidon | -3.705 | -2.476 | -0.725 | -0.022 |
| 2 | estimated | QC non-transparent | -4.418 | -2.592 | -0.012 | +0.094 |
| 2 | estimated | QC transparent | -3.933 | -2.098 | -0.498 | -0.400 |
| 2 | estimated | DFT baseline | -4.431 | -2.498 | 0 | 0 |
| 2 | ideal | Sidon | -7.281 | -6.042 | +0.589 | +1.153 |
| 2 | ideal | QC | -6.672 | -4.981 | -0.020 | +0.092 |
| 2 | ideal | DFT baseline | -6.691 | -4.889 | 0 | 0 |

QC transparent 相对 non-transparent 的接收代价为：AL1 在 10%/1% 分别为 `0.353/0.253 dB`；AL2 分别为 `0.486/0.494 dB`。

| AL | 波形 | estimated-to-ideal 10% / dB | estimated-to-ideal 1% / dB |
|---:|---|---:|---:|
| 1 | Sidon | 2.762 | 2.713 |
| 1 | QC | 1.714 | 1.965 |
| 1 | DFT baseline | 1.545 | 1.998 |
| 2 | Sidon | 3.576 | 3.566 |
| 2 | QC | 2.254 | 2.389 |
| 2 | DFT baseline | 2.261 | 2.390 |

estimated CE NMSE 在 `-10 dB → 5 dB` 的端点如下：

| AL | Sidon / dB | QC non-transparent / dB | QC transparent / dB | DFT / dB |
|---:|---:|---:|---:|---:|
| 1 | -1.152 → -8.813 | -2.265 → -10.875 | -1.610 → -10.187 | -2.551 → -11.156 |
| 2 | -1.313 → -8.426 | -2.762 → -11.657 | -2.009 → -10.611 | -2.825 → -11.581 |

4Rx 的所有曲线在 `[-10,5] dB` 内均有原始 `BLER<=0.01` 点，满足本次 1% 可见性验收。

### 2.1 4T2R 追加结果

4T2R AL1 在原 `-10:1:5 dB` 网格上部分曲线尚未达到 1%，因此按研究者要求只追加 6、7、8 dB；原 5 dB 数据不重跑。合并后 AL1 有 19 个 SNR 点和 133 行，AL2 保持 16 个 SNR 点和 112 行。

| AL | CSI | 方案 | 10% SNR / dB | 1% SNR / dB | 相对 DFT 10% / dB | 相对 DFT 1% / dB |
|---:|---|---|---:|---:|---:|---:|
| 1 | estimated | Sidon | 3.455 | 5.648 | +0.691 | +2.250 |
| 1 | estimated | QC non-transparent | 3.823 | 7.270 | +0.323 | +0.627 |
| 1 | estimated | QC transparent | 4.066 | 7.438 | +0.080 | +0.460 |
| 1 | estimated | DFT baseline | 4.146 | 7.898 | 0 | 0 |
| 1 | ideal | Sidon | 1.073 | 3.092 | +1.726 | +3.050 |
| 1 | ideal | QC | 2.268 | 5.410 | +0.531 | +0.732 |
| 1 | ideal | DFT baseline | 2.799 | 6.142 | 0 | 0 |
| 2 | estimated | Sidon | -0.765 | 0.935 | +0.002 | +0.989 |
| 2 | estimated | QC non-transparent | -0.869 | 1.745 | +0.105 | +0.179 |
| 2 | estimated | QC transparent | -0.527 | 2.036 | -0.236 | -0.112 |
| 2 | estimated | DFT baseline | -0.764 | 1.924 | 0 | 0 |
| 2 | ideal | Sidon | -3.762 | -2.194 | +1.067 | +2.008 |
| 2 | ideal | QC | -2.731 | -0.265 | +0.036 | +0.078 |
| 2 | ideal | DFT baseline | -2.695 | -0.186 | 0 | 0 |

4T2R 的 QC transparency 代价在 AL1 的 10%/1%点为 `0.243/0.168 dB`，在 AL2 为 `0.342/0.292 dB`。estimated-to-ideal 代价如下：

| AL | 波形 | 10% / dB | 1% / dB |
|---:|---|---:|---:|
| 1 | Sidon | 2.381 | 2.555 |
| 1 | QC | 1.555 | 1.860 |
| 1 | DFT baseline | 1.347 | 1.756 |
| 2 | Sidon | 2.997 | 3.129 |
| 2 | QC | 1.861 | 2.009 |
| 2 | DFT baseline | 1.931 | 2.110 |

4T2R AL1 在 8 dB 的四条 estimated 原始 BLER 分别为 Sidon `4/22000=0.000182`、QC non-transparent `126/22000=0.005727`、QC transparent `151/22000=0.006864`、DFT `204/22000=0.009273`，均小于或等于 1%。其中 DFT 的 Wilson 95%上限为 `0.010628`，略高于 1%；该点因达到 200 errors 停止，不能表述为“95%置信上限低于 1%”。

## 3. 分析与结论

事实：ideal CSI 下 Sidon 在 AL1 和 AL2 均明显优于 DFT，1%增益分别为 `1.851 dB` 和 `1.153 dB`。estimated CSI 下该优势大幅缩小：AL1 Sidon 在 1%仍有 `1.137 dB` 增益，但在 10%略差 `0.159 dB`；AL2 Sidon 在 10%差 `0.725 dB`，1%基本持平 `-0.022 dB`。这与 Sidon 的 estimated-to-ideal 代价和较差 CE NMSE一致。

事实：QC non-transparent 在 estimated CSI 下与 DFT 接近或略好；AL1 的 1%增益为 `0.532 dB`，AL2 为 `0.094 dB`。同一 QC 波形使用 physical-fullband 协方差后，目标 SNR 增加约 `0.25–0.49 dB`，说明透明接收机的协方差失配存在可测代价。

事实：减少到 2Rx 后，所有方案的目标 SNR整体升高，但 Sidon 相对 DFT 的 1%优势增大：AL1 estimated 从 4Rx 的 `1.137 dB` 增至 2Rx 的 `2.250 dB`，AL2 从近似持平的 `-0.022 dB` 增至 `0.989 dB`。这是本次有限配置下的观察，不应外推为接收天线减少必然提升 CDD 相对增益。

推断：本结果支持“严格 Sidon 的波形级分集增益可能被单-symbol comb-4 信道估计误差抵消”，尤其在 AL2。该解释由 ideal/estimated 目标差和 CE NMSE排序共同支持，但没有把估计误差分解为 bias、variance 或模型失配成分，因此不是机制的完整证明。

统计限制：表中目标差为原始点插值的点估计。本次尚未执行 plan-036 原定的 4,000 次 absolute-trial paired bootstrap，因此不提供目标差置信区间，也不能把接近零的差值（例如 AL2 QC/DFT 和 Sidon/DFT 的 1%差）解释为统计显著。单点 Wilson 区间保存在原始 CSV。

范围限制：只执行了 5 Hz。没有 1100 Hz 正式数据，不能据此报告 Doppler 不变性；同时当前平台不含 symbol 内时变或 ICI，本结果也不能支持“高移动性无损”的结论。结果仅适用于表列的 4Tx/4Rx 与 4Tx/2Rx、独立 Rx、TDL-C300、单层和接收机知识口径。

## 4. 未完成项

1. 未运行 1100 Hz；这是研究者明确缩减的范围。
2. 未执行 4,000 次 paired bootstrap，目标差区间仍缺失。
3. 图已由本地脚本生成并通过文件存在性和 `1368×917` 像素尺寸检查，但未由 Agent 加载；横纵坐标、图例遮挡和缩小到约 13 cm 后的可读性仍需研究者人工核验。
4. result 尚待研究者确认；确认前不更新全局结论或创建 Git checkpoint。
