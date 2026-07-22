# result-024：分段边界诊断与 Sidon/QC V-aware 链路终审

> 对应 `research/plan-024.md`。无图版：`research/result-024-text.md`。代码：`tools/run_experiment024_segment_sidon_qc.py`。

## 1. TL;DR

**E1：分段对齐估计不是 B3/B4 的通用修正。** 48 PRB 的 B3 改善 6/16、B4 改善 0/100；24 PRB 的 B3 改善 7/16、B4 改善 0/100。48 PRB 无通过点，24 PRB 只有 `B3_cc_nseg8_T6_seq` 在已知边界的半透明条件下通过。

**E2：Sidon 优势传递到有限码长 estimated-CSI BLER。** 在 48 PRB、平坦 8 分支信道、相同 DMRS 和双方 matched LMMSE 下，Sidon 的 10% BLER SNR 改善 0.33 dB，保守 95% 区间 [0.20,0.45] dB；1% 改善 0.85 dB，保守 95% 区间 [0.63,1.07] dB。matched CE NMSE 差最大绝对值 0.060 dB。

## 2. 配置回执

E1：8 Tx / 1 Rx、30 kHz、48/24 PRB、平坦模型、DMRS comb 24、16 dB NMSE 点、全部 B3/B4。R4 知道分段数和边界，不知道相位；每段只用段内导频；均匀 PDP 先验 [0,300 ns]。

E2：48 PRB、`K=576`、平坦 `h~CN(0,I_8)`。QC 索引 `[0,9,18,27,36,45,54,63]`，Sidon 索引 `[0,1,3,7,12,20,30,65]`。DMRS comb 24；双方各用真实 `R_g=VV^H` 的 matched LMMSE；16QAM、MCS 8、码率 553/1024、8 次 LDPC；共同随机数。粗扫 400 trials/点，10% 和 1% 区间加密为 3000 trials/点。

## 3. E1 结果

`gain_R4_vs_R1 = NMSE_R1(dB) - NMSE_R4(dB)`，正值表示 R4 更好。

| 带宽 | 家族 | 候选数 | 改善数 | 中位变化 | 最大改善 | 最差变化 |
|---|---:|---:|---:|---:|---:|---:|
| 48 PRB | B3 | 16 | 6 | -0.53 dB | +1.27 dB | -3.43 dB |
| 48 PRB | B4 | 100 | 0 | -2.46 dB | -2.02 dB | -2.93 dB |
| 24 PRB | B3 | 16 | 7 | -0.43 dB | +1.63 dB | -3.24 dB |
| 24 PRB | B4 | 100 | 0 | -4.41 dB | -3.06 dB | -5.97 dB |

![E1 分段对齐 NMSE 改善](../docs/figures/v_design_progress_h1_h3/segment_nmse_gain.png)

| 带宽 | 候选 | R1 NMSE | R4 NMSE | 改善 | H1 复判 |
|---|---|---:|---:|---:|---|
| 48 PRB | `B3_cc_nseg16_T6_bitrev` | -6.05 | -7.31 | +1.27 dB | 不通过 |
| 48 PRB | `B3_cc_nseg16_T6_seq` | -6.59 | -7.85 | +1.26 dB | 不通过 |
| 24 PRB | `B3_cc_nseg8_T6_seq` | -6.50 | -8.13 | +1.63 dB | 通过 |
| 24 PRB | `B3_cc_nseg8_T6_bitrev` | -5.97 | -7.57 | +1.59 dB | 不通过 |

E1 只支持“边界不对齐是部分 B3 损失来源”，不支持完全透明分段方案成功。

## 4. E2 结果

目标 SNR 使用加密区间全部点做二项 logit 拟合。差值区间未利用正配对协方差，因此属于保守近似。

| 目标 | QC SNR，95% 区间 | Sidon SNR，95% 区间 | 改善 | 改善的保守 95% 区间 |
|---|---:|---:|---:|---:|
| 10% BLER | 14.65，[14.56,14.74] dB | 14.32，[14.23,14.41] dB | 0.33 dB | [0.20,0.45] dB |
| 1% BLER | 17.23，[17.06,17.40] dB | 16.38，[16.24,16.52] dB | 0.85 dB | [0.63,1.07] dB |

![E2 Sidon 与 QC 的 estimated-CSI BLER](../docs/figures/v_design_progress_h1_h3/sidon_qc_bler_refined.png)

1% 区间直接计数：

| SNR | QC 错误/3000 | QC BLER | Sidon 错误/3000 | Sidon BLER |
|---:|---:|---:|---:|---:|
| 16.00 | 119 | 3.967% | 48 | 1.600% |
| 16.25 | 77 | 2.567% | 30 | 1.000% |
| 16.50 | 64 | 2.133% | 31 | 1.033% |
| 16.75 | 58 | 1.933% | 19 | 0.633% |
| 17.00 | 41 | 1.367% | 12 | 0.400% |
| 17.25 | 30 | 1.000% | 13 | 0.433% |
| 17.50 | 19 | 0.633% | 8 | 0.267% |

加密点 `NMSE_Sidon(dB)-NMSE_QC(dB)` 最大绝对值 0.060 dB，且正负均有。

## 5. 验收判定

| 判据 | 实际结果 | 判定 |
|---|---|---|
| E1 改善不低于 1 dB | 有 B3 点超过 1 dB | 第一条件通过 |
| E1 同时满足 H1 门槛 | 48 PRB 0；24 PRB 1 个半透明点 | 不支持完全透明方案 |
| E2 10% 改善不低于 0.15 dB | 0.33 dB；区间下界 0.20 dB | 通过 |
| E2 1% 附近至少 30 错误 | 两候选均有直接锚点 | 通过 |
| E2 CE 公平性 | 最大 NMSE 差 0.060 dB | 通过 |

## 6. 解释、边界和复现

result-023 outage 预测的 10%/1% 改善为 0.210/0.368 dB，本轮 BLER 为 0.33/0.85 dB。方向一致，深尾 BLER 改善更大。高阶频域相关结构与有限码长译码是待验证解释。

结论只适用于 48 PRB、平坦 8 分支信道、V-aware matched LMMSE。未验证 5–100 ns PDP、定时误差、协方差失配、其他带宽和移动性。

原始数据位于 `outputs/experiment024_segment_sidon_qc/20260716_main/`。复现：

```bash
python tools/run_experiment024_segment_sidon_qc.py --stage segment
python tools/run_experiment024_segment_sidon_qc.py --stage link \
  --snrs 16,16.25,16.5,16.75,17,17.25,17.5 \
  --trials 3000 --batch-size 20 \
  --out outputs/experiment024_segment_sidon_qc/20260716_main/refine_1pct
python tools/analyze_experiment024.py
```
