# result-021

> 回填记录：本轮实验早于新工作流建立，原始过程记录在 `docs/archive/v_design_piecewise_tradeoff_experiment.md` 与英文综合报告 §5。本文件是按新格式对其结论的结构化回填，作为 plan-022 的起点。

## TL;DR（一句话结论）

N-series 分段线性相位在 **ideal-CSI BLER** 上确实优于 CDD（48PRB 多数点 +0.36~+1.55 dB），但换到 **estimated-CSI**（同一组稀疏 DMRS，Algorithm 1 matched RMMSE）后 **12/12 个候选 NMSE 全倒扣**（20dB 处 +1.8~+8.0 dB），BLER 无一稳定超过 CDD。→ 分集优势被稀疏导频下的估计代价抵消，符合 [[KNOWLEDGE.md]] B5/B6。**下一步（plan-022）**：把真实导频下的 RMMSE NMSE 直接纳入 `V` 的筛选目标,而非只用矩阵指标。

**相关图表**（`docs/figures/`）：ideal-CSI BLER `experiment21_ideal_csi_bler_{48,36,24}prb_full.svg`；matched RMMSE NMSE `experiment21_mmse_nmse_{48,36,24}prb.svg`。

## 仿真目的与验收标准

本轮（回填）验证的假设：矩阵层面 Pareto 前沿更优的 N-series 分段线性相位候选，其 ideal-CSI 分集增益应能传导到 estimated-CSI（稀疏 DMRS + Algorithm 1 matched RMMSE）的 NMSE/BLER 增益。隐含验收标准：存在候选在相同 DMRS 开销下 estimated-CSI BLER 不劣于匹配的 CDD 参考。本轮无正式 plan 文件（早于新工作流），原始实验定义见 `docs/archive/v_design_piecewise_tradeoff_experiment.md` §5。

## 配置回执（实际使用的参数）

| 项目 | 取值 |
|---|---|
| 物理信道 | Static TDL，指数 PDP，delay spread 5 ns |
| 天线 | 8Tx / 1Rx，单层 |
| SCS / FFT | 30 kHz / 4096，采样率 122.88 MHz |
| PDSCH | 10 OFDM symbols |
| 带宽 | 48 / 36 / 24 PRB（576 / 432 / 288 active subcarriers） |
| Segment 数 | 8（每段 72 / 54 / 36 subcarriers） |
| DMRS | symbol [2,7]，频域间隔 24 subcarriers |
| MCS | NR 256QAM table MCS 8，实际 16QAM，R=553/1024 |
| LDPC | 8 iterations |
| Ideal-CSI trials | 400 TB / SNR |
| RMMSE NMSE trials | 300 / SNR，SNR ∈ {0,4,8,12,16,20,24} dB |
| 估计器 | Algorithm 1 matched full-covariance RMMSE（各候选 `V` 用各自完整非平稳协方差，非 shifted-PDP 近似） |
| 候选 `V` | N-series（分段线性相位，相位连续）N3/N6/N8/N23，对照 CDD 参考点 C0-C7（按各带宽最接近的 `B_0.5` 匹配） |

## 关键指标（表格，对照计划里的预期值）

预期（来自 N-series 矩阵层面 Pareto 前沿优于 CDD 的假设）：ideal-CSI 分集增益应该同时传导到 estimated-CSI BLER/NMSE 增益。

**Ideal-CSI BLER 增益（N-series 相对匹配 CDD 参考点，正值=N更好）**

| 带宽 | 候选/CDD参考 | 10% BLER 增益 | 1% BLER 增益 |
|---|---|---:|---:|
| 48 PRB | N3/C0 | +1.03 dB | +3.00 dB |
| 48 PRB | N6/C0 | +1.55 dB | +1.60 dB |
| 48 PRB | N8/C1 | +0.36 dB | +1.44 dB |
| 48 PRB | N23/C2 | +0.62 dB | +1.00 dB |
| 36 PRB | N3/C5 | -0.02 dB | +1.00 dB |
| 36 PRB | N6/C6 | -0.47 dB | -0.06 dB |
| 36 PRB | N8/C3 | +1.02 dB | +2.71 dB |
| 36 PRB | N23/C4 | +0.58 dB | +0.80 dB |
| 24 PRB | N3/C6 | +0.01 dB | +1.10 dB |
| 24 PRB | N6/C6 | +0.14 dB | 0.00 dB |
| 24 PRB | N8/C6 | 0.00 dB | +0.33 dB |
| 24 PRB | N23/C4 | +0.29 dB | +2.25 dB |

**Matched RMMSE NMSE 倒扣（20 dB 处，正值=N比CDD差）**

| 带宽 | 候选/CDD参考 | 倒扣 |
|---|---|---:|
| 48 PRB | N3/C0 | +2.57 dB |
| 48 PRB | N6/C0 | +2.49 dB |
| 48 PRB | N8/C1 | +3.26 dB |
| 48 PRB | N23/C2 | +2.79 dB |
| 36 PRB | N3/C5 | +1.92 dB |
| 36 PRB | N6/C6 | +1.82 dB |
| 36 PRB | N8/C3 | +2.62 dB |
| 36 PRB | N23/C4 | +1.81 dB |
| 24 PRB | N3/C6 | +2.22 dB |
| 24 PRB | N6/C6 | +1.90 dB |
| 24 PRB | N8/C6 | **+8.04 dB** |
| 24 PRB | N23/C4 | +2.43 dB |

## 与预期的偏差

- 预期：ideal-CSI 分集增益（表1，多数为正）应传导为 estimated-CSI BLER 增益。
- 实际：所有 12 个候选点的 matched RMMSE NMSE 全部倒扣（表2 全正值，N-series 更差），且**没有一个候选在 estimated-CSI BLER 上稳定超过对应 CDD 参考点**。24PRB N8 尤其严重（+8.04 dB 倒扣，对应约 11-15% BLER floor）。
- 结论：矩阵层面/ideal-CSI 层面的分集优势被稀疏导频（spacing=24）下的信道估计代价完全抵消，符合 [[KNOWLEDGE.md]] 第 5 条。这不是本轮"失败"，是明确了当前 N-series + Algorithm 1 组合在稀疏导频下不可行，为 plan-022 的联合目标设计提供了明确的失败边界。

## 异常现象

- 24PRB N8：NMSE 从 12dB 的 4.80e-2 缓降到 24dB 的 2.31e-2，插值残差很强，怀疑目标 RE 相对稀疏 pilot 集合可预测性太差（该 segment 内非平稳协方差观测不足）。
- 当把 CDD 参考扩大到大 delay 候选（不限于匹配 B_0.5 的点）时，`step=64 samples` 的 C7 在 1% BLER 有最佳 ideal-CSI 尾部（48/36/24 PRB 分别 16.00/16.25/16.20 dB），说明仅用矩阵指标匹配 CDD 参考点会低估 CDD 本身的最优表现，需要在 plan-022 里同时纳入大 delay CDD 作为参考基线。

## 原始数据路径

- `outputs/v_design_balanced_slope_scan/v_design_balanced_slope_scan_20260619_162622/`：N-series 矩阵扫描与 Pareto front。
- `outputs/v_design_new_front_link_ideal_24prb_main/`、`_merged/`、`_tail/`，`outputs/v_design_new_front_link_large_cdd_ideal_36/`、`_48/`、`_48_tail/`：ideal-CSI BLER 结果。
- `outputs/v_design_new_front_link_nmse/`、`_nmse_24prb/`，`outputs/v_design_new_front_link_large_cdd_nmse_36/`、`_48/`：matched RMMSE NMSE 结果。
- `outputs/nmse_same_overhead_search_unified/search_20260619_164606/unified_gain_map.csv`：Algorithm 2/3 统一 420 点网格（用于 KNOWLEDGE.md 第 2/3 条）。
