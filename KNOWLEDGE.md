# KNOWLEDGE：Agent 精简知识

## 1. 使用方法

本文件记录实验或明确推导支持的结论、已排除方向、开放问题和证据路径。制定 plan 时同时读取 `GOALS.md`、`DESIGN.md` 和最新无图版 result。历史推导按需读取 `docs/design/DESIGN_ANNOTATED.md` 或 `docs/archive/`。

## 2. 当前状态

在 48 PRB、8 Tx / 1 Rx、平坦分支信道、DMRS comb 24、V-aware matched LMMSE 条件下，Sidon delay 集 `[0,1,3,7,12,20,30,65]` 相对等差 QC CDD 的 10% BLER 改善 0.33 dB，1% 改善 0.85 dB，两个保守 95% 区间均不跨 0。下一阶段验证 5–100 ns PDP、失配、定时误差和带宽变化。

## 3. 平台与命名

- `run.py` 是通用入口；YAML 在 `configs/`；专题入口在 `tools/`；产物在 `outputs/`。
- 已实现 static TDL、单层 PDSCH、2/4/8 Tx、1/4 Rx、CDD/PRG/一般 `\mathbf V`、IDEAL/RMMSE/reconstruction 类估计和 bit-level BLER。
- Algorithm 1：直接等效信道 RMMSE。
- Algorithm 2B：多导频 delay-domain basis LMMSE。
- Algorithm 2C：local-window deterministic。
- Algorithm 3：non-CDD per-port DMRS。
- 符号以 `DESIGN.md` 为准。

## 4. 已确立结论

### K1. Algorithm 1 是共享同一组合导频观测时的主要基线

Algorithm 2B/2C 没有增加独立观测，在统一网格上没有稳定增益。证据：`docs/FINDINGS.md` B1–B2。

### K2. Algorithm 2B/2C 的限制是观测秩和条件数

已测重构矩阵条件数约 `10^{12}`–`10^{13}`，未知量零空间比例可超过 70%。缩短 delay support、增大 diagonal loading 和 pairwise 解耦均未解决。证据：`docs/FINDINGS.md` B2。

### K3. Algorithm 3 的增益来自参考信号改变

代价是额外 DMRS 开销或更低每端口导频密度，比较必须核算开销。证据：`docs/FINDINGS.md` B3。

### K4. 矩阵指标和 ideal-CSI 不能替代 estimated-CSI

N-series 在矩阵与 ideal-CSI 上可优于 CDD，但 result-021 的 12/12 个 estimated-CSI NMSE 点全部劣化，最大倒扣 8.04 dB。候选必须经过实际 DMRS 下的 NMSE 和编码 BLER。证据：`research/result-021.md`。

### K5. 大 CDD delay 的性能不是单调函数

性能由全带正交、导频折叠位置、物理 PDP 和接收机共同决定。证据：`docs/FINDINGS.md` B7、`research/result-023.md`。

### K6. 当前完全透明分段方向未达到验收标准

result-023 中 48 PRB 和 24 PRB 共 232 个 CC/CN 挑战者均未达到透明 CDD/PRG 基线支配门槛。证据：`research/result-023.md`。

### K7. 分段边界信息不是 B3/B4 的通用修正

result-024 E1：48 PRB B3 改善 6/16、B4 改善 0/100；24 PRB B3 改善 7/16、B4 改善 0/100。48 PRB 无通过点，24 PRB 只有 `B3_cc_nseg8_T6_seq` 半透明诊断点通过。证据：`research/result-024-text.md`。

### K8. 栅格等差 CDD 在平坦模型下同轨道等价

等差步长 2、4、9 相对步长 1 的 10%/1% outage 差在 0.0001 dB 内。证据：`research/result-023.md` H3。

### K9. Sidon delay 集改善平坦模型 outage

`[0,1,3,7,12,20,30,65]` 相对等差 CDD 的 10%/1% outage 改善为 0.210/0.368 dB。该集合 36 个无序整数成对和均不同。证据：`research/result-023.md`、`research/result-024.md`。

### K10. Sidon 改善传递到 V-aware estimated-CSI BLER

| 目标 | QC SNR | Sidon SNR | 改善 | 保守 95% 区间 |
|---|---:|---:|---:|---:|
| 10% BLER | 14.65 dB | 14.32 dB | 0.33 dB | [0.20, 0.45] dB |
| 1% BLER | 17.23 dB | 16.38 dB | 0.85 dB | [0.63, 1.07] dB |

1% 区间使用 3000 trials/点。证据：`research/result-024-text.md`。

### K11. result-024 的改善不是 CE NMSE 差导致

全部加密点中 `NMSE_{\rm Sidon}-NMSE_{\rm QC}` 最大绝对值 0.060 dB，且正负均有。高阶频域相关结构与有限码长译码是待验证解释。

### K12. Sidon 结论只在限定条件内成立

已验证：48 PRB、8 Tx / 1 Rx、单层、DMRS comb 24、平坦独立分支信道、V-aware matched LMMSE。未验证：5–100 ns PDP、定时误差、协方差失配、24/36 PRB 和移动性。

### K13. LLR 尚未显式加入信道估计误差项

当前 LLR 噪声方差没有加入 CE-error-aware 项。该问题尚未验证，不构成推翻 result-024 成对比较的证据。证据：`docs/FINDINGS.md` B8。

## 5. 已排除或暂停方向

1. 不再单独调整 delay support、diagonal loading 或 pairwise 解耦优化 Algorithm 2B/2C，除非改变导频观测。
2. 不使用矩阵代理量作为最终验收。
3. 不把等差 CDD 大 delay 步长本身解释为平坦模型分集增益。
4. 不继续把完全透明 CC/CN 家族作为当前主要方向。
5. 不把单个半透明通过点表述为完全透明方案成功。
6. 完成 Sidon 稳健性检查前，不启动新候选大规模搜索。

## 6. 当前开放问题

1. 5–100 ns PDP 下的导频矩阵条件、matched NMSE 和 BLER；
2. PDP 或协方差失配；
3. 定时误差；
4. 24 / 36 / 48 PRB 一致性；
5. 深尾优势大于 outage 预测的原因；
6. CE-error-aware LLR 的影响。

## 7. 证据索引

`GOALS.md`、`DESIGN.md`、`docs/design/DESIGN_ANNOTATED.md`、`research/result-021.md`、`research/result-023.md`、`research/result-024.md`、`research/result-024-text.md`。
