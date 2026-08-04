# KNOWLEDGE：Agent 精简知识

## 1. 使用方法

本文件记录实验或明确推导支持的结论、已排除方向、开放问题和证据路径。制定 plan 时同时读取 `GOALS.md`、`DESIGN.md` 和最新无图版 result。历史推导按需读取 `docs/design/DESIGN_ANNOTATED.md` 或 `docs/archive/`。

## 2. 当前状态

在 48 PRB、8 Tx / 1 Rx、DMRS comb 24 和 V-aware matched LMMSE 条件下，严格 Sidon 相对 QC 的 estimated-CSI BLER 优势已在平坦独立分支信道和近似平坦 TDL-A 中复现；TDL-A 5 ns 下只保留小幅 10% 优势，1% 排序不能确定。只使用 $T_\epsilon$ 的厚 Sidon 规则在 TDL-A 10 ns、TDL-C 5 ns 和 TDL-C 10 ns 中得到约 0.16–0.20 dB 的 10% ideal-CSI outage 增益，且 16 dB CE 变化不超过 0.152 dB。result-026 尚未隔离完整协方差知识与候选空间放宽的独立贡献；下一阶段应使用协方差加权 $M_2^{\rm eff}/M_4^{\rm eff}$、展宽感知等差基线和匹配的 T1/T2 对照回答该问题。

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

### K14. TDL-A 5 ns 下原 Sidon 只保持小幅 10% BLER 优势

result-025 的 delay-matched `/576` 补做保持 plan-024 的 CDD 相位矩阵、物理时延、资源、DMRS、功率、MCS 和统计口径，只把底层平坦分支信道替换为 Sionna TDL-A 5 ns，并使接收端协方差与各候选真实复合信道匹配。Sidon 相对 QC 的 10% BLER 改善为 0.156 dB，保守 95% 区间 [0.045, 0.266] dB；1% 改善为 0.090 dB，保守区间 [-0.210, 0.391] dB。TDL 加权后二阶代理为 QC 略低，四阶代理仍为 Sidon 明显较低；Sidon 的 data-RE CE NMSE 较差 0.754–1.271 dB。该结果说明严格 Sidon 四元条件在物理展宽信道下不足以保证平坦模型中的增益幅度，不证明任意 TDL 下 Sidon 或厚 Sidon 都无效。证据：`research/result-025-text.md`。

### K15. 严格 Sidon 的 estimated-CSI 优势在近似平坦 TDL-A 中恢复

result-026 在 TDL-A 0.1/1/5 ns 下得到严格 Sidon 相对 QC 的 10% estimated-CSI BLER 增益分别为 0.304/0.402/0.174 dB，保守 95% 区间分别为 [0.169, 0.439]、[0.256, 0.547]、[0.055, 0.292] dB。对应 1% 点估计为 1.091/1.242/-0.003 dB；0.1 ns 的 QC 目标明显外推，1 ns 的 QC 目标轻微外推，5 ns 区间 [-0.231, 0.224] dB 跨 0。结合 result-024/025，只能确认严格 Sidon 在平坦或近似平坦 NLoS 下有优势；物理展宽增加后，CE 与相关结构共同决定最终 BLER，无权 Sidon 条件不是充分保证。证据：`research/result-026-text.md` E3。

### K16. 厚 Sidon 在三个预定 NLoS 场景中表现出小幅 ideal-CSI outage 潜力

只使用 $T_\epsilon$ 和系统几何生成的硬厚 Sidon 候选，在 TDL-A 10 ns、TDL-C 5 ns、TDL-C 10 ns 下相对 T1 等差基线的 10% ideal-CSI outage 增益分别为 0.160/0.202/0.174 dB，1% 增益分别为 0.234/0.291/0.248 dB；16 dB matched CE NMSE 变化分别为 +0.152/+0.026/-0.065 dB。该结果支持厚 Sidon 作为后续 estimated-CSI 候选生成规则，但不等价于 BLER 增益。TDL-A 20/30 ns 虽存在硬几何候选，预冻结 outage 子集没有覆盖它们，因此不能据此判定硬厚 Sidon 在这两个场景中的性能。证据：`research/result-026-text.md` E2。

### K17. result-026 未证明完整协方差知识具有独立 outage 增益

TDL-A 5 ns 下，完整协方差搜索得到的 E1_0028 相对本轮等差 B1 的 10%/1% ideal-CSI outage 增益为 0.161/0.281 dB，16 dB matched CE 劣化 0.688 dB，但没有通过预定的近零噪声 CE 门槛。更重要的是，本轮没有在相同非等差候选空间内设置只知道 $T_\epsilon$ 的匹配选择器，因此该结果只能说明存在比该等差基线 outage 更好的非等差候选，不能把改善唯一归因于额外知道完整协方差。后续应分别比较 T1/T2 知识等级和等差/非等差候选空间。证据：`research/result-026-text.md` E1。

### K18. result-026 的当前 JCDD 数值核不能继续用于候选排序

正式结果后的数值复核发现，当前 10 阶四维 Gauss–Hermite 实现会在中高 SNR 下把大段中等相关系数对应的互信息协方差算成负值，随后经非负裁剪和单调化变成数值零；保存的 16/18 dB 核分别只有 7/201 和 1/201 个相关系数网格点为非零。Monte Carlo 复核显示这些区间的真实协方差并非零。因此，result-026 的 JCDD 绝对值、数值地板和基于该核的 Pareto 排序不能作为后续知识；Monte Carlo outage、matched CE 以及直接计算的 $M_2^{\rm eff}/M_4^{\rm eff}$ 不受该问题影响。修复并验证核之前，后续搜索应直接使用有效相关矩或 Monte Carlo 判据。证据：`cdd_lls/design/cdd_metrics.py`，以及 `outputs/experiment026_cdd_design/20260724_main/validation/` 下的 `jcdd_kernel_14db.npz`、`jcdd_kernel_16db.npz`、`jcdd_kernel_18db.npz`。

### K19. result-026 的 TDL-D/E 大增益不能推广为一般 LoS Sidon 结论

本轮 Sionna TDL-D/E 在未传入空间相关矩阵时，对不同 Tx 生成独立 diffuse 分量，但把同一个 specular 复系数广播到全部 Tx，等效于未加入阵列方向相位的全 1 LoS steering。该实现下得到的约 6.6 dB 10% BLER 分离是有价值的实现诊断，但没有经过 CDL、显式阵列 steering、K-factor 或无 CDD 基线消融，不能推广为一般 38.901 多天线 LoS 信道下的 Sidon 增益。证据：`cdd_lls/phy/channel_tdl.py`、本轮使用的 Sionna `tdl.py`、`research/result-026-text.md` E3。

## 5. 已排除或暂停方向

1. 不再单独调整 delay support、diagonal loading 或 pairwise 解耦优化 Algorithm 2B/2C，除非改变导频观测。
2. 不使用矩阵代理量作为最终验收。
3. 不把等差 CDD 大 delay 步长本身解释为平坦模型分集增益。
4. 不继续把完全透明 CC/CN 家族作为当前主要方向。
5. 不把单个半透明通过点表述为完全透明方案成功。
6. 不把原严格 Sidon 集直接推广为任意 TDL 下的最优设计；后续 CDD 搜索必须同时检查实际物理协方差、CE 和 outage/BLER。

## 6. 当前开放问题

1. 在 TDL-A 1/5/10/30/100 ns 下，完整协方差加权 $M_2^{\rm eff}/M_4^{\rm eff}$ 选择相对展宽感知等差基线的 10%/1% outage 增益；
2. 严格 Sidon 和厚 Sidon 的适用范围与 $T_\epsilon$、导频密度、发射天线数、子载波数及物理 PDP 的关系；
3. 更高导频密度能否在保留 Sidon/厚 Sidon 分集潜力时改善 matched CE，以及导频能量和开销的公平口径；
4. PDP 或协方差失配、定时误差和 24/36/48 PRB 一致性；
5. 深尾优势大于 outage 预测的原因；
6. CE-error-aware LLR 的影响；
7. JCDD 数值核的稳定计算，以及几何 LoS/CDL 条件下 Sidon 与 CDD 的重新评价。

## 7. 证据索引

`GOALS.md`、`DESIGN.md`、`docs/design/DESIGN_ANNOTATED.md`、`research/result-021.md`、`research/result-023.md`、`research/result-024.md`、`research/result-024-text.md`、`research/result-025.md`、`research/result-025-text.md`、`research/result-026.md`、`research/result-026-text.md`。
