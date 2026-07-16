# FINDINGS —— 已确立结论与证据（给人看的详尽版）

> **定位**：本文件面向研究者本人,是 `KNOWLEDGE.md` 中每条结论的**证据展开**——带仿真数字、条件、图表。
> 分析 agent 默认不读本文件（它读精简的 `KNOWLEDGE.md` 即可）;只有当我要复核某条结论的原始依据、或写论文/汇报时才看这里。
> 每当某条结论**固化或被推翻**,我在更新 `KNOWLEDGE.md` 的同时更新这里对应小节。小节编号与 `KNOWLEDGE.md` 的 B1~B8 一一对应。
> 图表路径相对本文件（`docs/`）,即 `figures/...`。

## 目录

- [B1. Algorithm 1 是最强通用 baseline](#b1)
- [B2. Algorithm 2B/2C 在统一网格上没有稳健增益](#b2)
- [B3. Algorithm 3 在稀疏导频+大 CDD delay 能拿大增益](#b3)
- [B4. N-series 分段线性相位在矩阵层面改善 Pareto 前沿](#b4)
- [B5. ideal-CSI 分集增益换不到 estimated-CSI BLER 增益（核心痛点）](#b5)
- [B6. 矩阵指标只能筛选,不能当最终目标](#b6)
- [B7. CDD delay 增大对 BLER 非单调](#b7)
- [B8. 接收机 LLR 未注入 CE 误差自噪声](#b8)

---

<a id="b1"></a>
## B1. Algorithm 1（matched direct RMMSE）是最强通用 baseline

**主张**：在所有估计器共享同一组 CDD-combined DMRS 观测的前提下,直接对等效信道 `g[k]` 做 matched full-covariance RMMSE 是最强的通用 baseline。

**为什么**：
- 它直接估计检测器真正需要的量（等效信道 `g[k]`）,不绕道去分离物理分支;
- 用上了已知 CDD 的 shifted-PDP 统计（matched covariance）,理论上就是等效信道的 LMMSE;
- 不引入分支分离带来的额外病态。

**证据**：统一 420 点算法对比网格上,没有任何 Algorithm 2B/2C 配置能稳健超过它（见 B2 数字）。这是后续所有"能不能超过 baseline"讨论的锚点。

**出处**：`docs/archive/QC_explicit_CDD_link_level_simulation_plan.md` §6;英文综合报告 §2.1 与 §6#1。实验编号 E1–E7。

---

<a id="b2"></a>
## B2. Algorithm 2B/2C 在统一网格上没有稳健 NMSE 增益

**主张**：物理分支重构类算法（2B basis-LMMSE / 2C local-window）相对 Algorithm 1 没有稳健增益。

**证据（数字）**：
- 统一 420 点网格上,**最佳观测到的 Algorithm 2B 增益仅 +0.014 dB**——小到无法与有限 trial 波动区分。
- 重构敏感矩阵 `Phi` 的病态诊断（实验 8）:

| Case | Delay spread | CDD delay | DMRS spacing | Pilots | 未知 taps | `Phi` cond | Nullity |
|---|---:|---:|---:|---:|---:|---:|---:|
| ds30_d16_sf6 | 30 ns | 16 | 6 | 96 | 34 | 2.25e13 | 13 |
| ds30_d32_sf6 | 30 ns | 32 | 6 | 96 | 34 | 8.45e12 | 10 |
| ds100_d16_sf6 | 100 ns | 16 | 6 | 96 | 114 | 1.94e13 | 84 |
| ds100_d32_sf6 | 100 ns | 32 | 6 | 96 | 114 | 1.55e13 | 81 |
| ds30_d64_sf2 | 30 ns | 64 | 2 | 288 | 34 | 7.59e12 | 6 |

**为什么**：这是**观测受限（observation-limited）**问题,不是估计器不够聪明。从单层 CDD-combined 标量导频观测反推多个物理分支的 delay-domain tap,矩阵条件数普遍 1e12~1e13、nullity 可达未知量的 70%+。重构编码了有用的物理结构,但**不增加独立观测数**,因此性能被敏感矩阵条件数/模型 support/局部平坦误差卡死。

**已排除的修补**（都试过、无效,别再动）：缩短 delay support 到 90/95/99% 能量、diagonal loading 扫 1e-6~1e-2、pairwise 解耦（Reconstruction-A,条件数同样爆炸）。要突破必须**改变导频观测方式本身**（这正是 Algorithm 3 有效的原因,见 B3）。

**图**：`figures/reconstruction_conditioning_summary.png`、`figures/experiment8_covariance_correlation_with_bc.png`。
**出处**：英文综合报告 §2.2、§3.4、§6#2;实验记录 §8。实验编号 E7、E8。

---

<a id="b3"></a>
## B3. Algorithm 3 在「稀疏导频 + 大 CDD delay」regime 能拿非常大的增益

**主张**：让 DMRS 不加 CDD、按端口正交发送并估计更平滑的物理端口信道,再用已知 delay 合成等效信道,在特定 regime 下大幅优于 Algorithm 1。

**证据（数字）**：在 CDD delay = 512 samples、DMRS spacing = 24 subcarriers 下,15 个 delay-spread × SNR 组合**全部正增益,最大 +18.237 dB**。

**为什么**：稀疏导频会对快变的 CDD 等效信道产生混叠(aliasing),但对更平滑的物理端口信道仍然够用。Algorithm 3 改的是参考信号设计本身,绕开了 B2 的观测受限,而不是在同一组受限观测上换估计器。

**代价与公平性**：Algorithm 3 需要额外 DMRS 开销,或在同总开销下降低每端口导频密度。任何系统级比较必须显式声明公平性模式（equal-total-overhead vs equal-per-port-density,定义见 `docs/archive/QC_explicit_CDD_link_level_simulation_plan.md` §4.3）,否则增益不可信。

**图**：`figures/experiment10_unified_alg3_gain_map_snr*.png`。
**出处**：英文综合报告 §2.4、§3.3、§6#3;实验记录 §17–§19。实验编号 E17–E19。

---

<a id="b4"></a>
## B4. N-series 分段线性相位在矩阵层面确实改善分集/相干带宽 Pareto 前沿

**主张**：相位连续的局部 slope 切换（N-series）在很大一部分分集/相干带宽平面上优于常规全带 CDD,且该优势在 ideal-CSI BLER 上可见。

**证据（数字,ideal-CSI BLER 增益,正值=N-series 更好）**：

| 带宽 | 候选/CDD参考 | 10% BLER 增益 | 1% BLER 增益 |
|---|---|---:|---:|
| 48 PRB | N3/C0 | +1.03 dB | +3.00 dB |
| 48 PRB | N6/C0 | +1.55 dB | +1.60 dB |
| 48 PRB | N8/C1 | +0.36 dB | +1.44 dB |
| 48 PRB | N23/C2 | +0.62 dB | +1.00 dB |

（36/24 PRB 结果变为候选相关、区域相关,48PRB 的排序不能直接照搬——见 B6。）

**为什么**：CDD 只有 `N` 个 delay 自由度,一般 `V` 有 `K×N` 个相位自由度;分段 slope 切换能让不同频段的 off-diagonal 相关项在全带累加时更好抵消,从而在相同 group-delay spread 下拿到更低 mutual coherence / 更大最小奇异值。

**图**：`figures/experiment18_v_design_tradeoff_scatter.svg`、`figures/experiment21_ideal_csi_bler_48prb_full.svg`。
**出处**：`docs/archive/V_design_diversity_CE_tradeoff.md`;英文综合报告 §4、§5.2、§6#4。实验编号 E21。

---

<a id="b5"></a>
## B5.【核心痛点】ideal-CSI 分集增益换不到 estimated-CSI BLER 增益

**主张**：B4 里 N-series 的 ideal-CSI 分集优势,在稀疏 DMRS + Algorithm 1 RMMSE 的真实估计链路上消失甚至倒扣。**这就是"结果不符合预期"的确切位置。**

**证据（数字,matched RMMSE NMSE 在 20dB 处的倒扣,正值=N-series 比 CDD 更差）**：

| 带宽 | 候选/CDD参考 | N NMSE @20dB | CDD NMSE @20dB | 倒扣 |
|---|---|---:|---:|---:|
| 48 PRB | N3/C0 | 3.985e-3 | 2.203e-3 | +2.57 dB |
| 48 PRB | N6/C0 | 3.907e-3 | 2.203e-3 | +2.49 dB |
| 48 PRB | N8/C1 | 4.179e-3 | 1.973e-3 | +3.26 dB |
| 36 PRB | N23/C4 | 4.350e-3 | 2.870e-3 | +1.81 dB |
| 24 PRB | **N8/C6** | **2.726e-2** | 4.280e-3 | **+8.04 dB** |

- 全部 12 个候选点 NMSE 倒扣,**没有一个候选在 estimated-CSI BLER 上稳定超过对应 CDD 参考**;最接近的只是在 400-trial 分辨率内打平。
- 24PRB N8 尤其病态:NMSE 从 12dB 的 4.80e-2 只缓降到 24dB 的 2.31e-2,对应约 **11–15% BLER floor**。

**为什么**：N-series 的 slope 切换制造了**非平稳等效协方差**;CDD 的等效协方差是平稳的 shifted-PDP（Toeplitz）,更容易被稀疏导频插值。稀疏 comb（spacing=24）下,非平稳协方差的目标 RE 相对采样导频集合可预测性差,即使用 matched covariance 也补不回来。这是分集增益与 CE 处理增益之间折中的真实体现,**不是实现 bug**。

**图**：`figures/experiment21_mmse_nmse_48prb.svg`、`figures/experiment21_mmse_nmse_24prb.svg`。
**出处**：英文综合报告 §5.3、§5.4、§6#5、§6#7;`research/result-021.md`。实验编号 E21。

---

<a id="b6"></a>
## B6. 矩阵指标是筛选工具,不是充分的最终目标

**主张**：`log det(V^H V)` 与 `B_0.5`(相干带宽阈值),单独或合起来,都不足以决定候选 `V` 的真实链路性能。

**证据**：
- B4/B5 直接反例:相似的 Gram 行列式 + 相似 `B_0.5` 的候选,estimated-CSI NMSE/BLER 可以差 1.8~8 dB。
- 当把 CDD 参考从"`B_0.5` 最接近"放宽到包含大 delay 候选时,`step=64 samples` 的 C7 在 ideal-CSI 1% BLER 尾部最好（48/36/24 PRB 分别 16.00/16.25/16.20 dB）——**仅用矩阵指标匹配 CDD 参考会系统性低估 CDD 的真实最优**。这也是为什么 plan-022 强制把大 delay CDD 纳入参考基线。

**为什么**：`B_0.5` 是对角平均阈值,RMMSE 依赖的是完整非平稳协方差 + 稀疏导频梳的联合可预测性;矩阵体积/相干宽度都看不到这层。

**结论落地**：任何候选必须过两关——matched full-cov RMMSE NMSE + estimated-CSI 编码 BLER。这直接催生了 plan-022 的导频感知联合目标 `J(V)=J_div - mu·L_CE`。

**出处**：英文综合报告 §5.2、§5.4、§6#5。

---

<a id="b7"></a>
## B7. CDD delay 增大对 BLER 非单调

**主张**：增大 CDD delay 不必然改善 BLER;存在与 DMRS spacing 交互决定的非平凡最优点。

**证据**：等效信道相干带宽随 delay 增大而单调压缩（实验 8 相干带宽表:ds30 原始 branch `B_0.5`=311 sc,CDD d=32 后压到 43 sc,d=64 后 `B_0.9` 仅 10 sc）。分集变好但可估计性变差,两者反向,故 BLER 关于 delay 非单调。d=512 之所以不是必然最优,正是这个折中。

**图**：`figures/experiment8_covariance_correlation_with_bc.png`。
**出处**：实验记录 §18。实验编号 E18。

---

<a id="b8"></a>
## B8. 接收机 LLR 未注入信道估计误差自噪声（潜在修正点,未验证）

**主张**：当前解调 LLR 的噪声方差只含 AWGN,没有把信道估计误差方差加进去。

**可能影响**：大 NMSE 候选（典型是 24PRB N8）出现 BLER floor 而非平滑退化,部分可能源于此——接收机"过度自信"于错误的信道估计。加入 CE-error-aware 噪声缩放后,floor 结论是否改变尚未验证。

**状态**：候选修正方向,未排期,不阻塞 plan-022。若 plan-022 出现"NMSE 差但 BLER floor 反常"的点,应优先回来验证本条。

**出处**：英文综合报告 §5.4。
