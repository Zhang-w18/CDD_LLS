# V 设计理论框架（无注释版）

## 1. 用途与研究问题

本文件用于制定实验计划。推导解释和历史修订见 `docs/design/DESIGN_ANNOTATED.md`。制定 plan 时同时读取 `GOALS.md`、`KNOWLEDGE.md` 和最新无图版 result。

频域恒模相位预编码矩阵定义为

$$
\mathbf V\in\mathbb C^{K\times N_t},\qquad V_{k,n}=e^{j\phi_{k,n}},
$$

其中 `K` 是有效子载波数，`N_t§ 是发射分支数，`k` 和 `n` 是索引，`\phi_{k,n}` 是相位。带内平坦模型为

$$
\mathbf g=\mathbf V\mathbf h,\qquad
\mathbf h\sim\mathcal{CN}(\mathbf0,\mathbf I_{N_t}),
$$

其中 `\mathbf h` 是底层分支信道，`\mathbf g` 是等效频域信道。目标是在相同 DMRS 开销下比较一般 `\mathbf V` 与 CDD 的分集、信道估计误差和 estimated-CSI BLER，并确定适用条件。

物理信道随子载波变化时，必须改用

$$
g_k=\sum_{n=0}^{N_t-1}V_{k,n}h_{k,n},
$$

不得继续把 `\mathbf g=\mathbf V\mathbf h` 当作精确模型。

## 2. 默认系统条件

| 参数 | 默认值 |
|---|---|
| 发射与接收 | 8 Tx / 1 Rx，单层 |
| 子载波间隔 `\Delta f` | 30 kHz |
| `K` | 576 / 432 / 288，对应 48 / 36 / 24 PRB |
| DMRS comb `S_f` | 24 子载波 |
| 导频数 `N_p=K/S_f` | 24 / 18 / 12 |
| DMRS 符号 | 2；静态信道下只降噪 |
| MCS | 16QAM，码率 553/1024 |
| 当前未完成验证 | Sidon 在 5–100 ns 频率选择性信道下的稳健性 |

plan 改变任一默认条件时必须明确列出。

## 3. 分集评价

$$
I(\mathbf h)=\frac1K\sum_{k=0}^{K-1}I_{\rm QAM}(\mathrm{snr}|g_k|^2),
\qquad
P_{\rm out}(\mathrm{snr},R)=\Pr[I(\mathbf h)<R],
$$

其中 `I_{\rm QAM}` 是 QAM/BICM 每 RE 互信息函数，`R` 是目标谱效率。至少报告 10% 和 1% outage 所需 SNR、样本、种子、插值和置信区间。

全带列正交条件与归一化相关系数为

$$
\mathbf V^H\mathbf V=K\mathbf I_{N_t},\qquad
\rho_{kl}=\frac{(\mathbf V\mathbf V^H)_{kl}}{N_t}.
$$

`\sum_{k\ne l}|\rho_{kl}|^4` 可诊断高阶相关结构。Gram、log-det、条件数和相干带宽只用于粗筛，最终判定使用 outage 和 BLER。

## 4. 信道估计评价

### V-aware

导频观测为

$$
\mathbf y_P=\mathbf V_P\mathbf h+\mathbf n,\qquad
\mathbf V_P^H\mathbf V_P=N_p\mathbf I_{N_t},
$$

其中 `P` 是导频集合，`\mathbf V_P` 是导频行子矩阵。必须报告秩、条件数和最小奇异值。

一般协方差下，数据集合 `D` 的 matched LMMSE NMSE 为

$$
\mathrm{NMSE}=
\frac{\operatorname{tr}\!\left[
\mathbf R_{DD}-\mathbf R_{DP}
(\mathbf R_{PP}+\sigma_{LS}^2\mathbf I)^{-1}
\mathbf R_{PD}\right]}
{\operatorname{tr}(\mathbf R_{DD})}.
$$

`\mathbf R` 下标表示协方差子块，`\sigma_{LS}^2` 是导频 LS 噪声方差。

### V-agnostic

固定通用先验、窗口、PDP、边界处理和可用导频必须写入 plan，且不得按候选单独优化。若接收机获得分段边界、候选家族或相位参数，该方案不再是完全 V-agnostic，必须单列。

## 5. 导频混叠与 CDD

$$
\tau_{\rm alias}=\frac1{S_f\Delta f}=1.389\ \mu\mathrm{s},\qquad
\tau_n=\frac{j_n}{K\Delta f}.
$$

`j_n` 是 DFT 栅格 delay 索引。设计必须检查：

1. 全带列正交；
2. `j_n\bmod N_p` 互异；
3. 物理 PDP 展宽后的可辨识性；
4. 接收机协方差与实际 `\mathbf V`、PDP 一致。

近平坦、V-aware matched LMMSE、`N_p\ge N_t` 且条件良好时，DFT 栅格 CDD 可同时达到全带和导频域列正交。一般 `\mathbf V` 的改进空间主要来自高阶相关结构。

## 6. 已验证结论

- result-021：N-series 在矩阵与 ideal-CSI 上可优于 CDD，但 12/12 个 estimated-CSI NMSE 点均劣化。
- result-023：固定 4RB V-agnostic 接收机下，232 个 CC/CN 挑战者均未达到透明基线支配标准。
- result-024 E1：已知分段边界仅改善部分 B3，B4 全部恶化；48 PRB 无通过点，24 PRB 只有一个半透明诊断点。

因此，不把分段边界对齐作为通用修正。

整数 delay 集 `A=\{j_0,\ldots,j_{N_t-1}\}` 若满足

$$
a+b=c+d\Longrightarrow\{a,b\}=\{c,d\},
$$

称为 Sidon 集。result-023 中 `[0,1,3,7,12,20,30,65]` 相对等差 CDD 的 10%/1% outage 改善为 0.210/0.368 dB。

result-024 E2 在 48 PRB 平坦模型、相同 DMRS、双方 matched LMMSE 下验证：

- 10% BLER 改善 0.33 dB，保守 95% 区间 [0.20, 0.45] dB；
- 1% BLER 改善 0.85 dB，保守 95% 区间 [0.63, 1.07] dB；
- matched CE NMSE 差最大绝对值 0.06 dB。

该结论尚未推广到频率选择性 TDL、定时误差、协方差失配、其他带宽或移动信道。

## 7. 后续实验顺序

1. 5–100 ns PDP 的导频矩阵条件、matched NMSE 和 BLER；
2. PDP 或 `\mathbf V` 协方差失配；
3. 定时误差；
4. 24 / 36 / 48 PRB；
5. 更稀疏导频、更多发射分支或移动信道。

每项使用相同 DMRS 开销、随机样本和 QC 基线。先做条件数与闭式 NMSE，再决定是否运行 LDPC BLER。

## 8. 规范符号

| 符号 | 含义 |
|---|---|
| `N_t` / `K` | 发射分支数 / 有效子载波数 |
| `\Delta f` | 子载波间隔 |
| `\mathbf V` | 恒模频域相位矩阵 |
| `\mathbf h` / `\mathbf g` | 底层分支信道 / 等效频域信道 |
| `S_f` / `N_p` | DMRS comb 间隔 / 导频数 |
| `P,D` | 导频 / 数据 RE 集合 |
| `\tau_n` / `j_n` | CDD delay / DFT 栅格索引 |
| `\tau_{\rm alias}` | 无混叠时延周期 |
| `\rho_{kl}` | 归一化相关系数 |
| `P_{\rm out}` | outage 概率 |

新增符号先加入本表。
