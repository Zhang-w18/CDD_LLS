# CDD 时延选择规则

## 1. 文档定位

本文整理项目截至 `result-028` 形成的四类结构化 CDD 时延生成规则：严格 Sidon、RMS 展宽感知等差、$99\%$ 有效支撑等差和厚 Sidon。目标是给出一套可复现的“系统参数到候选时延”推导过程，便于在改变带宽、DMRS comb、发射天线数或信道时延展宽后重新生成候选。

本文是 `DESIGN.md` 的操作性说明，不替代其中的理论定义和最终评价准则。若两者冲突，以 `DESIGN.md` 为准。实验是否支持某一候选，应查阅对应的 `research/result-NNN-text.md`；不能仅凭本文的几何条件宣布性能胜出。

`docs/FINDINGS.md` 记录实验 1--21 的历史证据，不适合继续承载当前设计规范，因此本规则单独成文。

## 2. 输入参数与公共派生量

### 2.1 输入参数

| 输入 | 含义 |
|---|---|
| RB 数 | 当前有效带宽包含的连续资源块数 |
| $N_t$ | 发射分支或发射天线数 |
| $\Delta f$ | 子载波间隔，单位 Hz |
| $S_f$ | DMRS 频域 comb 间隔，单位为子载波 |
| $T_{\rm RMS}$ | 物理信道 RMS delay spread |
| $\epsilon$ | 有效支撑允许忽略的 PDP 能量比例；当前主值为 $0.01$ |
| $T_\epsilon$ | 覆盖 $1-\epsilon$ PDP 能量的最短连续时延区间宽度 |
| $T_{\rm margin}$ | 厚 Sidon 二元和保护的额外余量 |
| $T_{\rm sync}$ | DMRS 折叠保护包含的同步或定时误差余量 |

### 2.2 有效带宽和时延栅格

对于当前 NR 资源映射，每个 RB 含 12 个子载波，因此先计算

$$
K=12\times\text{RB 数}.
$$

CDD 频域相位定义为

$$
V_{k,n}
=
\exp\left(-j2\pi k\Delta f\tau_n\right),
\qquad
k=0,\ldots,K-1.
$$

采用有效带宽 DFT 栅格时，

$$
\tau_n
=
\frac{j_n}{K\Delta f},
\qquad
j_n\in\mathbb Z_K.
$$

一个整数 delay index 对应的物理时延为

$$
\frac{1}{K\Delta f}.
$$

人工时延的循环周期为

$$
\frac{1}{\Delta f}.
$$

因此 $j_n$ 与 $j_n+K$ 产生相同的频域相位。候选应在模 $K$ 意义下规范化，并消除公共循环移位和普通天线排列产生的重复。

若实现使用 $N_{\rm FFT}$ 点整数 sample 循环移位，则同一物理时延对应

$$
d_n
=
j_n\frac{N_{\rm FFT}}{K}.
$$

只有 $d_n$ 为整数时才能直接用整数 sample buffer 精确实现；否则需要逐子载波相位旋转或分数时延实现。不得把 `/K` 相位模型与 `/N_{\rm FFT}` 相位模型混为一谈。

### 2.3 DMRS 折叠周期

本文的简化推导要求 active subcarrier 和 DMRS comb 对齐，且 $K/S_f$ 为整数。每个 DMRS symbol 的唯一频率导频数为

$$
N_p
=
\frac{K}{S_f}.
$$

无歧义时延周期为

$$
\tau_{\rm alias}
=
\frac{1}{S_f\Delta f}
=
\frac{N_p}{K\Delta f}.
$$

在导频子载波 $k=pS_f$ 上，

$$
V_{pS_f,n}
=
\exp\left(-j2\pi p\frac{j_n}{N_p}\right).
$$

所以整数栅格候选在 DMRS 上只由

$$
r_n
=
j_n\bmod N_p
$$

决定。纯平坦信道下，$r_n$ 互异可保证各 CDD 分支在导频域可辨识；均匀、最大圆周分离的 $r_n$ 通常具有更好的条件数。物理展宽信道还必须计算实际 $\mathbf R_{PP}$、$\mathbf R_{DP}$、CE NMSE 和零噪声误差地板。

如果实际 DMRS 位置不是严格等间隔或 $K/S_f$ 不是整数，应直接使用真实 pilot index 构造 $\mathbf V_P$，不能继续套用单一 $N_p$ 折叠圆公式。

## 3. 严格 Sidon

### 3.1 设计目标

选择

$$
\mathcal J
=
\{j_0,\ldots,j_{N_t-1}\}
\subset\mathbb Z_K.
$$

对全部无序天线对

$$
0\le a\le c<N_t,
$$

计算模 $K$ 二元和

$$
(j_a+j_c)\bmod K.
$$

严格 Sidon 要求 $N_t(N_t+1)/2$ 个无序二元和全部互不相同。等价地，

$$
j_a+j_c
\equiv
j_b+j_d
\pmod K
$$

只能在

$$
\{a,c\}=\{b,d\}
$$

时成立。

该准则消除平坦信道下的非平凡四阶加性共振。它直接依赖 $K$ 和 $N_t$，不直接依赖 $T_{\rm RMS}$、$T_\epsilon$ 或 $S_f$。

### 3.2 必要装填条件

模 $K$ 圆周必须容纳所有无序二元和，因此必要条件为

$$
K
\ge
\frac{N_t(N_t+1)}{2}.
$$

该条件仅是必要条件，不保证一定存在满足其他 DMRS 和实现约束的候选。

### 3.3 可复现构造

1. 利用公共循环移位等价性固定 $j_0=0$。
2. 规定确定的候选遍历顺序，例如从小到大尝试 $j_1,j_2,\ldots$。
3. 加入新 index $x$ 前，计算 $(x+j_a)\bmod K$ 和 $(x+x)\bmod K$。
4. 只有新二元和彼此不冲突、且不与已有二元和冲突时才接受 $x$。
5. 达到 $N_t$ 个元素后，按公共移位和天线排列规范化去重。
6. 多个严格 Sidon 均可行时，必须预先规定代表选择规则，例如最小最大时延、最大 DMRS fold 间距或固定字典序。

当前 $K=576$、$N_t=8$ 的固定历史参考为

$$
\mathbf j_{\rm S0}
=
[0,1,3,7,12,20,30,65].
$$

其 36 个无序二元和均不相同。

### 3.4 DMRS 可用性检查

严格 Sidon 本身不保证信道估计可用。生成后至少检查：

$$
N_p\ge N_t,
$$

以及

$$
j_a\bmod N_p
\ne
j_b\bmod N_p,
\qquad a\ne b.
$$

还应报告 fold 最小圆周间距、pilot rank、condition number 和实际协方差下的 CE NMSE。满足严格 Sidon 只说明平坦模型中的四阶几何成立，不保证物理展宽信道中的 outage 或 BLER 改善。

## 4. RMS 展宽感知等差

### 4.1 生成公式

直接把 RMS delay spread 作为相邻人工时延的公差：

$$
\tau_n^{\rm RMS}
=
nT_{\rm RMS},
\qquad
n=0,\ldots,N_t-1.
$$

对应的有效带宽 DFT 坐标为

$$
j_n^{\rm RMS}
=
K\Delta f\tau_n^{\rm RMS}
=
nK\Delta fT_{\rm RMS}.
$$

这里的 $j_n^{\rm RMS}$ 允许为实数。项目中的该基线直接按物理时延构造相位，不把它舍入到整数 DFT 栅格。

### 4.2 参数依赖

- $N_t$ 和 $T_{\rm RMS}$ 唯一决定物理时延数组；
- $K$ 和 $\Delta f$ 决定物理时延到 DFT 坐标的换算；
- $S_f$ 不参与时延生成，但决定 DMRS 折叠和信道估计可用性。

例如，48 PRB、$K=576$、$\Delta f=30$ kHz、$N_t=8$、$T_{\rm RMS}=100$ ns 时，

$$
\boldsymbol\tau^{\rm RMS}
=
[0,100,200,300,400,500,600,700]\ {\rm ns},
$$

$$
\mathbf j^{\rm RMS}
=
[0,1.728,3.456,5.184,6.912,8.64,10.368,12.096].
$$

### 4.3 生成后检查

应排除或明确标注满足以下任一关系的分支：

$$
\tau_a-\tau_b
=
m\frac1{\Delta f},
$$

因为它们的全带频域相位完全等价；还应检查

$$
\tau_a-\tau_b
\approx
m\tau_{\rm alias},
$$

因为它们可能在 DMRS 上重合或接近重合。对于连续坐标，必须直接构造实际 $\mathbf V_P$，不能只使用整数 residue 规则。

## 5. $99\%$ 有效支撑等差

### 5.1 有效支撑计算

设归一化连续 PDP 为 $p(\tau)$，当前取 $\epsilon=0.01$。有效支撑定义为

$$
T_\epsilon
=
\min_{a\le b}
\left\{
b-a:
\int_a^bp(\tau)\,d\tau
\ge
1-\epsilon
\right\}.
$$

对于按时延递增排列的离散抽头 $(\delta_i,P_i)$，先令 $\sum_iP_i=1$，再计算

$$
T_\epsilon
=
\min_{i\le l}
\left\{
\delta_l-\delta_i:
\sum_{m=i}^{l}P_m
\ge
1-\epsilon
\right\}.
$$

这里取覆盖目标能量的最短连续区间，而不是默认使用首末抽头总跨度。

### 5.2 生成公式

把有效支撑作为相邻人工时延的公差：

$$
\tau_n^\epsilon
=
nT_\epsilon,
\qquad
n=0,\ldots,N_t-1.
$$

对应连续 DFT 坐标为

$$
j_n^\epsilon
=
nK\Delta fT_\epsilon.
$$

同 RMS 等差一样，不进行整数栅格舍入。

例如，TDL-A 100 ns 在当前 PDP 离散化下有

$$
T_{0.01}=479.66\ {\rm ns},
$$

因此

$$
\boldsymbol\tau^epsilon
=
[0,479.66,959.32,1438.98,1918.64,2398.30,2877.96,3357.62]\ {\rm ns},
$$

$$
\mathbf j^epsilon
=
[0,8.2885248,16.5770496,24.8655744,33.1540992,
41.442624,49.7311488,58.0196736].
$$

### 5.3 参数依赖和解释

该方案依赖 $N_t$、$\epsilon$、$T_\epsilon$、$K$ 和 $\Delta f$。$S_f$ 只在生成后的 CE 检查中进入。

该构造的解释是让相邻人工移位按主要 PDP 支撑宽度分离。它是 T1 知识条件下的确定性基线，不是“所有物理 PDP 副本完全不重叠”的保证，也不以完整 PDP 协方差进行优化。

## 6. 厚 Sidon

### 6.1 物理保护条件

严格 Sidon 只排除二元和完全相等。物理信道具有有效宽度后，每个人工二元和会扩展成近似宽度 $2T_\epsilon$ 的簇，因此对不同无序天线对要求

$$
d_{\rm circ}(s_{ac},s_{bd})
>
2T_\epsilon+T_{\rm margin},
\qquad
\{a,c\}\ne\{b,d\},
$$

其中

$$
s_{ac}
=
\tau_a+\tau_c.
$$

CE 侧独立要求

$$
d_{\rm fold}(\tau_a,\tau_b)
>
T_\epsilon+T_{\rm sync},
\qquad a\ne b.
$$

$d_{\rm circ}$ 是在人工时延周期 $1/\Delta f$ 上的圆周距离；$d_{\rm fold}$ 是在 $\tau_{\rm alias}$ 上的圆周距离。这两个条件分别控制近四阶共振和 DMRS 可辨识性，不能互相替代。

### 6.2 整数栅格要求

一个整数 delay index 对应 $1/(K\Delta f)$ 秒。严格大于保护宽度所需的最小整数距离为

$$
g_\Sigma
=
\left\lfloor
(2T_\epsilon+T_{\rm margin})K\Delta f
\right\rfloor+1,
$$

$$
g_{\rm fold}
=
\left\lfloor
(T_\epsilon+T_{\rm sync})K\Delta f
\right\rfloor+1.
$$

因此硬厚 Sidon 必须同时满足：

1. 任意两个不同无序二元和在模 $K$ 圆周上的最小距离不小于 $g_\Sigma$；
2. 任意两个 $j_n\bmod N_p$ 在模 $N_p$ 圆周上的最小距离不小于 $g_{\rm fold}$。

### 6.3 residue/lift 构造

将整数时延写为

$$
j_n
=
r_{\pi(n)}+N_pq_n,
$$

其中 $r_n$ 是模 $N_p$ 的导频折叠余数，$q_n$ 是整数 lift，$\pi$ 是 residue 与 lift 的配对。

对规范化到 $0\le j_n<K$ 的候选，lift 可取

$$
q_n\in\{0,\ldots,S_f-1\}.
$$

构造步骤为：

1. 在模 $N_p$ 圆上选择 $N_t$ 个 residue，使 fold 最小圆周距离不小于 $g_{\rm fold}$；
2. 为每个 residue 分配 lift，形成 $j_n=r_{\pi(n)}+N_pq_n$；
3. 规范化并去除公共移位、普通天线排列的等价候选；
4. 计算全部 $N_t(N_t+1)/2$ 个模 $K$ 无序二元和；
5. 只保留 pair-sum 最小圆周距离不小于 $g_\Sigma$ 的候选；
6. 对多个硬可行候选，按预先规定的最大 pair gap、最大 fold gap或平衡几何规则冻结代表。

在独立同分布发射分支下，普通天线排列等价；若各天线 PDP 或空间相关不同，必须保留并搜索配对 $\pi$。

### 6.4 快速不可行判据

模 $K$ 圆周需要放置 $N_t(N_t+1)/2$ 个二元和，模 $N_p$ 圆周需要放置 $N_t$ 个 fold 中心，因此必要条件为

$$
g_\Sigma
\le
\left\lfloor
\frac{K}{N_t(N_t+1)/2}
\right\rfloor,
$$

$$
g_{\rm fold}
\le
\left\lfloor
\frac{N_p}{N_t}
\right\rfloor.
$$

任一条件不满足时，硬厚 Sidon 必然不存在；两项都满足仍不保证组合约束一定可同时实现。

例如，$K=576$、$N_t=8$ 时，pair-sum 装填上界为

$$
\left\lfloor\frac{576}{36}\right\rfloor=16.
$$

TDL-A 100 ns 的 $T_{0.01}=479.66$ ns 给出

$$
g_\Sigma
=
\left\lfloor
2(479.66\ {\rm ns})(576)(30\ {\rm kHz})
\right\rfloor+1
=17.
$$

因为 $17>16$，该带宽和天线数下，无论提高到哪一档当前规则内的 DMRS comb，硬厚 Sidon 都受 pair-sum 装填限制而不可行。

### 6.5 硬约束不可行时

不得把软候选标记为硬厚 Sidon。应计算每个候选的实际 pair-sum 最小距离和 fold 最小距离，并保留二者构成的 Pareto 前沿，包括：

- 最大 pair-gap 候选；
- 最大 fold-gap 候选；
- 最大化两类归一化间距最小值的平衡候选；
- 其他未被同时支配的软几何候选。

这些候选属于 `GEO_T1_CTRL` 一类，只使用 $T_\epsilon$ 和系统几何生成，不得在生成阶段读取完整 PDP 功率或协方差。

## 7. 四类规则的参数依赖

| 方案 | 直接生成参数 | DMRS 密度的作用 | 输出是否唯一 |
|---|---|---|---|
| 严格 Sidon | $K,N_t$ | 生成后检查 fold、pilot rank 和 CE | 否，通常有多个可行集合 |
| RMS 展宽感知等差 | $N_t,T_{\rm RMS}$ | 生成后检查 fold、pilot rank 和 CE | 是 |
| $99\%$ 有效支撑等差 | $N_t,\epsilon,T_\epsilon$ | 生成后检查 fold、pilot rank 和 CE | 是 |
| 厚 Sidon | $K,N_t,N_p,\Delta f,T_\epsilon,T_{\rm margin},T_{\rm sync}$ | 通过 $N_p=K/S_f$ 和 fold 保护直接进入构造 | 否，通常需要组合搜索 |

所有物理时延到相位坐标的换算还依赖 $K\Delta f$。FFT sample 实现额外依赖 $N_{\rm FFT}$。

## 8. 给定新参数时的执行清单

给定 RB 数、$N_t$、$\Delta f$、$S_f$、$T_{\rm RMS}$、$T_\epsilon$ 和保护余量后，按以下顺序执行：

1. 计算 $K=12\times\text{RB 数}$、$N_p=K/S_f$、一个 DFT delay index 的秒值和 $\tau_{\rm alias}$。
2. 直接生成 RMS 等差和 $T_\epsilon$ 等差的连续时延数组。
3. 在 $\mathbb Z_K$ 上通过二元和唯一性生成严格 Sidon 候选池。
4. 计算 $g_\Sigma$、$g_{\rm fold}$ 和两项装填上界。
5. 硬厚条件可行时搜索硬候选；不可行时生成并明确标记软几何 Pareto 候选。
6. 保存每个候选的 `delay_ns`、`delay_grid_coordinates`；仅对整数栅格候选保存 `delay_indices`、residue 和 lift。
7. 使用真实 DMRS index 构造 $\mathbf V_P$，计算 rank、condition number、实际 $\mathbf R_{PP}$、$\mathbf R_{DP}$、matched 或指定失配接收机下的 CE NMSE和零噪声误差地板。
8. 使用实际物理协方差计算相关代理或 Monte Carlo ideal-CSI outage，仅作为筛选。
9. 对入围候选运行实际信道估计、LLR、LDPC 编译码和译码器，最终报告达到 $10\%$ 和样本允许时 $1\%$ estimated-CSI BLER 所需 SNR、错误块计数和置信区间。

## 9. 适用边界

- 本文公式基于 active-band `/K` 相位定义；改变为物理 FFT-bin `/N_{\rm FFT}` 定义时必须重新推导和生成候选。
- 人工 delay 是数字循环移位或等价频域相位，不作为真实传播时延计入 CP；实际硬件仍需检查相位量化和分数时延实现。
- 严格或厚 Sidon 是候选生成条件，不是最终性能定理。
- RMS 和 $T_\epsilon$ 等差是知识受限的确定性基线，不代表在完整 PDP 已知时的最优等差解。
- DMRS 满秩不保证 data RE 可预测性；最终 CE 必须用实际 $\mathbf R_{PP}$ 和 $\mathbf R_{DP}$ 计算。
- 跨带宽、跨 DMRS comb 或跨 PDP 的候选必须重新生成或至少重新完成全部可行性与 CE 检查。

## 10. 依据

- 理论定义与最终评价：`DESIGN.md` §3.4、§3.8、§4.5、§5.3--5.4、§6。
- 严格和厚 Sidon 的可复现构造：`research/plan-027-有效矩-Sidon-DMRS.md` §5--6。
- 厚 Sidon 已验证候选与不可行区域：`research/result-026-text.md` §4、`research/result-027-有效矩-Sidon-DMRS-text.md` §3。
- A30/A100/A300 不同展宽和 DMRS comb 的最终链路表现：`research/result-027-有效矩-Sidon-DMRS-text.md` §6、`research/result-028-comb6三类CSI-TDL-A300ns-text.md` §3。
