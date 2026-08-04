# V 设计理论框架（无注释版）

## 1. 文档用途

本文件定义项目当前使用的总体目标、问题 framing、理论分析框架和研究地图。推导解释与历史修订见 `docs/design/DESIGN_ANNOTATED.md`；实验事实、已排除方向和开放问题分别记录在 `KNOWLEDGE.md`、`GOALS.md` 和对应的 `research/result-NNN-text.md`。

本文件不记录阶段性实验结论，也不规定实验执行顺序。制定 plan 时，应从本文件选定明确的信道模型、发射端知识、接收端知识、候选设计空间和评价指标。

## 2. 总体目标

频域恒模相位预编码矩阵定义为

$$
\mathbf V\in\mathbb C^{K\times N_t},
\qquad
V_{k,n}=e^{j\phi_{k,n}},
$$

其中 $K$ 是有效子载波数，$N_t$ 是发射分支数，$k$ 和 $n$ 分别是子载波和发射分支索引。

项目总体目标是：

> 设计具有明确结构、可解释机理和可实现参数化的 $\mathbf V$，在相同 DMRS 开销和公平接收机条件下，其分集、信道估计和 estimated-CSI BLER 的 Pareto 前沿优于一般性 CDD。

“优于一般性 CDD”不是只与单个等差时延基线比较，而是与经过合理优化的 CDD 设计空间比较，包括：

- DFT 栅格和非栅格时延；
- 等差、非等差和折叠位置受控的时延集合；
- 在给定物理信道知识下优化的 CDD；
- 满足相同功率、DMRS、循环移位定义和接收机知识约束的 CDD。

候选 $\mathbf V$ 必须具有可复现的生成规则。完全随机、仅依赖单次搜索种子的无结构矩阵只能作为对照或搜索初始化，不作为最终设计。

## 3. 问题 framing：理论分析框架

### 3.1 信道模型

#### 纯平坦分支信道

每个发射分支在带内只有一个复系数：

$$
\mathbf g=\mathbf V\mathbf h,
\qquad
\mathbf h\sim\mathcal{CN}(\mathbf0,\mathbf I_{N_t}),
$$

其中 $\mathbf h\in\mathbb C^{N_t}$ 是底层分支信道，$\mathbf g\in\mathbb C^K$ 是等效频域信道。

该模型用于隔离 $\mathbf V$ 本身的相位轨迹、相关结构和自平均能力。

#### 物理展宽信道

物理信道随子载波变化时，

$$
g_k=\sum_{n=0}^{N_t-1}V_{k,n}h_{k,n}.
$$

若各发射分支独立、具有相同的归一化物理频率协方差 $\mathbf R_{\rm phy}$，则

$$
\mathbf R_g
=
\mathbf R_{\rm phy}\odot
\frac{\mathbf V\mathbf V^H}{N_t},
$$

其中 $\odot$ 为 Hadamard 积，$\mathbf R_g$ 是归一化等效信道协方差。存在天线相关、分支相关或分支 PDP 不同时，应从完整的分支间协方差构造 $\mathbf R_g$，不得继续使用上述简化式。

### 3.2 分集目标函数

对一次信道实现，块平均互信息定义为

$$
I(\mathbf h)
=
\frac1K\sum_{k=0}^{K-1}
I_{\rm QAM}\!\left(\mathrm{snr}|g_k|^2\right),
$$

其中 $I_{\rm QAM}$ 是实际 QAM/BICM 每 RE 互信息函数。目标谱效率为 $R$ 时，outage 定义为

$$
P_{\rm out}(\mathrm{snr},R)
=
\Pr[I(\mathbf h)<R].
$$

分集侧的最终判据是指定 outage 概率所需 SNR，通常至少包含 10% 和 1% 两个尾部。Gram、log-det、相关矩和条件数只能作为筛选或机理诊断，不能替代 outage 或链路 BLER。

### 3.3 纯平坦信道下的模式分解

记 $\mathbf V$ 第 $k$ 行的相位向量为

$$
\boldsymbol\phi_k
=
(\phi_{k,0},\ldots,\phi_{k,N_t-1}),
$$

并定义相位环面函数

$$
f_{\mathbf h}(\boldsymbol\theta)
=
I_{\rm QAM}\!\left(
\mathrm{snr}
\left|
\sum_{n=0}^{N_t-1}h_ne^{j\theta_n}
\right|^2
\right).
$$

则

$$
I(\mathbf h)
=
\frac1K\sum_{k=0}^{K-1}
f_{\mathbf h}(\boldsymbol\phi_k).
$$

对 $f_{\mathbf h}$ 作 Fourier 展开：

$$
f_{\mathbf h}(\boldsymbol\theta)
=
\sum_{\mathbf m\in\mathbb Z^{N_t}}
c_{\mathbf m}(\mathbf h)
e^{j\mathbf m^T\boldsymbol\theta}.
$$

相位环面均值为零模系数

$$
I_\infty(\mathbf h)=c_{\mathbf0}(\mathbf h),
$$

而有限 $K$ 的偏差为

$$
I(\mathbf h)-I_\infty(\mathbf h)
=
\sum_{\mathbf m\ne\mathbf0}
c_{\mathbf m}(\mathbf h)W_K(\mathbf m),
$$

其中

$$
W_K(\mathbf m)
=
\frac1K\sum_{k=0}^{K-1}
e^{j\mathbf m^T\boldsymbol\phi_k}
$$

是由设计决定的 Weyl 和。

由于 $f_{\mathbf h}$ 对所有分支同时增加公共相位不变，只有满足

$$
\sum_n m_n=0
$$

的偶数阶模式可能非零。低阶模式通常比高阶模式具有更大的 Fourier 系数，因此设计应优先压低低阶非零模式。

二阶模式为 $\mathbf m=\mathbf e_a-\mathbf e_b$，对应

$$
W_K(\mathbf e_a-\mathbf e_b)
=
\frac1K(\mathbf V^H\mathbf V)_{b,a}.
$$

因此

$$
\mathbf V^H\mathbf V=K\mathbf I_{N_t}
$$

等价于消除所有二阶非平凡模式。不同全带列正交设计之间的区别来自四阶及更高阶模式。

### 3.4 DFT 栅格 CDD 与 Sidon 指导

CDD 的线性相位形式为

$$
V_{k,n}=e^{-j2\pi k\Delta f\tau_n}.
$$

取 DFT 栅格时延

$$
\tau_n=\frac{j_n}{K\Delta f},
\qquad j_n\in\mathbb Z,
$$

则

$$
W_K(\mathbf m)
=
\mathbb 1\!\left[
\sum_nm_nj_n\equiv0\pmod K
\right].
$$

不同 $j_n\bmod K$ 消除二阶共振。四阶共振满足

$$
j_a+j_c\equiv j_b+j_d\pmod K.
$$

若除 $\{a,c\}=\{b,d\}$ 外不存在其他解，则人工 delay 索引构成 Sidon 型集合。对纯平坦信道，四阶相关代理满足

$$
\sum_{k,l}|\rho_{kl}|^4
=
\frac{K^2}{N_t^4}
\#\left\{
(a,b,c,d):
j_a+j_c\equiv j_b+j_d\pmod K
\right\},
$$

其中

$$
\rho_{kl}
=
\frac{(\mathbf V\mathbf V^H)_{kl}}{N_t}.
$$

因此，在“纯平坦分支信道 + DFT 栅格 CDD + 接收端知道 $\mathbf V$”的假设下，无非平凡四元关系是压低四阶模式的直接设计指导。该等价关系不得直接推广为物理展宽信道下的全局设计定理。

### 3.5 物理展宽信道下的相关结构

物理展宽信道下，实际归一化相关系数为

$$
\rho_{\rm eff}(k,l)
=
R_{{\rm phy},kl}\rho_{kl}.
$$

对应的四阶相关代理为

$$
M_4(\mathbf V;\mathbf R_{\rm phy})
=
\sum_{k\ne l}
\left|
R_{{\rm phy},kl}\rho_{kl}
\right|^4.
$$

对 CDD delay 集，将其展开后可写成

$$
M_4
=
\frac1{N_t^4}
\sum_{a,b,c,d}
\mathcal K_{\rm PDP}
(j_a+j_c-j_b-j_d),
$$

其中 $\mathcal K_{\rm PDP}$ 是由物理频率协方差决定的加权核。纯平坦信道下，该核退化为模 $K$ 指示函数；物理展宽下，精确四元关系、近四元关系和不同频差位置均可能产生不同权重。

若物理 PDP 的抽头时延为 $\delta_p$，复合路径时延为

$$
\tau_{n,p}=\tau_n+\delta_p.
$$

即使人工中心时延严格满足 Sidon 条件，复合路径集合仍天然包含

$$
(\tau_a+\delta_p)+(\tau_c+\delta_q)
=
(\tau_a+\delta_q)+(\tau_c+\delta_p)
$$

等关系。因此，物理展宽信道下不应以“整个复合路径集合严格无四元关系”为目标，而应直接优化由 PDP 功率和频率协方差加权的相关结构。

### 3.6 统一的分集目标与块方差代理

在“各发射分支独立、具有相同归一化 PDP”的简化假设下，纯平坦和物理展宽信道统一写为

$$
\mathbf g
\sim
\mathcal{CN}(\mathbf0,\mathbf R_g),
\qquad
\mathbf R_g
=
\mathbf R_{\rm phy}
\odot
\frac{\mathbf V\mathbf V^H}{N_t},
\qquad
R_{g,kk}=1.
$$

纯平坦信道对应 $\mathbf R_{\rm phy}=\mathbf1\mathbf1^H$；物理展宽信道使用实际 PDP 构造的 $\mathbf R_{\rm phy}$。定义

$$
Y_k
=
I_{\rm QAM}\!\left(\mathrm{snr}|g_k|^2\right),
\qquad
I_{\rm blk}
=
\frac1K\sum_kY_k.
$$

最终分集判据仍是

$$
P_{\rm out}(\mathbf V;\mathrm{snr},R,\theta)
=
\Pr[I_{\rm blk}< R],
\qquad
\mathrm{snr}_p(\mathbf V;R,\theta)
=
\inf\left\{
\mathrm{snr}:
P_{\rm out}\le p
\right\},
$$

其中 $\theta\in\Theta$ 表示物理信道统计参数，至少取 $p=0.1$ 和 $p=0.01$。归一化条件下 $\mathbb E[I_{\rm blk}]$ 与设计无关，而

$$
\operatorname{Var}(I_{\rm blk})
=
\frac1{K^2}
\sum_{k,l}
\operatorname{Cov}(Y_k,Y_l).
$$

当 $R<\mathbb E[I_{\rm blk}]$ 时，Cantelli 不等式给出

$$
\Pr[I_{\rm blk}< R]
\le
\frac{
\operatorname{Var}(I_{\rm blk})
}{
\operatorname{Var}(I_{\rm blk})
+
\left(\mathbb E[I_{\rm blk}]-R\right)^2
}.
$$

因此块方差是 outage 的有依据代理，但不决定 1% 左尾。对归一化联合 proper complex Gaussian 子载波对，令

$$
I_{\rm QAM}(\mathrm{snr}x)
-
\mathbb E[I_{\rm QAM}(\mathrm{snr}X)]
=
\sum_{q\ge1}a_q(\mathrm{snr})L_q(x),
\qquad
X\sim\operatorname{Exp}(1),
$$

其中 $L_q$ 是指数权重下的 Laguerre 多项式。记

$$
\rho_{\rm eff}(k,l)
=
R_{{\rm phy},kl}
\frac{(\mathbf V\mathbf V^H)_{kl}}{N_t},
\qquad
C(r;\mathrm{snr})
=
\sum_{q\ge1}a_q^2(\mathrm{snr})r^{2q},
$$

并定义偶数阶有效相关矩

$$
M_{2q}^{\rm eff}
=
\sum_{k\ne l}
|\rho_{\rm eff}(k,l)|^{2q}.
$$

则设计相关的块方差代理为

$$
J_{\rm corr}(\mathbf V;\mathbf R_{\rm phy},\mathrm{snr})
=
\sum_{k\ne l}
C\!\left(
|\rho_{\rm eff}(k,l)|;
\mathrm{snr}
\right)
=
\sum_{q\ge1}
a_q^2(\mathrm{snr})M_{2q}^{\rm eff},
$$

且

$$
\operatorname{Var}(I_{\rm blk})
=
\frac1{K^2}
\left[
K C(1;\mathrm{snr})
+
J_{\rm corr}
\right].
$$

这个式子说明只有$J_{\rm corr}$是$\operatorname{Var}(I_{\rm blk})$中的可设计部分。低成本筛选可截断为 $a_1^2M_2^{\rm eff}+a_2^2M_4^{\rm eff}+a_3^2M_6^{\rm eff}$。纯平坦信道下，这些相关矩汇总对应阶数的 Weyl 和能量；物理展宽信道下，应使用 $\mathbf R_{\rm phy}$ 与设计联合决定的加权相关矩或 PDP 加权核。上述代理只用于筛选和机理分析，最终仍使用 Monte Carlo outage 或链路 BLER。

在接收端已知 $\mathbf V$ 和 PDP 的 R2 条件下，最终优化对象是实际信道估计、LLR、LDPC 编译码和译码器均启用时，10% 和 1% estimated-CSI BLER 所需的 SNR。$\mathcal R_\Theta[J_{\rm corr}]$、Monte Carlo outage 和 $\mathcal R_\Theta[L_{\rm CE}]$ 是不同成本和不同精度的筛选或解释指标，不得单独替代最终链路判据。

默认不要求候选的平均 CE NMSE 不差于等差 CDD 基线，也不得仅因平均 NMSE 略差而淘汰候选。plan 若因稳定性或实现要求设置 $L_{\rm CE,max}$、零噪声误差地板或其他 CE 门槛，必须在正式结果前预先规定；未规定时，$L_{\rm CE}$ 作为独立诊断量或 Pareto 坐标报告。$\mathbf V_P^H\mathbf V_P=N_p\mathbf I_{N_t}$ 且无导频混叠在纯平坦 DFT 栅格模型下可保证分支可辨识并避免导频噪声增强，但在物理展宽信道下不是充分条件；实际 $\mathbf R_g$ 下的 matched 或 mismatched LMMSE NMSE 仍必须计算。

### 3.7 V-aware 信道估计目标

纯平坦信道下，导频观测为

$$
\mathbf y_P=\mathbf V_P\mathbf h+\mathbf n,
$$

其中 $P$ 是导频集合，$\mathbf V_P$ 是 $\mathbf V$ 的导频行子矩阵。导频域列正交条件为

$$
\mathbf V_P^H\mathbf V_P
=
N_p\mathbf I_{N_t}.
$$

一般协方差下，数据集合 $D$ 的 matched LMMSE NMSE 为

$$
L_{\rm CE}(\mathbf V;\mathbf R_{\rm phy})
=
\frac{
\operatorname{tr}\!\left[
\mathbf R_{DD}
-
\mathbf R_{DP}
(\mathbf R_{PP}+\sigma_{LS}^2\mathbf I)^{-1}
\mathbf R_{PD}
\right]
}{
\operatorname{tr}(\mathbf R_{DD})
}.
$$

该式同时包含导频处的可观测性和从导频到数据 RE 的可预测性。$\operatorname{rank}(\mathbf V_P)$、条件数、最小奇异值和折叠位置只作为诊断；物理展宽信道下的最终 CE 判据应使用实际 $\mathbf R_g$ 计算的 NMSE。

### 3.8 导频混叠

频域 DMRS comb 间隔为 $S_f$ 时，无歧义时延周期为

$$
\tau_{\rm alias}
=
\frac1{S_f\Delta f}.
$$

相差整数倍 $\tau_{\rm alias}$ 的路径在导频上不可区分。对 DFT 栅格 CDD，令

$$
N_p=\frac K{S_f},
$$

则导频上的 CDD 相位只由

$$
j_n\bmod N_p
$$

决定。

纯平坦信道下，折叠索引互异可保证 CDD 导频列可辨识；均匀折叠位置可进一步改善条件数。物理展宽信道下，每个折叠中心被 PDP 展宽，必须同时考虑：

1. 折叠中心之间的圆周距离；
2. 有效 PDP 支撑和定时误差；
3. $\mathbf R_{PP}$ 的数值条件；
4. $\mathbf R_{DP}$ 对数据 RE 的预测能力；
5. 接收机使用的协方差是否与实际 $\mathbf V$ 和 PDP 匹配。

分集相关目标和 CE NMSE 是不同的矩阵泛函，但由同一个 $\mathbf V$ 和 $\mathbf R_{\rm phy}$ 决定，因此应作为两个独立目标做 Pareto 设计。

## 4. 研究地图

### 4.1 三组假设轴

#### 信道模型

| 编号 | 假设 | 理论重点 |
|---|---|---|
| C0 | 纯平坦分支信道 | 模式分解、Weyl 和、DFT 栅格、二阶正交和高阶加性结构 |
| C1 | 物理展宽信道 | $\mathbf R_{\rm phy}\odot\mathbf V\mathbf V^H$、PDP 加权相关结构、复合时延和 CE |

#### 发射端对物理 PDP 的知识

| 编号 | 假设 | 可使用的信息 |
|---|---|---|
| T0 | 不知道 PDP | 只能使用预先规定的运行包络或信道模型集合 |
| T1 | 知道时延展宽 | 知道 RMS delay spread、有效支撑或其上界，不知道完整抽头功率和位置 |
| T2 | 知道完整 PDP 参数 | 知道 profile、抽头时延、功率及构造 $\mathbf R_{\rm phy}$ 所需的统计参数 |

#### 接收端知识

| 编号 | 假设 | 接收机能力 |
|---|---|---|
| R0 | 不知道 $\mathbf V$，也不知道底层物理 PDP | 使用固定通用先验、窗口或局部插值规则 |
| R1 | 完美知道底层物理 PDP，但不知道 $\mathbf V$ | 能匹配物理信道，不能匹配人工预编码造成的等效协方差 |
| R2 | 知道 $\mathbf V$ 和底层物理 PDP | 可构造实际 $\mathbf R_g$ 并进行 matched MMSE 估计 |

“发射端知识”和“接收端知识”必须分别声明。发射端根据 PDP 选择 $\mathbf V$，不代表接收端自动知道所选 $\mathbf V$ 或相应等效协方差。

### 4.2 纯平坦信道地图

纯平坦信道没有需要适配的 PDP 形状，因此 T0、T1、T2 的区别退化。主要分支由接收端知识决定：

| 接收端 | 研究方向 | 理论指导 |
|---|---|---|
| R0 | 透明、可由固定接收机估计的一般 $\mathbf V$ | 限制局部 group delay、相位连续性和局部模型维度 |
| R1 | 与 R0 接近；已知平坦物理模型不能补偿未知 $\mathbf V$ | 评价人工相位造成的协方差失配 |
| R2 | V-aware 分集与 matched CE 联合设计 | 先消除二阶模式，再优化四阶及更高阶模式；DFT 栅格下产生 Sidon 型指导 |

纯平坦 C0+R2 是 Sidon 理论的直接适用分支。该分支用于理解一般 $\mathbf V$ 的高阶相关自由度，不代表物理展宽信道下的最优设计。

### 4.3 物理展宽信道的完整地图

| 发射端知识 | R0：不知道 $\mathbf V$/PDP | R1：知道 PDP、不知道 $\mathbf V$ | R2：知道 $\mathbf V$/PDP |
|---|---|---|---|
| T0：不知道 PDP | 在预定信道包络内设计透明、局部可估计的固定 $\mathbf V$ | 设计不能依赖具体 PDP；同时限制人工结构对物理匹配估计器的失配 | 在信道模型集合上优化通用、PDP 鲁棒的 V-aware 设计 |
| T1：知道时延展宽 | 用时延支撑约束局部 group delay，但接收端仍使用通用模型 | 以已知物理 PDP 接收机的失配 MSE 为目标，限制复合支撑 | 按展宽分档设计厚 Sidon、均匀折叠 lift 或低维一般 $\mathbf V$ |
| T2：知道完整 PDP | PDP 匹配的发射设计可能对未知 $\mathbf V$ 接收机造成额外失配，必须约束透明性 | 发射端和接收端都知道 PDP，但接收端不知道 $\mathbf V$；优化必须显式计算等效协方差失配 | 使用实际 $\mathbf R_{\rm phy}$ 联合优化分集相关结构与 matched NMSE |

当前主要关心 C1+R2，并分别研究 T0、T1、T2。该研究方向不要求整个复合路径集合严格满足 Sidon，而是寻找比等差 CDD 和无结构随机 delay 更好的结构化 delay 或一般 $\mathbf V$。

### 4.4 C1+T0+R2：发射端不知道 PDP

完全没有信道范围上界时，不存在对任意 PDP 都保证可辨识和分集改善的有限设计。因此必须定义运行包络 $\Theta_{\rm all}$，例如允许的 TDL profile、时延展宽、定时误差和移动性范围。

设计目标为

$$
\min_{\mathbf V}
\left(
\mathcal R_{\Theta_{\rm all}}[J_{\rm div}(\mathbf V)],
\mathcal R_{\Theta_{\rm all}}[L_{\rm CE}(\mathbf V)]
\right),
$$

其中 $\mathcal R$ 可以是最坏情形、均值加尾部惩罚或其他预先规定的风险泛函。

理论指导：

- 导频折叠位置优先保持均匀和最大分离；
- 在多个 PDP 上联合压低加权相关代理，而不是针对单一 $\mathbf R_{\rm phy}$；
- 保持全带列正交或对其偏离施加明确惩罚；
- 候选必须对未参与优化的 PDP 做外推验证；
- 目标是固定、无需 PDP 反馈的鲁棒设计。

### 4.5 C1+T1+R2：发射端知道时延展宽

只知道 RMS delay spread 时，必须进一步规定从 RMS 值到有效支撑的模型集合。定义 $T_\epsilon$ 为覆盖 $1-\epsilon$ PDP 能量的有效时延支撑。

人工二元和为

$$
s_{ac}=\tau_a+\tau_c.
$$

物理展宽后，每个二元和形成近似宽度 $2T_\epsilon$ 的簇。可将严格 Sidon 放宽为带保护距离的“厚 Sidon”指导：

$$
d_{\rm circ}(s_{ac},s_{bd})
>
2T_\epsilon+T_{\rm margin},
\qquad
\{a,c\}\ne\{b,d\}.
$$

CE 侧独立要求折叠中心满足

$$
d_{\rm fold}(\tau_a,\tau_b)
>
T_\epsilon+T_{\rm sync}.
$$

两个保护距离分别控制近四元共振和导频可辨识性，不得互相替代。硬约束不可行时，应最大化两类最小距离或使用软惩罚形成 Pareto 前沿。

适合的设计形式是按时延展宽分档的结构化码本，而不是为每个具体 PDP 单独优化。

### 4.6 C1+T2+R2：发射端知道完整 PDP

发射端可以从完整 PDP 参数构造 $\mathbf R_{\rm phy}$，直接优化

$$
J_{\rm corr}(\mathbf V;\mathbf R_{\rm phy})
=
\sum_{k\ne l}
C\!\left(
\left|
R_{{\rm phy},kl}
\frac{(\mathbf V\mathbf V^H)_{kl}}{N_t}
\right|;
\mathrm{snr}
\right),
$$

其中 $C(\cdot;\mathrm{snr})$ 是两个相关复高斯子载波对应的互信息协方差核。低成本筛选可用 $M_4$，最终使用 MC outage。

CE 目标使用第 3.7 节的实际 matched NMSE。完整 PDP 已知时，低功率远端抽头不应与主抽头按相同权重计数，因此不再优先使用无权“复合 Sidon”条件。

标称 PDP 的最优设计还应在其参数邻域内验证：

$$
\theta'\in\mathcal N(\hat\theta),
$$

防止设计对 PDP 估计误差、定时误差或 profile 偏差过度敏感。

## 5. 设计路线

本节规定候选的生成与筛选路线。三条路线都只产生候选或低成本排序，不直接证明 estimated-CSI BLER 改善。正式比较必须使用相同 DMRS、发射功率、信道样本、接收机知识和链路配置，并同时保留等差 CDD 和经过优化的 CDD 作为参考。

### 5.1 路线一：基函数参数化并优化 $J_{\rm corr}$

一般 $\mathbf V$ 可以按具体 $(k,l)$ 子载波对重新分配相关能量。推荐从一个基准 CDD 出发，使用低维附加相位

$$
V_{k,n}
=
e^{-j2\pi kj_n/K}
e^{j\psi_{k,n}},
\qquad
\psi_{k,n}
=
\sum_{b=1}^{B}c_{b,n}B_b(k),
$$

其中 $B$ 是基函数数量，$B_b(k)$ 是第 $b$ 个实值频域基函数，$c_{b,n}\in\mathbb R$ 是待优化系数。pilot-anchored 参数化进一步要求

$$
B_b(k_p)=0,
\qquad
k_p\in P,
$$

其中 $k_p$ 是导频集合 $P$ 中的子载波索引。因此导频行上的附加相位为零，$\mathbf V_P$ 与基准 CDD 完全相同。该约束保持 DMRS 开销、导频相位码以及纯平坦模型下的 pilot Gram 不变，但不固定物理展宽信道下的 $\mathbf R_{DP}$，所以不能保证 CE NMSE 不变。

给定发射端知识对应的信道集合 $\Theta$ 后，搜索目标为

$$
\min_{\{c_{b,n}\}}
\mathcal R_\Theta\!\left[
J_{\rm corr}
\bigl(
\mathbf V(\{c_{b,n}\});
\mathbf R_{\rm phy},
\mathrm{snr}
\bigr)
\right],
$$

其中使用的 SNR 集合、风险泛函和参数边界必须在 plan 中预先规定。可采用投影梯度、Riemannian 优化、增广拉格朗日、坐标下降或其他可复现的迭代方法。

基函数的主要作用是降维和规定结构，不保证单独使用某一类基函数就能降低 $J_{\rm corr}$。可选的低维结构包括：

- pilot 区间内零端点的二次相位或 chirp；
- pilot 处为零的连续样条；
- 少量连续分段斜率；
- 由 $J_{\rm corr}$ 在基准 CDD 附近的梯度、Hessian 或物理协方差加权算子导出的主要方向。

pilot anchoring 会缩小可行域，因此可能牺牲一部分理论分集收益。应把 pilot-anchored 方案作为主要结构化路线，并以允许导频行变化的 unanchored 优化作为数值上界或消融；两者的 outage 和链路差异用于量化 anchoring 代价。

最终候选应使用少量共享、对称、可量化或可写成明确规则的系数。若只能用大量独立浮点系数保持性能，该结果只能作为数值上界或提炼新结构的证据，不作为最终可解释设计。

发射端知识决定 $J_{\rm corr}$ 的具体风险对象：T2 使用完整 PDP 构造的实际 $\mathbf R_{\rm phy}$ 直接优化，并在其参数邻域内验证；T1 在与已知 RMS delay spread 或有效支撑相容的 PDP 场景集上做最坏值或尾部风险优化；T0 在预先规定的运行包络上优化一个固定参数集，并用未参与搜索的 PDP 检查外推。

### 5.2 路线二：优化 $J_{\rm CDD}$ 并搜索 residue/lift

对于频率平稳的物理协方差，记

$$
R_{{\rm phy},kl}
=
r_{\rm phy}(k-l),
$$

并定义 CDD delay 集的阵列因子

$$
A_\tau(d)
=
\frac1{N_t}
\sum_{n=0}^{N_t-1}
e^{-j2\pi d\Delta f\tau_n},
$$

其中 $d$ 是子载波频差索引。CDD 下的相关代理为

$$
J_{\rm CDD}
=
2\sum_{d=1}^{K-1}
(K-d)
C\!\left(
\left|
r_{\rm phy}(d)A_\tau(d)
\right|;
\mathrm{snr}
\right).
$$

完整 PDP 已知时直接使用对应的 $r_{\rm phy}(d)$；只知道时延展宽或不知道具体 PDP 时，在预先规定的场景集合上优化 $J_{\rm CDD}$ 的最坏值、均值加尾部惩罚或其他风险泛函。

对有效带宽 DFT 栅格 CDD，

$$
\tau_n
=
\frac{j_n}{K\Delta f},
\qquad
N_p
=
\frac K{S_f},
$$

可写成

$$
j_n
=
r_{\pi(n)}
+N_pq_n.
$$

其中 $r_m$ 是 $j_n\bmod N_p$ 的导频折叠余数，$q_n\in\mathbb Z$ 是整数 lift，$\pi$ 表示 residue 与 lift 的配对。对导频子载波 $k=pS_f$，其中 $p$ 是整数导频频域索引，lift 引入的附加相位满足

$$
\exp\left(
-j2\pi
\frac{pS_fN_pq_n}{K}
\right)
=
e^{-j2\pi pq_n}
=
1,
$$

所以 lift 不改变导频行相位，但会改变数据 RE 相位、全带相关旁瓣和加性结构。

当各天线独立且具有相同 PDP 时，普通天线标签排列不影响性能，可直接令每个 residue 对应一个 $q_r$，把 $\pi$ 吸收到 lift 的定义中。固定 lift 多重集合后搜索 residue/lift 配对，或各天线 PDP、空间相关不同而 delay 分配会影响性能时，必须显式保留 $\pi$。

该路线是离散组合优化。搜索变量至少包括 residue 集、整数 lift 及必要时的配对；可使用有限范围穷举、coordinate exchange、branch-and-bound、局部搜索或退火。T0/T1 下应优先从均匀、最大圆周分离的 residue 开始；T2 下可在实际 CE 和鲁棒性检查允许时联合搜索非均匀 residue。每个候选都必须重新计算实际 $\mathbf R_{PP}$、$\mathbf R_{DP}$、闭式 CE NMSE 和零噪声误差地板。

T2 直接使用完整 PDP 对应的 $r_{\rm phy}(d)$ 搜索 residue/lift；T1 在支撑或展宽约束生成的多 PDP 集合上优化鲁棒 $J_{\rm CDD}$，并同时保持导频折叠中心的保护距离；T0 固定均匀或最大分离的 residue，在预定信道包络上搜索通用 lift，不能把针对单一 PDP 的最优集合表述为无知识通用设计。

### 5.3 路线三：优化 $J_{\rm CDD}$ 的四阶部分并设计 weighted almost-Sidon

CDD 下有

$$
J_{\rm CDD}
=
a_1^2M_2^{\rm eff}
+a_2^2M_4^{\rm eff}
+\sum_{q\ge3}
a_q^2M_{2q}^{\rm eff}.
$$

其中四阶部分可写成

$$
M_4^{\rm eff}
=
\frac1{N_t^4}
\sum_{a,b,c,d}
\mathcal K_{\rm PDP}
\left(
j_a+j_c-j_b-j_d
\right).
$$

纯平坦 DFT 栅格信道下，$\mathcal K_{\rm PDP}$ 退化为模 $K$ 共振计数，严格 Sidon 条件通过消除非平凡二元和相等来降低 $M_4^{\rm eff}$。物理展宽信道下，精确相等、近相等和位于不同频差的二元和差具有不同权重，因此应最小化上述 PDP 加权加性能量，得到 weighted almost-Sidon，而不是只判断是否存在严格四元关系。

weighted almost-Sidon 不是独立于 $J_{\rm CDD}$ 的最终性能目标，而是 $J_{\rm CDD}$ 四阶分量的结构化代理。它适合：

- 根据二元和及其加权间距生成可解释候选；
- 对 residue/lift 搜索做剪枝或初始化；
- 解释候选相对等差 CDD 改变了哪些四阶共振；
- 在完整 $J_{\rm CDD}$ 优化前进行低成本筛选。

最终仍需检查完整 $J_{\rm CDD}$。仅降低 $M_4^{\rm eff}$、但使 $M_2^{\rm eff}$、更高阶矩、CE 或链路 BLER 恶化的候选，不能表述为更优设计。只知道有效支撑 $T_\epsilon$ 时，第 4.5 节的厚 Sidon 保护距离可作为 $\mathcal K_{\rm PDP}$ 未知时的几何近似。

具体而言，T2 在完整 PDP 已知时直接最小化对应的 $\mathcal K_{\rm PDP}$ 加权核；T1 只知道有效支撑或其上界时，使用厚 Sidon 二元和保护距离，或在与该支撑相容的核集合上优化最坏加性能量；T0 只能在预定 PDP 包络上搜索通用 weighted almost-Sidon，严格无权 Sidon 可作为平坦场景初始化，但不能作为任意展宽信道的保证。

### 5.4 推荐优化流程

1. **固定问题定义。** 明确物理信道集合、发射端和接收端知识、DMRS、功率归一化、SNR 集合、循环移位定义、基线、候选预算、随机种子和停止条件。
2. **确定性生成与粗筛。** 分别按路线一优化 $J_{\rm corr}$、按路线二优化完整 $J_{\rm CDD}$、按路线三优化或约束 $M_4^{\rm eff}$。同时计算 pilot rank、条件数、实际协方差下的 CE NMSE 和零噪声误差地板。默认不要求平均 CE NMSE 不差于等差 CDD，但应排除不可辨识、数值失稳或违反 plan 预设门槛的候选。
3. **Monte Carlo outage。** 对粗筛候选使用实际 QAM/BICM 互信息计算 10% 和 1% outage 所需 SNR。该阶段评价 ideal-CSI 分集潜力，输出是 outage，不是 BLER。
4. **estimated-CSI 链路预扫。** 对保留不同 outage/CE 取舍的少量 Pareto 候选运行真实信道估计、LLR、LDPC 编译码和译码器的小样本链路预扫，尽早识别代理指标与 BLER 排序不一致的候选。
5. **正式链路验收。** 对入围候选使用相同随机信道、payload、噪声、DMRS、MCS 和译码器运行成对正式仿真，报告 10% 和样本允许时的 1% estimated-CSI BLER 所需 SNR、错误块计数和置信区间。
6. **鲁棒性与解释。** 在未参与搜索的 PDP、协方差偏差、定时误差和带宽上验证；同时报告 $J_{\rm corr}$ 或 $J_{\rm CDD}$、各阶相关矩、CE 指标和生成参数，用于区分事实、机理推断和待验证解释。

最终胜出条件是实际 estimated-CSI BLER 相对预定基线改善，并满足 plan 中明确规定的实现和稳定性要求。理想 CSI outage、平均 CE NMSE、$M_4^{\rm eff}$、Gram 和条件数均不能单独作为胜出判据。

### 5.5 注意事项

1. **时延和相位单位必须统一。** 使用有效带宽 DFT 栅格时，

   $$
   V_{k,n}
   =
   e^{-j2\pi kj_n/K},
   \qquad
   \tau_n
   =
   \frac{j_n}{K\Delta f}.
   $$

   使用 FFT sample 循环移位时，

   $$
   V_{k,n}
   =
   e^{-j2\pi kd_n/N_{\rm FFT}},
   \qquad
   \tau_n
   =
   \frac{d_n}{N_{\rm FFT}\Delta f}.
   $$

   两种定义产生不同的相位矩阵。配置、输出和 result 必须记录相位分母、$k$ 是有效带宽局部索引还是物理 FFT bin、索引单位以及秒或采样点换算，不得把 `/K` 与 `/N_{\rm FFT}` 视为同一 delay 集。
2. **人工时延是循环移位。** 当前模型中的人工 delay 通过循环移位或等价频域线性相位实现，不占用真实传播时延的 CP 预算，因此不施加“人工时延跨度加物理 PDP 必须小于 CP”的约束。仍需规定循环周期、允许索引、整数或分数移位、相位量化和实际实现方式。
3. **消除等价候选。** 在独立同分布分支模型下，公共循环移位和普通天线标签排列不改变分集性能；搜索与结果汇总应使用规范化代表，避免重复计数。存在分支相关或不同 PDP 时不得使用该等价化。
4. **区分 pilot 几何与数据可预测性。** $\mathbf V_P$ 相同、满秩或条件数为 1 都不保证 $\mathbf R_{DP}$ 相同，也不保证 CE NMSE 或 BLER 相同。
5. **CE NMSE 是诊断量而非默认硬比较门槛。** 平均 NMSE 略差的候选仍可能获得更好的最终 BLER；但高 SNR 误差地板、不可辨识方向、协方差失配和局部极端误差必须报告。若 plan 设置硬门槛，应在看到正式结果前固定。
6. **保持可解释和可复现。** 一般 $\mathbf V$ 优先使用低维共享参数，CDD 优先报告 residue、lift、二元和及加权核贡献。完全随机矩阵或大量无法压缩的浮点系数只能作为初始化或数值上界。
7. **发射端和接收端知识分别记录。** 发射端按 PDP 选择 $\mathbf V$ 不代表接收端知道所选矩阵；R0/R1 必须使用实际固定或失配估计器的 MSE 和 BLER，不能引用 R2 matched 结果。

## 6. 参考系统与公平性

默认参考条件为：

| 参数 | 默认值 |
|---|---|
| 发射与接收 | 8 Tx / 1 Rx，单层 |
| 子载波间隔 $\Delta f$ | 30 kHz |
| 有效子载波数 $K$ | 576 / 432 / 288，对应 48 / 36 / 24 PRB |
| DMRS comb $S_f$ | 24 子载波 |
| 导频数 $N_p=K/S_f$ | 24 / 18 / 12 |
| DMRS 符号 | 2 |
| 调制与码率 | 16QAM，553/1024 |

比较不同 $\mathbf V$ 时必须固定：

- DMRS 图样和开销；
- 发射功率和 $\mathbf V$ 归一化；
- MCS、数据 RE 和解码器；
- 信道模型与发射端/接收端知识；
- 接收机是否 matched；
- 随机样本、统计口径和置信区间；
- 循环移位的相位分母、索引语义、周期、量化和实现约束。

## 7. 规范符号

| 符号 | 含义 |
|---|---|
| $N_t$ / $K$ | 发射分支数 / 有效子载波数 |
| $\Delta f$ | 子载波间隔 |
| $\mathbf V$ | 恒模频域相位矩阵 |
| $\phi_{k,n}$ | 第 $k$ 个子载波、第 $n$ 个分支的相位 |
| $\mathbf h$ / $\mathbf g$ | 底层分支信道 / 等效频域信道 |
| $\mathbf R_{\rm phy}$ / $\mathbf R_g$ | 物理信道 / 等效信道频率协方差 |
| $S_f$ / $N_p$ | DMRS comb 间隔 / 每个 DMRS 符号的导频数 |
| $P,D$ | 导频 / 数据 RE 集合 |
| $I_{\rm QAM}$ | 实际 QAM/BICM 每 RE 互信息函数 |
| $I(\mathbf h)$ / $I_\infty(\mathbf h)$ | 块平均互信息 / 相位环面平均 |
| $Y_k$ / $I_{\rm blk}$ | 第 $k$ 个子载波的互信息贡献 / 块平均互信息 |
| $R$ / $P_{\rm out}$ | 目标谱效率 / outage 概率 |
| $\mathrm{snr}_p$ | 达到 outage 概率 $p$ 所需的最小 SNR |
| $\boldsymbol\phi_k$ | $\mathbf V$ 第 $k$ 行的相位向量 |
| $\mathbf m$ / $c_{\mathbf m}$ | 相位环面 Fourier 模式 / 模式系数 |
| $W_K(\mathbf m)$ | 设计轨迹的 Weyl 和 |
| $\rho_{kl}$ / $\rho_{\rm eff}(k,l)$ | $\mathbf V$ 行相关 / 物理协方差加权后的相关 |
| $M_{2q}^{\rm eff}$ / $M_4^{\rm eff}$ | 偶数阶有效相关矩 / 其中的四阶相关代理；上下文明确时简写为 $M_4$ |
| $\mathcal K_{\rm PDP}$ | PDP 加权加性核 |
| $L_q$ / $a_q$ | Laguerre 多项式 / 互信息非线性的 Laguerre 系数 |
| $C(\cdot;\mathrm{snr})$ | 相关复高斯子载波对应的互信息协方差核 |
| $r_{\rm phy}(d)$ / $A_\tau(d)$ | 物理频率相关的 lag 表示 / CDD delay 集的阵列因子 |
| $\tau_n$ / $j_n$ | CDD 时延 / DFT 栅格索引 |
| $N_{\rm FFT}$ / $d_n$ | FFT 长度 / FFT sample 循环移位索引 |
| $\delta_p$ | 物理 PDP 第 $p$ 个抽头时延 |
| $\tau_{n,p}$ | 第 $n$ 个人工时延与第 $p$ 个物理抽头的复合路径时延 |
| $\tau_{\rm alias}$ | 导频无混叠时延周期 |
| $\epsilon$ / $T_\epsilon$ | 允许忽略的 PDP 能量比例 / 覆盖 $1-\epsilon$ 能量的有效时延支撑 |
| $s_{ac}$ | 两个人工时延的二元和 $\tau_a+\tau_c$ |
| $T_{\rm margin}$ / $T_{\rm sync}$ | 二元和保护余量 / 定时误差余量 |
| $J_{\rm div}$ / $J_{\rm corr}$ / $J_{\rm CDD}$ | 分集目标 / 一般 $\mathbf V$ 互信息相关代理 / CDD 特化相关代理 |
| $L_{\rm CE}$ / $L_{\rm CE,max}$ | 实际协方差下的信道估计 NMSE / 允许的 NMSE 上限 |
| $\mathcal V$ | 满足功率、DMRS 和实现约束的候选设计集合 |
| $\theta$ / $\Theta$ | 一组信道统计参数 / 信道参数集合 |
| $\mathcal R_\Theta$ / $\mathcal N(\hat\theta)$ | 信道集合上的风险泛函 / 标称信道参数的失配邻域 |
| $d_{\rm circ}$ / $d_{\rm fold}$ | 二元和圆周距离 / 导频混叠圆距离 |
| $r_n$ / $q_n$ / $\pi$ | 导频折叠余数 / 整数 lift / residue 与 lift 的配对排列 |
| $B$ / $B_b(k)$ / $c_{b,n}$ | 基函数数量 / 第 $b$ 个频域基函数 / 第 $n$ 个分支的基函数系数 |
| $\psi_{k,n}$ | pilot-anchored 一般 $\mathbf V$ 的附加相位 |

新增符号必须先加入本表。历史文档与本文件的当前符号定义冲突时，以本文件为准。
