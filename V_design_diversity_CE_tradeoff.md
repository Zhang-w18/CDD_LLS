# 一般性相位矩阵 V 设计：分集增益与信道估计处理增益之间的折中

## 1. 问题背景

在显式 CDD 或更一般的频域相位预编码中，不同发射分支在不同子载波上叠加不同相位，使得接收端看到的等效频域信道不再是简单的多天线同相叠加，而是

$$
\mathbf h_{\mathrm{eff}} = \mathbf V \mathbf h.
$$

其中：

- $\mathbf h \in \mathbb C^{N}$ 表示 $N$ 个发射分支上的底层独立信道；
- $\mathbf h_{\mathrm{eff}} \in \mathbb C^{K}$ 表示一个频域局部窗口内 $K$ 个子载波上的等效信道；
- $\mathbf V \in \mathbb C^{K\times N}$ 是发射端施加的频域相位矩阵；
- $K$ 可以理解为一个相干带宽或一个信道估计处理窗口内包含的子载波数；
- $N$ 是参与叠加的 i.i.d. 发射分支数。

若底层信道满足

$$
\mathbf h \sim \mathcal{CN}(0,\mathbf I_N),
$$

则等效频域信道协方差为

$$
\mathbf R_{\mathrm{eff}}
=\mathbb E[\mathbf h_{\mathrm{eff}}\mathbf h_{\mathrm{eff}}^H]
=\mathbf V\mathbf V^H.
$$

非零特征值等价于 $\mathbf V^H\mathbf V$ 的特征值，也就是 $\mathbf V$ 的奇异值模方。因此，$\mathbf V$ 的奇异值分布决定了不同发射分支在频域上被“展开”的程度。

本文关注的问题是：

> 能否设计比传统 CDD 更一般的 $\mathbf V$，在分集增益和信道估计处理增益之间获得更优的 Pareto 前沿？

传统 CDD 对应的是全带线性相位矩阵；本文考虑是否可以使用分段线性相位、chirp 相位、平滑随机相位码等更一般的结构。

---

## 2. 基本模型

### 2.1 底层发射分支信道

先考虑一个简化模型：在一个较窄频域窗口内，$N$ 个发射分支的底层信道近似为频率不变，记为

$$
\mathbf h = [h_0,h_1,\ldots,h_{N-1}]^T,
$$

并假设

$$
h_n \sim \mathcal{CN}(0,1),
\qquad
\mathbb E[h_nh_m^*]=\delta_{n,m}.
$$

更一般地，也可以令

$$
\mathbf h \sim \mathcal{CN}(0,\mathbf R_h),
$$

此时等效信道协方差为

$$
\mathbf R_{\mathrm{eff}}=\mathbf V\mathbf R_h\mathbf V^H.
$$

为了突出 $\mathbf V$ 本身的结构，后文主要以 $\mathbf R_h=\mathbf I_N$ 为例。

### 2.2 一般性频域相位矩阵

令

$$
V_{k,n}=e^{j\phi_{k,n}},
\qquad k=0,\ldots,K-1,
\quad n=0,\ldots,N-1.
$$

这里暂时只考虑 constant-modulus 约束：

$$
|V_{k,n}|=1.
$$

这对应每个发射分支在每个子载波上只做相位旋转，不改变功率分配。

传统 CDD 是其中一个特殊形式：

$$
V_{k,n}=e^{-j2\pi k\Delta f\tau_n},
$$

其中 $\tau_n$ 是第 $n$ 个发射分支的固定 cyclic delay 或 true time delay。

因此，传统 CDD 的相位是全带线性的：

$$
\phi_{k,n}=-2\pi k\Delta f\tau_n.
$$

而一般性 $\mathbf V$ 可以允许 $\phi_{k,n}$ 随频率非线性变化，例如分段线性、二次相位、平滑随机相位等。

---

## 3. 分集增益指标

### 3.1 奇异值与可观测维度

对 $\mathbf V$ 做 SVD：

$$
\mathbf V = \mathbf U\mathbf \Sigma\mathbf W^H.
$$

若 $K\ge N$，则 $\mathbf V^H\mathbf V$ 有 $N$ 个非零特征值：

$$
\lambda_i = \sigma_i^2,
\qquad i=1,\ldots,N.
$$

若 $\mathbf V$ 满列秩，则 $N$ 个独立发射分支在频域上都可以被等效信道观察到。若某些奇异值很小，则对应的发射分支组合虽然理论上存在，但在接收端很难被观测到，实际分集或模式增益会变差。

因此可以用以下指标衡量 $\mathbf V$ 的分集展开能力。

### 3.2 log-det 指标

$$
J_{\mathrm{div}}(\mathbf V)
=
\log\det\left(\frac{1}{K}\mathbf V^H\mathbf V+\epsilon\mathbf I\right).
$$

当 $\epsilon$ 很小时，该指标等价于最大化所有奇异值模方的乘积：

$$
\prod_{i=1}^{N}\sigma_i^2.
$$

它反映的是总体模式体积，也可以理解为一种 coding gain / conditioning gain。

### 3.3 最小奇异值指标

也可以考虑

$$
J_{\min}(\mathbf V)=\lambda_{\min}(\mathbf V^H\mathbf V).
$$

该指标强调最弱模式，避免某个发射分支组合几乎不可观测。

### 3.4 条件数指标

$$
\kappa(\mathbf V)
=
\frac{\sigma_{\max}(\mathbf V)}{\sigma_{\min}(\mathbf V)}.
$$

理想情况下希望

$$
\mathbf V^H\mathbf V \approx K\mathbf I_N,
$$

也就是所有奇异值尽量相等。

在 constant-modulus 约束下，每一列范数固定为

$$
\|\mathbf v_n\|^2=K.
$$

根据 Hadamard 不等式：

$$
\det(\mathbf V^H\mathbf V)
\le
\prod_{n=0}^{N-1}\|\mathbf v_n\|^2
=K^N.
$$

等号成立当且仅当列向量两两正交：

$$
\mathbf V^H\mathbf V=K\mathbf I_N.
$$

因此，对固定 $K,N$ 而言，最大分集展开的理想目标是构造一个近似 tight frame 或局部 DFT 子矩阵。

---

## 4. 信道估计处理增益指标

分集展开要求不同发射分支的相位序列在频域上足够不同；但信道估计处理增益要求等效信道在频域上足够平滑。二者天然存在冲突。

### 4.1 频域平滑性

令相邻子载波上的相位差为

$$
\Delta\phi_{k,n}=\phi_{k+1,n}-\phi_{k,n}.
$$

对应的局部 group delay 可定义为

$$
\tau_{k,n}^{\mathrm{group}}
=
-\frac{1}{2\pi\Delta f}\Delta\phi_{k,n}.
$$

若 $\tau_{k,n}^{\mathrm{group}}$ 较小且随 $k$ 变化平缓，则 $V_{k,n}$ 在频域上平滑，等效信道相干带宽较大。

可定义粗糙度代价：

$$
J_{\mathrm{rough},1}(\mathbf V)
=
\sum_{n=0}^{N-1}\sum_{k=0}^{K-2}
|\phi_{k+1,n}-\phi_{k,n}|^2.
$$

也可以使用二阶差分，约束 group delay 的变化速度：

$$
J_{\mathrm{rough},2}(\mathbf V)
=
\sum_{n=0}^{N-1}\sum_{k=0}^{K-3}
|\phi_{k+2,n}-2\phi_{k+1,n}+\phi_{k,n}|^2.
$$

### 4.2 有效人工时延扩展

对 CDD 而言，人工时延扩展近似为

$$
\tau_{\mathrm{art}}
=\max_n\tau_n-\min_n\tau_n.
$$

等效信道的总时延扩展可以粗略写为

$$
\tau_{\mathrm{eff}}
\approx
\tau_{\mathrm{phy}}+	au_{\mathrm{art}},
$$

因此有效相干带宽为

$$
B_{c,\mathrm{eff}}
\approx
\frac{1}{\tau_{\mathrm{eff}}}.
$$

对于一般性 $\mathbf V$，可以用 group delay spread 来定义人工时延扩展：

$$
\tau_{\mathrm{art}}
=
\max_{k,n}\tau_{k,n}^{\mathrm{group}}
-
\min_{k,n}\tau_{k,n}^{\mathrm{group}}.
$$

或者用 RMS 形式：

$$
\tau_{\mathrm{art,rms}}^2
=
\frac{1}{KN}\sum_{k,n}
\left(\tau_{k,n}^{\mathrm{group}}-ar\tau\right)^2.
$$

### 4.3 信道估计处理增益

频域信道估计处理增益来自多个 pilot / data RE 之间的相关性。若等效信道在一个窗口内高度相关，就可以通过插值、平滑或 LMMSE 获得处理增益。

若 $K_c$ 表示等效相干带宽内可联合处理的子载波数，则处理增益粗略正比于

$$
G_{\mathrm{CE}}\sim K_c.
$$

人工相位变化越快，等效信道频域相关性越弱，$K_c$ 越小，信道估计处理增益越低。

因此，一个合理的优化目标不是单独最大化分集，而是联合考虑：

$$
\max_{\mathbf V}
\quad
J_{\mathrm{div}}(\mathbf V)
-
\lambda J_{\mathrm{CE\ loss}}(\mathbf V).
$$

---

## 5. CDD 的本质与局限

### 5.1 CDD 是全带线性相位

CDD 的相位为

$$
\phi_{k,n}=-2\pi k\Delta f\tau_n.
$$

不同发射分支之间的列相关为

$$
\langle \mathbf v_m,\mathbf v_n\rangle
=
\sum_{k=0}^{K-1}
 e^{-j2\pi k\Delta f(\tau_n-\tau_m)}.
$$

这是一个 Dirichlet kernel。若希望两列近似正交，需要相对相位在 $K$ 个子载波内至少旋转一圈左右：

$$
K\Delta f |\tau_n-\tau_m|\gtrsim 1.
$$

等价于

$$
|\tau_n-\tau_m|\gtrsim \frac{1}{K\Delta f}.
$$

这说明：为了获得列正交性，CDD 必须引入足够大的相对时延；但这会增加等效时延扩展，降低相干带宽。

### 5.2 CDD 是受限最优，不是一般最优

在“每个发射分支只能使用一个全带固定 delay”的约束下，CDD 是自然且接近最优的选择。若 delay 取在 DFT grid 上，可以在给定 $K$ 个子载波内得到正交列。

但是 CDD 的自由度只有 $N$ 个 delay：

$$
\{\tau_0,\tau_1,\ldots,\tau_{N-1}\}.
$$

而一般性相位矩阵有 $K\times N$ 个相位自由度：

$$
\{\phi_{k,n}\\}.
$$

因此，在允许频率相关相位设计的前提下，一般性 $\mathbf V$ 的搜索空间严格包含 CDD：

$$
\mathcal V_{\mathrm{CDD}}
\subset
\mathcal V_{\mathrm{smooth\ phase}}.
$$

所以一般性 $\mathbf V$ 的 Pareto 前沿理论上不会差于 CDD，实际有限维系统中有可能严格优于 CDD。

---

## 6. 为什么一般性 V 可能获得更优 Pareto 前沿

分集指标主要由 $\mathbf V^H\mathbf V$ 决定；信道估计处理增益主要由 $\mathbf V$ 在频域上的局部平滑性决定。这两个属性虽然相关，但并不完全等价。

CDD 使用一个全带固定斜率来同时控制二者：

$$
\phi_{k,n}=a_n k+b_n.
$$

这会导致一个问题：为了让 $\mathbf V^H\mathbf V$ 全带接近对角，CDD 需要在全带持续引入相对相位旋转；这等价于全带人工 delay spread。

一般性 $\mathbf V$ 可以更灵活地分配相位变化。例如：

1. 在不同频段使用不同相位斜率，使 off-diagonal 项在全带累加时互相抵消；
2. 在局部频段内保持相对平滑，以保留局部信道估计处理增益；
3. 使用非线性相位，让相位差更均匀地覆盖单位圆，降低最坏相关性；
4. 使用优化得到的平滑相位码，在固定 roughness 下最大化 $\log\det(\mathbf V^H\mathbf V)$。

因此可能实现：

$$
\text{相同 } J_{\mathrm{div}}
\quad \text{但更小的 } J_{\mathrm{CE\ loss}},
$$

或

$$
\text{相同 } J_{\mathrm{CE\ loss}}
\quad \text{但更大的 } J_{\mathrm{div}}.
$$

这就是一般性 $\mathbf V$ 可能优于纯 CDD Pareto 前沿的原因。

---

## 7. 候选方向一：分段线性相位

### 7.1 基本形式

将 $K$ 个子载波划分为 $S$ 个 segment：

$$
\mathcal S_1,\mathcal S_2,\ldots,\mathcal S_S.
$$

在第 $s$ 个 segment 内：

$$
V_{k,n}
=
e^{-j2\pi (k-k_s)\Delta f\tau_{s,n}+j\theta_{s,n}},
\qquad k\in\mathcal S_s.
$$

其中：

- $\tau_{s,n}$ 是第 $s$ 个 segment 内第 $n$ 个发射分支的局部 delay；
- $\theta_{s,n}$ 是该 segment 的初始相位；
- $k_s$ 是 segment 起始子载波索引。

### 7.2 全带相关矩阵分解

总 Gram 矩阵可以写成各 segment 贡献之和：

$$
\mathbf V^H\mathbf V
=
\sum_{s=1}^{S}\mathbf V_s^H\mathbf V_s.
$$

即使每个 segment 内 $\mathbf V_s^H\mathbf V_s$ 不是完全对角，不同 segment 的 off-diagonal 项也可能因为相位方向不同而互相抵消。

这给了分段线性相位优于全带 CDD 的机会。

### 7.3 连续性约束

如果 segment 边界相位不连续，会在时域产生长尾，相当于引入额外人工 delay spread。因此应至少满足相位连续：

$$
\phi_{s,n}(k_{s+1}^{-})
=
\phi_{s+1,n}(k_{s+1}^{+}).
$$

更稳健的设计还应限制斜率突变：

$$
|\tau_{s+1,n}-\tau_{s,n}|
\le
\Delta\tau_{\max}.
$$

这样分段线性相位可以被理解为“缓慢变化的 group delay”，而不是硬切换的频域相位码。

### 7.4 适合验证的问题

分段线性相位可以用来验证：

> 在相同最大局部 delay 或相同 group delay spread 约束下，是否能够获得比全带 CDD 更好的 $\mathbf V^H\mathbf V$ 条件数或更高的 $\log\det$ 指标？

---

## 8. 候选方向二：chirp / 二次相位

### 8.1 基本形式

可以令

$$
V_{k,n}=e^{j(a_n k^2+b_n k+c_n)}.
$$

此时相位一阶差分为

$$
\phi_{k+1,n}-\phi_{k,n}
\approx
2a_n k+a_n+b_n.
$$

对应 group delay 随频率近似线性变化：

$$
\tau_{k,n}^{\mathrm{group}}
\propto
-(2a_n k+b_n).
$$

### 8.2 潜在优势

相比硬分段线性相位，chirp 相位没有 segment 边界突变，因此时域旁瓣可能更可控。

同时，chirp 相位可以让不同发射分支之间的相位差以非线性方式展开，可能降低 Dirichlet kernel 型旁瓣的峰值。

### 8.3 需要注意的问题

chirp 相位不再对应传统 CDD 的固定 physical delay，而是频率相关 group delay。它更接近频域数字预编码或频变相位网络。因此实现复杂度和标准化解释性都弱于 CDD。

---

## 9. 候选方向三：平滑随机相位码

### 9.1 基本形式

令

$$
V_{k,n}=e^{j\phi_{k,n}},
$$

其中 $\phi_{k,n}$ 由随机序列生成后经过低通滤波，使其满足频域平滑约束。

### 9.2 优化目标

可以求解如下优化问题：

$$
\max_{\{\phi_{k,n}\}}
\quad
\log\det\left(\frac{1}{K}\mathbf V^H\mathbf V+\epsilon\mathbf I\right)
-
\lambda
\sum_{k,n}|\phi_{k+1,n}-2\phi_{k,n}+\phi_{k-1,n}|^2,
$$

subject to

$$
|V_{k,n}|=1.
$$

其中第二项约束相位曲率，也就是约束 group delay 的变化速度。

### 9.3 潜在价值

平滑随机相位码的可解释性较弱，但它可以作为数值上探索 Pareto 上界的工具。

如果优化相位码明显优于 CDD，则说明 CDD 不是一般性 $\mathbf V$ 空间中的 Pareto 最优；随后可以尝试把优化结果结构化为分段线性或 chirp 相位，增强工程可实现性。

---

## 10. 候选方向四：frame / Welch-bound 视角下的相位矩阵设计

也可以把 $\mathbf V$ 设计成 constant-modulus frame。目标是让列之间互相关尽量小：

$$
\mu(\mathbf V)
=
\max_{m\neq n}
\frac{|\mathbf v_m^H\mathbf v_n|}
{\|\mathbf v_m\|\|\mathbf v_n\|}.
$$

理想目标是最小化 mutual coherence：

$$
\min_{\mathbf V}\quad \mu(\mathbf V).
$$

同时加入平滑性约束：

$$
|\phi_{k+1,n}-\phi_{k,n}|
\le
2\pi\Delta f\tau_{\max},
$$

或者

$$
\sum_{k,n}|\phi_{k+1,n}-2\phi_{k,n}+\phi_{k-1,n}|^2
\le C.
$$

这个角度有助于把问题从传统 CDD 推广到“平滑 constant-modulus frame”设计。

---

## 11. 与信道估计模型的关系

### 11.1 直接等效信道协方差

若 UE 知道 $\mathbf V$，则等效信道协方差为

$$
\mathbf R_{g}=\mathbf V\mathbf R_h\mathbf V^H.
$$

在 pilot 位置 $P$ 和 data 位置 $D$ 上，有

$$
\mathbf R_{g,PP}=\mathbf V_P\mathbf R_h\mathbf V_P^H,
$$

$$
\mathbf R_{g,DP}=\mathbf V_D\mathbf R_h\mathbf V_P^H.
$$

因此 LMMSE 信道估计器为

$$
\hat{\mathbf g}_D
=
\mathbf R_{g,DP}
\left(\mathbf R_{g,PP}+\sigma^2\mathbf I\right)^{-1}
\tilde{\mathbf g}_P.
$$

### 11.2 CDD 是特殊闭式形式

对 CDD，$\mathbf V$ 是线性相位，因此 $\mathbf R_g$ 可以等价地由 shifted PDP 计算，得到 Toeplitz 型频域协方差。

对一般性 $\mathbf V$，$\mathbf R_g$ 可能不再是简单 Toeplitz 结构，但仍可以直接由

$$
\mathbf R_g=\mathbf V\mathbf R_h\mathbf V^H
$$

构造。

这说明一般性 $\mathbf V$ 并不会破坏 LMMSE 估计理论，只是 UE 需要知道 $\mathbf V$，并且 estimator 需要使用 $\mathbf V$-aware covariance。

### 11.3 处理增益的变化

如果 $\mathbf V$ 频域变化很快，$\mathbf R_{g,PP}$ 的相关性会快速衰减，LMMSE 无法从远处 pilot 获得太多增益。

如果 $\mathbf V$ 频域平滑，$\mathbf R_{g,PP}$ 更强相关，可以支持更宽频域窗口内的联合估计。

因此评价 $\mathbf V$ 时不能只看 $\mathbf V^H\mathbf V$，还必须看基于 $\mathbf R_g$ 的 CE NMSE 或 BLER。

---

## 12. 建议仿真评价指标

建议同时输出以下指标。

### 12.1 矩阵侧指标

1. $\log\det\left(\frac{1}{K}\mathbf V^H\mathbf V+\epsilon I\right)$；
2. $\lambda_{\min}(\mathbf V^H\mathbf V)$；
3. condition number：

$$
\kappa(\mathbf V)=\frac{\sigma_{\max}}{\sigma_{\min}};
$$

4. mutual coherence：

$$
\mu(\mathbf V)=\max_{m\neq n}
\frac{|\mathbf v_m^H\mathbf v_n|}
{\|\mathbf v_m\|\|\mathbf v_n\|}.
$$

### 12.2 平滑性指标

1. 最大 group delay spread：

$$
\tau_{\mathrm{art}}
=
\max_{k,n}\tau_{k,n}^{\mathrm{group}}
-
\min_{k,n}\tau_{k,n}^{\mathrm{group}};
$$

2. RMS group delay spread；
3. 一阶相位 roughness；
4. 二阶相位 roughness；
5. 等效频域相关函数半功率宽度。

### 12.3 信道估计指标

1. CE NMSE：

$$
\mathrm{NMSE}
=
\frac{\mathbb E\|\hat{\mathbf g}-\mathbf g\|^2}
{\mathbb E\|\mathbf g\|^2};
$$

2. data RE 上的 interpolation / LMMSE MSE；
3. estimator condition number；
4. 对 covariance mismatch 的鲁棒性。

### 12.4 链路级指标

1. BLER；
2. throughput / goodput；
3. required SNR at target BLER；
4. 相同 DMRS density 下的性能；
5. 相同 CE NMSE 下的分集收益。

---

## 13. 建议仿真方案

建议至少比较以下 $\mathbf V$：

| 方案 | 结构 | 目的 |
|---|---|---|
| No CDD | $V_{k,n}=1$ | 无人工频率选择性基线 |
| 全带 CDD | $V_{k,n}=e^{-j2\pi k\Delta f\tau_n}$ | 传统线性相位基线 |
| DFT-grid CDD | delay 取 DFT grid | 最大化列正交性的 CDD 上界 |
| 分段线性相位 | 每个 segment 一组局部 delay | 测试局部正交 + 局部平滑折中 |
| 连续分段线性相位 | 相位连续、斜率分段变化 | 避免边界长尾 |
| Chirp 相位 | $e^{j(a_nk^2+b_nk+c_n)}$ | 测试平滑非线性相位 |
| 平滑优化相位 | 数值优化 $\log\det - \lambda\mathrm{roughness}$ | 探索 Pareto 上界 |

### 13.1 推荐的横纵轴

可以画 Pareto 曲线：

$$
x = \mathrm{CE\ NMSE}
\quad \text{或} \quad
x = \tau_{\mathrm{art,rms}},
$$

$$
y = \log\det\left(\frac{1}{K}\mathbf V^H\mathbf V+\epsilon I\right)
\quad \text{或} \quad
 y = \mathrm{BLER/throughput\ gain}.
$$

如果某个一般性 $\mathbf V$ 在相同 $x$ 下获得更高 $y$，则说明它优于 CDD Pareto 前沿。

### 13.2 推荐扫描参数

1. $K$：局部处理窗口大小；
2. $N$：发射分支数；
3. DMRS density；
4. physical delay spread；
5. maximum artificial group delay spread；
6. segment 长度；
7. 相位连续性约束；
8. covariance 是否匹配。

---

## 14. 可能的理论命题

可以尝试形成如下理论命题。

### 命题 1：CDD 是全带固定 delay 约束下的自然最优结构

在每个发射分支只能使用一个全带固定 delay 的约束下，$\mathbf V$ 是 Vandermonde 矩阵。若 delay 取在 DFT grid 上，则可在给定 $K$ 个子载波内实现列正交，从而最大化 $\det(\mathbf V^H\mathbf V)$。

### 命题 2：一般性平滑相位矩阵的 Pareto 前沿包含 CDD

由于 CDD 是平滑相位矩阵的特殊情形，若允许一般性 $\phi_{k,n}$ 并加入同等平滑性约束，则其可行域包含 CDD。因此其最优 Pareto 前沿不劣于 CDD。

### 命题 3：有限维系统中一般性 V 可能严格优于 CDD

在有限 $K,N$ 下，CDD 的列相关由 Dirichlet kernel 决定，旁瓣结构固定。分段线性、chirp 或优化相位可以改变 off-diagonal 项的相位分布，使其在全带累加时更好抵消。因此在相同 group delay spread 下，可能获得更低 mutual coherence 或更大的最小奇异值。

### 命题 4：分集增益与 CE 处理增益存在基本折中

若 $\mathbf V$ 的列要正交，则不同列的相对相位必须在频域内充分展开；这等价于引入人工 delay spread 或 group delay spread，从而降低等效相干带宽。因此不存在同时无限提高分集展开和 CE 处理增益的设计。

---

## 15. 主要结论

1. 传统 CDD 对应全带线性相位，是一个很自然、可解释、易实现的 $\mathbf V$ 设计。
2. 从固定 delay 约束看，CDD 可以接近最优，尤其当 delay 取在 DFT grid 上时，$\mathbf V^H\mathbf V$ 可以接近 $K\mathbf I$。
3. 但 CDD 不是一般性相位矩阵空间中的全局最优，因为它只有 $N$ 个 delay 自由度，而一般性 $\mathbf V$ 有 $K\times N$ 个相位自由度。
4. 分集增益和信道估计处理增益之间存在天然折中：列正交需要相位快速变化，频域平滑需要相位慢变。
5. 分段线性相位、chirp 相位和平滑优化相位有可能在有限维实际系统中获得比 CDD 更优的 Pareto 前沿。
6. 最值得优先尝试的是连续分段线性相位，因为它保留了 CDD 的可解释性，同时增加了频域局部设计自由度。
7. 判断新 $\mathbf V$ 是否真的优于 CDD，不能只看 $\log\det(\mathbf V^H\mathbf V)$，还必须同时看 CE NMSE、BLER、throughput、有效相干带宽和 covariance mismatch 鲁棒性。

---

## 16. 建议下一步工作

建议下一阶段按照如下路径推进：

1. 先在纯矩阵层面比较不同 $\mathbf V$ 的 $\log\det$、最小奇异值、condition number 和 mutual coherence；
2. 再把 $\mathbf V$ 放入等效信道模型，计算 $\mathbf R_g=\mathbf V\mathbf R_h\mathbf V^H$，比较 CE NMSE；
3. 然后接入当前 CDD LLS 平台，在相同 DMRS density 和相同 total transmit power 下比较 BLER；
4. 最后扫描 physical delay spread、人工 group delay spread、segment size，画出 Pareto 曲线。

最关键的验证图是：

$$
\text{CE NMSE} \quad \text{vs.} \quad \log\det(\mathbf V^H\mathbf V),
$$

以及

$$
\text{BLER/throughput} \quad \text{vs.} \quad \tau_{\mathrm{art,rms}}.
$$

如果连续分段线性相位或 chirp 相位在这些图上支配全带 CDD，则可以说明：

> 全带 CDD 不是分集增益和信道估计处理增益折中的一般最优方案；更一般的平滑频域相位矩阵 $\mathbf V$ 可以带来更优的有限维 Pareto 前沿。
