# QC 显式 CDD 链路级仿真复现与结构化信道估计算法对比计划书

版本：v0.1  
日期：2026-06-08  
目标读者：链路级仿真平台开发人员、算法研究人员  

---

## 0. 背景与总体目标

QC 提案关注 6GR/NR 演进中的单层 transmit diversity。核心背景是：NR 没有显式定义类似 LTE SFBC 的 transmit diversity，而是默认透明方案，例如 PRG-level precoder cycling 或 transparent small-delay CDD，可以提供一定分集增益。QC 提案进一步指出，PRG-level precoder cycling 在窄带资源分配下可用的 distinct precoder 数量有限，并且跨 PRG 不保证相位连续，因此 UE 通常只能做 PRG-level narrowband channel estimation；CDD 的循环时延在频域等价于线性相位斜坡，可以看成 phase-continuous RE-level precoder cycling，因此允许 UE 使用 wideband DMRS channel estimation。提案还强调，若 UE 知道 cyclic delay value，则可以正确建模 CDD 引起的 PDP shift，从而避免 RMMSE 信道估计的 PDP mismatch。

本计划书的目标是：

1. 在 QC 提案的主要链路级仿真条件下复现显式 CDD 相对 PRG-level precoder cycling 的 BLER 增益趋势；
2. 在相同 TDL 信道、相同 PDSCH/DMRS 资源配置下，对比三类基于已知 CDD delay 的信道估计算法：
   - **Direct equivalent-channel RMMSE**：直接估计 CDD 后的等效信道；
   - **Structural / deterministic reconstruction**：利用已知 CDD delay，先估计底层物理分支信道，再重构 CDD 等效信道；
   - **Non-CDD per-port DMRS physical-channel estimation**：DMRS 不施加 CDD，按 Tx/effective port 正交发送和估计底层端口信道，再用已知 CDD delay 合成数据 RE 上的 CDD 等效信道；
3. 对不同底层信道时延扩展、不同 CDD delay、不同 DMRS 频域密度进行公平对比；
4. 用 BLER-vs-SNR 曲线作为主性能指标，并加入接收端理想 CSI 作为性能上界；
5. 第一版暂不考虑 UE 移动性和时变信道；后续扩展到 CDL 信道、宽带预编码、以及在宽带预编码后的 effective ports 上施加 CDD。

---

## 1. 需要复现和扩展的 QC 仿真要点

### 1.1 QC 提案中的主要仿真条件

根据提案文本，QC CDD evaluation 的主要 LLS 设置如下：

| 项目 | QC 设置 |
|---|---|
| Carrier BW | 100 MHz |
| SCS | 30 kHz |
| Channel model | TDL |
| Delay spread | 30 ns / 100 ns |
| UE speed | 3 km/h / 60 km/h |
| gNB array | 2Tx / 4Tx |
| UE array | 4Rx |
| PDSCH F/TDRA | 8 RB 或 48 RB × 10 symbols，含 2 个 DMRS symbols |
| MCS | MCS 8，16QAM，R = 553/1024，256QAM table |
| DMRS channel estimation | RMMSE-based，PRG bundling size = 4 或 wideband |
| Baseline | 4-RB PRG-level QPSK-based precoder cycling |
| CDD delay | 选择使 10% BLER 所需 SNR 最小的 cyclic delay |

第一版实现中，移动性先关闭，即不引入 Doppler / time variation；后续再补齐 QC 的 3 km/h 和 60 km/h 场景。

### 1.2 第一版必须复现的 QC 现象

第一版需要复现以下趋势，不要求数值完全等同 QC 提案图，但曲线关系和增益来源应一致：

1. **8 RB 窄带分配**：PRG-level precoder cycling 只有 2 个 4-RB PRG，因此 distinct precoder 数量有限；CDD 即使只用 4-RB channel estimation，也应体现明显 diversity gain。
2. **48 RB 宽带分配**：PRG-level precoder cycling 可以让 4 个 precoder 在 48 RB 内循环多次，diversity 已经比较充分；CDD + 4-RB CE 相对 PRG cycling 的增益有限。
3. **CDD + wideband CE**：在 48 RB 场景中，CDD 由于相位连续，可以从 4-RB processing 扩展到 48-RB processing，RMMSE 信道估计获得额外 frequency-domain processing gain。
4. **known delay vs unknown delay**：若 UE 知道 CDD delay，可构造 CDD-shifted PDP/covariance；若不知道，则使用 non-CDD PDP 造成 mismatch，低 delay spread 场景下退化应更明显。
5. **ideal CSI 上界**：接收机直接使用真实等效信道时，应给出所有估计算法的 BLER 性能上界。

---

## 2. 系统模型

### 2.1 OFDM 与资源配置

设系统使用 OFDM，子载波间隔为

\[
\Delta f = 30\ \mathrm{kHz}.
\]

PDSCH 频域分配大小为

\[
N_{\mathrm{RB}} \in \{8,48\},
\]

每个 RB 含 12 个子载波，因此 active subcarrier 数为

\[
N_{\mathrm{sc}} = 12 N_{\mathrm{RB}}.
\]

PDSCH 时域长度为 10 个 OFDM symbols，其中 2 个 symbols 用于 DMRS。第一版中，时域信道近似静态，因此主要关注频域信道估计；2 个 DMRS symbols 可以作为重复观测来降低噪声，也可以在实现中先合并为等效频域 pilot observation。

### 2.2 天线与层数

第一版采用单层 PDSCH：

\[
N_{\mathrm{layer}} = 1.
\]

gNB transmit branches：

\[
N_t \in \{2,4\},
\]

UE receive antennas：

\[
N_r = 4.
\]

这里的 \(N_t=2/4\) 不必解释为真实基站只有 2 或 4 根物理天线。在后续 massive MIMO 扩展中，可以解释为大规模阵列经宽带 beam/basis 预编码降维后的 effective transmit ports，例如 panel、polarization、beam 或 basis ports。

### 2.3 TDL 信道模型

第一版使用 TDL 信道。每个 Tx-Rx 分支的时域离散信道表示为

\[
h_{r,m}[n] = \sum_{\ell=0}^{L_h-1} h_{r,m,\ell}\,\delta[n-\ell],
\]

其中：

- \(r=0,\ldots,N_r-1\)：接收天线索引；
- \(m=0,\ldots,N_t-1\)：发射分支索引；
- \(\ell\)：离散 delay tap 索引；
- \(L_h\)：有效 delay support；
- \(p[\ell]=\mathbb E[|h_{r,m,\ell}|^2]\)：TDL PDP。

频域信道为

\[
H_{r,m}[k]
=
\sum_{\ell=0}^{L_h-1}h_{r,m,\ell}e^{-j2\pi k\ell/N_{\mathrm{FFT}}}.
\]

底层 delay spread 至少包含 QC 的 30 ns 和 100 ns；建议扩展为

\[
\tau_{\mathrm{DS}} \in \{10,30,100,300\}\ \mathrm{ns},
\]

用于观察 reconstruction 方法在不同物理信道平滑程度下的优势或失效条件。

---

## 3. 发射方案

### 3.1 单层 CDD 发射

CDD 对不同 transmit branches 施加不同 cyclic delay：

\[
d_m \in \mathbb Z,
\quad m=0,\ldots,N_t-1,
\]

单位为 OFDM waveform sampling period。频域中，第 \(m\) 个分支的 CDD 相位为

\[
\alpha_m[k]
=
e^{-j2\pi k d_m/N_{\mathrm{FFT}}}.
\]

单层发射向量为

\[
\mathbf x[k]
= \mathbf c_{\mathrm{CDD}}[k]s[k],
\]

其中

\[
\mathbf c_{\mathrm{CDD}}[k]
=\frac{1}{\sqrt{N_t}}
\begin{bmatrix}
 e^{-j2\pi k d_0/N_{\mathrm{FFT}}}\\
 e^{-j2\pi k d_1/N_{\mathrm{FFT}}}\\
 \vdots\\
 e^{-j2\pi k d_{N_t-1}/N_{\mathrm{FFT}}}
\end{bmatrix}.
\]

接收端第 \(r\) 根天线看到的 CDD 后等效信道为

\[
g_r[k]
=
\mathbf H_r[k]\mathbf c_{\mathrm{CDD}}[k]
=
\frac{1}{\sqrt{N_t}}
\sum_{m=0}^{N_t-1}H_{r,m}[k]e^{-j2\pi k d_m/N_{\mathrm{FFT}}},
\]

其中

\[
\mathbf H_r[k]=[H_{r,0}[k],\ldots,H_{r,N_t-1}[k]].
\]

接收信号为

\[
y_r[k] = g_r[k]s[k] + w_r[k].
\]

向量形式为

\[
\mathbf y[k] = \mathbf g[k]s[k] + \mathbf w[k],
\]

其中

\[
\mathbf g[k]=[g_0[k],\ldots,g_{N_r-1}[k]]^T.
\]

### 3.2 PRG-level QPSK-based precoder cycling baseline

QC 提案只明确说明 baseline 是 4-RB PRG-level QPSK-based precoder cycling，并未给出具体 precoder 矩阵。因此第一版需要实现为可配置 codebook。默认实现建议如下。

PRG size：

\[
N_{\mathrm{PRG}}=4\ \mathrm{RB}.
\]

第 \(b\) 个 PRG 内使用固定 precoder：

\[
\mathbf c[k]=\mathbf w_{b\bmod N_p},\quad k\in\mathrm{PRG}\ b.
\]

#### 2Tx 默认 QPSK relative-phase codebook

\[
\mathcal W_{2\mathrm{Tx}}
=\left\{
\frac{1}{\sqrt 2}
\begin{bmatrix}1\\1\end{bmatrix},
\frac{1}{\sqrt 2}
\begin{bmatrix}1\\j\end{bmatrix},
\frac{1}{\sqrt 2}
\begin{bmatrix}1\\-1\end{bmatrix},
\frac{1}{\sqrt 2}
\begin{bmatrix}1\\-j\end{bmatrix}
\right\}.
\]

#### 4Tx 默认 QPSK DFT-like codebook

\[
\mathcal W_{4\mathrm{Tx}}
=\left\{
\frac{1}{2}
\begin{bmatrix}1\\1\\1\\1\end{bmatrix},
\frac{1}{2}
\begin{bmatrix}1\\j\\-1\\-j\end{bmatrix},
\frac{1}{2}
\begin{bmatrix}1\\-1\\1\\-1\end{bmatrix},
\frac{1}{2}
\begin{bmatrix}1\\-j\\-1\\j\end{bmatrix}
\right\}.
\]

48 RB 场景下共有

\[
48/4=12
\]

个 PRG。如果使用 4 个 precoder，则 4 个 precoder 在 48 RB 内循环 3 次。

> 开发要求：baseline codebook、PRG size、cycling order 必须作为配置项输出到结果 JSON/CSV 中，避免后续与 QC 精确复现时出现不可追踪差异。

---

## 4. DMRS 设计与开销公平性

### 4.1 DMRS 时域配置

与 QC 对齐，PDSCH 时域长度为 10 symbols，2 个 DMRS symbols。第一版中不考虑时变信道，可将两个 DMRS symbols 的 LS 估计做平均，得到每个 DMRS frequency position 上的等效 LS 观测：

\[
\tilde g_r[p] = g_r[p] + n_r[p].
\]

若两个 DMRS symbols 独立噪声平均，则等效噪声方差为

\[
\sigma_{\mathrm{LS,eff}}^2 = \frac{\sigma_w^2}{N_{\mathrm{DMRS,sym}} |x_{\mathrm{DMRS}}|^2}.
\]

### 4.2 DMRS 频域密度

为了公平比较 DMRS overhead，必须支持可配置 DMRS frequency density。定义 DMRS 频域间隔：

\[
S_f \in \{2,3,4,6,12\}
\]

单位为 subcarrier。则频域 pilot set 为

\[
\mathcal P = \{p: p\equiv p_0 \pmod {S_f}\}.
\]

对应每 RB 的 pilot 数约为

\[
N_{\mathrm{pilot/RB}}=12/S_f.
\]

总 DMRS RE 数为

\[
N_{\mathrm{DMRS,RE}} = N_{\mathrm{DMRS,sym}}\cdot |\mathcal P|.
\]

DMRS overhead 为

\[
\eta_{\mathrm{DMRS}}
=\frac{N_{\mathrm{DMRS,RE}}}{N_{\mathrm{RB}}\cdot 12\cdot N_{\mathrm{PDSCH,sym}}}.
\]

### 4.3 开销公平比较原则

主 BLER 对比建议采用以下原则：

1. 对 Direct RMMSE 和 CDD-combined structural reconstruction 使用**相同 CDD-DMRS pattern**；
2. 对 Non-CDD per-port DMRS 方法，必须单独标明公平性模式：
   - equal-total-overhead：总 DMRS RE 与 CDD-DMRS 方法相同，但每端口 pilot density 降低；
   - equal-per-port-density：每端口 pilot density 与 CDD-DMRS 方法相同，但总 DMRS overhead 增加；
3. 对某一个 DMRS density / overhead 配置，所有方法使用相同数据 RE 集合或明确扣除不同 DMRS overhead 后比较 net goodput；
4. 固定 MCS index、modulation order 和目标码率，例如 MCS 8、16QAM、R=553/1024；
5. TBS 依据可用 data RE 数重新计算，使不同 DMRS overhead 下保持近似相同 nominal code rate，而不是固定 payload 后被动抬高实际码率；
6. 若要研究“相同 payload 下增加 DMRS overhead 的代价”，另设 secondary experiment。

主结果使用 BLER；建议同时输出 goodput / net spectral efficiency，避免仅看 BLER 而忽略 DMRS overhead 代价：

\[
\mathrm{Goodput} = (1-\mathrm{BLER})\cdot N_{\mathrm{info,bits}}/T_{\mathrm{slot}}.
\]

---

## 5. 接收机与解调

### 5.1 Ideal CSI 上界

Ideal CSI 接收机直接使用真实 CDD 等效信道：

\[
\mathbf g[k] = \mathbf H[k]\mathbf c[k].
\]

单层 LMMSE/MRC equalizer 可写为

\[
\hat s[k]
=\frac{\hat{\mathbf g}^H[k]}{\|\hat{\mathbf g}[k]\|^2+\sigma_w^2/E_s}\mathbf y[k].
\]

Ideal CSI 中

\[
\hat{\mathbf g}[k]=\mathbf g[k].
\]

该曲线作为所有非理想信道估计算法的 BLER 性能上界。

### 5.2 非理想 CSI 接收机

对 RMMSE 和 reconstruction 方法，先由 DMRS 得到

\[
\hat{\mathbf g}[k],\quad k\in\mathcal D,
\]

再用相同 equalizer、LLR 计算、LDPC 解码流程评估 BLER。所有 CE 方法的后端 equalization / LLR / decoding 必须完全相同，只替换 channel estimator。

---

## 6. 信道估计算法一：Direct equivalent-channel RMMSE

### 6.1 估计对象

Direct RMMSE 直接估计 CDD 后的等效信道 \(g_r[k]\)，不尝试恢复每个发射分支的物理信道 \(H_{r,m}[k]\)。对每根 RX 天线独立执行相同估计。

DMRS LS 观测为

\[
\tilde{\mathbf g}_{r,\mathcal P}
=
\mathbf g_{r,\mathcal P}+\mathbf n_{r,\mathcal P}.
\]

目标数据频点集合为 \(\mathcal D\)。RMMSE 估计器为

\[
\hat{\mathbf g}_{r,\mathcal D}
=
\mathbf R_{\mathcal D\mathcal P}^{\mathrm{CDD}}
\left(
\mathbf R_{\mathcal P\mathcal P}^{\mathrm{CDD}}
+
\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1}
\tilde{\mathbf g}_{r,\mathcal P}.
\]

其中

\[
\mathbf R_{\mathcal D\mathcal P}^{\mathrm{CDD}}
=\mathbb E[\mathbf g_{\mathcal D}\mathbf g_{\mathcal P}^H],
\]

\[
\mathbf R_{\mathcal P\mathcal P}^{\mathrm{CDD}}
=\mathbb E[\mathbf g_{\mathcal P}\mathbf g_{\mathcal P}^H].
\]

这里的 \(R\) 表示 correlation/covariance matrix。若实现称为 RMMSE，R 可理解为 robust / regularized / reduced-complexity 的工程实现变体，但数学核心是基于频域相关矩阵的 MMSE 信道估计。

### 6.2 CDD-shifted PDP 与频域相关矩阵

若不同 Tx branches 的底层信道相互独立，且其 PDP 为 \(p_m[n]\)，则 CDD 后等效 PDP 近似为

\[
p_g[n]
=\frac{1}{N_t}\sum_{m=0}^{N_t-1}p_m[n-d_m].
\]

若各 Tx branches 具有相同 PDP：

\[
p_m[n]=p[n],
\]

则

\[
p_g[n]
=\frac{1}{N_t}\sum_{m=0}^{N_t-1}p[n-d_m].
\]

频域相关函数为

\[
R_G[k,k']
=
\mathbb E[G[k]G^*[k']]
=
\sum_n p_g[n]e^{-j2\pi(k-k')n/N_{\mathrm{FFT}}}.
\]

如果信道在频域统计平稳，可写成

\[
R_G[\Delta k]
=\sum_n p_g[n]e^{-j2\pi\Delta k n/N_{\mathrm{FFT}}}.
\]

然后根据 \(\mathcal P\)、\(\mathcal D\) 抽取子矩阵 \(\mathbf R_{\mathcal D\mathcal P}\)、\(\mathbf R_{\mathcal P\mathcal P}\)。

### 6.3 4-RB RMMSE 与 wideband RMMSE

#### 4-RB RMMSE

将 PDSCH 频域分配切分为多个 4-RB bundles：

\[
\mathcal B=\bigcup_i \mathcal B_i,
\quad |\mathcal B_i|=4\ \mathrm{RB}.
\]

对每个 bundle 独立估计：

\[
\hat{\mathbf g}_{r,\mathcal D_i}
=
\mathbf R_{\mathcal D_i\mathcal P_i}
\left(
\mathbf R_{\mathcal P_i\mathcal P_i}+\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1}
\tilde{\mathbf g}_{r,\mathcal P_i}.
\]

这对应 QC 曲线命名中的 CDD_4rb 或 narrowband CE。

#### Wideband RMMSE

对整个 PDSCH 分配内的 DMRS 统一处理：

\[
\mathcal P=\bigcup_i\mathcal P_i,
\quad
\mathcal D=\bigcup_i\mathcal D_i.
\]

估计器为

\[
\hat{\mathbf g}_{r,\mathcal D}
=
\mathbf R_{\mathcal D\mathcal P}
\left(
\mathbf R_{\mathcal P\mathcal P}+\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1}
\tilde{\mathbf g}_{r,\mathcal P}.
\]

这对应 QC 的 CDD_wideband。48 RB 场景下，它表示 channel-estimation processing window 从 4 RB 扩展到 48 RB，并不是 PDSCH 分配从 4 RB 变成 48 RB。

### 6.4 Matched 与 mismatched RMMSE

#### Known CDD delay, matched covariance

若 UE 知道 \(d_m\)，则使用 CDD-shifted PDP 构造

\[
\mathbf R^{\mathrm{CDD}}.
\]

该方法称为 matched RMMSE。

#### Unknown CDD delay, mismatched covariance

若 UE 不知道 \(d_m\)，则只能使用 TRS / non-CDD PDP 构造

\[
\mathbf R^{\mathrm{TRS}}.
\]

但真实 DMRS/PDSCH 等效信道服从 CDD-shifted covariance。此时估计器为

\[
\hat{\mathbf g}_{r,\mathcal D}^{\mathrm{mismatch}}
=
\mathbf R_{\mathcal D\mathcal P}^{\mathrm{TRS}}
\left(
\mathbf R_{\mathcal P\mathcal P}^{\mathrm{TRS}}+\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1}
\tilde{\mathbf g}_{r,\mathcal P}.
\]

### 6.5 MSE 量化表达式

Matched LMMSE/RMMSE 的误差协方差为

\[
\mathbf C_e
=
\mathbf R_{\mathcal D\mathcal D}^{\mathrm{CDD}}
-
\mathbf R_{\mathcal D\mathcal P}^{\mathrm{CDD}}
\left(
\mathbf R_{\mathcal P\mathcal P}^{\mathrm{CDD}}+\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1}
\mathbf R_{\mathcal P\mathcal D}^{\mathrm{CDD}}.
\]

总 MSE 定义为目标频点集合上的估计误差平方和均值：

\[
\mathrm{MSE}_{\mathrm{total}}
=
\mathbb E[\|\mathbf g_{\mathcal D}-\hat{\mathbf g}_{\mathcal D}\|^2]
=\mathrm{tr}(\mathbf C_e).
\]

每频点平均 MSE 为

\[
\mathrm{MSE}_{\mathrm{avg}}
=\frac{1}{|\mathcal D|}\mathrm{tr}(\mathbf C_e).
\]

归一化 NMSE 为

\[
\mathrm{NMSE}
=\frac{\mathrm{tr}(\mathbf C_e)}{\mathrm{tr}(\mathbf R_{\mathcal D\mathcal D})}.
\]

单个目标频点 \(k\) 的 MSE 为

\[
\mathrm{MSE}(k)
=R[k,k]
-\mathbf r_{k\mathcal P}
\left(\mathbf R_{\mathcal P\mathcal P}+\sigma_{\mathrm{LS}}^2\mathbf I\right)^{-1}
\mathbf r_{\mathcal P k}.
\]

当 CDD delay 过大导致 \(\mathbf r_{k\mathcal P}\approx 0\) 时，

\[
\mathrm{MSE}(k)\approx R[k,k],
\]

说明导频对该目标频点提供的信息很少。这也是需要评估 DMRS frequency density 的原因。

---

## 7. 信道估计算法二：Structural / deterministic physical-channel reconstruction

### 7.1 基本思想

Direct RMMSE 直接估计 CDD 后的等效信道：

\[
g_r[k]
=\frac{1}{\sqrt{N_t}}\sum_m H_{r,m}[k]e^{-j2\pi k d_m/N_{\mathrm{FFT}}}.
\]

CDD delay 增大时，\(g_r[k]\) 的频域选择性可能很强，频域相关性下降，RMMSE 插值误差变大。Structural reconstruction 的思想是：CDD 造成的快速频域起伏是已知确定性相位项，因此不应把 \(g_r[k]\) 当作普通未知频率选择性信道直接插值，而应利用已知 \(d_m\) 先恢复较平滑的底层物理分支信道 \(H_{r,m}[k]\)，再重构等效信道。

该方法的目标链路为

\[
\tilde g_r[p]
\rightarrow
\{\hat H_{r,m}[k]\}_{m=0}^{N_t-1}
\rightarrow
\hat g_r[k].
\]

### 7.2 Reconstruction-A：pilot-pair algebraic decoupling，适合作为 v0 baseline

以 2Tx 为例，设

\[
d_0=0,\quad d_1=d.
\]

两个 DMRS 子载波 \(p_1,p_2\) 上的 LS 观测为

\[
\begin{bmatrix}
\tilde g_r[p_1]\\
\tilde g_r[p_2]
\end{bmatrix}
=\frac{1}{\sqrt 2}
\begin{bmatrix}
1 & e^{-j2\pi p_1 d/N_{\mathrm{FFT}}}\\
1 & e^{-j2\pi p_2 d/N_{\mathrm{FFT}}}
\end{bmatrix}
\begin{bmatrix}
H_{r,0}\\
H_{r,1}
\end{bmatrix}
+
\mathbf z_r.
\]

这里假设在 \(p_1,p_2\) 之间底层信道近似不变：

\[
H_{r,0}[p_1]\approx H_{r,0}[p_2]\triangleq H_{r,0},
\]

\[
H_{r,1}[p_1]\approx H_{r,1}[p_2]\triangleq H_{r,1}.
\]

记

\[
\mathbf A(p_1,p_2)
=\frac{1}{\sqrt 2}
\begin{bmatrix}
1 & e^{-j2\pi p_1 d/N_{\mathrm{FFT}}}\\
1 & e^{-j2\pi p_2 d/N_{\mathrm{FFT}}}
\end{bmatrix}.
\]

ZF 解耦为

\[
\hat{\mathbf h}_r
=\mathbf A^{-1}\tilde{\mathbf g}_r.
\]

更推荐 regularized LS：

\[
\hat{\mathbf h}_r
=\left(\mathbf A^H\mathbf A+\lambda\mathbf I\right)^{-1}\mathbf A^H\tilde{\mathbf g}_r,
\]

或局部 MMSE：

\[
\hat{\mathbf h}_r
=\mathbf R_h\mathbf A^H
\left(\mathbf A\mathbf R_h\mathbf A^H+\sigma_{\mathrm{LS}}^2\mathbf I\right)^{-1}
\tilde{\mathbf g}_r.
\]

得到稀疏位置上的 \(\hat H_{r,0}\)、\(\hat H_{r,1}\) 后，对两个分支独立做插值/平滑：

\[
\hat H_{r,m}[k],\quad k\in\mathcal D.
\]

最后重构 CDD 等效信道：

\[
\hat g_r[k]
=\frac{1}{\sqrt 2}
\left(
\hat H_{r,0}[k]
+
\hat H_{r,1}[k]e^{-j2\pi kd/N_{\mathrm{FFT}}}
\right).
\]

#### 条件数诊断

该方法必须输出每个 pilot pair 的条件数。设

\[
\delta=2\pi(p_2-p_1)d/N_{\mathrm{FFT}}.
\]

则 2Tx 解耦矩阵的条件数近似为

\[
\kappa(\mathbf A)
=
\sqrt{\frac{1+|\cos(\delta/2)|}{1-|\cos(\delta/2)|}}.
\]

当 \(\delta\to 0\) 时，\(\kappa\to\infty\)，噪声会被严重放大。因此 adjacent DMRS pair 未必最优；应支持可配置 pair spacing，并记录 flatness error 与 condition number 的 tradeoff。

### 7.3 Reconstruction-B：multi-pilot delay-domain basis LMMSE，推荐作为主重构算法

为了避免 pairwise 解耦的病态问题，推荐第一版实现一个更稳健的 multi-pilot basis reconstruction。

底层物理分支信道写成 delay-domain taps：

\[
H_{r,m}[k]
=\sum_{\ell\in\mathcal L}h_{r,m}[\ell]e^{-j2\pi k\ell/N_{\mathrm{FFT}}}.
\]

CDD 后观测为

\[
\tilde g_r[p]
=\frac{1}{\sqrt{N_t}}
\sum_{m=0}^{N_t-1}\sum_{\ell\in\mathcal L}
h_{r,m}[\ell]
e^{-j2\pi p(\ell+d_m)/N_{\mathrm{FFT}}}
+n_r[p].
\]

将所有 DMRS 位置堆叠为

\[
\tilde{\mathbf g}_{r,\mathcal P}
=\mathbf \Phi\mathbf h_r+\mathbf n_r.
\]

其中

\[
[\mathbf\Phi]_{p,(m,\ell)}
=\frac{1}{\sqrt{N_t}}e^{-j2\pi p(\ell+d_m)/N_{\mathrm{FFT}}},
\]

\[
\mathbf h_r
= [h_{r,0}[\ell_1],\ldots,h_{r,0}[\ell_L],h_{r,1}[\ell_1],\ldots,h_{r,N_t-1}[\ell_L]]^T.
\]

若使用 TDL PDP 作为 tap prior，则

\[
\mathbf R_h=\mathbb E[\mathbf h_r\mathbf h_r^H]
\]

通常可设为 block diagonal：

\[
\mathbf R_h=\mathrm{diag}(p[\ell_1],\ldots,p[\ell_L],\ldots,p[\ell_1],\ldots,p[\ell_L]).
\]

LMMSE tap estimator 为

\[
\hat{\mathbf h}_r
=
\mathbf R_h\mathbf\Phi^H
\left(
\mathbf\Phi\mathbf R_h\mathbf\Phi^H+\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1}
\tilde{\mathbf g}_{r,\mathcal P}.
\]

若不使用精确 PDP，也可使用 ridge estimator：

\[
\hat{\mathbf h}_r
=\left(\mathbf\Phi^H\mathbf\Phi+\lambda\mathbf I\right)^{-1}
\mathbf\Phi^H\tilde{\mathbf g}_{r,\mathcal P}.
\]

重构底层分支信道：

\[
\hat H_{r,m}[k]
=\sum_{\ell\in\mathcal L}\hat h_{r,m}[\ell]e^{-j2\pi k\ell/N_{\mathrm{FFT}}}.
\]

再重构等效信道：

\[
\hat g_r[k]
=\frac{1}{\sqrt{N_t}}\sum_{m=0}^{N_t-1}\hat H_{r,m}[k]e^{-j2\pi kd_m/N_{\mathrm{FFT}}}.
\]

#### Reduced-support basis variants

实验 8 的诊断显示，全 99% PDP support 下的 \(\mathbf\Phi\) 可能严重病态。因此算法二需要显式评估 reduced-support 版本，而不是默认把所有有效 tap 都作为未知量。

设完整候选 tap 集为

\[
\mathcal L_{\mathrm{full}}=\{0,1,\ldots,L_{\mathrm{full}}-1\}.
\]

Reduced-support reconstruction 使用子集

\[
\mathcal L_q\subseteq \mathcal L_{\mathrm{full}},
\quad |\mathcal L_q| \ll L_{\mathrm{full}}.
\]

常用选择包括累计能量阈值：

\[
\mathcal L_q
=\left\{0,\ldots,L_q-1\right\},
\quad
\sum_{\ell=0}^{L_q-1}p[\ell]\ge q,
\]

其中

\[
q\in\{0.90,0.95,0.99\}.
\]

也可以使用固定 tap 数：

\[
\mathcal L_K=\{0,1,\ldots,K-1\},
\quad K\in\{K_1,K_2,\ldots\}.
\]

或者使用最强 tap 位置：

\[
\mathcal L_K^{\mathrm{top}}
=
\underset{\mathcal S:|\mathcal S|=K}{\arg\max}
\sum_{\ell\in\mathcal S}p[\ell].
\]

对任意 reduced support \(\mathcal L_q\)，仍使用

\[
\tilde{\mathbf g}_{r,\mathcal P}
=\mathbf\Phi_q\mathbf h_{r,q}+\mathbf n_r,
\]

其中

\[
[\mathbf\Phi_q]_{p,(m,\ell)}
=
\frac{1}{\sqrt{N_t}}
e^{-j2\pi p(\ell+d_m)/N_{\mathrm{FFT}}},
\quad \ell\in\mathcal L_q.
\]

对应 LMMSE estimator 为

\[
\hat{\mathbf h}_{r,q}
=
\mathbf R_{h,q}\mathbf\Phi_q^H
\left(
\mathbf\Phi_q\mathbf R_{h,q}\mathbf\Phi_q^H
+\sigma_{\mathrm{LS}}^2\mathbf I
+\epsilon\mathbf I
\right)^{-1}
\tilde{\mathbf g}_{r,\mathcal P}.
\]

Reduced-support 的目的不是声称远端 tap 不存在，而是在“模型截断误差”和“矩阵病态/噪声放大”之间寻找更好的折中。实验必须同时输出：

\[
\mathrm{NMSE}_H,\quad
\mathrm{NMSE}_g,\quad
\kappa(\mathbf\Phi_q),\quad
\mathrm{rank}_{\mathrm{eff}}(\mathbf\Phi_q),
\quad |\mathcal L_q|.
\]

### 7.4 Reconstruction-C：coherence-band blockwise physical-channel reconstruction

这是算法二的一个新的局部结构化分支，用来避免 Reconstruction-B 把全带宽所有 DMRS 拼成一个大 \(\mathbf\Phi\) 后产生严重病态。其基本假设是：底层 physical branch channel 的相干带宽明显大于 CDD equivalent channel 的相干带宽。因此在一个不超过底层相干带宽的局部频域块内，可近似认为每个 Tx/effective port 的底层信道为常数。

先由原始 PDP 计算底层物理信道的频域相关函数：

\[
\rho_H[\Delta k]
=
\frac{R_H[\Delta k]}{R_H[0]}
=
\frac{
\sum_\ell p[\ell]e^{-j2\pi\Delta k\ell/N_{\mathrm{FFT}}}
}{
\sum_\ell p[\ell]
}.
\]

定义底层信道相干带宽：

\[
B_{c,H}(\gamma)
=
\min\left\{\Delta k:|\rho_H[\Delta k]|\le \gamma\right\},
\]

例如 \(\gamma=0.9\) 或 \(\gamma=0.5\)。将 active allocation 划分为多个局部块：

\[
\mathcal K
=
\bigcup_b \mathcal K_b,
\quad
|\mathcal K_b|\le B_{c,H}(\gamma).
\]

在第 \(b\) 个块内，对每根 RX 天线 \(r\)，近似：

\[
H_{r,m}[k]\approx h_{r,m}^{(b)},
\quad k\in\mathcal K_b,
\quad m=0,\ldots,N_t-1.
\]

令该块内 DMRS pilot 集为

\[
\mathcal P_b=\mathcal P\cap\mathcal K_b.
\]

CDD-combined pilot observation 写成局部线性模型：

\[
\tilde g_r[p]
=
\frac{1}{\sqrt{N_t}}
\sum_{m=0}^{N_t-1}
h_{r,m}^{(b)}
e^{-j2\pi p d_m/N_{\mathrm{FFT}}}
+n_r[p],
\quad p\in\mathcal P_b.
\]

堆叠后：

\[
\tilde{\mathbf g}_{r,\mathcal P_b}
=
\mathbf A_b\mathbf h_r^{(b)}
+\mathbf n_{r,b},
\]

其中

\[
\mathbf h_r^{(b)}
=
\left[
h_{r,0}^{(b)},\ldots,h_{r,N_t-1}^{(b)}
\right]^T,
\]

\[
[\mathbf A_b]_{p,m}
=
\frac{1}{\sqrt{N_t}}
e^{-j2\pi p d_m/N_{\mathrm{FFT}}},
\quad p\in\mathcal P_b.
\]

若采用 ridge / regularized LS：

\[
\hat{\mathbf h}_r^{(b)}
=
\left(
\mathbf A_b^H\mathbf A_b+\lambda\mathbf I
\right)^{-1}
\mathbf A_b^H
\tilde{\mathbf g}_{r,\mathcal P_b}.
\]

若采用局部 LMMSE，并假设块内各端口 prior 为

\[
\mathbf R_{h,b}=\sigma_h^2\mathbf I,
\]

则

\[
\hat{\mathbf h}_r^{(b)}
=
\mathbf R_{h,b}\mathbf A_b^H
\left(
\mathbf A_b\mathbf R_{h,b}\mathbf A_b^H
+\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1}
\tilde{\mathbf g}_{r,\mathcal P_b}.
\]

得到局部端口信道后，对块内所有 data / target subcarriers 合成 CDD 等效信道：

\[
\hat g_r[k]
=
\frac{1}{\sqrt{N_t}}
\sum_{m=0}^{N_t-1}
\hat h_{r,m}^{(b)}
e^{-j2\pi k d_m/N_{\mathrm{FFT}}},
\quad k\in\mathcal K_b.
\]

该方法的关键点是：每个 block 只估计 \(N_t\) 个未知量，而不是估计 \(N_t|\mathcal L|\) 个 delay-domain taps。它利用了底层信道相干带宽的先验，但不依赖全带宽大矩阵求逆。

该分支必须记录：

1. block 宽度 \(|\mathcal K_b|\)，以及其来源是 \(B_{c,H}(0.9)\)、\(B_{c,H}(0.5)\) 还是人工指定；
2. 每个 block 内 pilot 数 \(|\mathcal P_b|\)；
3. 每个 block 的 \(\kappa(\mathbf A_b)\)；
4. 若 \(|\mathcal P_b|<N_t\) 或 \(\kappa(\mathbf A_b)\) 过大，应标记该 block 不可辨识或依赖强 regularization；
5. 输出 blockwise branch NMSE 和最终 equivalent-channel NMSE。

### 7.5 Reconstruction 方法的实现要求

1. 第一版必须至少实现 Reconstruction-A；推荐同时实现 Reconstruction-B 作为主方法。
2. 对每个 DMRS density、每个 CDD delay，输出 reconstruction matrix 的 condition number 或 effective rank。
3. 输出底层分支信道 NMSE：

\[
\mathrm{NMSE}_{H_m}
=\frac{\sum_k|H_m[k]-\hat H_m[k]|^2}{\sum_k|H_m[k]|^2}.
\]

4. 输出最终等效信道 NMSE：

\[
\mathrm{NMSE}_{g}
=\frac{\sum_k|g[k]-\hat g[k]|^2}{\sum_k|g[k]|^2}.
\]

5. Reconstruction-B 的 delay support \(\mathcal L\) 应可配置：
   - ideal support：使用 TDL 真实 taps；
   - truncated support：只保留累计能量达到 90%、95%、99% 的 taps；
   - fixed-tap-count support：固定使用前 \(K\) 个 tap 或最强 \(K\) 个 tap；
   - mismatched support：用于鲁棒性测试。
6. Reconstruction-C 应可配置 block 宽度来源：
   - `bc_0p9`：使用原始底层物理信道 \(B_{c,H}(0.9)\)；
   - `bc_0p5`：使用原始底层物理信道 \(B_{c,H}(0.5)\)；
   - `manual_sc`：手动指定 block subcarrier 数。
7. Reconstruction-C 不应把全带宽所有 pilot 拼成一个大 \(\mathbf\Phi\)，而应逐 block 独立估计局部端口信道。

---

## 8. 信道估计算法三：Non-CDD per-port DMRS physical-channel estimation

### 8.1 基本思想

前两类算法默认 DMRS 与 PDSCH 数据使用相同的 CDD precoding，因此 UE 在 DMRS 上直接观测到 CDD 后的等效信道：

\[
\tilde g_r[p] = g_r[p] + n_r[p].
\]

Direct RMMSE 直接估计 \(g_r[k]\)。Structural reconstruction 则试图从这个单层 CDD-combined observation 中反推出多个底层物理分支信道 \(H_{r,m}[k]\)。上一轮诊断显示，这个反问题可能严重病态，因为一个 CDD-combined scalar observation 同时混合了多个 Tx/effective ports。

算法三改变 DMRS 设计：**数据 RE 仍然使用 CDD 发射，但 DMRS RE 不施加 CDD，并且每个 Tx/effective port 独立发送正交 DMRS**。UE 在 DMRS 上直接估计底层 physical/effective port channel：

\[
H_{r,m}[k],
\quad m=0,\ldots,N_t-1,
\]

然后在数据 RE 上利用已知 CDD delay 合成真实用于均衡的 CDD 等效信道：

\[
\hat g_r[k]
=
\frac{1}{\sqrt{N_t}}
\sum_{m=0}^{N_t-1}
\hat H_{r,m}[k]
e^{-j2\pi k d_m/N_{\mathrm{FFT}}}.
\]

该算法的核心区别是：DMRS 阶段不再从 CDD-combined observation 中解耦端口，而是通过 port-orthogonal DMRS 直接观测每个端口。因此它应避免 Reconstruction-B 中 \(\mathbf\Phi\) 矩阵病态的问题。它的代价是 DMRS 设计更复杂，并且在相同总 DMRS overhead 下，每个端口可用的 pilot density 可能降低；若保持每端口 pilot density 不变，总 DMRS overhead 会随端口数增加。

### 8.2 DMRS 发射模型：无 CDD 的端口正交导频

设第 \(m\) 个 Tx/effective port 的 DMRS pilot set 为

\[
\mathcal P_m,
\quad m=0,\ldots,N_t-1.
\]

端口正交可以通过 FDM、CDM/OCC、TDM 或它们的组合实现。第一版建议先实现 FDM port-orthogonal DMRS，即不同端口占用不同 pilot RE：

\[
\mathcal P_m \cap \mathcal P_{m'} = \varnothing,
\quad m\ne m'.
\]

若使用 CDM/OCC 正交，则不同端口可以共享相同 RE，但在解扩后得到等效的端口独立 LS 观测。为了统一表达，下文假设接收端已经完成 port de-orthogonalization，得到每个端口的 LS observation。

DMRS 不施加 CDD，因此第 \(m\) 个端口在 DMRS RE \(p\in\mathcal P_m\) 上的接收信号为

\[
y_{r,m}^{\mathrm{DMRS}}[p]
=
H_{r,m}[p]x_m^{\mathrm{DMRS}}[p]
+w_r[p].
\]

若 DMRS 符号满足 \(|x_m^{\mathrm{DMRS}}[p]|=1\)，则 LS 估计为

\[
\tilde H_{r,m}[p]
=
\frac{y_{r,m}^{\mathrm{DMRS}}[p]}{x_m^{\mathrm{DMRS}}[p]}
=
H_{r,m}[p] + n_{r,m}[p].
\]

若两个 DMRS symbols 在静态信道下提供重复观测，则可以对同一端口、同一频点的 LS observation 求平均，等效噪声方差为

\[
\sigma_{\mathrm{LS},m}^2
=
\frac{\sigma_w^2}
{N_{\mathrm{DMRS,sym}} |x_m^{\mathrm{DMRS}}|^2}.
\]

若使用 CDM/OCC 且多个端口同 RE 复用，则应在 \(\sigma_{\mathrm{LS},m}^2\) 中包含解扩后的噪声增强或 CDM orthogonality loss。

### 8.3 Per-port physical-channel RMMSE

对每根 RX 天线 \(r\) 和每个 Tx/effective port \(m\)，独立估计底层信道 \(H_{r,m}[k]\)。目标数据频点集合仍记为 \(\mathcal D\)。

底层物理端口信道的频域相关矩阵由 non-CDD PDP 构造：

\[
R_H[k,k']
=
\mathbb E[H_{r,m}[k]H_{r,m}^*[k']]
=
\sum_n p[n]e^{-j2\pi(k-k')n/N_{\mathrm{FFT}}}.
\]

对端口 \(m\)，RMMSE 估计器为

\[
\hat{\mathbf H}_{r,m,\mathcal D}
=
\mathbf R_{\mathcal D\mathcal P_m}
\left(
\mathbf R_{\mathcal P_m\mathcal P_m}
+\sigma_{\mathrm{LS},m}^2\mathbf I
\right)^{-1}
\tilde{\mathbf H}_{r,m,\mathcal P_m}.
\]

其中

\[
\tilde{\mathbf H}_{r,m,\mathcal P_m}
=
\left[
\tilde H_{r,m}[p]
\right]_{p\in\mathcal P_m}.
\]

如果不同端口具有不同 PDP 或不同 beam-dependent PDP，可使用端口相关的 covariance：

\[
R_{H,m}[k,k']
=
\sum_n p_m[n]e^{-j2\pi(k-k')n/N_{\mathrm{FFT}}}.
\]

第一版 TDL static 场景中，可先假设各端口共享相同 PDP：

\[
p_m[n] = p[n].
\]

### 8.4 数据 RE 上的 CDD 等效信道合成

数据 RE 仍使用 CDD transmit diversity。第 \(m\) 个端口在数据频点 \(k\) 上的 CDD 相位为

\[
\alpha_m[k]
=
e^{-j2\pi k d_m/N_{\mathrm{FFT}}}.
\]

单层 CDD 数据等效信道估计为

\[
\hat g_r[k]
=
\frac{1}{\sqrt{N_t}}
\sum_{m=0}^{N_t-1}
\hat H_{r,m}[k]\alpha_m[k],
\quad k\in\mathcal D.
\]

接收端后续 equalizer、LLR、LDPC decoding 与其他 CE 方法完全相同：

\[
\hat s[k]
=
\frac{\hat{\mathbf g}^H[k]}
{\|\hat{\mathbf g}[k]\|^2+\sigma_w^2/E_s}
\mathbf y[k].
\]

因此算法三与 Direct RMMSE / Structural reconstruction 的对比仍然只改变 channel estimator，不改变后端 detector 和 decoder。

### 8.5 与算法二的区别

算法二的 multi-pilot basis reconstruction 使用 CDD 后的混合观测：

\[
\tilde g_r[p]
=
\frac{1}{\sqrt{N_t}}
\sum_{m=0}^{N_t-1}\sum_{\ell\in\mathcal L}
h_{r,m}[\ell]
e^{-j2\pi p(\ell+d_m)/N_{\mathrm{FFT}}}
+n_r[p],
\]

并试图从 \(\tilde g_r[p]\) 反推出所有 \(h_{r,m}[\ell]\)。这会形成

\[
\tilde{\mathbf g}_{r,\mathcal P}=\mathbf\Phi\mathbf h_r+\mathbf n_r,
\]

其中 \(\mathbf\Phi\) 可能 rank deficient 或 condition number 极大。

算法三则在 DMRS 上直接获得端口分离后的观测：

\[
\tilde H_{r,m}[p] = H_{r,m}[p] + n_{r,m}[p].
\]

因此它不需要从一个 CDD-combined scalar observation 中解耦多个端口，避免了 \(\mathbf\Phi\) 逆问题。它更像是“DMRS 端口正交 + 数据端口 CDD”的显式信令/显式参考信号方案。

### 8.6 DMRS overhead 与公平性

算法三引入了一个重要的公平性问题：port-specific non-CDD DMRS 的资源开销如何与 CDD-combined DMRS 对齐。

建议同时评估两种模式。

#### 模式 A：每端口 pilot density 保持不变

每个端口使用与算法一/二相同的 pilot density：

\[
|\mathcal P_m| = |\mathcal P|,
\quad m=0,\ldots,N_t-1.
\]

总 DMRS RE 数近似为

\[
N_{\mathrm{DMRS,RE}}^{\mathrm{port}}
\approx
N_t
N_{\mathrm{DMRS,RE}}^{\mathrm{CDD-combined}}.
\]

该模式估计质量最好，但 overhead 更高。输出 goodput / net spectral efficiency 时必须扣除额外 DMRS RE。

#### 模式 B：总 DMRS overhead 保持不变

总 DMRS RE 数与算法一/二相同：

\[
\sum_{m=0}^{N_t-1}|\mathcal P_m|
=
|\mathcal P|.
\]

若各端口平均分配 pilot，则

\[
|\mathcal P_m|
\approx
\frac{|\mathcal P|}{N_t}.
\]

该模式 overhead 公平，但每个端口的 pilot density 降低，per-port RMMSE 插值误差可能增大。

两种模式都应输出：

1. `dmrs_overhead_total`：总 DMRS overhead；
2. `dmrs_overhead_per_port`：每端口 pilot overhead；
3. `pilot_density_per_port`：每端口频域 pilot density；
4. `data_re`：扣除 DMRS 后实际 data RE；
5. `goodput` / `net spectral efficiency`：考虑 DMRS overhead 后的净效率。

### 8.7 MSE 与 NMSE 度量

端口信道估计误差定义为

\[
\mathrm{NMSE}_{H}
=
\frac{
\sum_{r,m,k\in\mathcal D}
|H_{r,m}[k]-\hat H_{r,m}[k]|^2
}{
\sum_{r,m,k\in\mathcal D}
|H_{r,m}[k]|^2
}.
\]

合成后的 CDD 等效信道误差定义为

\[
\mathrm{NMSE}_{g}
=
\frac{
\sum_{r,k\in\mathcal D}
|g_r[k]-\hat g_r[k]|^2
}{
\sum_{r,k\in\mathcal D}
|g_r[k]|^2
}.
\]

其中真实 CDD 等效信道为

\[
g_r[k]
=
\frac{1}{\sqrt{N_t}}
\sum_{m=0}^{N_t-1}
H_{r,m}[k]e^{-j2\pi k d_m/N_{\mathrm{FFT}}}.
\]

若各端口估计误差相互独立，等效信道误差方差可近似为

\[
\mathbb E[|g_r[k]-\hat g_r[k]|^2]
\approx
\frac{1}{N_t}
\sum_{m=0}^{N_t-1}
\mathbb E[|H_{r,m}[k]-\hat H_{r,m}[k]|^2].
\]

若端口 DMRS 使用 CDM/OCC 或共享 RE，端口间估计误差可能相关，应在仿真中直接统计 empirical NMSE。

### 8.8 实现要求

该算法建议命名为：

| 模式名 | 描述 |
|---|---|
| `PORT_DMRS_RECON_RMMSE` | DMRS 不施加 CDD，port-specific DMRS 估计底层端口信道，再按已知 CDD delay 合成数据等效信道 |

第一版实现要求：

1. 支持 2Tx/4Tx port-specific non-CDD DMRS；
2. 支持 FDM orthogonal pilot allocation，后续可扩展 CDM/OCC；
3. 支持 overhead 模式 A：每端口 pilot density 不变；
4. 支持 overhead 模式 B：总 DMRS overhead 不变；
5. 每个端口独立 RMMSE / LMMSE interpolation，使用 non-CDD PDP；
6. 数据 RE 上按 CDD delay 合成 \(\hat g_r[k]\)；
7. 输出 `ce_nmse_branch` 和 `ce_nmse_eff`；
8. 输出 per-port pilot density、总 DMRS overhead、data RE 和 goodput；
9. 与 Direct RMMSE / Structural reconstruction 使用相同 LDPC、equalizer、LLR、decoder；
10. 对比时必须明确写明 overhead 模式，避免把更高 DMRS overhead 的结果误解为纯算法增益。

---

## 9. 需要实现的仿真模式

### 9.1 发射方案模式

| 模式名 | 描述 |
|---|---|
| PRG_CYCLING_4RB | 4-RB PRG-level QPSK-based precoder cycling baseline |
| CDD | 显式 CDD，UE 可配置 known / unknown delay |
| NO_CDD | 无 CDD，用于 sanity check |

### 9.2 信道估计模式

| 模式名 | 发射方案适用 | 描述 |
|---|---|---|
| IDEAL_CSI | 所有 | 使用真实等效信道 |
| LS_LINEAR | CDD / NO_CDD | DMRS LS + 线性插值 baseline |
| RMMSE_4RB_KNOWN | CDD | 使用 CDD-shifted PDP，4-RB processing |
| RMMSE_WB_KNOWN | CDD | 使用 CDD-shifted PDP，whole allocation processing |
| RMMSE_4RB_UNKNOWN | CDD | 使用 non-CDD PDP，4-RB processing |
| RMMSE_WB_UNKNOWN | CDD | 使用 non-CDD PDP，whole allocation processing |
| PRG_RMMSE_4RB | PRG_CYCLING_4RB | PRG 内 RMMSE，不跨 PRG |
| RECON_PAIRWISE | CDD | pairwise algebraic decoupling + interpolation + recombination |
| RECON_BASIS_LMMSE | CDD | multi-pilot delay-domain basis LMMSE + recombination |
| PORT_DMRS_RECON_RMMSE | CDD data + non-CDD port DMRS | 每端口正交 DMRS 估计底层端口信道，再按已知 CDD delay 合成数据等效信道 |

---

## 10. 实验矩阵

### 10.1 Phase 0：基础验证

目的：确认 OFDM、TDL、CDD、DMRS LS、equalizer、LDPC/BLER 流程无误。

配置：

| 项目 | 值 |
|---|---|
| Channel | AWGN / flat fading / simple 2-tap TDL |
| Tx/Rx | 2Tx, 4Rx |
| PDSCH | 8 RB, 10 symbols, 2 DMRS symbols |
| MCS | MCS 8 |
| CE | IDEAL_CSI, LS_LINEAR, RMMSE |

验收标准：

1. AWGN + ideal CSI BLER 与理论趋势一致；
2. flat fading 下 CDD 不应产生额外 frequency selectivity；
3. RMMSE 的 empirical NMSE 与理论 MSE 公式趋势一致；
4. reconstruction 在无 CDD 或简单 CDD 下能够重构等效信道。

### 10.2 Phase 1：QC 静态 TDL 趋势复现

配置：

| 参数 | 值 |
|---|---|
| Channel | TDL |
| Delay spread | 30 ns, 100 ns |
| Mobility | static / no Doppler |
| Tx/Rx | 2Tx/4Rx, 4Tx/4Rx |
| PDSCH RB | 8, 48 |
| PDSCH symbols | 10 |
| DMRS symbols | 2 |
| MCS | MCS 8, 16QAM, R=553/1024 |
| SNR sweep | 例如 -5:1:20 dB，按结果调整 |
| BLER target | 10% |

对比曲线：

1. PRG_CYCLING_4RB + PRG_RMMSE_4RB；
2. CDD + RMMSE_4RB_KNOWN；
3. CDD + RMMSE_WB_KNOWN；
4. CDD + RMMSE_4RB_UNKNOWN；
5. CDD + RMMSE_WB_UNKNOWN；
6. CDD + IDEAL_CSI。

输出：

- BLER vs SNR；
- 10% BLER 所需 SNR；
- 相对 PRG_CYCLING_4RB 的 SNR gain；
- CE NMSE vs SNR；
- CDD delay sweep 结果。

预期趋势：

1. 8 RB 下，CDD_4RB_KNOWN 应明显优于 PRG_CYCLING_4RB；
2. 48 RB 下，CDD_4RB_KNOWN 相对 PRG_CYCLING_4RB 增益有限；
3. 48 RB 下，CDD_WB_KNOWN 显著优于 CDD_4RB_KNOWN；
4. UNKNOWN delay 情况低于 KNOWN delay；
5. 30 ns delay spread 下 unknown delay 的退化应比 100 ns 更明显；
6. IDEAL_CSI 为性能上界。

### 10.3 Phase 2：RMMSE vs Reconstruction 对比

配置：

| 参数 | 值 |
|---|---|
| Channel | TDL |
| Delay spread | 10, 30, 100, 300 ns |
| Tx/Rx | 2Tx/4Rx first，后续 4Tx/4Rx |
| PDSCH RB | 48 RB 为主，8 RB 辅助 |
| DMRS spacing | S_f = 2, 3, 4, 6, 12 subcarriers |
| CDD delay | sweep，单位为 samples |
| CE | RMMSE_WB_KNOWN, RECON_PAIRWISE, RECON_BASIS_LMMSE, PORT_DMRS_RECON_RMMSE, IDEAL_CSI |

核心问题：

1. 在底层物理信道较平滑、CDD 造成等效信道强频率选择性时，reconstruction 是否优于 direct RMMSE？
2. 当 DMRS 频域密度降低时，reconstruction 是否比 direct RMMSE 更能利用已知 CDD delay？
3. 当 CDD delay 过大导致 CDD-combined reconstruction 解耦矩阵病态时，port-specific non-CDD DMRS 是否能通过直接观测底层端口信道获得收益？
4. 当 CDD delay 过大导致 pairwise 解耦矩阵病态时，reconstruction 是否退化？
5. Reconstruction-B 是否比 Reconstruction-A 更稳健？
6. 增加 DMRS density 对 direct RMMSE、CDD-combined reconstruction 和 port-specific DMRS reconstruction 的收益是否不同？

输出：

- BLER vs SNR；
- 10% BLER SNR gain；
- equivalent-channel NMSE；
- physical-branch NMSE；
- reconstruction condition number；
- DMRS overhead vs 10% BLER SNR；
- DMRS overhead vs goodput。

### 10.4 Phase 3：完整 QC 移动性补齐

在 Phase 1/2 成熟后加入 UE speed：

\[
v\in\{3,60\}\ \mathrm{km/h}.
\]

需要处理：

1. DMRS symbols 间信道时变；
2. time-domain interpolation / filtering；
3. Doppler covariance；
4. CDD delay 与 Doppler 对 RMMSE covariance 的联合影响。

### 10.5 Phase 4：CDL 与宽带预编码扩展

后续扩展模型：

\[
\mathbf H_{\mathrm{phys}}[k]
\in\mathbb C^{N_r\times N_{\mathrm{ant}}}.
\]

宽带预编码：

\[
\mathbf W_{\mathrm{WB}}
\in\mathbb C^{N_{\mathrm{ant}}\times P}.
\]

降维后的 effective channel：

\[
\mathbf H_{\mathrm{eff}}[k]
=\mathbf H_{\mathrm{phys}}[k]\mathbf W_{\mathrm{WB}}.
\]

在 effective ports 上施加 CDD：

\[
\mathbf y[k]
=\mathbf H_{\mathrm{eff}}[k]
\mathbf D_{\mathrm{CDD}}[k]
\mathbf v s[k]+\mathbf n[k].
\]

这里

\[
\mathbf D_{\mathrm{CDD}}[k]
=\mathrm{diag}(e^{-j2\pi kd_0/N_{\mathrm{FFT}}},\ldots,e^{-j2\pi kd_{P-1}/N_{\mathrm{FFT}}}).
\]

该阶段研究 CDL 下 beam-dependent PDP、port correlation、wideband beam squint 对 CDD/RMMSE/reconstruction 的影响。

---

## 11. CDD delay 选择方法

QC 提案中 CDD delay 选择为使 10% BLER 所需 SNR 最小。第一版建议实现两种模式。

### 11.1 Delay sweep 模式

对每组场景，枚举

\[
d_m\in\mathcal D_{\mathrm{CDD}}.
\]

2Tx 默认：

\[
[d_0,d_1]=[0,d].
\]

4Tx 默认：

\[
[d_0,d_1,d_2,d_3]=[0,d,2d,3d]
\]

或 configurable delay vector。

对每个 delay vector 运行 BLER-vs-SNR，记录：

\[
\mathrm{SNR}_{10\%}(d).
\]

选取：

\[
d^*=\arg\min_d \mathrm{SNR}_{10\%}(d).
\]

### 11.2 Fixed delay 模式

用于算法对比时固定 \(d\)，避免每种 CE 算法都选择不同最优 delay 导致不公平。推荐流程：

1. 先用 RMMSE_WB_KNOWN 选择 QC-style \(d^*\)；
2. 然后所有 CE 算法使用同一个 \(d^*\)；
3. 另做 per-algorithm optimal delay 作为辅助结果。

---

## 12. 结果统计与绘图要求

### 12.1 主图

1. **BLER vs SNR**：每个 scenario 一张图；
2. **10% BLER SNR gain bar chart**：相对 PRG_CYCLING_4RB；
3. **CE NMSE vs SNR**：对比 RMMSE、reconstruction、LS；
4. **DMRS density sweep**：横轴 DMRS overhead，纵轴 10% BLER SNR 或 goodput；
5. **CDD delay sweep**：横轴 cyclic delay samples，纵轴 10% BLER SNR；
6. **Condition number heatmap**：横轴 CDD delay，纵轴 DMRS pair spacing / DMRS density。

### 12.2 表格

每个 scenario 输出 CSV：

| 字段 | 含义 |
|---|---|
| scenario_id | 唯一场景 ID |
| channel_model | TDL/CDL |
| delay_spread_ns | delay spread |
| speed_kmh | UE speed |
| n_tx | Tx branches/effective ports |
| n_rx | Rx antennas |
| pdsch_rb | 8/48 |
| dmrs_spacing_sc | DMRS 频域间隔 |
| dmrs_overhead | DMRS overhead |
| cdd_delay_vector | CDD delay samples |
| tx_scheme | CDD/PRG/NO_CDD |
| ce_method | IDEAL/RMMSE/RECON |
| snr_db | SNR |
| n_trials | trials 数 |
| bler | BLER |
| ce_nmse_eff | 等效信道 NMSE |
| ce_nmse_branch | 分支信道 NMSE，可选 |
| cond_number | reconstruction condition number，可选 |
| tbs_bits | TBS |
| data_re | data RE 数 |
| goodput | goodput |

---

## 13. 信道估计算法复杂度对比

本节比较三类核心信道估计算法的计算复杂度、存储复杂度和实现风险。复杂度分为两类：

1. **预计算复杂度**：由 resource allocation、DMRS pattern、PDP/covariance、CDD delay、SNR/noise variance 决定，通常每个 scenario / SNR / delay 配置计算一次；
2. **每 trial 在线复杂度**：对每个 channel realization、每根 RX 天线的 LS observation 执行矩阵乘法、插值、recombination 等操作。

### 13.1 符号定义

| 符号 | 含义 |
|---|---|
| \(N_{\mathrm{sc}}\) | active subcarrier 数，例如 48 RB 时为 576 |
| \(N_D\) | 需要估计的目标频点数，NMSE-only 频域实验中通常 \(N_D=N_{\mathrm{sc}}\) |
| \(N_P\) | CDD-combined DMRS pilot 子载波数，多个 DMRS symbols 平均后按频域 pilot 数计 |
| \(N_{P,m}\) | 算法三中第 \(m\) 个 port 的 pilot 子载波数 |
| \(N_t\) | Tx branches / effective ports 数 |
| \(N_r\) | RX antennas 数 |
| \(L_g\) | CDD equivalent PDP 的有效 support 长度 |
| \(L\) | 算法二 physical-channel delay-domain basis support 长度 |
| \(B\) | 算法二 blockwise reconstruction 的频域 block 数 |
| \(K_b\) | 第 \(b\) 个 block 的目标子载波数 |
| \(P_b\) | 第 \(b\) 个 block 内 pilot 数 |

下文假设矩阵求解使用 Cholesky / Hermitian solver 或 `solve()`，不显式求逆。若为了在线加速预先形成 interpolation matrix，则会增加存储但降低 per-trial 计算。

### 13.2 算法一：Direct equivalent-channel RMMSE

算法一直接估计 CDD 后的 equivalent channel：

\[
\hat{\mathbf g}_{\mathcal D}
=
\mathbf R_{\mathcal D\mathcal P}
\left(
\mathbf R_{\mathcal P\mathcal P}
+\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1}
\tilde{\mathbf g}_{\mathcal P}.
\]

#### Wideband RMMSE

预计算：

1. 构造 covariance：

\[
\mathbf R_{\mathcal P\mathcal P}\in\mathbb C^{N_P\times N_P},
\quad
\mathbf R_{\mathcal D\mathcal P}\in\mathbb C^{N_D\times N_P}.
\]

若直接由 PDP 求和，复杂度约为

\[
O(N_P^2L_g+N_DN_PL_g).
\]

2. 分解 / 求解 pilot covariance：

\[
O(N_P^3).
\]

3. 若显式形成 weight matrix

\[
\mathbf W_{\mathrm{RMMSE}}
=
\mathbf R_{\mathcal D\mathcal P}
\left(
\mathbf R_{\mathcal P\mathcal P}
+\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1},
\]

额外复杂度约为

\[
O(N_DN_P^2).
\]

每 trial 在线复杂度：

\[
\hat{\mathbf g}_{r,\mathcal D}
=
\mathbf W_{\mathrm{RMMSE}}\tilde{\mathbf g}_{r,\mathcal P},
\quad r=0,\ldots,N_r-1.
\]

因此

\[
O(N_rN_DN_P).
\]

存储复杂度：

\[
O(N_P^2+N_DN_P).
\]

#### 4-RB / bundle RMMSE

若按 bundle 独立处理，设第 \(q\) 个 bundle 有 \(N_{P,q}\) 个 pilot、\(N_{D,q}\) 个目标频点，则预计算复杂度为

\[
\sum_q O(N_{P,q}^3+N_{D,q}N_{P,q}^2),
\]

每 trial 在线复杂度为

\[
\sum_q O(N_rN_{D,q}N_{P,q}).
\]

4-RB RMMSE 的复杂度和存储低于 wideband RMMSE，但不能利用跨 bundle 的 pilot correlation。

#### 特点

算法一的复杂度与 \(N_t\) 不直接成正比，因为 CDD 后多个 Tx branch 已经被折叠为一个 equivalent PDP / covariance。只要 UE 知道 CDD delay 和 PDP，算法一是三类算法中实现最直接、数值最稳健的 matched LMMSE baseline。

### 13.3 算法二：Structural / deterministic reconstruction

算法二有多个实现分支，其复杂度差异很大。

#### Reconstruction-A：pairwise / local algebraic decoupling

每个 local pilot group 解

\[
\tilde{\mathbf g}_{r,\mathcal P_i}
=
\mathbf A_i\mathbf h_{r,i}+\mathbf n_i,
\quad
\mathbf A_i\in\mathbb C^{N_t\times N_t}
\]

或 overdetermined local LS。若共有 \(N_{\mathrm{grp}}\) 个 local groups，则预计算每个 group 的 regularized inverse：

\[
O(N_{\mathrm{grp}}N_t^3).
\]

每 trial 在线复杂度：

\[
O(N_rN_{\mathrm{grp}}N_t^2)
\]

加上 branch interpolation 和 CDD recombination：

\[
O(N_rN_tN_D).
\]

该方法矩阵小、复杂度低，但对 pilot pair spacing、CDD phase difference 和底层信道局部平坦假设非常敏感。

#### Reconstruction-B：delay-domain basis LMMSE

观测模型：

\[
\tilde{\mathbf g}_{r,\mathcal P}
=
\mathbf\Phi\mathbf h_r+\mathbf n_r,
\quad
\mathbf\Phi\in\mathbb C^{N_P\times N_tL}.
\]

使用 observation-domain LMMSE：

\[
\hat{\mathbf h}_r
=
\mathbf R_h\mathbf\Phi^H
\left(
\mathbf\Phi\mathbf R_h\mathbf\Phi^H
+\sigma_{\mathrm{LS}}^2\mathbf I
\right)^{-1}
\tilde{\mathbf g}_{r,\mathcal P}.
\]

预计算：

1. 构造 \(\mathbf\Phi\)：

\[
O(N_PN_tL).
\]

2. 构造 observation covariance：

\[
\mathbf S
=
\mathbf\Phi\mathbf R_h\mathbf\Phi^H
+\sigma_{\mathrm{LS}}^2\mathbf I,
\quad
\mathbf S\in\mathbb C^{N_P\times N_P}.
\]

复杂度约为

\[
O(N_P^2N_tL).
\]

3. 分解 / 求解 \(\mathbf S\)：

\[
O(N_P^3).
\]

若使用 tap-domain normal equation，则可能出现

\[
O((N_tL)^3)
\]

的求解复杂度；当 \(N_tL>N_P\) 时该形式还容易 rank deficient。

4. 输出 condition number / effective rank 需要 SVD：

\[
O\left(\min(N_P,N_tL)^2\max(N_P,N_tL)\right).
\]

每 trial 在线复杂度：

1. 由 pilot observation 估计 taps：

\[
O(N_rN_tLN_P).
\]

2. 由 taps 重构 branch frequency response：

\[
O(N_rN_tLN_D).
\]

3. CDD recombination：

\[
O(N_rN_tN_D).
\]

存储复杂度：

\[
O(N_PN_tL+N_P^2+N_tLN_D).
\]

该方法的复杂度和数值风险随 \(L\) 线性到三次增长。若 99% PDP support 很长，例如 100 ns delay spread 下 \(L\) 很大，则 \(N_tL\) 可能超过 \(N_P\)，反问题欠定，condition number 也可能非常大。

#### Reduced-support basis

若使用 reduced support \(\mathcal L_q\)，令

\[
L_q=|\mathcal L_q|<L.
\]

则上面的复杂度中 \(L\) 替换为 \(L_q\)。它可以显著降低

\[
O(N_P^2N_tL_q),
\quad
O(N_rN_tL_qN_D)
\]

等项，也可能改善 condition number。但代价是 support truncation bias，若丢弃 tap 能量过多，NMSE 会变差。

#### Reconstruction-C：coherence-band blockwise reconstruction

每个 block 内只估计 \(N_t\) 个 port coefficient：

\[
\tilde{\mathbf g}_{r,\mathcal P_b}
=
\mathbf A_b\mathbf h_r^{(b)}+\mathbf n_{r,b},
\quad
\mathbf A_b\in\mathbb C^{P_b\times N_t}.
\]

预计算复杂度：

\[
\sum_{b=1}^{B}
O(P_bN_t^2+N_t^3)
\]

若使用 observation-domain solve，也可写成

\[
\sum_{b=1}^{B}O(P_b^3),
\]

但通常 \(P_b\) 和 \(N_t\) 都很小。

每 trial 在线复杂度：

\[
\sum_{b=1}^{B}
O(N_rP_bN_t+N_rK_bN_t).
\]

因为

\[
\sum_b K_b=N_D,
\quad
\sum_b P_b\approx N_P,
\]

所以总量级约为

\[
O(N_rN_t(N_P+N_D)).
\]

该方法复杂度最低、存储最少，但依赖 block 内底层 physical channel 近似常数。若 block 宽度选择过大，会产生 flatness error；若 block 宽度过小，pilot 数不足，局部矩阵不可辨识。

### 13.4 算法三：Non-CDD per-port DMRS physical-channel RMMSE

算法三先对每个 port 独立估计底层 physical/effective port channel：

\[
\hat{\mathbf H}_{r,m,\mathcal D}
=
\mathbf R_{\mathcal D\mathcal P_m}
\left(
\mathbf R_{\mathcal P_m\mathcal P_m}
+\sigma_{\mathrm{LS},m}^2\mathbf I
\right)^{-1}
\tilde{\mathbf H}_{r,m,\mathcal P_m},
\]

再合成 CDD equivalent channel：

\[
\hat g_r[k]
=
\frac{1}{\sqrt{N_t}}
\sum_{m=0}^{N_t-1}
\hat H_{r,m}[k]
e^{-j2\pi kd_m/N_{\mathrm{FFT}}}.
\]

预计算复杂度：

\[
\sum_{m=0}^{N_t-1}
O(N_{P,m}^2L+N_DN_{P,m}L+N_{P,m}^3+N_DN_{P,m}^2).
\]

每 trial 在线复杂度：

\[
O\left(
N_rN_D\sum_{m=0}^{N_t-1}N_{P,m}
\right)
\]

加上 CDD synthesis：

\[
O(N_rN_tN_D).
\]

#### Equal-total-overhead

同总 DMRS overhead 下：

\[
\sum_m N_{P,m}=N_P.
\]

若各端口平均分配：

\[
N_{P,m}\approx \frac{N_P}{N_t}.
\]

此时求解 cubic 项约为

\[
\sum_m N_{P,m}^3
\approx
N_t\left(\frac{N_P}{N_t}\right)^3
=
\frac{N_P^3}{N_t^2}.
\]

相比算法一的 \(O(N_P^3)\)，预计算矩阵求解可能更低；但在线应用仍约为

\[
O(N_rN_DN_P)+O(N_rN_tN_D),
\]

与算法一同量级。代价是每端口 pilot density 下降，per-port interpolation error 可能增加。

#### Equal-per-port-density

若每个 port 使用与算法一相同的 pilot density：

\[
N_{P,m}=N_P,
\]

则总 pilot RE 约为算法一的 \(N_t\) 倍。复杂度变为

\[
O(N_tN_P^3)
\]

级别的预计算求解，以及

\[
O(N_rN_tN_DN_P)
\]

级别的在线应用。该模式 NMSE 往往更好，但计算复杂度、DMRS overhead 和 data RE loss 都显著增加，不能与算法一做同开销直接比较。

### 13.5 典型 48 RB 场景数量级

以 48 RB、\(N_{\mathrm{sc}}=N_D=576\)、2Tx、4Rx 为例：

| DMRS spacing | \(N_P\) |
|---:|---:|
| 2 sc | 288 |
| 4 sc | 144 |
| 6 sc | 96 |
| 12 sc | 48 |

若 \(S_f=6\)，则 \(N_P=96\)。

| 算法 | 典型矩阵规模 | 预计算主项 | 每 trial 主项 | 备注 |
|---|---|---:|---:|---|
| Alg1 WB RMMSE | \(96\times96\) pilot covariance | \(O(96^3)\) | \(O(4\cdot576\cdot96)\) | matched 时数值稳健 |
| Alg2 Basis, DS30 E99 | \(N_tL=34\), \(N_P=96\) | \(O(96^2\cdot34)+O(96^3)\) | \(O(4\cdot2\cdot17\cdot576)\) | 仍可能病态 |
| Alg2 Basis, DS100 E99 | \(N_tL=114\), \(N_P=96\) | \(O(96^2\cdot114)+O(96^3)\) | \(O(4\cdot2\cdot57\cdot576)\) | 未知数多于 pilot，欠定风险高 |
| Alg2 Blockwise | per block \(P_b\times2\) | 小矩阵逐 block | \(O(4\cdot2\cdot(N_P+576))\) | 低复杂度但近似误差大 |
| Alg3 equal-total | 两个 \(48\times48\) per-port covariance | \(2O(48^3)\) | \(O(4\cdot576\cdot96)+O(4\cdot2\cdot576)\) | 同开销，复杂度同量级 |
| Alg3 equal-per-port | 两个 \(96\times96\) per-port covariance | \(2O(96^3)\) | \(O(4\cdot2\cdot576\cdot96)\) | 约 2 倍 DMRS overhead |

这说明：

1. 算法一不是复杂度最高的方案；它的主要成本是 \(N_P^3\) 和 \(N_rN_DN_P\)，实现成熟、数值稳定。
2. 算法二 basis reconstruction 的计算和存储会随 \(L\) 增长，且有明显病态风险；即便复杂度可以接受，数值稳定性也可能是主要瓶颈。
3. 算法二 blockwise reconstruction 复杂度最低，但性能受 block flatness assumption 限制。
4. 算法三 equal-total 的计算复杂度与算法一同量级，甚至 per-port 小矩阵求解更便宜；但它需要不同 DMRS 设计，并且对 CDD delay synthesis error 更敏感。
5. 算法三 equal-per-port 的 NMSE 通常更好，但资源开销和在线复杂度都随 \(N_t\) 增加，不是同 overhead 比较。

### 13.6 综合比较

| 维度 | 算法一 Direct RMMSE | 算法二 Structural reconstruction | 算法三 Port-DMRS reconstruction |
|---|---|---|---|
| 观测 | CDD-combined DMRS | CDD-combined DMRS | non-CDD per-port DMRS |
| 估计对象 | CDD equivalent channel | physical branch/taps，再合成 equivalent channel | per-port physical channel，再合成 equivalent channel |
| 预计算复杂度 | 中等，\(O(N_P^3)\) | basis 较高且病态；blockwise 低 | equal-total 中等；equal-per-port 较高 |
| 每 trial 复杂度 | \(O(N_rN_DN_P)\) | basis \(O(N_rN_tL(N_P+N_D))\)；blockwise 低 | equal-total 与 Alg1 同量级，外加 synthesis |
| 存储 | \(O(N_DN_P)\) | basis \(O(N_PN_tL+N_tLN_D)\) | equal-total \(O(N_D\sum_mN_{P,m})\) |
| 数值稳定性 | matched covariance 下最好 | basis 容易 rank deficient / ill-conditioned | per-port RMMSE 稳定，但依赖 port-DMRS 设计 |
| DMRS overhead | baseline | 与算法一相同 | equal-total 相同；equal-per-port 增加 |
| CDD delay 误差敏感性 | 主要是 covariance mismatch | \(\Phi\) 和 recombination 都受影响 | synthesis phase error 直接影响，通常最敏感 |
| 适合场景 | 已知 CDD delay / PDP，matched LMMSE baseline | 有额外 probing 或强物理低维 prior | 底层信道平滑、CDD equivalent 快变、delay signaling 准确 |

当前实现和实验结果建议：

- 若只使用同一组 CDD-combined DMRS，算法一是最重要的 matched baseline，复杂度和稳定性都好。
- 算法二要获得实际收益，不能只靠更复杂的 inversion；更需要改善 observation design 或引入更强物理先验。
- 算法三 equal-total 在复杂度上可接受，且同开销下存在胜出区域；但其系统代价是 DMRS pattern 变化和 CDD delay 精度要求提高。

---

## 14. 关键判据与预期结论

### 14.1 QC 复现判据

第一版静态 TDL 结果应满足：

1. CDD_4RB_KNOWN 在 8 RB 下明显优于 PRG_CYCLING_4RB；
2. CDD_4RB_KNOWN 在 48 RB 下相对 PRG_CYCLING_4RB 只有有限增益；
3. CDD_WB_KNOWN 在 48 RB 下显著优于 CDD_4RB_KNOWN；
4. CDD_WB_UNKNOWN 低于 CDD_WB_KNOWN，尤其 30 ns delay spread；
5. IDEAL_CSI 曲线最好。

### 14.2 RMMSE vs Reconstruction 预期

预计存在以下 tradeoff：

1. 当底层 TDL delay spread 小、物理分支信道较平滑、CDD delay 较大、且 DMRS density 足够支持解耦时，reconstruction 可能优于 direct RMMSE；
2. 当 SNR 较低、DMRS 稀疏、CDD 相位差不足、解耦矩阵病态时，pairwise reconstruction 可能差于 RMMSE；
3. Reconstruction-B multi-pilot basis LMMSE 应比 Reconstruction-A pairwise decoupling 更稳健；
4. 增加 DMRS frequency density 可以同时改善 RMMSE 和 reconstruction，但二者收益曲线可能不同；
5. 若 CDD delay 过大，direct RMMSE 的频域相关性下降，远距离 pilot 对估计帮助有限；此时 either 需要更密 DMRS，或者需要 reconstruction 利用 CDD 确定性结构。
6. PORT_DMRS_RECON_RMMSE 预期能避开 CDD-combined reconstruction 的解耦矩阵病态问题；但若采用 equal-total-overhead，每端口 pilot density 下降可能抵消一部分收益，若采用 equal-per-port-density，则必须用 net goodput 体现额外 DMRS overhead 的代价。

---

## 15. 开发实现建议

### 15.1 模块划分

建议开发为以下模块：

1. `config`: 场景配置、MCS、DMRS、CDD delay、CE method；
2. `resource_grid`: PDSCH/DMRS RE 映射，支持 CDD-combined DMRS 和 port-specific non-CDD DMRS 两类 pilot pattern；
3. `channel_tdl`: TDL 信道生成，支持 PDP、delay spread、Tx/Rx branch；
4. `precoder`: CDD、PRG cycling、NO_CDD；
5. `tx_chain`: TB、LDPC、rate matching、QAM、OFDM resource mapping；
6. `rx_ls`: CDD-combined DMRS LS 估计，以及 port de-orthogonalization 后的 per-port LS 估计；
7. `rx_ce_rmmse`: direct RMMSE 估计；
8. `rx_ce_reconstruction`: pairwise 和 basis reconstruction；
9. `rx_ce_port_dmrs`: non-CDD per-port physical-channel RMMSE，以及按已知 CDD delay 合成数据等效信道；
10. `rx_equalizer`: 单层 LMMSE/MRC equalizer；
11. `rx_decode`: LLR、LDPC decode、CRC/BLER；
12. `metrics`: BLER、NMSE、condition number、goodput；
13. `plots`: BLER/NMSE/gain/overhead/delay sweep 图。

### 15.2 可重复性要求

每个 trial 必须记录：

- random seed；
- channel realization ID；
- DMRS pattern；
- CDD delay vector；
- CE method；
- PDP/covariance assumption；
- SNR；
- decoding result；
- CE NMSE。

### 15.3 数值稳定性要求

1. RMMSE 矩阵求逆使用 Cholesky 或 Hermitian solver，不直接 `inv()`；
2. Reconstruction-A 输出 condition number，若超过阈值，例如 \(10^3\)，标记 unstable；
3. Reconstruction-B 支持 diagonal loading：

\[
\mathbf\Phi\mathbf R_h\mathbf\Phi^H+\sigma^2\mathbf I+\epsilon\mathbf I.
\]

4. 所有 CE 输出应检查 NaN/Inf；
5. 对所有算法统计 empirical NMSE，避免 BLER 异常时无法定位原因。

---

## 16. 第一版最小可交付范围

第一版不需要实现所有扩展，最低要求如下：

1. TDL static channel；
2. 2Tx/4Rx 和 4Tx/4Rx；
3. 8 RB 和 48 RB PDSCH；
4. 2 DMRS symbols，DMRS frequency density 可配置；
5. PRG_CYCLING_4RB baseline；
6. CDD transmit diversity；
7. IDEAL_CSI；
8. RMMSE_4RB_KNOWN 和 RMMSE_WB_KNOWN；
9. RMMSE_WB_UNKNOWN；
10. RECON_PAIRWISE；
11. 推荐实现 RECON_BASIS_LMMSE；
12. Phase 2 推荐新增 PORT_DMRS_RECON_RMMSE 的 FDM v0，至少支持 2Tx、equal-total-overhead 和 equal-per-port-density 两种模式；
13. BLER vs SNR 曲线；
14. CE NMSE vs SNR；
15. 10% BLER SNR extraction；
16. CDD delay sweep；
17. DMRS density sweep。

---

## 17. 后续扩展计划

### 17.1 加入 UE mobility

加入 Doppler 后，RMMSE covariance 需要扩展到 time-frequency 相关矩阵：

\[
R[(k,t),(k',t')]
=R_f[k-k']R_t[t-t'].
\]

需要研究 2 个 DMRS symbols 在时间上如何联合处理，以及 60 km/h 高速下 CDD 的收益是否仍成立。

### 17.2 CDL 信道

CDL 下信道为

\[
\mathbf H[k]
=\sum_\ell \alpha_\ell e^{-j2\pi k\tau_\ell\Delta f}
\mathbf a_{\mathrm{rx}}(\theta_\ell^{\mathrm{rx}})
\mathbf a_{\mathrm{tx}}^H(\theta_\ell^{\mathrm{tx}}).
\]

需要研究：

- beam-dependent PDP；
- effective port correlation；
- CDL angle-delay coupling；
- CDD delay 与 beamformed PDP 的关系。

### 17.3 Massive MIMO 宽带预编码 + CDD

大规模阵列先做宽带预编码：

\[
\mathbf H_{\mathrm{eff}}[k]
=\mathbf H_{\mathrm{phys}}[k]\mathbf W_{\mathrm{WB}}.
\]

再在 effective ports 上做 CDD：

\[
\mathbf y[k]
=\mathbf H_{\mathrm{eff}}[k]
\mathbf D_{\mathrm{CDD}}[k]
\mathbf v s[k]+\mathbf n[k].
\]

该阶段验证 2Tx/4Tx effective-port 抽象在 128T/512T 阵列下是否仍能代表实际性能。

---

## 18. 最终希望回答的研究问题

1. QC 显式 CDD 的 BLER 增益是否可以在独立链路级平台中复现？
2. CDD 增益中 diversity gain 与 frequency-domain processing gain 的相对贡献如何随 PDSCH 带宽变化？
3. RMMSE direct equivalent-channel estimation 是否会在大 CDD delay / 强频率选择性下受限？
4. 基于 known CDD delay 的 structural reconstruction 是否能在相同 DMRS overhead 下优于 RMMSE？
5. 若要弥补 CDD 后强频选导致的信道估计误差，应优先增加 DMRS density、减小 CDD delay，还是采用 reconstruction estimator？
6. 在低 delay spread 和高 delay spread 下，最优 CDD delay 和最优信道估计算法是否不同？
7. DMRS 不做 CDD、每端口独立估计底层物理信道后再合成 CDD 数据等效信道，是否能在解耦矩阵病态场景下优于 CDD-combined direct RMMSE 和 structural reconstruction？
8. 第一版 TDL 结果是否能为后续 CDL + wideband precoding + CDD 的系统设计提供明确方向？

---

## 19. 推荐的第一版实验交付物清单

1. `config_qc_static_tdl.yaml`：QC 静态 TDL 复现配置；
2. `config_recon_vs_rmmse.yaml`：RMMSE vs reconstruction 对比配置；
3. `bler_curves_qc_static_tdl/`：QC 复现 BLER 图；
4. `bler_curves_recon_vs_rmmse/`：重构算法对比 BLER 图；
5. `ce_nmse_curves/`：信道估计 NMSE 图；
6. `cdd_delay_sweep/`：CDD delay sweep 图和表；
7. `dmrs_density_sweep/`：DMRS overhead 公平性图和表；
8. `port_dmrs_overhead_sweep/`：算法三 equal-total-overhead / equal-per-port-density 对比图和表；
9. `summary_10pct_bler_snr.csv`：所有场景的 10% BLER SNR 汇总；
10. `trial_metrics.csv`：trial-level 或 batch-level 统计结果；
11. `README_results.md`：结果解释和与 QC 趋势的对应关系。

---

## 20. 简短结论

本仿真计划把 QC 显式 CDD 的复现与进一步算法创新分为两层。第一层复现 QC 的核心机制：CDD 通过 RE-level phase-continuous precoder cycling 获得分集，并通过 wideband RMMSE 信道估计获得频域处理增益；显式 delay signaling 使 UE 能正确构造 CDD-shifted PDP/covariance。第二层引入两类进一步信道估计方案：其一是 structural / deterministic reconstruction，利用已知 CDD delay 从 CDD-combined DMRS 中恢复底层物理分支信道，再重构 CDD 等效信道；其二是 non-CDD per-port DMRS，在导频阶段直接估计每个端口的底层信道，再在数据 RE 上按已知 CDD delay 合成等效信道。二者共同用于验证在强 CDD 频率选择性、稀疏 DMRS、不同 delay spread 下是否能优于 direct RMMSE，并量化算法收益与 DMRS overhead 的权衡。

第一版采用 static TDL，可以快速隔离 delay spread、CDD delay、DMRS density 和信道估计算法本身的作用；后续再加入 mobility、CDL、大规模阵列宽带预编码和 effective-port CDD。

---

## 21. 第一版 CDD LLS 平台实现说明

本节记录当前仓库中已经实现的第一版平台，便于后续复现实验、扩展算法和解读输出结果。实现目标是先形成一个本地可运行、结果可追踪的 bit-level 闭环：Sionna 负责真实 LDPC 编码/解码，CDD/PRG/TDL/DMRS 信道估计和均衡检测由平台显式实现。这样既保留了 Sionna 的信道编码可信度，又能把 CDD delay、PDP covariance、wideband/4RB processing、reconstruction matrix 等研究对象暴露出来。

### 21.1 项目结构

当前平台结构如下：

```text
CDD_LLS/
  run.py
  README.md
  configs/
    smoke.yaml
    smoke_variants.yaml
    config_qc_static_tdl.yaml
    config_recon_vs_rmmse.yaml
  cdd_lls/
    core/
      config.py
      mcs.py
    phy/
      channel_tdl.py
      estimators.py
      ldpc.py
      precoding.py
      qam.py
      resource_grid.py
    sim/
      orchestrator.py
      stats.py
    utils/
      plotting.py
  tests/
    test_core.py
```

各部分职责：

- `run.py`：统一运行入口，读取 YAML，启动仿真 orchestrator。
- `core/config.py`：定义 dataclass 配置、YAML 合并、配置校验、resolved config 保存。
- `core/mcs.py`：保存第一版需要的 NR MCS 表条目，生成单层单 CW 的 TBS、coded bits 和 CB 切分。
- `phy/resource_grid.py`：生成 PDSCH 子载波、DMRS pilot 子载波、data RE 序列和 DMRS overhead。
- `phy/channel_tdl.py`：生成 static TDL tap 信道，默认指数 PDP，并计算每个 active subcarrier 上的频域 MIMO 信道。
- `phy/precoding.py`：实现 `CDD`、`NO_CDD`、`PRG_CYCLING_4RB`，并由 branch channel 生成单层等效信道。
- `phy/estimators.py`：实现 `IDEAL_CSI`、`LS_LINEAR`、direct RMMSE、pairwise reconstruction、basis LMMSE reconstruction，并输出 NMSE/condition number。
- `phy/qam.py`：NumPy QAM 调制和 max-log soft demapper。LLR 符号约定与 Sionna decoder 使用保持一致：正 LLR 倾向 bit=1。
- `phy/ldpc.py`：Sionna `LDPC5GEncoder` / `LDPC5GDecoder` 适配层，每个 CB 独立编码/解码。
- `sim/orchestrator.py`：仿真主循环，展开 scenarios/variants/sweeps，运行 SNR sweep，保存统计结果。
- `sim/stats.py`：CSV/JSON 保存、SNR range 展开、10% BLER SNR 插值。
- `utils/plotting.py`：生成 BLER 和 CE NMSE 曲线。
- `tests/test_core.py`：标准库 `unittest` 基础测试，不依赖 pytest。

### 21.2 配置展开逻辑

YAML 顶层配置由 base config、`scenarios`、`variants` 和 `sweeps` 组合而成。

- base config 定义公共参数，例如天线数、资源网格、MCS、默认 SNR sweep。
- `scenarios` 通常用于改变信道和资源，例如 `tdl30_8rb`、`tdl100_48rb`。
- `variants` 通常用于改变发射方案和 CE 方法，例如 PRG baseline、CDD wideband RMMSE、CDD unknown delay、ideal CSI。
- `sweeps.cdd_base_delays` 会把 `cdd_delay_vector` 设为 `[0,d]` 或 `[0,d,2d,3d]` 形式。
- `sweeps.dmrs_spacing_sc` 会覆盖 DMRS 频域间隔，用于 DMRS density sweep。

每个 scenario/variant/sweep 组合都会独立运行全部 SNR 点，并在输出表中用 `scenario_id` 和 `variant_id` 标识。这样同一份 `summary.csv` 可以同时包含 QC 复现场景、delay sweep 和 DMRS density sweep。

### 21.3 发射链路实现

第一版固定为单层 PDSCH、单 CW。每个 trial 的流程是：

1. 根据 `resource` 生成 active subcarriers、DMRS pilots 和 data RE。
2. 根据 MCS 计算 coded-bit 容量、TBS 和 CB 切分。
3. 随机生成每个 CB 的 payload bits。
4. 使用 Sionna LDPC5GEncoder 对每个 CB 编码。
5. 串接 CB coded bits，按 data RE 顺序做 QAM 调制。
6. 生成 static TDL branch channel `H[r,m,k]`。
7. 根据发射方案生成 precoder `c[k]`，得到等效信道 `g[r,k]`。
8. data RE 上发送 `s[k]`，接收端看到 `y[r,k]=g[r,k]s[k]+w[r,k]`。

当前没有显式 OFDM IFFT/CP waveform。CDD 在频域以确定性相位斜坡实现，和计划书中的频域表达一致。TDL tap 到频域信道的变换使用 active subcarrier index 和 `n_fft`，因此 CDD delay sample、TDL tap delay sample、DMRS frequency spacing 在同一个频域相位模型中定义。

### 21.4 信道和 PDP 实现

第一版 TDL 是静态 tapped-delay line：

- 默认 PDP 为指数衰落，由 `channel.delay_spread_ns` 和 OFDM sample period 生成。
- sample period 由 `n_fft * scs_khz` 得到。
- `channel.max_delay_factor` 控制 PDP 截断长度。
- 每个 Tx-Rx branch 独立生成复高斯 taps。
- taps 能量默认归一化到 1。

这不是完整 3GPP TDL-A/B/C profile，而是用于隔离 CDD delay、DMRS density、RMMSE covariance 和 reconstruction 行为的第一版静态 TDL。后续如需精确对齐 QC 提案图，可在 `channel_tdl.py` 中替换为 38.901 TDL profile。

### 21.5 发射方案实现

`CDD`：

- 2Tx 默认 delay vector 为 `[0,d]`。
- 4Tx 默认 delay vector 为 `[0,d,2d,3d]`。
- 也可以在 YAML 中直接写 `cdd_delay_vector`。
- 每个 Tx branch 的频域权重为 `exp(-j 2 pi k d_m / N_fft) / sqrt(N_tx)`。

`PRG_CYCLING_4RB`：

- 默认 codebook 使用计划书中的 2Tx QPSK relative-phase codebook 和 4Tx DFT-like QPSK codebook。
- 每个 PRG 内 precoder 固定。
- `resource.prg_size_rb` 控制 PRG 大小，默认 4 RB。
- `transmission.prg_cycling_order` 可配置 precoder 循环顺序。

`NO_CDD`：

- 所有 delay 置零，可用于 sanity check。

输出表会记录 `tx_scheme`、`cdd_delay_vector`、`prg_size_rb`、`prg_codebook`，便于后续追踪 QC baseline 差异。

### 21.6 DMRS 和 LS 观测实现

DMRS frequency pattern 由 `resource.dmrs_spacing_sc` 和 `resource.dmrs_offset_sc` 决定：

```text
pilot subcarriers = offset, offset + spacing, offset + 2*spacing, ...
```

DMRS symbols 由 `resource.dmrs_symbol_indices` 给出。第一版假设信道在 PDSCH 时间内静态，因此多个 DMRS symbols 只用于降低 LS 观测噪声。若有 2 个 DMRS symbols，则 LS noise variance 使用 `noise_var / 2`。

当前已实现的 CE 方法使用同一组 CDD-combined pilot RE 和同一组 data RE，保证 Direct RMMSE 与 CDD-combined reconstruction 的 overhead 公平。算法三 `PORT_DMRS_RECON_RMMSE` 需要新增 port-specific non-CDD pilot pattern：equal-total-overhead 模式保持总 pilot RE 不变但降低每端口密度，equal-per-port-density 模式保持每端口密度但增加总 pilot RE。输出表中的 `dmrs_overhead`、`n_dmrs_re`、`data_re` 可用于检查不同 DMRS density 下的资源开销；实现算法三后还应增加 `dmrs_overhead_per_port` 和 `pilot_density_per_port`。

### 21.7 信道估计算法实现

`IDEAL_CSI`：

- 直接使用真实等效信道 `g[r,k]`，作为性能上界。

`LS_LINEAR`：

- 在 DMRS pilot 上做 LS，频域线性插值得到全部 active subcarriers 的等效信道。

`RMMSE_4RB_KNOWN` / `RMMSE_WB_KNOWN`：

- 直接估计 CDD 后的等效信道。
- known delay 时，先把底层 PDP 按 CDD delay shift 后平均，构造 CDD-shifted PDP。
- 由 PDP 得到频域 covariance，再做 LMMSE/RMMSE 插值。
- `4RB` 按 `channel_estimation.rmmse_bundle_rb` 独立处理；`WB` 在整个 allocation 内联合处理。

`RMMSE_4RB_UNKNOWN` / `RMMSE_WB_UNKNOWN`：

- 真实观测仍来自 CDD 等效信道。
- 估计器假设 non-CDD PDP，因此 covariance mismatch 会反映在 `ce_nmse_eff` 和 BLER 中。

`PRG_RMMSE_4RB`：

- 用于 PRG baseline。
- 每个 4RB PRG 内独立估计，不跨 PRG 做 wideband processing。

`RECON_PAIRWISE`：

- 使用 pilot group 解耦底层 Tx branch channel。
- 对每个 pilot group 构造 CDD phase matrix，做 regularized LS。
- 对解出的 branch channel anchor 做频域插值。
- 最后按已知 CDD delay 重构等效信道。
- 输出 `cond_number`，用于判断 pairwise 解耦是否病态。

`RECON_BASIS_LMMSE`：

- 把底层 branch channel 展开到 delay-domain basis。
- pilot 观测矩阵显式包含 `ell + d_m` 的 CDD shifted phase。
- 使用 PDP 作为 tap prior 做 LMMSE tap estimation。
- 支持 `basis_support: ideal` 或 `truncated`；`truncated` 默认保留累计能量达到 99% 的 taps。
- 输出 `ce_nmse_branch`、`ce_nmse_eff`、`cond_number` 和 `effective_rank`。

`PORT_DMRS_RECON_RMMSE`，计划扩展：

- DMRS RE 不施加 CDD，按 Tx/effective port 分配正交 pilot。
- 接收端先对每个端口做 LS 和 per-port RMMSE，估计底层物理端口信道 `H[r,m,k]`。
- 数据 RE 上按已知 `cdd_delay_vector` 乘回 CDD 相位斜坡，并合成等效信道 `g[r,k]`。
- 需要同时记录 `ce_nmse_branch`、`ce_nmse_eff`、总 DMRS overhead、每端口 pilot density 和 net goodput。
- 第一版建议先实现 FDM port DMRS，后续再扩展 CDM/OCC。

所有矩阵求解使用 `np.linalg.solve`，并加 `channel_estimation.diagonal_loading`，避免直接求逆。

### 21.8 均衡、LLR 和 LDPC 解码

接收端使用单层 MRC/ZF 风格均衡：

```text
z[k] = g_hat[:,k]^H y[:,k] / ||g_hat[:,k]||^2
```

解调噪声方差近似为：

```text
noise_eff[k] = noise_var / ||g_hat[:,k]||^2
```

然后用 max-log QAM demapper 产生 coded-bit LLR，按 CB 的 rate-matching 长度切片，送入 Sionna LDPC5GDecoder。BLER 定义为任意 CB 解码失败则 TB 失败；CB-BLER 也单独输出。Goodput 统计成功 CB 的 payload bits。

### 21.9 运行方法

推荐使用本机已验证的 Sionna 环境：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python run.py --config configs/smoke.yaml
```

多算法 smoke：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python run.py --config configs/smoke_variants.yaml
```

QC 静态 TDL 复现：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python run.py --config configs/config_qc_static_tdl.yaml
```

RMMSE vs reconstruction 对比：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python run.py --config configs/config_recon_vs_rmmse.yaml
```

基础测试：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python -m unittest discover -s tests
```

### 21.10 输出文件解读

每次运行会在 `simulation.output_dir` 下创建 `sim_YYYYMMDD_HHMMSS` 目录。

`resolved_config.yaml`：

- 保存本次运行的完整配置。
- 用于确认 YAML 默认值、scenario/variant 基础设置和运行参数。

`summary.csv` / `summary.json`：

- 每行对应一个 scenario/variant/SNR。
- 主要性能字段：
  - `bler`：TB BLER，主性能指标。
  - `cb_bler`：CB 维度错误率。
  - `goodput_bits_per_slot`：平均成功 payload bits。
  - `goodput_se_per_re`：按 data RE 归一化的 goodput。
  - `ce_nmse_eff`：最终等效信道 `g` 的 NMSE。
  - `ce_nmse_branch`：reconstruction 方法的底层 branch channel NMSE。
  - `cond_number`：reconstruction matrix 或 weighted basis matrix 条件数。
  - `effective_rank`：basis reconstruction 的有效秩。
- 主要配置字段：
  - `scenario_id`、`variant_id`：曲线分组标识。
  - `tx_scheme`、`ce_method`：发射和信道估计方法。
  - `n_tx`、`n_rx`、`pdsch_rb`、`delay_spread_ns`。
  - `dmrs_spacing_sc`、`dmrs_overhead`、`data_re`。
  - `cdd_delay_vector`、`prg_size_rb`、`prg_codebook`。
  - `mcs_table`、`mcs_index`、`tbs_bits`、`coded_bits`。

`summary_10pct_bler_snr.csv`：

- 按 `simulation.bler_target`，默认 10% BLER，对每条曲线做 SNR 插值。
- `snr_at_target_bler_db` 是该算法达到目标 BLER 所需 SNR。
- 若同一 scenario 中存在 PRG baseline，`snr_gain_vs_prg_db` 给出相对 PRG 的 SNR gain，正值表示该算法更好。
- 如果曲线没有穿过目标 BLER，则该字段为 `nan`，需要扩大 SNR 范围或增加 trial 数。

`trial_metrics.csv`：

- 只有 `simulation.save_trial_metrics: true` 时输出。
- 用于 debug 单 trial 的 CE NMSE、条件数和解码结果。

`bler_curves.png`：

- 按 `scenario_id | variant_id` 分组画 BLER vs SNR。
- 正式论文/报告出图建议从 `summary.csv` 重新按需要筛选和美化。

`ce_nmse_curves.png`：

- 按同样分组画 equivalent-channel NMSE vs SNR。
- 当 BLER 曲线异常时，优先检查该图和 `cond_number`。

### 21.11 如何看第一版结果

QC 静态 TDL 复现时，重点比较：

1. `PRG_CYCLING_4RB + PRG_RMMSE_4RB`：baseline。
2. `CDD + RMMSE_4RB_KNOWN`：CDD diversity gain，但 CE processing window 与 PRG 对齐。
3. `CDD + RMMSE_WB_KNOWN`：CDD phase-continuous 带来的 wideband CE processing gain。
4. `CDD + RMMSE_WB_UNKNOWN`：UE 不知道 CDD delay 的 covariance mismatch 损失。
5. `CDD + IDEAL_CSI`：接收端 CSI 上界。

如果 8 RB 场景下 CDD_4RB_KNOWN 优于 PRG baseline，而 48 RB 场景下 CDD_WB_KNOWN 明显优于 CDD_4RB_KNOWN，就说明第一版已经复现了计划书中的主要机制趋势。

RMMSE vs reconstruction 对比时，重点看：

- `ce_nmse_eff`：最终等效信道估计是否更准。
- `ce_nmse_branch`：reconstruction 是否真的恢复了底层 branch channel。
- `cond_number`：pairwise 或 basis matrix 是否病态。
- `summary_10pct_bler_snr.csv`：同一 DMRS overhead 下 reconstruction 是否带来 SNR gain。

如果 reconstruction 的 `ce_nmse_eff` 较低但 BLER 没有改善，通常说明当前均衡/LLR noise model 低估了估计误差引入的自噪声；这属于后续接收机建模需要继续精细化的部分。

### 21.12 当前限制和后续扩展点

当前第一版限制：

- 仅实现 static TDL，不含 Doppler 和 DMRS time-domain interpolation。
- TDL PDP 是指数模型，不是完整 38.901 TDL profile。
- 单层、单 CW，不含 rank>1 layer mapping。
- CDD 在频域等效实现，不生成完整 OFDM waveform。
- 均衡后 LLR 的 noise variance 使用 MRC 近似，没有显式把 CE error 当作额外 self-noise。
- PORT_DMRS_RECON_RMMSE 尚未在当前代码中实现；当前代码只覆盖 CDD-combined DMRS 下的 Direct RMMSE 与 structural reconstruction。
- CDL、massive-MIMO wideband precoding、effective-port CDD 尚未实现。

后续优先扩展建议：

1. 加入 38.901 TDL profile，替换当前指数 PDP。
2. 加入 mobility，扩展 time-frequency RMMSE covariance。
3. 把 CE error variance 注入 LLR noise model，改善非理想 CSI BLER 可信度。
4. 实现 PORT_DMRS_RECON_RMMSE，包括 FDM port-specific DMRS、两种 overhead 公平模式、per-port physical-channel RMMSE 和 CDD synthesis。
5. 加入 CDL channel builder 和 wideband beam/effective-port CDD。
6. 增加 batch-first LDPC 编解码，以提升大规模 sweep 速度。
