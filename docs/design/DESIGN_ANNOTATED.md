# V 设计理论框架（v2）：分集增益与信道估计增益的折中

> **文档角色**：这是本项目的**顶层理论与设计空间文档**。它给人看，也是分析 Claude 做**战略性规划**时的思想地基——每次要开辟新的 `V` 设计方向、或判断某条实验路线是否值得走时，对照这里的度量框架、塌缩条件和两条赛道。它**不是**每轮迭代都读的三个文件之一（那是 `GOALS.md` + `KNOWLEDGE.md` + 最新 `result`）；常规迭代不必读它。
> **版本沿革**：v1（原 `V_design_diversity_CE_tradeoff.md`，Gram/平滑性框架）已归档至 `docs/archive/DESIGN_v1.md`。v2 基于实验 1~21 的经验教训（尤其 result-021："矩阵指标赢 ≠ estimated-CSI BLER 赢"）与 2026-07 的原理讨论重写，**两条度量轴都换掉了**，见 §0。v1 中候选结构的原始细节（分段线性构造、chirp、frame/Welch 视角）仍可在归档版查阅。
> **符号规范**：全项目符号以本文**附录 A** 为唯一权威（发射分支数用 $N_t$、等效信道用 $\mathbf g$）。本文每个符号第一次出现处均附一句话解释。

---

## 0. v2 相对 v1 的修订摘要

| 维度 | v1 的坐标轴 | v2 的坐标轴 | 为什么改 |
|---|---|---|---|
| 分集增益度量 | Gram 谱指标：$\log\det(\mathbf V^H\mathbf V)$、$\lambda_{\min}$、条件数、互相干 | **有限 SNR 块互信息 outage**（§2） | Gram 不是 BLER 的充分统计量（§2.4），且 log-det 只在深尾渐近区有意义，10%/1% BLER 工作点不在那里 |
| CE 代价度量 | 相位平滑性（roughness、group delay spread） | **按接收机类型分开**：V-aware 看导频子采样可辨识性（§3.2–3.3）；V-agnostic 看时延窗代价（§3.4） | 平滑性只对"不知道 V、靠插值"的接收机成立；对 V-aware matched LMMSE，平滑性无关紧要，混叠/条件数才是硬约束 |
| 折中命题 | 分集 ↔ CE 处理增益普遍存在基本折中（v1 命题 4） | **分 regime**：V-aware + 近平坦物理信道 + 导频充足时折中**退化**（CDD 同时最优，§4.1）；V-agnostic 滑窗接收机下折中**真实存在**（§3.4） | 塌缩命题一举解释了 result-021 的全部现象 |
| 研究结构 | 单一探索："一般 V 能否优于 CDD" | **两条赛道**：Track A 非透明（收尾定基准）+ Track B 透明（开放前沿）（§5） | 研究者判据（2026-07）：透明设计若能明显优于透明基线、且接近/超过 QC 非透明方案，则值得做 |

---

## 1. 问题背景与系统模型

### 1.1 问题

在显式 CDD 或更一般的频域相位预编码中，$N_t$ 个发射分支（$N_t$：发射分支数，当前 8）在每个子载波上叠加不同相位，接收端看到的等效频域信道为

$$
\mathbf g = \mathbf V\,\mathbf h .
$$

其中：

- $\mathbf h \in \mathbb C^{N_t}$：底层各发射分支的信道，基准假设 $\mathbf h\sim\mathcal{CN}(\mathbf 0,\mathbf I_{N_t})$（$\mathbf I_{N_t}$ 为 $N_t\times N_t$ 单位矩阵）；
- $\mathbf g \in \mathbb C^{K}$：$K$ 个有效子载波上的等效信道（$K$：有效子载波数）；
- $\mathbf V \in \mathbb C^{K\times N_t}$：频域相位矩阵，$V_{k,n}=e^{j\phi_{k,n}}$，恒模约束 $|V_{k,n}|=1$（每分支每子载波只旋相位、不改功率）。

传统 CDD 是线性相位特例：$\phi_{k,n}=-2\pi k\,\Delta f\,\tau_n$，其中 $\Delta f$ 为子载波间隔、$\tau_n$ 为第 $n$ 分支的固定 cyclic delay。

北极星问题（同 `GOALS.md`）：**能否设计比 CDD 更一般的 $\mathbf V$ 及配套估计算法，在相同 DMRS 开销下同时获得更好的分集增益和更低的 estimated-CSI BLER？还是 CDD 已是实践最优？**

**建模注记（$\mathbf g=\mathbf V\mathbf h$ 隐含的平坦假设）**：此式成立的前提是**每个发射分支的底层信道在频域恒定**。$\mathbf h\in\mathbb C^{N_t}$ 是每分支一个标量、不带子载波下标 $k$；展开即 $g_k=\sum_{n=1}^{N_t}V_{k,n}\,h_n$，其中 $k$ 只进入 $V_{k,n}$，$h_n$ 与 $k$ 无关。因此 $\mathbf g$ 在频域的全部起伏都由 $\mathbf V$ 人工制造，分支信道被当作跨 $K$ 个子载波的常数——这正是 §1.3 近平坦假设（$\tau_{\rm phy}=5$ ns，相干带宽 $\sim$200 MHz $\gg$ 载波带宽）在系统模型层的体现。若分支信道有频率选择性，须写成 $h_{k,n}$，则 $g_k=\sum_n V_{k,n}h_{k,n}$ 是逐子载波的 $\mathbf V_k$ 与 $\mathbf h_k$ 配对，**无法坍缩成"$K\times N_t$ 矩阵乘 $N_t$ 维常向量"** 的形式。即：本式的可写性与近平坦假设等价。

### 1.2 仿真体制的关键参数（理论要对着它们讲）

| 参数 | 取值 | 理论上的角色 |
|---|---|---|
| $\Delta f$ | 30 kHz | 频域采样基本单位 |
| $K$ | 576 / 432 / 288（48/36/24 PRB） | 分集展开的频域维度 |
| DMRS comb 间隔 $S_f$ | 24 子载波 | 决定导频混叠周期 $\tau_{\rm alias}$（§3.3） |
| 每 DMRS 符号导频数 $N_p=K/S_f$ | 24 / 18 / 12 | CE 侧的"观测预算" |
| DMRS 符号数 | 2（symbol 2,7，时域平均） | **只降噪（$\sigma_{LS}^2=\sigma^2/2$），不增加频域采样点** |
| $N_t$ | 8 | CE 侧的"未知数预算"（平坦近似下） |
| MCS | 16QAM，码率 553/1024 | 谱效率 $R\approx 2.16$ bit/RE（outage 的门限） |
| 物理信道 | static TDL，delay spread $\tau_{\rm phy}=5$ ns | 近平坦（§1.3） |

### 1.3 一个决定理论走向的观察：物理信道近平坦

$\tau_{\rm phy}=5$ ns 对应相干带宽 $\sim$200 MHz，远大于 48 PRB 的 17.3 MHz。两个推论贯穿全文：

1. **频率选择性几乎全部是 $\mathbf V$ 人工制造的**——分集的原料不是天上掉的，是设计出来的；
2. **估计问题的本质是"未知数计数"**：平坦近似下等效信道由 $N_t=8$ 个复数完全参数化，估计就是"从 $N_p$ 个导频观测里解 $N_t$ 个未知数"。当前 $N_p/N_t = 3\times/2.25\times/1.5\times$，观测是充裕的。

$\mathbf R_h=\mathbf I_{N_t}$ 时等效信道协方差 $\mathbf R_g=\mathbb E[\mathbf g\mathbf g^H]=\mathbf V\mathbf V^H$；计入物理选择性时 $\mathbf R_g=\mathbf R_{\rm phy}\odot(\mathbf V\mathbf V^H)$（$\mathbf R_{\rm phy}$：物理信道频域协方差，$\odot$：Hadamard 积），效果是给每个人工时延分量做 $\tau_{\rm phy}$ 量级的展宽——设计时留足余量即可（§4.2 第 4 条）。

---

## 2. 分集侧理论：有限 SNR 的块互信息 outage

### 2.1 从实际链路推导

实际链路顺序是：LDPC 编码 → 16QAM 调制 → 映射到全带 RE → 逐子载波乘 $\mathbf V$ → 静态信道。outage 框架恰好适配这个顺序，推导只用三步：

**第一步（quasi-static）**：每个 TB 抽一次 $\mathbf h$，码字期间不变。所以一个码字看到的是"信道的一次随机抽样"，而不是遍历平均。

**第二步（并行子信道）**：理想 CSI 下，码字覆盖的 $K$ 个子载波 ×（数据符号数）个 RE 是一组并行 AWGN 子信道，第 $k$ 个子载波的接收 SNR 为 $\mathrm{snr}\cdot|g_k|^2$（$\mathrm{snr}$：发射符号信噪比）。静态信道下时间维只是复制，不添新随机性。$\mathbf V$ 的全部作用被压缩为：**决定 $\{|g_k|^2\}_{k=1}^K$ 的联合分布**。

**第三步（编码定理）**：对一次给定的信道实现，这组并行信道能可靠承载的最大速率是块平均互信息

$$
I(\mathbf h)\;=\;\frac1K\sum_{k=0}^{K-1} I_{\rm QAM}\!\big(\mathrm{snr}\,|g_k|^2\big),
$$

其中 $I_{\rm QAM}(\cdot)$ 是实际调制（16QAM）的 BICM 每 RE 互信息，一条一维查找曲线，封顶 4 bit。理想码：谱效率 $R<I(\mathbf h)$ 必然译对，$R>I(\mathbf h)$ 必然译错。于是

$$
\mathrm{BLER}(\mathrm{snr})\;\approx\;\Pr_{\mathbf h}\big[\,I(\mathbf h)<R\,\big]\;\equiv\;P_{\rm out}(\mathrm{snr},R).
$$

**物理含义**：误块的主导机制是"这次抽到的信道太差撑不起码率"（outage 事件），不是"噪声偶然太大"。所以 BLER 由随机变量 $I(\mathbf h)$ 的**左尾分布**决定，"分集增益"在有限 SNR 下的准确含义就是：**让 $I(\mathbf h)$ 的分布更集中、左尾更薄**。

**〔释疑 · 补注 2026-07-07〕为什么是"$K$ 个子载波容量的平均"，而不是"平均 SINR 的容量"？** 二者差一个**频率选择性惩罚**，看懂它就看懂了整根分集轴。

- **块平均互信息**（本文用的这个）$I(\mathbf h)=\frac1K\sum_k I_{\rm QAM}(\mathrm{snr}|g_k|^2)$：把 $K$ 个子载波当作 $K$ 条**互相独立的并行 AWGN 子信道**，第 $k$ 条按自己的接收 SNR $\mathrm{snr}|g_k|^2$ 贡献自己那份互信息。信息论基本事实——一组独立并行信道（各喂独立输入）的总互信息**等于各子信道互信息之和** $I(X_1..X_K;Y_1..Y_K)=\sum_k I(X_k;Y_k)$，除以 RE 数就是这里的平均。它才是物理正确的量：一个码字的比特被摊到全部 $K$ 个 RE，每个 RE 独立地过一个标量信道 $g_k$，RE 之间**不做跨子载波联合合并**（§2.3 前提 3）。
- **平均 SINR 的容量**（想当然但错的量）$I_{\rm QAM}\!\big(\mathrm{snr}\cdot\frac1K\sum_k|g_k|^2\big)$：先把各子载波的功率**汇合**成一条 SNR 取均值的单信道，再取容量。

因为 $I_{\rm QAM}(\cdot)$ 关于线性 SNR 是**凹函数**（I-MMSE 定理：$\frac{dI_{\rm QAM}}{d\rho}=\mathrm{mmse}(\rho)$，而 $\mathrm{mmse}(\rho)$ 随 $\rho$ 单调下降，故导数递减、函数凹；对任意输入分布都成立），由 Jensen 不等式

$$
\frac1K\sum_k I_{\rm QAM}(\mathrm{snr}|g_k|^2)\;\le\;I_{\rm QAM}\!\Big(\mathrm{snr}\cdot\tfrac1K\textstyle\sum_k|g_k|^2\Big).
$$

即**"平均 SINR 的容量"永远高估**，二者之差正是频率选择性的代价：强子载波多出来的 SNR **搬不去**补深衰子载波（每个 RE 固定调制、独立承载自己那段编码比特，深衰 RE 上的比特就是丢了），凹性精确记下了这份"不可汇合"。只有当信道平坦、或允许逐 RE 注水 / 自适应调制时两者才相等。

把它接到 §2.4.4 更清楚：全带列正交（$\mathbf V^H\mathbf V=K\mathbf I_{N_t}$）时 $\frac1K\sum_k|g_k|^2=\frac1K\mathbf h^H\mathbf V^H\mathbf V\mathbf h=\|\mathbf h\|^2$，于是"平均 SINR 容量" $=I_{\rm QAM}(\mathrm{snr}\|\mathbf h\|^2)$ 恰是 **MRC 相干合并上界**——一个恒模、只旋相位的设计根本达不到它（那要求接收机逐 RE 把 $N_t$ 个分支相干合并）。真实可达的天花板是自平均后的 $I_\infty(\mathbf h)\le I_{\rm QAM}(\mathrm{snr}\|\mathbf h\|^2)$，Jensen 的这道 gap 就是"无相干合并"的代价。

**〔释疑 · 补注 2026-07-08〕"总互信息 = 各子信道互信息之和"要求发送比特独立吗？——不要求。** 疑点：单流下一个 LDPC 码字铺到各 RE，coded bits 有奇偶校验约束、**并不独立**，那"独立并行输入"的前提还成立吗？答：成立，因为该式是**最大可达速率（容量）**，不是对实际发送符号的统计描述。互信息/容量由**输入分布**定义，取的是"让并行信道总 MI 最大"的那个分布（各 RE 独立、均匀 QAM）；这是**信道的性质**（= 随机编码定理里生成码本的 ensemble），不是"某个具体码字符号的经验分布"。三点：

1. **数学**：信道跨子载波无记忆（$p(y^K|x^K)=\prod_k p(y_k|x_k)$，各 RE 噪声独立），故对**任何**输入（相关或独立）都有 $I(X^K;Y^K)\le\sum_k I(X_k;Y_k)$，独立输入时取等。所以"和"是上确界、是可达目标；把输入弄成相关只会让联合 MI **变小**，不可能超过它——即输入相关顶多让你"够不满"这个和，不会在高估之外再引入别的错。
2. **编码定理**：只要码率 $R<\frac1K\sum_k I_{\rm QAM}$ 就存在能可靠译码的码。LDPC 的奇偶校验相关**正是它逼近该界的手段（受控冗余），不是违反**——好码必须有相关才有纠错能力；而每个 coded bit 的**边缘**分布 $\approx$ 均匀，与容量所用的输入边缘一致，这就够了。
3. **BICM 交织器专门干这个**：把码字比特打散映射到各 RE，使译码器看到的"每 RE 比特子信道"近似**独立并行**——这正是"**比特交织**编码调制"得名之由，也是 §2.3 前提 3（逐 RE 出 LLR、当独立并行信道喂译码器）成立的机制。工程上 RBIR/MIESM 就是拿**真实 LDPC** 大规模校准的、且好用，即经验佐证。

一句话：**"独立输入"指的是定义容量的那个输入集合，不是实际发出去的码字；$\frac1K\sum_k I_{\rm QAM}$ 是天花板，LDPC 的比特相关是够到它的梯子，方向相反、并不打架。**

**〔释疑 · 补注 2026-07-07〕$I_{\rm QAM}(\cdot)$ 是理论曲线吗？能算吗？** 是**确定性理论曲线，能精确计算，且只需算一次**。它是 16QAM 星座在 AWGN 下的 BICM 互信息，只是**线性/对数 SNR $\rho$ 的一维函数**，与信道随机性无关：

$$
I_{\rm QAM}(\rho)=\sum_{i=1}^{m} I(b_i;Y),\qquad Y=\sqrt\rho\,X+N,\ \ N\sim\mathcal{CN}(0,1),\ \ X\ \text{均匀取自单位能量 16QAM},\ \ m=4,
$$

每个 bit 子信道的 $I(b_i;Y)$ 由对噪声 $N$ 求期望给出（被积函数是若干 $\log_2$-of-exp 混合的闭式表达），用 Gauss–Hermite 求积或对 $N$ 数值积分即可算到任意精度。性质：单调增、低 SNR 趋 0、高 SNR 饱和于 $m=4$ bit。工程上就是在 SNR 网格上**预先算好存成查找表**（即 3GPP 链路抽象里的 RBIR / MI 映射曲线），仿真时查表——这正是 §2.5 说的"一维查找曲线"。注意它是**理论/确定性**曲线（由"星座 + AWGN"解析给出），**不是**从链路仿真拟合出来的经验曲线；§2.2 第 2 条的 $f_{\rm code}$ 才是要在 AWGN 上实测的经验曲线，两者切勿混淆。

> 澄清一个容易混淆的点：经典分集分析里的 $C/\det\cdot\mathrm{snr}^{-N_t}$ 公式是对 $P_{\rm out}$ 再做**高 SNR 渐近展开**的结果，只在深尾（BLER $\ll$ 1%）成立。在 10%/1% 工作点不适用的是那个渐近展开，不是 outage 本身——$P_{\rm out}$ 在任意 SNR 下都是定义良好的精确对象。

### 2.2 两个工程修正（从"近似适用"到"标准做法"）

1. **每 RE 互信息必须用 $I_{\rm QAM}$，不能用高斯输入的 $\log_2(1+\mathrm{snr})$**。16QAM 在高 SNR 处饱和于 4 bit；用 $\log(1+\mathrm{snr})$ 会给"把能量集中到少数超强 RE"的设计虚报收益。
2. **有限码长修正**：真实 LDPC 的 BLER 不是在 $I=R$ 处的阶跃，而是 $\mathrm{BLER}\approx\mathbb E_{\mathbf h}[f_{\rm code}(I(\mathbf h))]$，$f_{\rm code}$ 是该 MCS 在 AWGN 上的 BLER–互信息曲线（平滑 sigmoid）。所有候选共用同一 MCS/码长 ⇒ $f_{\rm code}$ 相同 ⇒ 全体候选的曲线平移同一个 coding gap，**排序保序**。

这套"链路性能压缩为每 RE 互信息再平均"的方法就是 3GPP 系统级仿真的标准 PHY abstraction（MIESM/RBIR），经过大规模校准验证，不是本项目自创的近似。

**〔释疑 · 补注 2026-07-07〕有限码长修正到底在修什么？需要它吗？**

*要修的对象*：§2.1 的理想模型把 BLER 写成硬阶跃 $\mathbb 1[I(\mathbf h)<R]$——块互信息一过码率线就 100% 译对、差一点就 100% 译错。这是"无限码长理想码"的行为。真实有限长 LDPC 没这么利落的门限：$I(\mathbf h)$ 略高于 $R$ 时仍有残余误块率，略低时偶尔也能救回。

*怎么修*：把阶跃换成一条平滑的 $f_{\rm code}(I)$——含义是"当本块有效互信息为 $I$ 时，这套具体 MCS/码长/译码器的译码失败概率"。这条曲线**在 AWGN 上实测一次**得到：拿真实的 16QAM + 这套 LDPC + 这个码长，在若干 SNR 上跑出 BLER；因为 AWGN 上块互信息就是确定值 $I_{\rm QAM}(\mathrm{snr})$，于是可把 BLER 重新参数化成"BLER 关于 $I$"的函数，得到一条平滑 sigmoid（码越长越陡、越逼近理想阶跃）。这就是所谓 **"BLER–互信息曲线"**。修正后 $\mathrm{BLER}\approx\mathbb E_{\mathbf h}[f_{\rm code}(I(\mathbf h))]$：先把频选信道压缩成一个有效互信息 $I(\mathbf h)$，再去这条 AWGN 参考曲线上读 BLER——就是 MIESM/RBIR 的标准做法。

*"候选"指什么*：指参与比较的各个 $\mathbf V$ 设计（CDD 各变体 C1–C7、N-series、chirp、一般 $\mathbf V$……）。整篇文档的任务就是在这些候选里挑最优；此处提"候选"是要说明修正对它们的**相对排序**有无影响。

*只比不同 $\mathbf V$ 的话——不需要做*（这正是第 2 条最后半句的意思）。所有候选共用同一 MCS/码长/译码器 ⇒ 共用同一条 $f_{\rm code}$ ⇒ 有限码长带来的只是一个对所有候选**共模**的水平平移（coding gap，约 1–2 dB）。共模平移不改变排序，所以**为筛选/排序 $\mathbf V$**，直接用便宜的理想 $P_{\rm out}$（§2.5 第 1 条 MC outage）即可，跳过修正。

*什么时候才真需要它*：① 要预测**绝对** BLER / 具体 dB 工作点（去与真实链路仿真对齐）——coding gap 是真实偏移，不能忽略；② 怀疑两候选的 $I(\mathbf h)$ 左尾**交叉**（而非整体平移）——那时 sigmoid 平滑可能在 gap 内翻转排序（见 §2.3"失效方式"），需直接画两条 $P_{\rm out}$ 曲线查有无交叉。一句话：这一条与其说是"必须多跑一步"，不如说是"论证了我们为什么可以省掉这一步"。

### 2.3 适用前提清单

1. 块内信道恒定（或至少：一个码字对应一次信道抽样；上 Doppler 后公式形式不变，但需重新审视）；
2. 码字均匀铺满被评估的频带（若码字只占子带，$I(\mathbf h)$ 只在该子带上平均）;
3. 接收机逐 RE 解调出 LLR，无跨 RE 联合检测（单层 OFDM 满足）；
4. 用实际调制的 BICM 互信息（§2.2 第 1 条）；
5. 所有候选共享同一 MCS/码长/译码器（保证 coding gap 相同）；
6. **理想 CSI**。这是最重要的边界：估计误差进来后 LLR 失配，PHY abstraction 不再干净。所以此指标**只筛分集轴**；CE 轴由 §3 的指标单独把关，两轴互不替代。

失效方式：若两个候选的 $I(\mathbf h)$ 分布左尾交叉（而非整体平移），有限码长平滑可能在 gap 内翻转排序——直接画两条 $P_{\rm out}$–snr 曲线即可检查有无交叉。

### 2.4 为什么 Gram 指标不够

#### 2.4.1 Hadamard 界到底界住了什么

恒模下每列范数固定 $\|\mathbf v_n\|^2=K$，Hadamard 不等式给出

$$
\det(\mathbf V^H\mathbf V)\;\le\;\prod_n\|\mathbf v_n\|^2 = K^{N_t},
$$

等号当且仅当列两两正交。**这是 Gram 行列式的代数上界，不是 BLER 的上界。** $\det$ 与 BLER 的联系只出现在深尾渐近式 $P_{\rm out}\simeq c\cdot\mathrm{snr}^{-N_t}/\det(\mathbf V^H\mathbf V/K)$ 里：$\mathrm{snr}^{-N_t}$ 是分集阶数（斜率），$\det$ 是 coding gain（平移量）。10%/1% 工作点不在渐近区，所以 log-det 只是"深尾正确、有限 SNR 未必"的代理。

#### 2.4.2 充分统计量是 $\mathbf V\mathbf V^H$，不是 $\mathbf V^H\mathbf V$

$\mathbf g=\mathbf V\mathbf h$ 是零均值高斯向量，其全部统计由 $K\times K$ 协方差 $\mathbf V\mathbf V^H$ 完全决定，因此 $I(\mathbf h)$ 的分布也由它决定。逻辑链：

$$
\text{BLER}\;\Leftarrow\;P_{\rm out}\;\Leftarrow\;I(\mathbf h)\text{ 的分布}\;\Leftarrow\;\mathbf V\mathbf V^H\;(K\times K)\;\Leftarrow\;\mathbf V .
$$

Gram $\mathbf V^H\mathbf V$（$N_t\times N_t$）与 $\mathbf V\mathbf V^H$ 共享非零特征值，但**不含特征向量**——即不含"这 $N_t$ 个模式沿频率怎么排布"的信息。log-det、$\lambda_{\min}$、条件数、互相干全是 Gram 的函数，共享这个盲区。

#### 2.4.3 反例：同 Gram，不同分布

$N_t=2$，行向量取 $(1, e^{j\theta_k})$，则 $g_k=h_1+e^{j\theta_k}h_2$，列正交 $\Leftrightarrow \sum_k e^{j\theta_k}=0$。两个设计：

- **A**：$\theta_k$ 沿频率均匀扫满 $[0,2\pi)$；
- **B**：$\theta_k\in\{0,\pi\}$ 各占一半。

两者 Gram 都精确等于 $K\,\mathbf I_2$，一切谱指标无差。但 B 的等效信道全带只取两个值 $|h_1\pm h_2|^2$（恰好是两个独立的指数随机变量），$I(\mathbf h)$ 是两点平均；A 则连续扫过 $h_1,h_2$ 的所有相对相位。两者 $I(\mathbf h)$ 分布不同 ⇒ **Gram 不是充分统计量**。（至于哪个更好需要算——这正是必须用 outage 直接算、而不是用矩阵指标推断的理由。）

#### 2.4.4 Gram 之外的第二根轴：自平均

把第 $k$ 行的相位 $(\phi_{k,1},\dots,\phi_{k,N_t})$ 看作 $N_t$ 维相位环面上一条随 $k$ 前进的轨迹，则 $I(\mathbf h)$ 是沿这条轨迹对函数 $f(\cdot)=I_{\rm QAM}(\mathrm{snr}|\cdot|^2)$ 取平均。轨迹扫得越快越密（**自平均**越好），$I(\mathbf h)$ 越接近其环面平均

$$
I_\infty(\mathbf h)=\mathbb E_{\boldsymbol\theta\sim U}\Big[I_{\rm QAM}\big(\mathrm{snr}\,\big|\textstyle\sum_n h_n e^{j\theta_n}\big|^2\big)\Big],
$$

这是一切恒模设计的集中上限（由 Jensen，$I_\infty(\mathbf h)$ 不超过 MRC 的 $I_{\rm QAM}(\mathrm{snr}\|\mathbf h\|^2)$）。方差越小左尾越薄，10%/1% BLER 越好。CDD 家族内"轨迹速度"就等于 delay 大小。

**项目内证据**：C 系列 grid-delay CDD 的 Gram 全部精确等于 $K\mathbf I_{N_t}$、log-det 无差，但 C7（step=64 samples，最大 delay）的 ideal-CSI 1% 尾部最好（result-021 异常现象节）。差别正是自平均，Gram 完全看不见。这也是 KNOWLEDGE B7（delay 对 BLER 非单调）中"增益来源"的那一半机理；"回落"的那一半在 §3.3。〔2026-07-09：下方补注把推导做严之后，本段归因需要修正——平坦模型下 C 系列各步长应打平，见补注第 5 步与预测 P7。〕

**〔释疑 · 补注 2026-07-09〕推导：$I(\mathbf h)$ 是怎么逼近 $I_\infty(\mathbf h)$ 的（Weyl 和分解）**

下面把"扫得越快越密 ⇒ 越接近环面平均"从比喻变成可计算的定理，并顺带修正上一段"C 系列证据"的归因。全程固定一次信道实现 $\mathbf h$（推导对每个 $\mathbf h$ 逐点成立），使用 §1.1 平坦模型 $\mathbf g=\mathbf V\mathbf h$。

**第 1 步（$I(\mathbf h)$ = 环面函数沿轨迹的经验平均）**。记第 $k$ 行相位向量 $\boldsymbol\phi_k=(\phi_{k,1},\dots,\phi_{k,N_t})$，定义相位环面函数 $f_{\mathbf h}(\boldsymbol\theta):=I_{\rm QAM}\big(\mathrm{snr}\,\big|\sum_n h_ne^{j\theta_n}\big|^2\big)$，则

$$
I(\mathbf h)=\frac1K\sum_{k=0}^{K-1}f_{\mathbf h}(\boldsymbol\phi_k),\qquad
I_\infty(\mathbf h)=\int_{[0,2\pi)^{N_t}}f_{\mathbf h}(\boldsymbol\theta)\,\frac{d\boldsymbol\theta}{(2\pi)^{N_t}} .
$$

问题化为经典**等分布**问题：设计点列 $\{\boldsymbol\phi_k\}$ 上的经验平均何时逼近整个环面的均匀平均。

**第 2 步（Fourier 展开 ⇒ 主恒等式）**。$f_{\mathbf h}$ 光滑，作环面 Fourier 展开 $f_{\mathbf h}(\boldsymbol\theta)=\sum_{\mathbf m\in\mathbb Z^{N_t}}c_{\mathbf m}(\mathbf h)\,e^{j\mathbf m^{\!\top}\boldsymbol\theta}$；零模系数恰是环面平均，$c_{\mathbf 0}=I_\infty(\mathbf h)$。代入第 1 步、交换求和（系数绝对可和，合法）：

$$
\boxed{\;I(\mathbf h)-I_\infty(\mathbf h)\;=\;\sum_{\mathbf m\ne\mathbf 0}c_{\mathbf m}(\mathbf h)\,W_K(\mathbf m),\qquad
W_K(\mathbf m):=\frac1K\sum_{k=0}^{K-1}e^{j\,\mathbf m^{\!\top}\boldsymbol\phi_k}\;}
$$

误差被**精确分离**成两个互不掺和的因子：Fourier 系数 $c_{\mathbf m}(\mathbf h)$ 只依赖（$\mathbf h$、snr、调制），与设计无关；**Weyl 和** $W_K(\mathbf m)$ 只依赖设计轨迹，与（$\mathbf h$、snr）无关。三条结构事实：

1. **模式筛选**：$f_{\mathbf h}$ 对全体分支同加常相位不变（$|g|^2$ 不变）⇒ 只有 $\sum_n m_n=0$ 的模式系数非零；最低非零阶是二阶成对模式 $\mathbf m=\mathbf e_n-\mathbf e_{n'}$（来自 $|g|^2$ 的交叉项 $h_nh_{n'}^*$）。
2. **系数衰减**：$I_{\rm QAM}$ 光滑有界 ⇒ $|c_{\mathbf m}(\mathbf h)|$ 随阶数 $\|\mathbf m\|_1$ 快速衰减 ⇒ 误差由**低阶模式主导**。
3. **Weyl 和范围**：$|W_K(\mathbf m)|\le1$，取 1 当且仅当该模式沿轨迹恒定（完全没被平均掉）。

所以"自平均好"的精确含义是：**让所有低阶非零模式的 $|W_K(\mathbf m)|$ 都小**——这就是 $I(\mathbf h)$ 逼近 $I_\infty(\mathbf h)$ 的全部机制。顺带得到 Gram 的准确地位：二阶 Weyl 和恰是归一化 Gram 非对角元，

$$
W_K(\mathbf e_n-\mathbf e_{n'})=\tfrac1K(\mathbf V^{H}\mathbf V)_{n'n},
$$

即 **Gram 正交 ⟺ 二阶模式全灭；Gram 是 Weyl 和族的"二阶切片"**。本节说"Gram 看不见自平均"的精确含义 = 它看不见 $\|\mathbf m\|_1\ge4$ 的模式（与 §2.6 第一层/第二层的划分一一对应）。

**第 3 步（CDD 特化：Dirichlet 核与"带内圈数"）**。CDD 轨迹是直线 $\boldsymbol\phi_k=k\boldsymbol\omega$（$\omega_n=-2\pi\Delta f\,\tau_n$），Weyl 和成为几何级数：

$$
|W_K(\mathbf m)|=\frac{1}{K}\left|\frac{\sin\!\big(\pi K\Delta f\,\tau_{\mathbf m}\big)}{\sin\!\big(\pi\Delta f\,\tau_{\mathbf m}\big)}\right|,
\qquad \tau_{\mathbf m}:=\sum_n m_n\tau_n\ \ (\text{模式 }\mathbf m\text{ 的复合时延}).
$$

$K\Delta f\,|\tau_{\mathbf m}|$ 就是该模式相位横跨整个带宽转过的**圈数**；远离整数共振时 $|W_K(\mathbf m)|\approx1/(\pi\times\text{圈数})$。这给出"扫得快"的定量版：**delay 越大 ⇒ 低阶模式的复合时延越大 ⇒ 带内圈数越多 ⇒ Weyl 和按 $O\!\big(1/(K\Delta f|\tau_{\mathbf m}|)\big)$ 衰减**。（共振判据是 $\Delta f\tau_{\mathbf m}$ 到最近整数的距离，即复合时延模 $1/\Delta f$——这正是 §3.3 导频混叠在数据栅格 $S_f{=}1$ 上的同构版本。）

**第 4 步（DFT 栅格特化：精确归零 + 共振残差）**。栅格 delay $\tau_n=j_n/(K\Delta f)$（$j_n$ 整数）时 $K\Delta f\,\tau_{\mathbf m}=j_{\mathbf m}:=\sum_nm_nj_n\in\mathbb Z$，Dirichlet 核塌成指示函数：

$$
W_K(\mathbf m)=\mathbb 1\big[\,j_{\mathbf m}\equiv0\ (\mathrm{mod}\ K)\,\big]
\;\Longrightarrow\;
I(\mathbf h)=I_\infty(\mathbf h)+\!\!\sum_{\mathbf m\ne\mathbf 0,\ j_{\mathbf m}\equiv0\,(K)}\!\!c_{\mathbf m}(\mathbf h).
$$

非共振模式**精确归零**（不是趋零），残差只剩**共振模式**——沿轨迹恒定、永远平均不掉的方向。逐阶清点：

- **二阶**：$j_n$ 两两不同 ⇒ 无共振（= 全带列正交的又一读法）。
- **四阶**：共振 ⟺ 非平凡**加性四元组** $j_a+j_c\equiv j_b+j_d\ (\mathrm{mod}\ K)$（含退化型 $2j_a=j_b+j_c$）。**等差集 $j_n=s\cdot n$ 因 $j_{n-1}+j_{n+1}=2j_n$ 携带一整批与步长 $s$ 无关的永久共振**。几何图像：换步长只是让同一条闭合曲线 $t\mapsto(1,e^{jt},e^{j2t},\dots,e^{j(N_t-1)t})$ 走快走慢、采样重数不同，**轨道本身不变** ⇒ 平坦模型下 $I(\mathbf h)$ 的分布与 $s$ 无关（步长间差异只剩 $\|\mathbf m\|_1\gtrsim10$ 的高阶碎屑，可忽略）。
- **与 §2.6 的 ρ-语言严格对接**：$\sum_{k,l}|\rho_{kl}|^4=\frac{K^2}{N_t^4}\times\#\{(a,b,c,d):j_a+j_c\equiv j_b+j_d\ (\mathrm{mod}\ K)\}$——**四阶矩 = 加性四元组计数**。$N_t=8$：等差集 344 个（其中平凡解 $2N_t^2-N_t=120$），**Sidon（$B_2$）集恰好只剩平凡的 120 个** ⇒ 把最低存活共振推到 $\ge6$ 阶。

**第 5 步（两级图景；对本节与 §4.2 的修正）**。推导把"扫得快、扫得密"拆成两个不同的箭头：

$$
I(\mathbf h)\ \xrightarrow{\ \text{速度：带内圈数}\ \uparrow\ }\ \bar I_{\rm orbit}(\mathbf h)\ \xrightarrow{\ \text{结构:低阶共振}\ \to\ \varnothing\ }\ I_\infty(\mathbf h),
$$

其中 $\bar I_{\rm orbit}(\mathbf h)$ 是 $f_{\mathbf h}$ 沿轨迹闭包（轨道）的平均。"快"只负责第一个箭头——收敛到**轨道平均**（栅格上一步到位、精确达成）；能否继续逼近环面平均由**轨道的低阶共振结构**决定，与速度无关。三条推论：

1. 本节上文"CDD 家族内'轨迹速度'就等于 delay 大小"只对 **off-grid** delay 成立（大步长加速第一个箭头的收敛）；栅格等差家族内改步长什么都不改变（第二个箭头卡在同一批 AP 共振上）。
2. **C 系列（全为栅格等差）在平坦模型下 $I(\mathbf h)$ 分布应彼此重合**，故"C7 尾部最好 = 自平均"的归因在平坦模型内不成立，须重检（候选解释：$\tau_{\rm phy}$ 展宽破坏其 64 子载波周期性、有限码长效应、trial 噪声）——登记为预测 P7。
3. 真正能压低共振残差的旋钮是 delay 集的**加性结构**：无非平凡四元关系的 **Sidon 型栅格集**预期尾部不劣于任何等差集（系数衰减启发，待 P7 检验）。可行性示例（同时满足栅格正交 + 模 $N_p{=}24$ 折叠互异，未优化）：$j_n\in\{0,1,3,7,12,20,30,65\}$；其混叠圆折叠最小间距 1 格 = 57.9 ns，仍 $\gg\tau_{\rm phy}{=}5$ ns，但劣于 C7 的 173.6 ns——分集残差与 CE 折叠余量的联合优化留给 plan。

### 2.5 分集指标的三层选择

1. **主指标：MC outage**。采样 $\mathbf h\sim\mathcal{CN}(0,\mathbf I_{N_t})$ 共 $10^5\sim10^6$ 次，对每个样本算 $I(\mathbf h)$，输出 $P_{\rm out}$–snr 曲线及 10%/1% 插值 SNR。8 维采样 + 矩阵乘，**每候选亚秒级，无需链路仿真、无需噪声、无需 LDPC**。
2. **解析近似（机理分析用）**：恒模 ⇒ 每个 $|g_k|^2$ 边缘分布相同（均值 $N_t$ 的指数）⇒ $\mathbb E[I]$ 对所有候选相同，差异全在方差与高阶矩。$\mathrm{Var}(I)=\frac1{K^2}\sum_{k,l}\mathrm{Cov}_{kl}$，每个协方差项只依赖相关系数 $|\rho_{kl}|$（$\rho_{kl}=(\mathbf V\mathbf V^H)_{kl}/N_t$），一维函数族可预表；$P_{\rm out}\approx Q((\mathbb E[I]-R)/\sigma_I)$。
3. **谱指标（log-det、$\lambda_{\min}$、条件数）只做粗筛**：谱严重退化（$\lambda_{\min}\approx0$）的候选直接扔；谱打满的候选之间它们无分辨力（§2.4.4）。**不得作为验收指标**（KNOWLEDGE B6 的强化版）。

**〔释疑 · 补注 2026-07-07〕为什么均值相同？"候选"？"一维函数族可预表"？**

*为什么 $\mathbb E[I]$ 对所有候选相同*：$\mathbb E_{\mathbf h}[I]=\frac1K\sum_k\mathbb E[I_{\rm QAM}(\mathrm{snr}|g_k|^2)]$ 只依赖每个 $|g_k|^2$ 的**边缘分布**。而 $g_k=\sum_n V_{k,n}h_n$ 是 $N_t$ 个独立 $\mathcal{CN}(0,1)$ 各乘一个单位模 $V_{k,n}$ 再相加——相位只旋转不改方差，故 $g_k\sim\mathcal{CN}(0,N_t)$、$|g_k|^2\sim\mathrm{Exp}(\text{均值 }N_t)$，**对每个 $k$、每个恒模候选都一模一样**。所以"均值"这根轴对所有恒模 $\mathbf V$ 完全无分辨力，差异全部落到**方差与高阶矩**，即落到子载波间的**联合**结构（$\mathbf V\mathbf V^H$ 的非对角 $\rho_{kl}$）上。这就是"分集 = 让 $I(\mathbf h)$ 更集中、左尾更薄"的解析版陈述。

*"候选"*：同上，指各个待比较的 $\mathbf V$ 设计。

*"一维函数族可预表"*：对圆对称复高斯，$\mathrm{Cov}(|g_k|^2,|g_l|^2)=|\mathbb E[g_kg_l^*]|^2=N_t^2|\rho_{kl}|^2$，**精确**只依赖单个标量 $|\rho_{kl}|\in[0,1]$。再过一层 $I_{\rm QAM}$ 后仍然如此：

$$
\mathrm{Cov}\big(I_{\rm QAM}(\mathrm{snr}|g_k|^2),\,I_{\rm QAM}(\mathrm{snr}|g_l|^2)\big)=C(|\rho_{kl}|;\mathrm{snr}),\qquad C(r;\mathrm{snr})=\sum_{m\ge1}a_m(\mathrm{snr})^2\,r^{2m},
$$

（$a_m$ 是把 $I_{\rm QAM}(\mathrm{snr}\cdot)$ 按指数变量的 Laguerre 多项式展开的系数，只由 $I_{\rm QAM}$ 与 snr 定）——这是一条**单调增、凸**的一维曲线。"一维函数族" = 每个 snr 对应一条以 $r=|\rho_{kl}|$ 为自变量的曲线（snr 是族参数）；"可预表" = 只需在 $r\in[0,1]$ 上把 $C$ **预先打表一次**，此后任一候选的 $\mathrm{Var}(I)=\frac1{K^2}\sum_{k,l}C(|\rho_{kl}|)$ 只靠查表求和即得，无需对每个候选重新积分或跑 MC。机理一眼可见：$\mathrm{Var}(I)$ 小 ⟺ 各 $|\rho_{kl}|$ 小（子载波去相关）⟺ 自平均好 ⟺ 左尾薄。末尾 $P_{\rm out}\approx Q((\mathbb E[I]-R)/\sigma_I)$ 是对 outage 的高斯（CLT）近似，只在分布主体准、深尾不保真，故标"机理分析用"；真正验收仍以第 1 条 MC outage 为准。

### 2.6 从指标反推设计：理论上什么样的 V 使分集最大化

> **〔补注 2026-07-07〕** 本节回答"能否在仿真/筛选之前，先从分集指标本身从理论上探明 $\mathbf V$ 的设计方向"。**结论前置**：分集轴上的最优结构本质已被"恒模 + rank $N_t$"两条约束锁死；一般 $\mathbf V$ 相对最优 CDD 的可赢空间是**高阶、微小**的。这从理论侧独立佐证了 §4.1 塌缩命题，也解释了为什么全局把一般 $\mathbf V$ 的机会让给了 CE 轴 / 透明滑窗 regime（§3.4、§5.3）。

目标函数（固定均值下）：最小化 $\mathrm{Var}(I(\mathbf h))=\frac1{K^2}\sum_{k,l}C(|\rho_{kl}|)$。对角项（$k=l$，共 $K$ 个，$|\rho|=1$）对所有候选相同，可分离出去；真正能设计的是**非对角相关** $\{|\rho_{kl}|\}_{k\ne l}$。$C$ 单调增且凸（见上条释疑）。于是"分集最优化" = **在恒模、rank $N_t$ 约束下，把非对角相关既压小、又铺平**。分三层看：

**第一层（二阶 → 全带列正交是最优，且这是个守恒律）**。非对角相关的总能量不是自由的，被恒模钉死：

$$
\sum_{k,l}|\rho_{kl}|^2=\frac1{N_t^2}\|\mathbf V\mathbf V^H\|_F^2=\frac1{N_t^2}\|\mathbf V^H\mathbf V\|_F^2=\frac1{N_t^2}\sum_i\lambda_i^2,
$$

而恒模给出 $\sum_i\lambda_i=\operatorname{tr}(\mathbf V^H\mathbf V)=KN_t$ 固定。$\sum_i\lambda_i^2$ 在 $\sum_i\lambda_i$ 固定下、当且仅当**所有 $\lambda_i=K$（即 $\mathbf V^H\mathbf V=K\mathbf I_{N_t}$，全带列正交）时取最小** $N_tK^2$。代入得非对角总相关能量的下确界

$$
\sum_{k\ne l}|\rho_{kl}|^2\;\ge\;\frac{K^2}{N_t}-K=\frac{K(K-N_t)}{N_t}\quad(\text{正交时取等}).
$$

由 $C(r)=a_1^2r^2+O(r^4)$，**到二阶 $\mathrm{Var}(I)$ 就正比于这份守恒的总相关能量** ⇒ 任何非正交恒模设计在二阶上严格更差，而**所有全带列正交设计到二阶完全打平**。这给了 §2.4.1 Hadamard 上界一个新读法：全带列正交不只是 coding-gain 的谱上界，它**同时是分集方差的二阶最小点**；二阶方差大致对应 **10% outage 点**。CDD 取 DFT 栅格 delay 即落此点（§4.2 第 1 条）。

**第二层（高阶 → 旁瓣铺平 / 自平均 / 等角紧框架）**。二阶打平之后，正交设计之间的胜负全在 $C$ 偏离纯二次的高阶项 $\sum_{m\ge2}a_m^2|\rho_{kl}|^{2m}$，大致对应更深的 **1% 尾部**。$C$ 凸 ⇒ 在 $\sum_{k\ne l}|\rho_{kl}|^2$ 固定下，$\sum_{k\ne l}C(|\rho_{kl}|)$ 由 Jensen 在**各 $|\rho_{kl}|$ 尽量相等**时最小——即把守恒的相关能量**均匀摊到所有子载波对**、不让任何一个 lag 扛大相关。设计语言：$\mathbf V$ 的行相关**旁瓣要低而平**；极限理想是各行构成**等角紧框架**（tight = 全带正交那一半，equiangular = 所有 $|\rho_{kl}|$ 等幅那一半）。这正是 §2.4.4 自平均的代数对偶：旁瓣越平、相位环面扫得越均匀，$I(\mathbf h)$ 越集中、左尾越薄。（本层的可计算代理由 §2.4.4 补注第 4 步给出：四阶矩 $\sum_{k\ne l}|\rho_{kl}|^4$ = 非平凡加性四元组计数，Sidon 型 delay 集使其达到理论最小；注意**同轨道**设计的各阶矩全部相同，故逐设计比较要用该矩而非旁瓣峰值，见 P6 的度量修正。）

**第三层（$N_t$ 决定的地板，不可逾越）**。把 $\mathbf V$ 扫到无穷快、旁瓣压到极平，$I(\mathbf h)\to I_\infty(\mathbf h)$（§2.4.4 的相位环面平均），它仍是 $\mathbf h\in\mathbb C^{N_t}$ 的随机函数，$\mathrm{Var}(I_\infty(\mathbf h))>0$ 是任何恒模 $\mathbf V$ 都下不去的**分集地板**，只由 $N_t$（与 snr、$R$）定。要击穿地板只能破恒模（不等功率注水，但需 CSIT、且改变了问题定义）或加大 $N_t$。

**落到设计方向的三条结论**：

1. **必要一阶条件 = 全带列正交**（$\mathbf V^H\mathbf V=K\mathbf I_{N_t}$）。任何违背它的恒模设计都在二阶方差上白白吃亏；这是筛 $\mathbf V$ 的第一道硬门槛（也解释了 §2.5 第 3 条为何谱严重退化的候选可直接扔）。
2. **正交之上 = 压平行相关旁瓣**（等价：相位环面上均匀快扫 = 自平均最优）。CDD 因相关只依赖 lag（平稳），**结构上无法把旁瓣压成全平**（有限 delay 集的自相关必有起伏），只能用 incommensurate / folded-grid delay 去逼近（C7 即此，§4.2 第 2–3 条）；**一般（非平稳）$\mathbf V$ 唯一的分集抓手正在这里——它能做出比任何平稳 CDD 更平的旁瓣**。但这只值二阶以上的高阶收益。
3. **因此：分集轴单独并不构成上一般 $\mathbf V$ 的理由。** 二阶已被正交锁死、地板已被 $N_t$ 锁死，一般 $\mathbf V$ 能赢的只是"有限 $K$ 下把环面铺得比直线轨迹更匀"这一点高阶差，量级微小（与 §4.1"没有可赢空间"、result-021"C7 仅在 1% 尾部略胜"一致）。**一般 $\mathbf V$ 真正的自由度不在分集轴，而在 CE 轴 / 透明滑窗 regime**——这就是本文把探索重心放在 §3.4、§5.3 的理论依据。

可证伪的具体预测（供后续 plan 采用，登记为 §7 P6）：**在所有全带列正交的恒模候选内部，其 10%/1% outage SNR 应落在一个很窄的带内（二阶被正交锁死打平），且带内排序与"行相关旁瓣峰值 $\max_{k\ne l}|\rho_{kl}|$（越低越好）"单调**。若实测排序与旁瓣峰值无关，则第二层机理需修正。

---

## 3. 信道估计侧理论：接收机决定坐标轴

### 3.1 接收机两类，坐标轴两套

- **V-aware（非透明）**：UE 知道 $\mathbf V$，用 matched LMMSE（即 Algorithm 1，协方差 $\mathbf R_g=\mathbf V\mathbf R_h\mathbf V^H$）。
- **V-agnostic（透明）**：UE 不知道 $\mathbf V$，把等效信道当"普通物理信道"处理。实际 UE 的典型形态：
  - LS + 频域线性/样条插值；
  - 鲁棒 Wiener/LMMSE：用**假设的**通用协方差（典型：$[0,\tau_w]$ 上均匀 PDP，$\tau_w$ 为接收机假设的时延窗）；
  - DFT 加窗：导频 LS → IFFT 到时延域 → 只留窗内抽头 → FFT 回来；
  - 作用域分**宽带**（一个窗覆盖全带）与**滑窗**（如 4RB 局部窗滑动）两种，结论截然不同（§3.4）。

v1 的平滑性指标只在 V-agnostic 下有意义：一阶 roughness $J_{\rm rough,1}=\sum_{n,k}|\phi_{k+1,n}-\phi_{k,n}|^2$ 惩罚 group delay 总量（相邻子载波相位差 $\propto$ 局部 group delay，$\tau^{\rm group}_{k,n}=-\frac{\Delta\phi_{k,n}}{2\pi\Delta f}$）；二阶 roughness $J_{\rm rough,2}=\sum_{n,k}|\phi_{k+2,n}-2\phi_{k+1,n}+\phi_{k,n}|^2$ 惩罚 group delay 的变率（任何 CDD 恒为 0，故它实际度量"偏离 CDD 结构的程度"）。**对 V-aware 接收机，这两个量与估计质量无关。**

### 3.2 V-aware：估计 = 参数估计 + 合成，一切由 $\mathbf V_P$ 决定

近平坦物理信道下，导频观测为

$$
\mathbf y_P=\mathbf V_P\,\mathbf h+\mathbf n,\qquad \mathbf V_P\in\mathbb C^{N_p\times N_t},
$$

$\mathbf V_P$ 是 $\mathbf V$ 在导频子载波行上的子采样。数据 RE 上的"插值"其实是合成：$\hat{\mathbf g}_D=\mathbf V_D\hat{\mathbf h}$（$\mathbf V_D$：数据 RE 行的子矩阵）。**因此估计质量由 $\mathbf V_P$ 的可辨识性/条件数决定，与 $\mathbf V$ 平不平滑无关**——$\mathbf V_D$ 变化再快，只要 $\hat{\mathbf h}$ 准，合成就准。

$\mathbf h$ 的 LMMSE 误差协方差为 $\mathbf E=(\mathbf R_h^{-1}+\mathbf V_P^H\mathbf V_P/\sigma_{LS}^2)^{-1}$。恒模 ⇒ $\operatorname{tr}(\mathbf V_P^H\mathbf V_P)=N_pN_t$ 固定；$\operatorname{tr}(\mathbf E)=\sum_i(1+\lambda_i/\sigma_{LS}^2)^{-1}$ 是特征值的对称凸（Schur 凸）函数 ⇒ **在迹约束下当且仅当特征值全相等、即导频域列正交 $\mathbf V_P^H\mathbf V_P=N_p\mathbf I_{N_t}$ 时最小**。这就是 CE 侧的"Hadamard 时刻"，与分集侧的全带列正交平行。

一般情形（不必平坦）有精确闭式 NMSE，**不需要 Monte Carlo**：

$$
\mathrm{NMSE}
=\frac{\operatorname{tr}\!\big(\mathbf R_{DD}-\mathbf R_{DP}(\mathbf R_{PP}+\sigma_{LS}^2\mathbf I)^{-1}\mathbf R_{PD}\big)}
{\operatorname{tr}(\mathbf R_{DD})},
$$

其中 $\mathbf R_{PP},\mathbf R_{DP},\mathbf R_{DD}$ 是 $\mathbf R_g$ 在导频/数据行列上的子块。逐候选逐 SNR 确定性算出，比 trial 仿真快约两个数量级、零方差，且可对相位求导做梯度优化。

### 3.3 导频混叠：V-aware 下唯一的硬约束

**推导**：一条时延为 $\tau$ 的分量在导频栅格（子载波 $k=pS_f$）上的相位序列是 $e^{-j2\pi(pS_f)\Delta f\,\tau}$。把 $\tau$ 换成 $\tau+m/(S_f\Delta f)$，多出的因子是 $e^{-j2\pi pm}=1$——**在导频上完全不可分辨**。这就是 Nyquist 采样定理用在频域：把 $H(f)$ 当"信号"，其"带宽"是时延扩展，间隔 $S_f\Delta f$ 的频域采样能无歧义支撑的时延范围是

$$
\tau_{\rm alias}=\frac{1}{S_f\,\Delta f}
=\frac{1}{24\times30\ \text{kHz}}=1.389\ \mu s\;(\approx170.7\ \text{samples @122.88 MHz}).
$$

混叠不是估计算法的缺陷，是**观测算子的零空间**——任何估计器（含 matched LMMSE）都无法从导频数据分开相差 $m\cdot\tau_{\rm alias}$ 的分量。影响分三档：

1. **精确混叠 → 不可辨识**：$\mathbf V_P$ 对应列共线、掉秩；LMMSE 靠先验兜底不发散，但歧义方向纯靠猜，且"导频上一致的两个假设在数据 RE 上不一致"——零噪声也有合成误差。
2. **接近混叠 → 条件数病态**：折叠后 delay 靠近但不重合，噪声沿弱方向放大。**特征信号：高 SNR 下 NMSE 出现 floor 而非随 SNR 线性下降**（误差主项是结构歧义不是噪声）。例：24PRB N8（result-021，12→24 dB 仅从 4.8e-2 缓降至 2.3e-2）。
3. **折叠但折得开 → 仍可辨识**：delay span 超过 $\tau_{\rm alias}$ 不必然完蛋，只要 delay 集在模 $\tau_{\rm alias}$ 的"**混叠圆**"上仍彼此分开。

**实例（C7，值得反复看）**：step=64 samples=520.8 ns，8 分支 span 3.65 µs $\gg\tau_{\rm alias}$，表面"混叠了"。但 $520.8/1388.9=3/8$，折叠位置为 $173.6\ \text{ns}\times(3n\bmod 8)$——3 与 8 互素，是个置换 ⇒ 折叠后恰为混叠圆上**均匀 8 点栅格**。换到导频域：相位 $e^{-j2\pi\cdot3pn/8}$，24 个导频走 3 个整周期 ⇒ **导频域精确列正交**；同时全带相位 $e^{-j2\pi kn/64}$，$576=9\times64$ 整周期 ⇒ **全带也精确列正交**。C7 同时做到 §2 的自平均最优（扫相最快）与 §3 的 CE 无代价——见预测 P1。

### 3.4 V-agnostic：按时延跨度付费

通用估计器不知道"$\mathbf g$ 其实只有 $N_t$ 个自由度"，它必须假设一个时延窗并估计窗内全部抽头，未知数 $W\approx$（复合时延跨度 $\times$ 处理带宽）个时延域抽头。处理增益 $\approx N_p/W$。与 V-aware 对比，这是全文档最重要的一句话：

> **V-aware 按"分支数"付费**：未知数 $=N_t$，CE 增益 $\approx N_p/N_t$，与 delay 摆多大无关（只要混叠圆上折得开）⇒ **分集与 CE 无折中**。
> **V-agnostic 按"时延跨度"付费**：未知数 $=W\propto$ 人工 delay 跨度，而分集偏要求跨度大 ⇒ **折中真实存在**。v1 命题 4 的正确适用域在此。

**〔释疑 · 补注 2026-07-09〕$W$ 是什么？"复合时延跨度 × 处理带宽"怎么来的？为什么增益是 $N_p/W$？**

*$W$ = 估计器要解的未知数个数。* V-aware 接收机手里有结构模型 $\mathbf y_P=\mathbf V_P\mathbf h+\mathbf n$，未知数就是 $N_t$ 个分支系数。V-agnostic 接收机没有这个模型，只能用"万能模型"：把 $\mathbf g$ 当一条普通频率选择性信道，其全部参数是时延域冲激响应，并假设它支撑在某个时延窗内。这个万能模型的参数就是窗内的**时延域抽头**，$W$ = 抽头个数 = 未知数个数。

*两个名词。* **处理带宽 $B_{\rm proc}$**：估计器一次联合处理的频域孔径——宽带窗接收机 $B_{\rm proc}=K\Delta f$（48PRB 为 17.28 MHz），滑窗接收机 $B_{\rm proc}=$ 窗宽（4RB 为 1.44 MHz）。Fourier 对偶给出**时延分辨率** $=1/B_{\rm proc}$（分别为 57.9 ns / 694 ns）：孔径内的观测只能把时延域分辨到这么细的格子，比格子更细的结构在该孔径内不可分辨、也不必参数化。**复合时延跨度 $\tau_{\rm span}$**：等效信道 $\mathbf g$ 的时延域支撑长度。"复合"指人工与物理的叠加——每个分支的物理冲激响应被 $\mathbf V$ 平移到自己的 $\tau_n$ 处，总支撑 $\approx(\max_n\tau_n-\min_n\tau_n)+\tau_{\rm phy}$（工程上再加定时误差余量）。于是

$$
W\;\approx\;\frac{\text{支撑长度}}{\text{分辨率}}\;=\;\tau_{\rm span}\times B_{\rm proc}\ \text{个分辨格，每格一个复未知数}.
$$

（与 §3.3 拼成同一套采样理论的两半：**导频间隔定范围**——无混叠范围 $\tau_{\rm alias}=1/(S_f\Delta f)$；**孔径定分辨率**——$1/B_{\rm proc}$。两者相除 = 该孔径内导频观测的总自由度 $=\tau_{\rm alias}B_{\rm proc}=$ 窗内导频数，自洽。）

*为什么处理增益 $\approx N_p/W$。* 处理增益指估计噪声相对单导频 LS（方差 $\sigma_{LS}^2$）的压低倍数。一句话推导：$N_p$ 个导频观测里的噪声均匀散布在 $N_p$ 维观测空间；线性估计把观测投影到 $W$ 维信号子空间（窗内抽头张成），只有落进子空间的那份噪声留下来，占比 $W/N_p$ ⇒ 每个数据 RE 的估计噪声方差 $\approx(W/N_p)\sigma_{LS}^2$，增益 $=N_p/W=$ **观测数/未知数**。DFT 加窗估计器把这句话执行得最字面：$N_p$ 个导频 LS → IFFT 得 $N_p$ 个时延 bin（分辨率 $1/B_{\rm proc}$、范围 $\tau_{\rm alias}$）→ 只留窗内 $W$ 个 bin（其余置零，那里只有噪声）→ FFT 回频域。宽带窗下还有个等价写法：$N_p/W=\tau_{\rm alias}/\tau_w$——**增益 = 假设时延窗占无混叠范围的比例的倒数**。

*更准确的计费表*（对引文的精化）：V-aware 的未知数其实是 $\min(N_t,\ \tau_{\rm span}B_{\rm proc})$——分支数**封顶**（$\mathbf V$ 已知时多余的时延精细结构不必估）；V-agnostic 的未知数是 $\tau_{\rm span}B_{\rm proc}$——**不封顶**。两张计费表只在大跨度 regime 分道，而大跨度恰是分集想要的——折中因此真实存在。

*当前体制下的两个硬事实*：① $S_f=24\Rightarrow\tau_{\rm alias}=1.389\ \mu s$。透明接收机若保守地把 $\tau_w$ 设成 CP 量级（30 kHz SCS 下 $\approx2.34\ \mu s$），则 $\tau_w>\tau_{\rm alias}$，窗盖满整个无混叠范围，$W\to N_p$、增益 $\to1$——稀疏 comb 下透明接收机**必须敢假设短窗**才有增益，而人工 delay 跨度每加大一分，就直接吃掉一分"敢假设短窗"的本钱。② C7 折叠后均匀铺满混叠圆（§3.3 实例）——对 V-aware 是导频正交的最优排布，对 V-agnostic 宽带窗却是**最坏情形**（折叠支撑 = 整个 $\tau_{\rm alias}$，$W\to N_p$，增益 $\approx1$）：同一个 $\mathbf V$，V-aware 拿 $24/8=3$（4.8 dB），V-agnostic 拿 $\approx1$（0 dB）——两种计费方式的极端演示。

作用域再分两种，结论不同：

- **宽带窗接收机**：代价由**全局** delay 跨度决定 ⇒ 凡 group delay 都落在窗内的设计代价相同，窗内摆 $N_t$ 个分开的 delay 就是 CDD——塌缩论证在透明 regime 重演（预测 P4）。
- **滑窗接收机**（真实 UE 常见形态，如 4RB 窗）：代价由**每个局部窗口内**的 delay 扩展决定 ⇒ 出现 CDD 结构性做不到的自由度（§5.3 的"局部聚簇 + 全局轮转"）。数字感：$S_f=24$ 时 4RB 窗内每 DMRS 符号只有 2 个导频，局部模型只有约 2 个自由度（近似"局部单一 delay + 常数"）——**局部 delay 聚簇是硬要求**。

**〔释疑 · 补注 2026-07-09〕宽带窗/滑窗到底指什么？窗大小怎么影响估计？工程上怎么取舍？**

*先分清两个"窗"，本节同时在用，别混*：① **时延窗 $\tau_w$**（模型假设）——估计器假设信道冲激响应支撑的时延范围，决定未知数 $W$；② **频域处理窗 $B_{\rm proc}$**（作用域）——一次联合处理的频谱孔径，决定时延分辨率与窗内可用导频数。"宽带窗 / 滑窗"分类指的是②。

*宽带窗接收机*：把整个分配带宽（全部 $N_p$ 个导频）喂给同一个估计器，一次输出全带估计。典型形态是宽带鲁棒 Wiener（假设 $[0,\tau_w]$ 均匀 PDP 的 LMMSE）或全带 DFT 加窗。特征：分辨率最细（$1/K\Delta f$）、可用导频最多，但隐含要求**全带共用一个模型**（信道统计沿频率恒定）。付费项是**全局跨度**：凡 group delay 都落进 $\tau_w$ 的设计代价全相同 ⇒ 透明设计空间在此塌缩回 CDD（P4）。

*滑窗接收机*：用一个固定宽度的局部窗（如 4RB）沿频率滑动，每个位置只用**窗内导频**估计窗中心附近的数据 RE，滑完拼出全带。这是真实 UE 的主流形态（原因见下面取舍第 4 条）。付费项是**每个窗口各自的局部跨度**。窗内账本：导频只有 $N_p^{\rm win}=B_{\rm proc}/(S_f\Delta f)$ 个（4RB、$S_f{=}24$：每 DMRS 符号 2 个），可辨识未知数至多 2 个 ⇒ 局部模型被钉死在"一个复幅度 + 一个公共斜率（局部单一 group delay）"；窗内分辨率也只有 $1/B_{\rm proc}=694$ ns ⇒ 各分支 delay 若聚在远小于此的范围内，窗内看就是**一个**抽头（2 参数模型恰好够用，估得准）；若散开几百 ns~µs，窗内既分辨不开也参数化不起（2 个观测撑不起多抽头模型），差额直接变成模型失配误差、高 SNR 出 floor。这就是"局部聚簇是硬要求"的定量根据，也是 CDD（delay 全带常数，每个窗都看到全跨度、都付全额）在滑窗下必然吃亏、而"局部聚簇 + 全局轮转"有生存空间的原因。

*窗大小如何影响估计（两个旋钮各自的 bias–variance）*：

| 旋钮 | 取小了 | 取大了 |
|---|---|---|
| 时延窗 $\tau_w$ | 截掉窗外真实信道能量 ⇒ 泄漏偏差，高 SNR error floor（与 SNR 无关的地板） | $W\uparrow$ ⇒ 增益 $N_p/W\downarrow$；$\tau_w\to\tau_{\rm alias}$ 时增益 $\to1$，超过则混叠折入、模型本身出错 |
| 频域处理窗 $B_{\rm proc}$ | 窗内导频少 ⇒ 降噪弱；分辨率粗 ⇒ 局部模型维度低，只容得下"局部近单抽头"的设计 | 违反"窗内一个模型"：透明场景撞 PRG 边界（预编码跳变被当成信道估进去）、非平稳 $\mathbf V$ 的局部斜率漂移；复杂度 $\uparrow$（MMSE 矩阵 $\mathcal O(N_p^{\rm win\,3})$ / 滤波器变长）；流水线延迟 $\uparrow$ |

*工程取舍的常用规则*：

1. **$\tau_w$ 下界** = 预期复合支撑（物理 delay spread + 人工 delay 跨度 + 定时同步误差余量），常取"PDP 能量 ~99% 点 + 余量"；**上界** = 明显小于 $\tau_{\rm alias}$，否则无增益。零知识的静态保守选择是 CP 量级——如上一条补注所示，$S_f{=}24$ 下这等于放弃全部处理增益，所以实际 UE 常做**自适应窗**：先粗估 PDP，再按能量收窗。
2. **$B_{\rm proc}$ 上界**：透明场景的硬顶是 **PRG**——预编码只保证在 PRG 内连续，跨界的窗把预编码跳变当信道估，必错（这正是 PRG bundling 存在的意义：向 UE 承诺一个可安全联合处理的最大孔径）；次约束是窗内统计近平稳（对我们的设计即"局部斜率近似恒定"）。
3. **$B_{\rm proc}$ 下界**：窗内导频数 $\ge$ 局部模型维度 × 条件数裕量（经验 $\ge2\times$）。$S_f=24$ 时 4RB = 每符号 2 导频已是"还能滑得动"的最小窗，再窄连局部斜率都估不出。
4. **实现偏好**：滑窗只需小矩阵（2×2~8×8）、定点友好、可流水线，且对统计失配天然更鲁棒（错也只错一个窗）；宽带 MMSE 是 $N_p\times N_p$ 的全带求逆。这两点是真实 UE 偏好滑窗的主因——**也是 Track B 把主参考定为 4RB 滑窗、次参考定为宽带鲁棒 Wiener 的依据**（§3.5、§5.3 第 1 条）。

*数字感（把两条计费表钉在同一把尺上，48PRB、$S_f{=}24$、$N_p{=}24$）*：V-aware matched 增益 $24/8=3$（4.8 dB），与跨度无关；透明小 delay CDD（复合跨度 ~200 ns）：$W\approx0.2\,\mu s\times17.28\,\mathrm{MHz}\approx3.5$，增益 $\approx24/3.5\approx6.9$（8.4 dB）——小跨度设计下透明接收机反而降噪更狠（V-aware 此时同样受益：未知数封顶公式 $\min(N_t,\tau_{\rm span}B_{\rm proc})$ 两边一致），代价是频选弱、分集差；跨度往上加，透明增益按 $\tau_{\rm alias}/\tau_w$ 反比缩水，直到 C7 类满圆折叠时 $\to1$（0 dB），而 V-aware 始终握住 3。**从 8.4 dB 滑向 0 dB 的这条连续拨盘，就是 Track B 基线"扫跨度与 $N_{\rm eff}$"要画出的 Pareto 曲线**（§5.3 第 2 条）。

### 3.5 CE 指标选择

| 轨道 | 验收指标 | 辅助诊断 |
|---|---|---|
| V-aware（Track A） | 闭式 matched NMSE（§3.2，实际 DMRS pattern） | $\mathrm{cond}(\mathbf V_P)$、折叠 delay 分布图 |
| V-agnostic（Track B） | **固定参考接收机**下的闭式失配 MSE：估计器矩阵 $\mathbf A$ 由假设先验定死，$\mathrm{MSE}=\operatorname{tr}\big(\mathbf A(\mathbf R_{PP}+\sigma^2\mathbf I)\mathbf A^H-\mathbf A\mathbf R_{PD}-\mathbf R_{DP}\mathbf A^H+\mathbf R_{DD}\big)/\operatorname{tr}(\mathbf R_{DD})$，同样是确定量 | roughness / group delay 分布（解释用，非验收） |

参考接收机是 Track B **问题定义的一部分**，必须在 plan 里写死（§5.3 第 1 条）。

---

## 4. CDD 分析

### 4.1 CDD 何时最优：折中塌缩命题

**命题（塌缩）**：若 (i) 物理信道带内近平坦（$\tau_{\rm phy}\ll 1/(K\Delta f)$），(ii) 接收机为 V-aware matched LMMSE，(iii) 导频充足 $N_p\ge N_t$ 有余量，则分集与 CE 两个目标可**同时**达到各自上界，Pareto 前沿塌缩为一点，且由 DFT-grid CDD 达到：

- **全带**：取 $\tau_n=n/(K\Delta f)$，则 $\mathbf V^H\mathbf V=K\mathbf I_{N_t}$（Hadamard 等号，分集谱上界）；
- **导频域**：导频上相位为 $e^{-j2\pi pn/N_p}$，$n=0..N_t{-}1<N_p$ 是 $N_p$ 点 DFT 的不同列 ⇒ $\mathbf V_P^H\mathbf V_P=N_p\mathbf I_{N_t}$（§3.2 Schur 凸性的最优点）。

⇒ 在这个 regime 里，一般性 $\mathbf V$ 在这两根轴上**没有任何可赢的空间**；CDD 家族内剩余的自由轴只有自平均（delay 大小），而折叠栅格技巧（§4.2 第 2 条）可以免费拿到它。

**这条命题一举解释 result-021 的全部现象**：N-series 在 ideal-CSI 赢（排布/自平均不差）而 estimated-CSI 全倒扣（$\mathbf V_P$ 可辨识性差，随机局部斜率在 $S_f=24$ 子采样下条件数恶化）；"矩阵指标匹配 CDD 参考点"低估了 CDD（Gram 打满的 CDD 之间还有自平均轴）。即 KNOWLEDGE B5/B6 的机理版（待预测 P2 定量确认后回写 FINDINGS）。

### 4.2 CDD 的设计原则（即使在它的最优 regime，delay 也不能乱选）

1. **delay 取有效带宽 DFT 栅格的整数倍** $m/(K\Delta f)$ ⇒ 全带精确列正交。栅格步长：57.9 / 77.2 / 115.7 ns（48/36/24 PRB，即 7.11/9.48/14.22 samples）。
2. **折叠位置设计**：delay 集在模 $\tau_{\rm alias}$ 的混叠圆上尽量均匀/最大分离。系统构造：步长 $\delta=(a/N_t)\,\tau_{\rm alias}$，$\gcd(a,N_t)=1$ ⇒ 折叠后为均匀 $N_t$ 点栅格 ⇒ 导频域精确正交（C7 即 $a=3$ 的实例）。
3. **满足 1–2 后，delay 步长尽量大** ⇒ 相位轨迹扫得快 ⇒ 自平均好 ⇒ 1% 尾部好。B7 的"非单调"由此化解：增大 delay 的收益来自自平均，回落来自折叠位置恶化；折叠栅格构造把回落避开了。〔修正 2026-07-09，据 §2.4.4 补注推导〕**"步长尽量大"只在 off-grid 情形加速自平均收敛**；满足第 1 条（栅格）后，等差家族的平坦模型自平均与步长无关（同一条轨道，走快走慢而已），残差由 delay 集的**加性结构**决定——Sidon 型集预期优于任何等差集（待 P7 检验）。本条的"越大越好"相应降级为 off-grid 经验规则。
4. **物理时延余量**：每个折叠格点被 $\tau_{\rm phy}$ 展宽，格点间距（$\tau_{\rm alias}/N_t=173.6$ ns）须 $\gg\tau_{\rm phy}$。当前 5 ns，余量 ~35×；$\tau_{\rm phy}$ 到 100 ns 量级时此原则开始吃紧（§4.3 第 2 条）。

### 4.3 CDD 的缺陷条件（= 一般性 V 还活着的 regime）

1. **导频稀缺（$N_p\lesssim N_t$）**：更疏的 DMRS（$S_f\ge48$）、更多分支（$N_t=16$）、或更窄带（8~12 PRB）。$\mathbf V_P$ 无法对所有列正交，设计变成"欠采样下分配可辨识性"的压缩感知式问题，CDD 的均匀分配未必最优。
2. **大物理时延扩展（100~300 ns TDL）**：折叠格点展宽吃掉混叠余量，人工 delay 摆放必须与物理 PDP 联合优化；非均匀 delay 分配、频变结构开始有内容。当前 5 ns 是对 CDD 最友好的设定。
3. **V-agnostic 滑窗接收机**：CDD 的 delay 是全带常数，每个局部窗口都看到全部 delay 跨度、都付全额 CE 代价；"局部聚簇"机制 CDD 结构性做不到（§5.3）。
4. **失配鲁棒性**：matched LMMSE 依赖准确的 $\mathbf R_g$（即准确的 delay/V 知识）；对失配的敏感度是一根从未测过的轴。

---

## 5. 一般性 V 的探索：两条赛道

### 5.1 竞争格局与判定不等式

| | 透明（UE 不知 $\mathbf V$） | 非透明（UE 知 $\mathbf V$） |
|---|---|---|
| **基线** | 透明小 delay CDD、PRG precoder cycling | QC 显式 CDD + V-aware wideband CE |
| **挑战者** | 一般 $\mathbf V$（§5.3，机制存在，开放） | 一般 $\mathbf V$（§4.1，塌缩 regime 内无空间） |

**研究者判据（2026-07-04）**：若存在透明设计，① 同 DMRS 开销下明显优于透明基线（透明 CDD、PRG cycling），且 ② 性能超过或接近 QC 非透明方案（接近时以"标准实现更简单/零 spec 影响"取胜），则透明方向值得做；两个方向都可探索。

诚实预期：透明接收机信息严格更少，纯性能反超 QC 不太可能；Track B 现实的赢法是"接近 + 零 spec 改动"。反过来，若 Track A 基准远超一切透明方案，等于用数据论证了 QC 非透明提案的必要性——**两种结局都是有价值的结论**。

### 5.2 Track A（非透明）：收尾定基准

- **目标**：确认塌缩命题（§4.1），产出"最优非透明 = folded-grid CDD 家族"的 10%/1% BLER 基准线，作为全项目的尺子。
- **手段**：plan-022 修正三条后照跑——(a) $L_{CE}$ 用闭式 NMSE（§3.2）替代 Monte Carlo；(b) 显式加入 DFT-grid / folded-grid CDD 参考点（§4.2 构造）；(c) 分集侧补 MC outage 指标（§2.5），与 log-det 并列输出以便对照淘汰。
- **预期**：最优一般候选逼近但不超过最强 CDD（±0.2 dB 内打平）。走 `GOALS.md` "假设不成立"分支即视为本轨道**封闭**，非失败。

### 5.3 Track B（透明）：开放前沿

1. **参考接收机（问题定义的一部分，写进 plan 前须研究者确认）**。默认建议：主参考 = **4RB 滑窗 MMSE**（假设局部均匀 PDP、不知 $\mathbf V$）；次参考 = **宽带鲁棒 Wiener**（均匀 PDP $[0,\tau_w]$）——后者专门用于检验预测 P4（宽带窗下透明空间同样塌缩），两个一起跑结论才说得圆。
2. **基线**：透明小 delay CDD（扫 delay 跨度与**有效分支数** $N_{\rm eff}$——多根天线共享同一 delay 即合并为一支，是"分集 ↔ CE"的离散折中拨盘）；PRG precoder cycling（扫 PRG 大小）。
3. **核心机制（CDD 做不到的自由度）：局部聚簇 + 全局轮转**。让各分支的局部 group delay 在每个接收窗口内聚在一起（局部近似单一 delay ⇒ 2 个导频就估得准，CE 几乎无代价），而聚簇的相位模式沿频率轮转（全局把 $\mathbf h$ 的 $N_t$ 维都扫到 ⇒ 分集保留）。CDD 的 delay 全带固定，每个局部窗口都看到全部跨度——结构上做不到这件事。
4. **候选家族（极限形态 = 相位连续的 precoder cycling）**：每段内所有分支共享同一斜率、只差每段一组常相位 $\theta_{s,n}$，段间相位连续拼接。它把两个透明基线连成连续统：段数 $N_{\rm seg}\to$ PRG 数且去掉连续性 = PRG cycling；$N_{\rm seg}=1$ = 透明 CDD。chirp 是它的"连续轮转"版本；N-series 的重生条件是斜率字母表按局部窗口 delay 预算封顶。（v1 §7–10 的候选结构细节可在归档版查阅，但评价一律换用 v2 指标。）
5. **历史重估**：N-series 死在错误的 regime——result-021 用 V-aware wideband matched LMMSE 考核它们，那是 CDD 不可战胜的主场（§4.1）；它们从未在滑窗透明接收机（自己的设计初衷所在）下被公平测过。
6. **先导轮 = 纯数值扫描，不跑链路**：在（固定参考接收机闭式 MSE，10%/1% outage SNR）平面上画三个家族的 Pareto——透明 CDD（扫跨度、$N_{\rm eff}$）、连续 cycling（扫 $N_{\rm seg}$、拼接宽度）、chirp/封顶 N-series。两轴指标都是确定性计算，判断路径短。有空间再上 BLER 链路验证。
7. **验收（写 plan 时量化）**：同 DMRS 开销下，① estimated-CSI（参考接收机）性能优于两个透明基线（幅度阈值待定）；② 距 Track A 基准 $\le\delta_{\rm bench}$（建议 0.3~0.5 dB，**待研究者确认**）。

### 5.4 与既有结论（KNOWLEDGE B4~B7）的关系

- **B4**（N-series 矩阵层面赢）：所用指标本身已降级（Gram 不充分，§2.4），结论保留但意义弱化。
- **B5**（ideal 赢、estimated 倒扣）：机理更新为"$\mathbf V_P$ 可辨识性/混叠圆条件数"，非"非平稳协方差难插值"；待 P2 定量确认后改写 FINDINGS §B5。
- **B6**（矩阵指标只能筛不能当目标）：强化并具体化——正确目标 = MC outage（分集轴）+ 闭式 NMSE（CE 轴）。
- **B7**（delay 非单调）：机理 = 自平均收益 vs 折叠位置恶化，folded-grid 构造（§4.2）化解。

---

## 6. 修订后的理论命题

- **命题 1（保留）**：在"每分支一个全带固定 delay"约束下，CDD 是自然最优结构；delay 取 DFT 栅格时全带列正交。
- **命题 2（保留但降格）**：一般 $\mathbf V$ 可行域包含 CDD，故其 Pareto 前沿不劣于 CDD——集合论上平凡，不指明 regime 就没有内容。
- **命题 3（修订）**：一般 $\mathbf V$ 严格优于 CDD **只可能**发生在 §4.3 列出的缺陷 regime；在塌缩 regime（§4.1 条件 i–iii）内不可能。
- **命题 4′（替换 v1 命题 4）**：折中的存在性取决于接收机的"计费方式"——V-aware 按分支数付费（增益 $N_p/N_t$，与 delay 无关，**无折中**）；V-agnostic 按时延跨度付费（增益 $N_p/W$，**折中真实存在**）。
- **命题 5（新）**：Gram $\mathbf V^H\mathbf V$ 不是有限 SNR BLER 的充分统计量（§2.4.3 构造反例 + C 系列同 Gram 不同尾部的实测证据）。
- **命题 6（新，塌缩）**：§4.1 全文；折叠栅格 CDD 在塌缩 regime 同时达到分集谱上界与 CE 误差下界。

## 7. 可证伪预测清单

理论要接受检验，以下预测任何一条被证伪，都要回到本文档修正对应小节，不得绕过：

- **P1**：C7（folded-grid 构造的近似实例）的 estimated-CSI BLER 应为全场最强或并列最强。——plan-022 直接检验。
- **P2**：result-021 的 12 个 NMSE 倒扣点可被闭式 NMSE / $\mathrm{cond}(\mathbf V_P)$ 完全复现排序（把 B5 从经验观察升级为被解释的机制）。——便宜的回算诊断。
- **P3**：MC outage 曲线与 ideal-CSI BLER 曲线近似恒定水平间距（coding gap）。——用 result-021 已有曲线校准。
- **P4**：宽带窗 V-agnostic 接收机下，透明设计空间同样塌缩到 CDD（凡 group delay 在窗内代价相同）。——Track B 次参考接收机检验。
- **P5**：滑窗 V-agnostic 接收机下，"局部聚簇 + 全局轮转"家族在（固定接收机 MSE，outage）平面上严格支配透明 CDD 与 PRG cycling 的连线。——Track B 先导扫描检验；这是 Track B 存亡的判定实验。
- **P6**：在所有全带列正交（$\mathbf V^H\mathbf V=K\mathbf I_{N_t}$）的恒模候选内部，10%/1% outage SNR 落在很窄带内（二阶方差被正交锁死打平），且带内排序与行相关旁瓣结构单调。〔度量修正 2026-07-09〕排序量改用**四阶矩 $\sum_{k\ne l}|\rho_{kl}|^4$ 超出平凡值的部分**（= 非平凡加性四元组计数，§2.4.4 补注第 4 步），不用旁瓣峰值 $\max|\rho_{kl}|$——后者会把同轨道等价设计误判为不同（C7 在 lag 64 处 $|\rho|=1$，峰值最差，却应与小步长打平）。——§2.6 第二层机理的检验；证伪则 §2.6 第二层需修正。
- **P7**：平坦模型 MC outage 下，**栅格等差 CDD 家族（C 系列各步长，含 C7）的 $I(\mathbf h)$ 分布应在数值精度内彼此重合**（同轨道论证，§2.4.4 补注第 4–5 步）；且满足同等栅格正交与折叠约束的 Sidon 型 delay 集，其 outage 尾部不劣于任何等差集。若 C7 在 MC outage 中复现出显著尾部优势 ⇒ 补注推导有误，回改 §2.4.4；若不复现 ⇒ result-021 的 C7 归因改为平坦模型外效应（$\tau_{\rm phy}$ 展宽/有限码长/trial 噪声），同步改写 §2.4.4 证据段与 §4.2 第 3 条。——与 P3 同批用 MC outage 便宜检验。

---

## 附录 A：规范符号表（全项目唯一权威）

| 符号 | 含义 |
|---|---|
| $N_t$ | 发射分支数（当前 8）。**不用** $N$ |
| $K$ | 有效子载波数（576/432/288） |
| $k,\ n$ | 子载波索引 / 发射分支索引 |
| $\Delta f$ | 子载波间隔（30 kHz） |
| $\mathbf V,\ V_{k,n}=e^{j\phi_{k,n}}$ | $K\times N_t$ 频域相位矩阵，恒模 |
| $\mathbf h$ | 底层分支信道 $\in\mathbb C^{N_t}$，基准 $\mathcal{CN}(0,\mathbf I_{N_t})$ |
| $\mathbf g=\mathbf V\mathbf h$ | 等效频域信道 $\in\mathbb C^K$。**不用** $\mathbf h_{\rm eff}$ |
| $\mathbf R_h,\ \mathbf R_g$ | $\mathbf h$ / $\mathbf g$ 的协方差；$\mathbf R_g=\mathbf R_{\rm phy}\odot(\mathbf V\mathbf V^H)$ |
| $\mathbf R_{\rm phy}$ | 物理信道频域协方差；$\odot$ 为 Hadamard 积 |
| $S_f$ | DMRS comb 频域间隔（子载波数，当前 24） |
| $N_p=K/S_f$ | 每 DMRS 符号导频数（24/18/12） |
| $P,\ D$（下标） | 导频 / 数据 RE 集合；$\mathbf V_P,\mathbf V_D,\mathbf R_{PP},\mathbf R_{DP}$ 等为对应行/列子阵 |
| $\tau_n$ | 第 $n$ 分支 CDD delay；$\tau_{\rm phy}$：物理 delay spread |
| $\tau_{\rm alias}=1/(S_f\Delta f)$ | 导频无混叠时延周期（当前 1.389 µs） |
| $\tau^{\rm group}_{k,n}$ | 局部 group delay $=-\frac{\phi_{k+1,n}-\phi_{k,n}}{2\pi\Delta f}$ |
| $\tau_w$ | V-agnostic 接收机假设的时延窗长（模型假设，决定 $W$） |
| $\tau_{\rm span}$ | 等效信道的复合时延跨度 $\approx(\max_n\tau_n-\min_n\tau_n)+\tau_{\rm phy}$（§3.4 补注） |
| $B_{\rm proc}$ | 估计器一次联合处理的频域孔径（处理带宽）：宽带窗 $=K\Delta f$，滑窗 $=$ 窗宽；时延分辨率 $=1/B_{\rm proc}$ |
| $N_p^{\rm win}=B_{\rm proc}/(S_f\Delta f)$ | 频域处理窗内的导频数（滑窗账本的观测预算） |
| $W\approx\tau_{\rm span}B_{\rm proc}$ | V-agnostic 接收机须估计的时延域抽头数 = 未知数个数；处理增益 $\approx N_p/W$（§3.4 补注） |
| $\sigma^2,\ \sigma_{LS}^2$ | 噪声方差 / 导频 LS 估计噪声方差（双 DMRS 符号平均后 $=\sigma^2/2$） |
| $\mathrm{snr}$ | 发射符号信噪比 |
| $I_{\rm QAM}(\cdot)$ | 实际调制（16QAM）的 BICM 每 RE 互信息曲线 |
| $I(\mathbf h)$ | 块平均互信息 $=\frac1K\sum_k I_{\rm QAM}(\mathrm{snr}|g_k|^2)$ |
| $I_\infty(\mathbf h)$ | $I(\mathbf h)$ 的相位环面平均（自平均集中上限） |
| $R$ | 码率对应谱效率（当前 $\approx2.16$ bit/RE） |
| $P_{\rm out}(\mathrm{snr},R)$ | 块 MI outage 概率 $=\Pr[I(\mathbf h)<R]$ |
| $C(|\rho_{kl}|;\mathrm{snr})$ | 子载波对互信息协方差核 $=\mathrm{Cov}(I_{\rm QAM}(\mathrm{snr}|g_k|^2),I_{\rm QAM}(\mathrm{snr}|g_l|^2))$；单调增凸、只依赖 $|\rho_{kl}|$ 的一维函数（§2.5、§2.6） |
| $\rho_{kl}=(\mathbf V\mathbf V^H)_{kl}/N_t$ | 子载波 $k,l$ 等效信道的归一化复相关系数 |
| $\boldsymbol\phi_k$ | 第 $k$ 行相位向量 $(\phi_{k,1},\dots,\phi_{k,N_t})$，设计在相位环面上的轨迹点 |
| $f_{\mathbf h}(\boldsymbol\theta)$ | 相位环面函数 $I_{\rm QAM}(\mathrm{snr}\,|\sum_nh_ne^{j\theta_n}|^2)$；$I(\mathbf h)$ 是它沿轨迹 $\{\boldsymbol\phi_k\}$ 的平均（§2.4.4 补注） |
| $c_{\mathbf m}(\mathbf h)$ | $f_{\mathbf h}$ 的环面 Fourier 系数（$\mathbf m\in\mathbb Z^{N_t}$）；$c_{\mathbf 0}=I_\infty(\mathbf h)$ |
| $W_K(\mathbf m)$ | 设计轨迹的 Weyl 和 $\frac1K\sum_ke^{j\mathbf m^{\!\top}\boldsymbol\phi_k}$；二阶模式时 $=\frac1K(\mathbf V^H\mathbf V)_{n'n}$ |
| $j_n$ | 第 $n$ 分支 delay 的 DFT 栅格整数索引（$\tau_n=j_n/(K\Delta f)$） |
| $\tau_{\mathbf m},\ j_{\mathbf m}$ | 模式 $\mathbf m$ 的复合时延 $\sum_nm_n\tau_n$ / 复合栅格索引 $\sum_nm_nj_n$（§2.4.4 补注） |
| $\bar I_{\rm orbit}(\mathbf h)$ | $f_{\mathbf h}$ 沿轨迹闭包（轨道）的平均；栅格 CDD 下 $I(\mathbf h)$ 精确等于它 |
| $f_{\rm code}(\cdot)$ | 该 MCS 在 AWGN 上的 BLER–互信息曲线（coding gap 来源） |
| $\mathbf E$ | $\mathbf h$ 的 LMMSE 误差协方差 $=(\mathbf R_h^{-1}+\mathbf V_P^H\mathbf V_P/\sigma_{LS}^2)^{-1}$ |
| $\lambda_i,\ \sigma_i(\mathbf V)$ | Gram 特征值 / $\mathbf V$ 奇异值（$\lambda_i=\sigma_i^2$） |
| $\kappa(\mathbf V),\ \mu(\mathbf V)$ | 条件数 / 列互相干（仅粗筛用） |
| $J_{\rm div},\ L_{CE}$ | 分集目标 / CE 代价（v2 下分别指 outage 指标与闭式 NMSE） |
| $J_{\rm rough,1},\ J_{\rm rough,2}$ | 一阶/二阶相位 roughness（仅 V-agnostic 解释用） |
| $N_{\rm eff}$ | 有效分支数（多天线共享同一 delay 合并计数） |
| $N_{\rm seg},\ \theta_{s,n}$ | 分段数 / 第 $s$ 段第 $n$ 分支的段常相位（连续 cycling 家族） |
| $\delta_{\rm bench}$ | Track B "接近非透明基准"的容差（建议 0.3~0.5 dB，待定） |

> 引入新符号时：先加进本表（含一句话含义），再在别处使用。历史文档与本表冲突时以本表为准。
