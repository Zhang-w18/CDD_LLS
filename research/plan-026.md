# plan-026：CDD residue/lift 协方差适配、厚 Sidon 潜力地图与近似平坦 TDL 复核（已确认）

> 状态：已确认。研究者已确认研究问题、物理路线、搜索预算、正式 trial 预算和数值验收阈值，可以按本 plan 启动实现、验证和正式运行。
>
> 执行回执：`20260724_main` 已完成，结果见 `research/result-026.md` 与 `research/result-026-text.md`；结论待研究者确认。

## 1. 研究问题、范围与可证伪假设

### 1.1 背景与目的

result-024 在 48 PRB、8 Tx / 1 Rx、平坦独立分支信道、DMRS comb 24 和 V-aware matched LMMSE 条件下，确认严格 Sidon delay 集相对等差 QC CDD 的 10%/1% estimated-CSI BLER 改善为 0.33/0.85 dB。result-025 的 delay-matched `/576` 补做只把底层信道替换为 Sionna TDL-A 5 ns，10% 改善缩小为 0.156 dB，1% 改善点估计为 0.090 dB且保守区间跨 0；同时，TDL 加权后二阶相关代理变为 QC 略低，Sidon 四阶代理仍明显较低，Sidon 的 data-RE CE NMSE 较差。

本轮不把严格 Sidon 条件继续当作物理展宽信道的充分设计定理，而是回答三个问题：

1. **E1，T2+R2：** 发射端和接收端都知道 TDL-A 5 ns 的完整物理协方差时，固定或放松导频折叠 residue 并优化 lift，能否找到在 $J_{\rm CDD}$、matched CE 和 QAM/BICM outage 上优于经过优化的等差 CDD 的结构化 delay 集？
2. **E2，T1+R2：** 发射端只知道 $T_\epsilon$，接收端仍知道实际 PDP 和候选 $\mathbf V$ 时，厚 Sidon 在哪些预定 TDL profile 和 RMS delay spread 下具有低成本指标优势？
3. **E3，固定设计+R2：** 当 TDL-A 的 RMS delay spread 逐步趋近零时，原 Sidon 相对 QC 的 estimated-CSI BLER 优势是否恢复；在近似平坦 TDL-D/E LoS 条件下，排序是否相同？

E1/E2 的 QAM/BICM 互信息阶段只输出 ideal-CSI outage，不称为 BLER 预测。E3 运行真实信道估计、LLR、LDPC 和译码器，输出 estimated-CSI BLER。

### 1.2 可证伪假设

#### H1：完整协方差下存在有潜力的 covariance-adapted CDD

在 TDL-A 5 ns 的完整协方差、相同 DMRS 和 matched LMMSE 条件下，分级 residue/lift 搜索至少找到一个非等差候选，相对本轮优化等差基线满足：

- 10% ideal-CSI outage 所需 SNR 改善不低于 0.15 dB；
- 16 dB matched data-RE CE NMSE 劣化不超过 1.0 dB；
- 80 dB 参考 SNR 下的近零噪声 CE NMSE相对基线劣化不超过 1.0 dB；若双方线性 NMSE 均低于 `1e-12`，按数值零地板并列处理；
- pilot rank 为 8，估计器数值有限且无预定门槛外失稳。

若没有候选满足全部条件，则 H1 不成立；只改善 $J_{\rm CDD}$ 而 outage、CE 或数值稳定性未通过，不算 H1 成立。

#### H2：存在厚 Sidon有潜力的预定物理区域

在 E2 全部预定场景均完整报告的前提下，至少一个只根据 $T_\epsilon$ 生成的厚 Sidon候选，在至少两个预定场景中相对同知识等级的等差基线满足：

- 10% ideal-CSI outage 所需 SNR改善不低于 0.10 dB；
- 16 dB matched data-RE CE NMSE 劣化不超过 1.0 dB；
- 同时满足预定的 pair-sum 和 pilot-fold 几何约束。

若只有单个场景通过，结论为“单点潜力”；若没有场景通过，H2 不成立。该判定只说明值得进入后续 estimated-CSI 链路，不证明 BLER 改善。

#### H3：近似平坦 TDL 下原 Sidon优势恢复

H3 分为两个物理问题：

- **H3a，NLoS 连续性主假设：** TDL-A 0.1 ns 下，原 Sidon 相对 QC 的 10% estimated-CSI BLER 所需 SNR改善不低于 0.15 dB，且保守 95% 区间不支持反向排序。TDL-A 0.1/1/5 ns 的全部结果用于报告增益随 RMS delay spread 的变化，不以事后单调性作为验收条件。
- **H3b，LoS 诊断：** TDL-D/E 0.1/1 ns 下分别报告 10% BLER 排序。改善区间全为正记为支持，区间全为负记为反向，区间跨 0 记为不可判定。H3b 不用于单独归因“时延扩展趋零”，因为 D/E 同时改变 LoS/NLoS 统计。

1% BLER 在目标附近每候选累计错误块不少于 30 时报告确定性结果，否则只报告先导点和区间。

## 2. 系统模型、知识条件与物理定义

### 2.1 固定链路条件

除本 plan 明确改变的 TDL profile、RMS delay spread 和 CDD delay 外，沿用 result-025 delay-matched 补做条件：

| 项目 | 固定值 |
|---|---|
| 有效带宽 | 48 PRB，$K=576$ active subcarriers |
| FFT / CP | 4096 / 288 samples |
| 子载波间隔 | $\Delta f=30$ kHz |
| 时域资源 | 10 OFDM symbols |
| Tx / Rx / layer | 8 / 1 / 1 |
| DMRS | symbols `[2,7]`，comb 24，offset 0 |
| pilot / data RE | 48 / 5712 |
| 速度 / 载频 | 0 km/h / 3.5 GHz |
| 调制编码 | 16QAM，MCS 8，码率 553/1024，LDPC 最多 8 次迭代 |
| 接收机 | known-$\mathbf V$、known-PDP matched 全带频域 LMMSE |
| DMRS 处理 | 两个静态 DMRS 等效平均，平均后 LS 噪声方差 $N_0/2$ |
| CDD 幅度 | 每分支恒模 1；数据噪声 $N_0=8/\mathrm{SNR}$ |
| 共同随机数 | 同一场景、SNR、trial 内候选共享 TDL、payload、平均 LS noise 和 data noise |

本轮 authoritative delay 单位为 ns。所有候选必须保存 8 个物理人工时延

$$
\boldsymbol\tau=[\tau_0,\ldots,\tau_7]\ {\rm ns},
$$

并由

$$
V_{k,n}
=
\exp(-j2\pi k\Delta f\tau_n)
$$

生成预编码矩阵。搜索限制在 result-025 delay-matched 的 576 点 DFT 栅格：

$$
\tau_n=\frac{j_n}{576\Delta f},
\qquad
j_n\in\{0,\ldots,575\}.
$$

一个栅格为 57.870370 ns。输出以 ns 为主，同时保存辅助字段 `j_n`、`j_n mod 24`、lift、`K=576`、`Delta_f_hz=30000` 和 `/576` 相位定义，防止复现时把同一数字误解释为 `/4096` FFT 整数采样延迟。

### 2.2 residue/lift 参数化和等价候选

定义

$$
j_n=r_n+24q_n,
\qquad
r_n=j_n\bmod24,
\qquad
q_n\in\{0,\ldots,23\}.
$$

导频行相位只由 residue 决定；lift 不改变导频行，但会改变数据 RE 相位、全带阵列因子、$J_{\rm CDD}$ 和 $\mathbf R_{DP}$。

当前各发射分支独立且具有相同 PDP，普通天线标签排列不影响性能。搜索和结果使用以下规范代表去重：

1. 最终 delay 集按圆周公共移位等价类规范化；
2. 8 个 delay 排序后取字典序最小代表；
3. 天线排列不重复计数；
4. 若未来引入不同分支 PDP 或空间相关，不得继续使用该等价化。

### 2.3 发射端与接收端知识

| 实验 | 发射端知识 | 接收端知识 |
|---|---|---|
| E1 | T2：知道 TDL-A 5 ns 完整 PDP/协方差 | R2：知道实际 PDP 和候选 $\mathbf V$ |
| E2 | T1：只知道 $\epsilon$、$T_\epsilon$、零同步误差和预定运行包络，不知道实际 profile | R2：知道评价场景实际 PDP 和候选 $\mathbf V$ |
| E3 | 固定 QC/Sidon，不根据信道适配 | R2：知道实际 profile、delay spread 和候选 $\mathbf V$ |

E2 候选生成代码不得接收实际 TDL profile、抽头功率或完整协方差。实际协方差只允许在生成结束后的评价阶段进入 $J_{\rm CDD}$、CE 和 outage 计算。

## 3. 基线、候选和评价指标

### 3.1 固定历史基线 B0

`B0_QC025` 使用

$$
\mathbf j=[0,9,18,27,36,45,54,63],
$$

即物理人工时延

$$
[0,\ 520.833333,\ 1041.666667,\ 1562.5,\ 2083.333333,\ 2604.166667,\ 3125,\ 3645.833333]\ {\rm ns}.
$$

`S0_Sidon025` 使用

$$
\mathbf j=[0,1,3,7,12,20,30,65],
$$

即

$$
[0,\ 57.870370,\ 173.611111,\ 405.092593,\ 694.444444,\ 1157.407407,\ 1736.111111,\ 3761.574074]\ {\rm ns}.
$$

S0 是 E2/E3 历史候选，不作为 E1 unrestricted 搜索的独立基线。

### 3.2 优化等差基线 B1

等差 CDD 定义为

$$
j_n=na\bmod576,
\qquad
a\in\{1,\ldots,575\}.
$$

消除公共移位、天线排列和产生相同无序集合的重复步长后：

- E1 使用完整 TDL-A 5 ns 协方差，在全部等价类上计算 $J_{\rm CDD}$ 和 CE，保留等差 Pareto 前沿；outage 最优且满足 CE 门槛的点定义为 `B1_T2_best_arithmetic`。
- E2 的公平基线 `B1_T1_arithmetic` 只使用 $T_\epsilon$ 和候选自身的无权几何量选择：先排除不满足 fold 约束的等差集合，再依次最小化无权 $M_2$、无权 $M_4$，最后用规范化 ns delay 的字典序打破完全并列。该规则不读取实际 TDL profile 或协方差。实际 profile 下按 $J_{\rm CDD}$/outage 事后选出的最优等差结果另标记为 `B1_oracle_arithmetic`，只作更强数值参考，不冒充 T1 公平基线。

不得在看到 E3 BLER 后重新选择 B1。

### 3.3 E1 候选

E1 采用两级搜索：

1. `E1_U`：固定最大分离 residue

   $$
   \mathbf r=[0,3,6,9,12,15,18,21],
   $$

   对 lift 做多起点 coordinate exchange/local search。

2. `E1_R`：从 `E1_U` 的预定 top Pareto 起点出发，允许 residue 局部移动 $\pm1,\pm2$；保持 residue 互异、pilot rank 8，并重新计算实际 CE。

3. `E1_O`：同一搜索预算下允许一般 residue/lift 的多起点数值参考。它是 unrestricted optimized CDD 上界，不是结构化候选基线。

搜索使用若干预定的 $J_{\rm CDD}$–CE 标量化权重生成不同 Pareto 区域。每个权重使用相同数量起点、最大 sweep 和停止条件；不得只增加某一家族预算。建议初值为每个权重 32 个确定性起点、最多 20 个完整 coordinate sweeps、连续 2 个 sweep 无改善停止。该预算须在正式执行前确认。

### 3.4 E2 候选与场景

固定

$$
\epsilon=0.01,\qquad
T_{\rm margin}=0,\qquad
T_{\rm sync}=0.
$$

零 margin/sync 与 result-025 的零定时误差条件一致；本轮不把结果外推到同步误差非零。

预定场景为：

| 场景 | Profile | RMS delay spread |
|---|---|---:|
| E2-A10 | TDL-A | 10 ns |
| E2-A20 | TDL-A | 20 ns |
| E2-A30 | TDL-A | 30 ns |
| E2-C5 | TDL-C | 5 ns |
| E2-C10 | TDL-C | 10 ns |

对每个场景先从标准化平均 PDP 确定性计算覆盖 99% 能量的最短连续区间宽度 $T_\epsilon$。候选生成器只接收 $T_\epsilon$ 和系统几何，要求

$$
d_{\rm circ}(\tau_a+\tau_c,\tau_b+\tau_d)
>
2T_\epsilon,
\qquad
\{a,c\}\ne\{b,d\},
$$

以及

$$
d_{\rm fold}(\tau_a,\tau_b)>T_\epsilon.
$$

若硬约束不可行，记录不可行证明或穷尽边界，转为最大化两类最小距离的软约束候选；软约束候选不得标记为满足厚 Sidon。

每个 $T_\epsilon$ 档位至少保留：

- 原 Sidon S0；
- QC B0；
- T1 可行等差基线；
- 最大 pair-sum 距离候选；
- 最大 fold 距离候选；
- 两类距离的 Pareto 折中候选。

候选评价时使用实际 TDL 协方差，但不得据此回改生成阶段。全部五个预定场景均报告，不删除负结果。

### 3.5 E3 候选和信道

E3 只比较 `B0_QC025` 与 `S0_Sidon025`：

- E3a：TDL-A RMS delay spread `[0.1,1,5]` ns；
- E3b：TDL-D RMS `[0.1,1]` ns；
- E3b：TDL-E RMS `[0.1,1]` ns。

TDL-A 是 NLoS 近似平坦连续性主实验。TDL-D/E 是 LoS 诊断，不把其结果单独归因于 delay spread。

### 3.6 低成本指标

#### 完整 $J_{\rm CDD}$

对频率平稳物理协方差

$$
R_{{\rm phy},kl}=r_{\rm phy}(k-l)
$$

计算

$$
A_\tau(d)
=
\frac18\sum_{n=0}^{7}
e^{-j2\pi d\Delta f\tau_n},
$$

$$
J_{\rm CDD}
=
2\sum_{d=1}^{575}
(576-d)
C\left(
\left|r_{\rm phy}(d)A_\tau(d)\right|;
\mathrm{snr}
\right).
$$

主参考 SNR 为 16 dB，并额外计算 `[14,16,18]` dB 以检查排序稳健性。核 $C(\cdot;\mathrm{snr})$ 必须用确定性积分、Gauss-Hermite 或预先固定的大样本 Monte Carlo 实现，并先通过独立数值验证；不得用未校准的多项式替代。

同时输出 $M_2^{\rm eff}$、$M_4^{\rm eff}$、最小 pair-sum 圆周距离和最小 fold 距离，仅作机理解释。

#### matched CE

每个候选、场景和 `[14,16,18]` dB 计算：

- 实际 $\mathbf R_{PP}$、$\mathbf R_{DP}$；
- matched data-RE 平均 NMSE及 dB 值；
- pilot rank、condition number、最小奇异值；
- 80 dB 参考 SNR 下的近零噪声 CE NMSE/error floor；线性 NMSE 低于 `1e-12` 时同时保存原值并标记为数值零地板；
- 数值非有限或半正定容差异常标志。

#### E1 二维散点图

E1 必须生成主图 `e1_jcdd_ce_scatter_16db.png`：

- 横轴：16 dB matched data-RE CE NMSE，dB，越左越好；
- 纵轴：$10\log_{10}(J_{\rm CDD}/J_{\rm CDD,B0})$，dB，越低越好；
- 颜色区分 B0/B1、固定 residue、relaxed residue 和 unrestricted reference；
- 叠加全体候选 Pareto 前沿；
- 标注 B0、B1、S0 和进入 outage 的候选；
- 图中不使用 BLER 或 outage 标签。

另保存散点源数据 CSV，包含所有点，不只保存 Pareto 点。

### 3.7 QAM/BICM outage

复用并验证 plan-023 的单位能量 Gray 16QAM BICM 互信息表。谱效率固定为

$$
R=4\times\frac{553}{1024}=2.1602\ {\rm bit/RE}.
$$

对每个入围候选使用实际 TDL 复合信道和共同随机数计算

$$
I_{\rm blk}
=
\frac1K\sum_k
I_{\rm QAM}\left(
\mathrm{snr}\frac{|g_k|^2}{8}
\right),
$$

以及

$$
P_{\rm out}=\Pr[I_{\rm blk}<R].
$$

建议使用 $2\times10^5$ 个共同信道样本、SNR `0:0.5:24` dB，并以对数域插值报告 10%/1% outage SNR。若 1% 尾部有效样本不足，增加共同样本而不改变候选集合或阈值。

## 4. 执行阶段与停止条件

### 4.1 Phase 0：实现与验证

预计新增：

- `cdd_lls/design/cdd_search.py`：CDD 等价类、residue/lift、几何距离、候选去重和搜索；
- `cdd_lls/design/cdd_metrics.py`：$J_{\rm CDD}$、相关矩和厚 Sidon几何指标；
- `tools/run_plan026_cdd_design.py`：E1/E2 确定性搜索和 outage 入口；
- `tools/run_plan026_near_flat_link.py`：E3 BLER 入口；
- `tools/analyze_plan026.py`：Pareto、散点图、目标 SNR、区间和汇总；
- `tests/test_cdd_design.py`、`tests/test_cdd_metrics.py`：新增单元测试。

若已有实现可安全复用，可以减少文件，但不得把可复用搜索逻辑复制到多个专题脚本。

必须验证：

1. ns、$j_n$、相位矩阵三者转换逐元素一致；
2. 公共循环移位和天线排列规范化不会重复计数；
3. residue/lift 重构回原 $j_n$；
4. $d_{\rm circ}$、$d_{\rm fold}$ 与直接枚举一致；
5. B0/S0 的 pair-sum 统计复现 result-025；
6. TDL 复合协方差 Hermitian、半正定且对角功率正确；
7. matched 闭式 NMSE 与小样本 Monte Carlo 差不超过 0.2 dB；
8. $C$ 核和 $J_{\rm CDD}$ 在 $\rho=0$、$\rho\rightarrow1$ 及代表相关系数上通过独立数值检查；
9. 16QAM MI 表和 plan-023 保存结果在容差内一致；
10. 同 seed 重放候选、搜索轨迹、outage 和 E3 smoke 哈希一致。

任一物理定义、转换或协方差验证失败时停止，不进入搜索。

### 4.2 Phase 1：E1 搜索和散点

执行顺序：

1. 穷举并去重全部等差 B1；
2. 固定最大分离 residue 搜索 lift；
3. 对预定 top Pareto 起点做 residue 局部松弛；
4. 运行同预算 unrestricted reference；
5. 保存全部候选、访问次数、起点、停止原因和搜索轨迹；
6. 计算 14/16/18 dB 的 $J_{\rm CDD}$ 和 CE；
7. 生成二维散点、Pareto 前沿和候选表；
8. 按预定规则保留每个家族最多 8 个 Pareto 代表，总 outage 候选不超过 32 个。

若搜索未在预算内收敛，报告当前 best-so-far 和覆盖率，不临时增加单一家族预算。若 unrestricted reference 与固定 residue 候选完全相同，也必须报告。

### 4.3 Phase 2：E1 outage

对 Phase 1 入围候选使用共同 TDL 抽头样本计算 10%/1% outage。输出相对 B0、B1 的成对目标 SNR差。

通过 H1 的候选最多保留 5 个，作为下一轮 estimated-CSI BLER 候选；本 plan 不运行 E1 新候选的正式 LDPC BLER。

### 4.4 Phase 3：E2 厚 Sidon潜力地图

对五个预定场景依次执行：

1. 确定性计算 99% 最短能量区间及 $T_\epsilon$；
2. 只把 $T_\epsilon$ 传入候选生成器；
3. 记录硬约束可行性、几何 Pareto 和精确 ns delay；
4. 冻结候选后加载实际 TDL 协方差；
5. 计算 $J_{\rm CDD}$、CE、相关矩和误差地板；
6. 每场景最多选择 8 个几何/CE Pareto 代表进入 outage；
7. 计算共同随机数 10%/1% outage；
8. 生成 profile、RMS、$T_\epsilon$、所需间距、可行性和性能的潜力地图。

不得只展示通过场景。若全部硬厚 Sidon不可行，H2 不成立，但仍报告软约束候选用于解释。

### 4.5 Phase 4：E3 smoke、粗扫与加密

#### Smoke

每个 E3 场景在 15 dB 运行 20 paired trials，只验证：

- TDL profile/delay spread 实际展开值；
- QC/Sidon $\mathbf V$ 与 result-025 一致；
- matched 协方差随 profile 和 delay spread 更新；
- 两个 DMRS 等效平均；
- 共同随机数和 seed 重放；
- CE、LLR、LDPC 和输出字段有限。

Smoke 不用于性能结论。

#### 粗扫

全部 7 个 E3 场景使用 SNR `13:0.5:18.5` dB、每点 400 paired trials。若任一候选未跨越 10% BLER，先报告并按预定规则仅向相邻方向扩展最多 1.0 dB；不得改变步长或删除场景。

#### 10% 加密

对每个场景，根据粗扫中 QC/Sidon 各自跨越 10% 的相邻区间取并集，以 0.25 dB 步长生成加密网格，每点 3000 paired trials。目标 SNR 使用预定区间全部点的二项 logit 拟合，单候选区间使用 delta method；改善区间采用不利用正配对协方差的保守近似，并保存逐点四格计数和 McNemar 精确检验作为方向性佐证。

#### 1% 加密

若粗扫显示 1% 交叉且预计预算可使目标附近每候选累计错误块不少于 30，则按同样规则运行 0.25 dB、3000 paired trials 加密；否则只报告粗扫先导结果，不追加事后选择的 SNR 点。

### 4.6 正式停止条件

出现以下任一情况时暂停对应阶段：

- 展开配置与本 plan 的时延、归一化、DMRS 或知识条件不一致；
- 候选设计阶段意外读取 E2 实际 profile/协方差；
- matched 协方差或估计器出现无法解释的非有限值；
- 共同随机数或 seed 重放失败；
- outage/BLER 目标不在预定或允许扩展网格内；
- 运行中需要改变基线、候选集合、物理相位定义或验收阈值。

停止后先更新 plan 并说明已有输出是否仍可比较，不用后验修改覆盖原结果。

## 5. 输出、统计与复现

### 5.1 输出目录

正式输出写入：

```text
outputs/experiment026_cdd_design/<run_id>/
```

建议结构：

```text
validation/
e1_search/
e1_outage/
e2_thick_sidon/
e3_smoke/
e3_prescan/
e3_refine_10pct/
e3_refine_1pct/
final/
```

每个阶段保存 `resolved_experiment.json`、`commands.json`、`environment.json`、日志和代码版本/工作区变更标识。

### 5.2 必需数据

至少输出：

- `candidate_catalog.csv`：候选 ID、family、8 个 ns delay、辅助 $j_n$、residue、lift、规范代表和生成知识等级；
- `search_trace.csv`：起点、权重、迭代、访问次数、目标、接受动作和停止原因；
- `e1_metrics.csv`：14/16/18 dB 的 $J_{\rm CDD}$、CE NMSE、误差地板、rank/condition、$M_2^{\rm eff}$、$M_4^{\rm eff}$ 和几何距离；
- `e1_pareto.csv`；
- `e1_jcdd_ce_scatter_16db.png`；
- `e1_outage_curves.csv` 和 `e1_outage_targets.csv`；
- `e2_scenarios.csv`：profile、RMS、$\epsilon$、$T_\epsilon$ 和两类保护距离；
- `e2_geometry.csv`、`e2_metrics.csv`、`e2_outage_targets.csv`；
- `e2_potential_map.png`；
- `e3_bler.csv`、`e3_paired_error_counts.csv`、`e3_ce_nmse.csv`；
- `e3_target_snr.csv` 和 `e3_gain_vs_delay_spread.png`；
- `final_summary.json`。

CSV 中 ns 字段名称必须带 `_ns`；SNR、NMSE、改善和区间字段必须带 `_db`；概率使用线性值。

### 5.3 result-026 必须回答

`research/result-026.md` 与 `research/result-026-text.md` 必须成对生成，并逐项回答：

1. 025 状态、delay-matched 物理定义和链路公平性是否准确继承；
2. E1 实际搜索空间、预算、覆盖、去重和停止条件是什么；
3. E1 二维散点中 B0、B1、固定 residue、relaxed residue 和 unrestricted reference 的 Pareto 关系；
4. 是否存在通过 H1 的候选，其 8 个时延 ns、辅助 $j_n$、residue/lift 和全部关键指标是什么；
5. 非均匀 residue 相对最大分离 residue 是否带来收益，代价是否出现在 CE；
6. E2 各场景的 $T_\epsilon$、硬厚 Sidon可行性和几何边界；
7. 是否存在通过 H2 的场景和候选；负场景是否完整报告；
8. E1/E2 的互信息结果是否始终称为 outage，而非 BLER；
9. TDL-A 0.1/1/5 ns 下 QC/Sidon 的 10% 和可用时的 1% BLER目标、错误数、区间与增益趋势；
10. TDL-D/E 0.1/1 ns 下排序是支持、反向还是不可判定；
11. CE NMSE、相关代理、outage 与 BLER 排序不一致的所有情况；
12. 哪些候选值得进入下一轮正式 estimated-CSI BLER，哪些方向应停止。

无图版必须包含关键小表、精确数据路径和所有验收判定，不得依赖图片才能核验结论。

## 6. 预计复现入口

以下命令为计划接口，实际脚本完成后必须在结果中保存原样命令：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe -m unittest tests.test_cdd_design tests.test_cdd_metrics tests.test_rmmse_time_frequency
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_cdd_design.py --stage validate --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_cdd_design.py --stage e1-search --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_cdd_design.py --stage e1-outage --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_cdd_design.py --stage e2 --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage smoke --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage prescan --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage refine-10pct --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan026_near_flat_link.py --stage refine-1pct --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_plan026.py --run-id <run_id>
```

## 7. 执行顺序

1. 研究者确认本草案的搜索预算、H1/H2 阈值、E3 场景和 trial 预算；
2. 实现可复用 CDD 搜索与指标模块；
3. 完成单元测试、协方差/CE/MI 校准和确定性重放；
4. 运行 E1 搜索并生成 $J_{\rm CDD}$–CE 二维散点；
5. 冻结 E1 Pareto 候选，运行 E1 outage；
6. 运行 E2 厚 Sidon几何、实际协方差评价和 outage；
7. 审核 E1/E2 低成本结果，但不得据此改变 E3 固定场景；
8. 运行 E3 smoke、粗扫、10% 加密和条件满足时的 1% 加密；
9. 生成 final 汇总和两版 result-026；
10. 研究者确认 result 后再更新 `GOALS.md`、`KNOWLEDGE.md` 和 `research/README.md`，并决定下一轮 estimated-CSI BLER 候选。
