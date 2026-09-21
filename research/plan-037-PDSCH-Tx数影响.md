# plan-037：PDSCH 1Rx、60 km/h 下 Tx 数对 Sidon、B0QC 与 PRG6 性能的影响

> 状态：**草案，待研究者确认后执行**。研究者于 2026-09-20 指定本轮比较
> `4Tx/16Tx/32Tx`，信道和接收口径参考 plan-035，运行 estimated CSI、ideal CSI 与
> CE NMSE；本 plan 已在正式仿真前冻结第 4 节的时延和第 5 节的 PRG cycling 映射。
> 时延只经过解析推导、组合构造和矩阵审计，没有使用 outage、NMSE、BLER 或链路预扫结果筛选。

## 1. 目的与解释边界

本轮在 plan-035 的 `1Rx/1 layer`、TDL-A 100 ns、60 km/h PDSCH 链路上改变发射天线数，回答：

1. 在 4Tx、16Tx 和 32Tx 下，B0QC 与严格 Sidon 相对同 Tx 数 `PRG6 precoder cycling`
   的 estimated-CSI 和 ideal-CSI 10%/1% BLER 目标 SNR 增益是多少；
2. 上述增益怎样随 Tx 数变化，变化来自 ideal-CSI 分集还是同时伴随 CE NMSE 和
   ideal/estimated 间隔变化；
3. 16Tx 严格 Sidon 是否仍可在当前 $K=576$、comb-6 条件下构造并保持导频满秩；
4. 32Tx 严格 Sidon 是否可行；不可行时只报告不可行证明，不以 almost-Sidon、截断或重复
   delay 冒充严格 Sidon。

主比较定义为同一 Tx 场景内相对 PRG6 的目标 SNR 差。令

$$
G_{X\leftarrow\mathrm{PRG6}}(p)
=
\mathrm{SNR}_{\mathrm{PRG6}}(p)
-\mathrm{SNR}_{X}(p),
\qquad p\in\{10\%,1\%\},
$$

正值表示方案 $X$ 比 PRG6 节省 SNR。estimated CSI 是正式性能判据；ideal CSI 和 CE NMSE
是机理诊断。32Tx 没有严格 Sidon，因此不产生 32Tx Sidon 曲线或虚构增益值。

本轮回答的是“按各 Tx 数预先规定的构造规则得到的端到端性能”，不是保持同一个
$\mathbf V$ 不变时 $N_t$ 的单因素因果效应：4Tx、16Tx 的 Sidon 集不同，PRG 使用对应维数的
DFT 码本，B0QC 在 16Tx/32Tx 使用分组重复。跨 Tx 数的差值不声明 paired 方差缩减。

该任务对应 `GOALS.md` 第 4 节关于 Sidon 适用范围的开放问题，并扩展其发射天线数维度。
移动性和大于 8Tx 当前不在旧版 `GOALS.md` 第 6 节默认范围内；本轮以研究者的明确要求作为
定向扩展，不据此修改全局阶段验收状态。理论与口径遵循 `DESIGN.md` 第 3.7--3.9、5.4--5.5、
第 6 节和 `docs/design/CDD_DELAY_SELECTION_RULES.md`。

## 2. 冻结系统条件

| 参数 | 冻结取值 |
|---|---|
| 场景 | `A100_NT4_NR1_V60`、`A100_NT16_NR1_V60`、`A100_NT32_NR1_V60` |
| 天线与层数 | 4Tx/1Rx、16Tx/1Rx、32Tx/1Rx；均为 1 layer |
| 波形 | 48 PRB，576 个连续有效子载波，30 kHz SCS，FFT 4096，CP 288 samples，10 个 PDSCH symbols |
| 信道 | Sionna 1.0.2 TDL-A，RMS delay spread 100 ns，3.5 GHz，60 km/h，20 sinusoids；Tx--Rx 分支独立、同 PDP、无空间相关 |
| 时变模型 | OFDM-symbol 采样的时变频域乘法；定时与载波同步；不模拟 CFO、ICI 或 ISI |
| DMRS | symbols `[2,7]`，comb-6；每 symbol 96 pilot RE，共 192 个独立 LS 观测；5568 data RE |
| estimated-CSI 接收机 | 两个 DMRS 不平均；使用真实 60 km/h 时间协方差和各方案指定的频域协方差做二维时频 LMMSE；单根 Rx 独立估计，data RE 上做 1Rx MRC |
| ideal-CSI 接收机 | 每个 data RE 使用当前真实等效信道做 1Rx MRC；不使用 DMRS 估计值，不统计 CE NMSE |
| 调制编码 | 16QAM；NR 256QAM MCS table 的 MCS 8；目标码率 553/1024；Sionna LDPC 最多 8 次迭代 |
| 发射端知识 | 不使用当前或过时瞬时 CSI；三类方案均为固定预编码规则 |
| 随机性 | base seed `20260727`；absolute-trial 可重放；同一 Tx 场景内全部方案和 CSI 模式成对 |
| 执行设备 | CPU-only；保存 CPU placement、峰值 RSS、batch size 和墙钟耗时 |

本轮不运行 aged-CSI MRT、small-delay CDD、CDL、空间相关、多层、反馈量化或协方差失配。
不把 plan-035 的 8Tx 数据并入本轮主统计；可在 result 中作为明确标注的历史背景引用，但不得
与本轮三个场景冒充 paired comparison。

## 3. 功率、SNR、接收机与 NMSE

CDD 相位和单位总功率预编码为

$$
V_{k,n}=\exp\left(-j2\pi k\frac{j_n}{576}\right),
\qquad
\mathbf W=\frac{\mathbf V}{\sqrt{N_t}},
\qquad
\lVert\mathbf W_k\rVert_2^2=1,
$$

其中 $k=0,\ldots,575$ 是从首个 active subcarrier 起算的局部索引。一个 delay-grid 单位为

$$
\frac{1}{576\times30\ \mathrm{kHz}}
=57.870370370\ \mathrm{ns},
$$

FFT-sample 等效值为 $j_n\times4096/576=j_n\times64/9$。人工 delay 是数字循环相位，不加入
真实传播时延或 CP 预算。

每根 Rx 分支的噪声方差和横轴 SNR 沿用 plan-035：

$$
\sigma_n^2=10^{-\mathrm{SNR}_{\rm dB}/10}.
$$

estimated 和 ideal MRC、LLR 以及 data-RE CE NMSE 完全沿用 plan-035 第 3 节。每个 trial 的
NMSE 先在单根 Rx 的全部 data RE 上按线性能量比聚合，再跨 trial 在线性域求均值并转 dB；
保存 trial 级线性 NMSE、和、平方和、均值及 95% Monte Carlo 区间。LLR 不额外加入
CE-error-aware 项。

CDD 的 estimated 接收机使用本方案真实等效协方差
$\mathbf R_g=\mathbf R_{\rm phy}\odot(\mathbf W\mathbf W^H)$。PRG6 接收机不知道 DFT 波束索引，
只在每个 PRG 内使用 physical covariance。重复 B0QC delay 的 pilot matrix 对独立 Tx 分支
不满列秩，但本轮接收机直接估计单层复合等效信道，不分别恢复每根 Tx 分支；仍必须把 rank
退化作为诊断完整报告，不得把它写成满秩设计。

## 4. 时延的解析推导与正式冻结表

### 4.1 公共几何

当前 $K=576$、DMRS comb $S_f=6$，所以每个 DMRS symbol 的唯一导频数和折叠周期为

$$
N_p=K/S_f=96,
\qquad r_n=j_n\bmod96.
$$

严格 Sidon 要求全部 $N_t(N_t+1)/2$ 个无序二元和
$(j_a+j_c)\bmod576$ 互不相同。所有候选固定 $j_0=0$，并按公共循环移位和普通天线排列规范化。

### 4.2 4Tx 严格 Sidon：同时最优化导频可辨识性

4 个 residue 在长度 96 的圆周上的最小间距上界为 $\lfloor96/4\rfloor=24$。达到上界时 residue
只能等间隔分布为 `[0,24,48,72]`（允许公共旋转）。固定 residue 0 的 lift 为 0 后，穷举其余
三个 residue 的 $q\in\{0,1,2,3,4,5\}$，只保留模 576 严格 Sidon，并按规范 key 去重；这是一项
$6^3=216$ 状态的确定性组合枚举，不运行信道或链路仿真。得到 12 个规范可行代表，按字典序
冻结第一个：

$$
\mathbf j_{\rm S,4}=[0,24,72,240].
$$

其 residue 排序为 `[0,24,48,72]`，fold 最小圆周间距为 24；10 个无序二元和全异，pair-sum
最小圆周间距为 24；真实 comb-6 pilot matrix 的 rank 为 4、2-norm condition number 为 1。
因此它同时达到当前 comb-6 下 4Tx 导频折叠间距的理论上界。该结论只说明几何和纯平坦
pilot Gram 最优，不预告 TDL-A 下的 CE 或 BLER 排序。

### 4.3 4Tx B0QC：从原 8Tx 优化集抽取的临时代表

原 8Tx B0QC 为 `[0,9,18,27,36,45,54,63]`。本轮 4Tx 不在全部等差时延中做性能筛选，
临时取原集合的偶数位置，得到 `[0,18,36,54]`；奇数位置 `[9,27,45,63]` 与其只差公共移位，
故冻结规范代表 `[0,18,36,54]`。该选择保持等差结构并在原 8 个值内取每隔一个 delay，
不是“4Tx 全等差空间最优”的结论。后续若要筛选全部等差值，必须另立 plan 并在正式链路前冻结。

### 4.4 16Tx 严格 Sidon：几何构造，不做性能筛选

16Tx 使用固定 seed `20260924` 的 geometry-only 逐点构造：从 `[0]` 开始，每次只接受同时满足
模 576 新增无序二元和不冲突、comb-6 residue 与已有 residue 圆周距离不小于 4 的 index；
候选 index 按升序形成列表，再由固定 seed 选择，失败则重启。冻结首个达到 16 个元素的规范代表：

$$
\mathbf j_{\rm S,16}=
[0,9,56,61,71,174,243,280,303,307,388,407,419,427,451,509].
$$

其 136 个无序二元和模 576 全异，pair-sum 最小圆周间距为 1；residue 排序为
`[0,4,9,15,19,23,29,35,43,51,56,61,67,71,78,88]`，fold 最小圆周间距为 4；
pilot rank 为 16、condition number 为 1。16 个点的 fold 间距理论上界为 6；若达到 6，
所有 index 必须同余于某个模 6 类，136 个 pair sum 只能落入 96 个模 576 位置，因而不可能严格
Sidon，所以 fold 间距 6 已解析排除。本 plan 不声称当前 fold 间距 4 是全局最优，也不以
outage、CE 或 BLER 对它做后验替换。

### 4.5 32Tx 严格 Sidon：解析不可行

严格 Sidon 会使所有 $N_t(N_t-1)$ 个非零有向差 $j_a-j_b$（$a\ne b$）互不相同；否则
$j_a-j_b=j_c-j_d$ 会产生一个非平凡无序二元和相等。模 576 只有 575 个非零元素，所以必须满足

$$
N_t(N_t-1)\le575.
$$

32Tx 要求 $32\times31=992>575$，因此当前定义下不存在 32 元严格 Sidon 集。32Tx Sidon
正式标记为 `NOT_FOUND_PROVEN_INFEASIBLE`，不运行 BLER/NMSE；必要条件
$N_t(N_t+1)/2\le576$ 对 32Tx 虽成立，但不是充分条件，不能据此误判可行。

### 4.6 冻结汇总表

| Tx | 方案 | 冻结 delay grid $\mathbf j$ / 分组规则 | 物理 delay（ns） | pair sums | fold gap | pilot rank / cond | 状态与选择规则 |
|---:|---|---|---|---:|---:|---|---|
| 4 | `B0_QC` | `[0,18,36,54]` | `[0,1041.666667,2083.333333,3125]` | 7/10 unique | 18 | 4 / 1 | 原 8T B0 的隔点规范子集；临时代表 |
| 4 | `S0_SIDON` | `[0,24,72,240]` | `[0,1388.888889,4166.666667,13888.888889]` | 10/10 unique | 24 | 4 / 1 | 严格 Sidon；fold 间距达到理论上界 |
| 16 | `B0_QC_GROUPED` | `[0,9,18,27,36,45,54,63] × 2`；$j_{8g+m}=9m$ | 基本组 `[0,520.833333,1041.666667,1562.5,2083.333333,2604.166667,3125,3645.833333] × 2` | 15/136 unique | 0 | 8 / rank deficient | 两个 8 天线组重复原 8T B0 |
| 16 | `S0_SIDON` | `[0,9,56,61,71,174,243,280,303,307,388,407,419,427,451,509]` | `[0,520.833333,3240.740741,3530.092593,4108.796296,10069.444444,14062.5,16203.703704,17534.722222,17766.203704,22453.703704,23553.240741,24247.685185,24710.648148,26099.537037,29456.018519]` | 136/136 unique | 4 | 16 / 1 | 固定 seed 的几何构造；不做性能筛选 |
| 32 | `B0_QC_GROUPED` | `[0,9,18,27,36,45,54,63] × 4`；$j_{8g+m}=9m$ | 上述基本组 `× 4` | 15/528 unique | 0 | 8 / rank deficient | 四个 8 天线组重复原 8T B0 |
| 32 | `S0_SIDON` | — | — | — | — | — | `NOT_FOUND_PROVEN_INFEASIBLE`；不运行 |

表中的 condition number 只对满列秩 pilot matrix 给有限值；重复 B0 的 full-column condition
number 记为 `inf`，可另报 8 维非零子空间 condition number 1。正式 delay manifest 必须保存完整
展开数组，不能只保存 `×2/×4` 记号；同时保存 ns、FFT-sample 等效值、residue、lift、全部无序
二元和、pair/fold gap、rank、奇异值和 condition number。

## 5. PRG6 precoder cycling 冻结规则

48 PRB 按 6 RB 划分为 8 个 PRG。$N_t$ 维单位范数空间 DFT 码本第 $m$ 个向量为

$$
w_m[n]=\frac1{\sqrt{N_t}}\exp\left(-j2\pi\frac{nm}{N_t}\right).
$$

按研究者要求：码本波束数小于 PRG 数时继续轮询，波束数大于 PRG 数时固定选取 8 个并各扫一次。
本轮不把“任意”留给运行时随机数，而是冻结如下映射：

| Tx | DFT 码本大小 | 8 个 PRG 的 vector index | 规则 |
|---:|---:|---|---|
| 4 | 4 | `[0,1,2,3,0,1,2,3]` | 两轮顺序循环全部 DFT4 波束 |
| 16 | 16 | `[0,1,2,3,4,5,6,7]` | 固定使用前 8 个 DFT16 波束，各一次 |
| 32 | 32 | `[0,1,2,3,4,5,6,7]` | 固定使用前 8 个 DFT32 波束，各一次 |

PRG 从首个 active RB 起编号；不得跨 PRG 使用波束，也不得在 trial 间随机重抽。正式 manifest
保存 DFT 大小、相位符号、归一化、vector index 和每个 RB/subcarrier 到 PRG 的映射。

## 6. 场景内配对与公平性

4Tx、16Tx 场景各有 3 个发射方案 × 2 个接收模式，共 6 条 BLER 曲线和 3 条 estimated-CSI
NMSE 曲线；32Tx 有 2 个发射方案 × 2 个接收模式，共 4 条 BLER 曲线和 2 条 NMSE 曲线。
每个 `scenario + SNR + absolute trial` 内必须共享：

- 当前时变 TDL realization、payload、编码比特和原始 data AWGN；
- 全部 estimated 方案的原始 DMRS AWGN；
- 同一发射方案 ideal/estimated 的发射波形、当前信道、payload 和 data AWGN。

不同方案只改变预编码和与之匹配或透明的接收协方差。所有方案使用相同 DMRS、数据 RE、MCS、
译码器、单位总发射功率和 trial 区间。同一 Tx 场景内的 gain 使用 paired bootstrap；跨 Tx 数
的 gain 变化使用独立场景 bootstrap 组合，不把相同 seed 标签解释为同一维度的配对信道。

## 7. 实现范围、测试与 smoke

现有 plan-035 核心把 Tx 限定为 4/8，并把六个候选写死；本轮不得复制完整链路主循环。实现应：

1. 将 plan-033/035 的移动 TDL、二维 LMMSE、MRC、LDPC、absolute-trial seed、interval 合并和
   adaptive 逻辑抽成接受候选 manifest 的公共核心，保持 033/035 的已有 schema、seed 和输出兼容；
2. 新增 `tools/run_plan037_pdsch_tx_scaling.py`，只接受 `n_tx in {4,16,32}`、`n_rx=1`、
   `speed_kmh=60`、`receiver_modes=[estimated,ideal]` 和第 4--5 节候选；32Tx 配置必须拒绝 Sidon；
3. 新增确定性 delay 审计/复现入口，逐项验证第 4 节数组、严格 Sidon、32Tx 不可行上界、
   residue、pair sum、rank、condition 和单位换算；运行时禁止重新搜索或替换候选；
4. 扩展 Tx 维度校验和 Sionna TDL shape 检查；只在本链路需要的公共接口开放 16/32Tx，
   不顺带改变 PDCCH 等未验证入口的支持范围；
5. 新增 `tools/analyze_result037_pdsch_tx_scaling.py`，从本地 CSV/NPY 生成 crossing、gain、
   paired bootstrap、independent cross-Tx comparison、ideal/estimated gap 和 NMSE 汇总。

定向测试至少覆盖：

1. 第 4 节所有展开 delay 与第 5 节 DFT 映射完全一致，预编码逐 RE 功率为 1；
2. 4Tx Sidon 的最优 fold 上界、216 状态枚举和规范代表可重放；
3. 16Tx Sidon 有 136 个唯一模 576 无序二元和、16 个唯一 residue、rank 16、condition 1；
4. 32Tx 的有向差上界证明和配置拒绝路径；
5. 16/32Tx grouped B0 的重复次数、完整列数、rank 8、matched effective covariance 与
   直接从展开 $\mathbf W$ 计算的结果一致；
6. DFT4 循环两轮，DFT16/32 只用 index 0--7，PRG 不跨边界；
7. 1Rx estimated/ideal MRC、NMSE 归属、共享 seed、interval identity 和合并键正确；
8. plan-033/035 的现有定向测试与至少一份已有配置 `--stage validate` 保持通过。

三个场景先各运行 `[6,16] dB`、每点 20 个共同 paired trials 的 smoke。初始 batch size 固定为
4Tx `25`、16Tx `8`、32Tx `4`；若 smoke 出现内存不足或峰值 RSS 超出执行机可用内存，
只允许依次降低为 16Tx `4/2/1`、32Tx `2/1`，不改变 trial、seed 或物理定义，并把最终 batch、
RSS 和理由写入展开配置。任一场景在 batch 1 下仍不能完成，停止正式流程并报告资源不可行。

## 8. 预扫描、正式预算与停止条件（待确认）

### 8.1 预扫描

每个场景先使用公共网格 `[6,8,10,12,14,16] dB`，每点 400 个共同 paired trials。该范围依据
result-035 的 4/8Tx 相关三方案 crossing 约位于 8.2--11.3 dB，只用于确定 bracket，不作为
本轮性能先验。若任一正式存在的 estimated/ideal 曲线没有同时形成 10%和1%双侧 bracket，
只向缺失方向按 2 dB 扩展，硬边界为 `[-2,22] dB`；到边界仍缺失时暂停并报告，不外推。

### 8.2 正式网格和预算

预扫后为每个 Tx 场景分别冻结一套该场景全部曲线共用的 SNR 网格：10%和1% crossing 邻域
最大间隔 0.25 dB，过渡区最大 0.5 dB，只保留 bracket 所需范围及两侧最多一个保护点。
formal trial 1 前把三个网格、配置 SHA-256、batch、预算上界和预计运行量回填到第 12 节并由
研究者确认。

正式预算沿用 plan-035 的预声明自适应规则：

1. 每个冻结 SNR 点先运行 1,000 个共同 paired trials；以 1,000 trials 为固定检查和可恢复区间，
   单点上限 50,000 trials；
2. 只从真实相邻采样点识别每条 estimated/ideal 曲线的 10%和1%双侧 bracket；缺失时停止；
3. 10% bracket 端点的观测 BLER 位于 5%--20%时至少需要 200 errors；1%端点位于
   0.5%--2%时至少需要 200 errors；区间外保护点不追错误数；
4. 同一 SNR 只要任一正式曲线未满足端点要求，该场景全部曲线共同追加 1,000 trials；
5. 到 50,000 trials 仍不足时标记 capped，保留原始点，不补伪计数；
6. crossing 只在真实双侧 bracket 内按 log-BLER 线性插值；不做外推、PAVA、PCHIP、平滑或
   单调修正；
7. CE NMSE 与 estimated 曲线使用相同 trials，不单独追加 NMSE-only trials。

单点 BLER 报 Wilson 95%区间。同场景方案 gain 和 ideal/estimated gap 使用至少 1,000 次
absolute-trial paired bootstrap；跨 Tx 的 gain 差使用独立 bootstrap。若有效 bootstrap 重复少于
80%，只报告点估计和失败原因，不给伪区间。CE NMSE 报线性域均值的 95% Monte Carlo 区间。

## 9. 分析问题与判定口径

每个存在的 `Tx × scheme × receiver × target` 必须报告 crossing、真实 bracket、两端
errors/trials 和区间。主表分别给出 B0QC、Sidon 相对 PRG6 的 $G$：

- $G$ 的 95%区间完全大于 0：该 Tx/CSI/目标下支持正增益；
- 区间跨 0：排序不能确定；
- 区间完全小于 0：支持相对 PRG6 劣化；
- 缺 bracket、达到样本上限或 32Tx Sidon 不可行：分别标记，不外推、不插补。

同时报告 estimated minus ideal 的目标 SNR penalty、各目标附近的 CE NMSE 及
`NMSE_scheme - NMSE_PRG6`。NMSE 不设置“不得差于基线”的硬门槛，也不单独决定胜负；若出现
高 SNR CE 地板、非有限滤波器或零噪声误差地板异常，则该曲线停止并作为实现/可辨识性失败报告。

跨 Tx 趋势至少回答：4Tx 到 16Tx 的 B0/PRG、Sidon/PRG gain 是否改变；16Tx 到 32Tx 只比较
B0/PRG，因为 32Tx Sidon 不存在。不同 Tx 下候选结构不同，result 必须把观察到的变化表述为
“当前冻结构造下的联合变化”，不得归因于纯 $N_t$ 效应。

## 10. 输出、图和 result 要求

正式产物写入

`outputs/experiment037_pdsch_tx_scaling/<run_id>/<stage>/<scenario_id>/`。

至少保存原 YAML、展开配置和 SHA-256、代码版本/工作区变更标识、完整 candidate/receiver
manifest 和 hash、delay/PRG audit、interval CSV、estimated/ideal error flags、estimated CE
trial arrays、合并 BLER/NMSE CSV、seed/pairing audit、滤波器与功率诊断、峰值 RSS、batch 和日志。

每个 Tx 场景至少生成 estimated-CSI BLER、ideal-CSI BLER、estimated-CSI CE NMSE 三张曲线图；
另生成两张跨 Tx 汇总图或等价表：10%/1% 的 B0/PRG gain，以及 4Tx/16Tx 的 Sidon/PRG gain，
32Tx Sidon 位置明确显示 `N/A: strict Sidon infeasible`。BLER 图使用对数纵轴；零误块 CSV 保留
真实零值，绘图才可用 `0.5/trials` 下界。所有图由本地脚本读取已保存数据生成，Agent 不加载
仿真结果图片。

正式完成后成对生成：

- `research/result-037-PDSCH-Tx数影响.md`；
- `research/result-037-PDSCH-Tx数影响-text.md`。

两版 result 必须逐项回答第 1、9 节，完整列出第 4--5 节冻结构造及 32Tx 不可行证明，并区分
事实、机理推断和待验证解释。结果未经研究者确认前不更新 `KNOWLEDGE.md`/`GOALS.md`，不创建
Git checkpoint。

## 11. 执行顺序与停止条件

1. 研究者确认本 plan，尤其是 4Tx 临时 B0、16Tx Sidon、32Tx 不运行 Sidon 以及第 8 节预算；
2. 实现公共核心、037 入口、几何审计、分析入口和定向测试；
3. 运行新增测试、033/035 回归 validate 和三个 smoke；
4. smoke 全部通过后运行三个 prescan；
5. 回填并确认第 12 节正式网格、配置 hash、batch 和预算后开始 formal trial 1；
6. 完成初始预算及自适应端点追加，运行数据一致性、配对、功率、delay 和 DFT 映射审计；
7. 从本地原始数据生成图、成对 result 和完整复现命令，交研究者确认。

以下任一情况立即停止受影响场景，不进入或不继续正式比较：候选/DFT manifest 与第 4--5 节不符；
单位总功率失败；absolute-trial 不连续或曲线间失配；32Tx 被错误配置为严格 Sidon；滤波矩阵非有限；
batch 1 仍资源不可行；prescan 在边界内缺 bracket。已经完成的其他场景保留，但 result 必须标明
部分完成，不能把缺失曲线视为零增益。

## 12. 正式冻结区

以下项目必须在 formal trial 1 前回填并由研究者确认：

- 4Tx、16Tx、32Tx 各自正式 SNR 网格；
- 三份 formal YAML 路径和 SHA-256；
- smoke 后的最终 batch size、峰值 RSS 与预计总运行量；
- 初始共同 trial 数、端点追加上限和严格预算上界复核；
- 研究者确认日期。

在上述项目确认前，本 plan 保持草案状态，只允许实现、测试、smoke 和 prescan，不启动正式 trial。
