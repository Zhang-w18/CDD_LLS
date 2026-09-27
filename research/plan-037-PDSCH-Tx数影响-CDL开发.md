# plan-037：PDSCH Tx 数影响与 CDL 平台开发

> 状态：**草案，待研究者确认后执行**。研究者于 2026-09-20 指定本轮比较
> `4Tx/16Tx/32Tx`，信道和接收口径参考 plan-035，运行 estimated CSI、ideal CSI 与
> CE NMSE；本 plan 已在正式仿真前冻结第 4 节的时延和第 5 节的 PRG cycling 映射。
> 时延只经过解析推导、组合构造和矩阵审计，没有使用 outage、NMSE、BLER 或链路预扫结果筛选。
> 2026-09-24 重新规划第 13 节独立增补：为 32T2R、CDL-E 下的 8 波束降维与
> 波束域 CDD 开发可复用平台能力。本次交付只实现、测试和 smoke 平台能力，不运行
> BLER prescan/formal，也不产生 CDL 性能结论。该增补不改写前述 TDL 实验的冻结定义、
> 已完成结果或原有 BLER 流程；原 CDL BLER 草案保留在附录 A，仅供追溯且不得执行。

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

## 13. 增补：32T2R CDL-E 下的 8 波束降维与波束域 CDD 平台开发

### 13.1 目标、范围与完成边界

> 状态：**已完成（2026-09-24）**。本增补使用模式 C `fixed_cdl_statistics`，已交付可由后续
> 实验 plan 调用的码本生成、波束权值生成与缓存、波束/时延域功率可视化、波束域 CDD
> 预编码和接收机协方差构造能力。本次不运行 BLER prescan/formal，不以 BLER、CE NMSE、
> outage 或目标 SNR 选择波束，也不向 `result-037` 追加性能结论。

平台首先在 32T2R、CDL-E 100 ns 上完成开发和 smoke；正规 DFT 码本、PDP、协方差和 CDD
原语不得写死为 32T、CDL-E 或固定 DFT16。第 13.3.2 节的两级 SSB 方法按研究者要求固定为
第一级 8 波束、第二级 8 波束。需要回答的是实现正确性和可复用性：

1. 能否根据天面、TXRU 映射、极化和目标分支数，自适应生成正规 DFT 码本并选择每个
   时延分支的预编码；
2. 能否用 8 个 SSB 宽波束覆盖指定 CDL 角域，按 CDL 长期接收功率选择最强 SSB 波束，
   再在该波束对应角域内生成 8 个均匀覆盖的 DFT 窄波束；
3. 能否把原始 ray、正规 DFT 波束、宽波束、子波束、原始 PDP 和逐波束 PDP 放在统一
   坐标与功率口径下计算、保存和绘图；
4. 能否依据仓库根目录 `CDL模式下PDP计算方法.md` 构造信道估计所需的频域协方差；
5. 对相同 CDL 定义和天面配置，能否用严格 key 命中已保存的权值与统计产物，避免重复优化。

本次允许运行的数值工作只有确定性单元测试、解析核验和小规模 smoke。smoke 可生成少量 CDL
样本以核对解析量和形状，但不得编码、译码、累计误块或扫描 SNR。

### 13.2 基准场景、配置接口与局部符号

首个开发基准为：Sionna 1.0.2 CDL-E、RMS delay spread 100 ns、下行 32T2R；发射端为
$(M,N,P)=(8,8,2)$ 天面，每极化 $(M_p,N_p)=(2,8)$ TXRU，共 32 TXRU；接收端为
$(1,1,2)$ 双极化 2R。AE-to-TXRU 映射、阵元间距、方向图、orientation、载频、CDL 原始
profile 和所有角度变换必须进入 resolved config，不得靠场景名称推断。

平台配置至少包含：

| 字段 | 含义与默认值 |
|---|---|
| `beam_design.method` | `regular_dft`、`wide_beam_split` 或同时生成两者的 `both`；未指定时默认 `regular_dft` |
| `beam_design.num_branches` | 波束域 CDD 分支数 $K_b$；基准为 8 |
| `beam_design.ssb_grid` | 第一级 SSB 码本角域与网格；基准为 AoD $[-60^\circ,60^\circ]$、ZoD $[90^\circ,110^\circ]$、垂直 2 × 水平 4 |
| `beam_design.narrow_grid` | 第二级 DFT 窄波束网格；基准为在最强 SSB 波束角域内垂直 2 × 水平 4，共 8 个波束 |
| `beam_design.dft_oversampling` | 每个有效空间维的 DFT 过采样倍数；默认 1，可按维配置 |
| `beam_design.weight_constraint` | `unit_norm` 或 `constant_modulus`；默认 `unit_norm` |
| `beam_design.cache_dir` | 权值与统计 manifest 的仓库内输出目录 |
| `visualization.enabled` | 是否生成第 13.4 节全部数值表和本地图 |

令 $\mathbf H_s[k]\in\mathbb C^{N_r\times N_t}$ 为 TXRU 域信道，$\mathbf T$ 为
AE-to-TXRU 映射，$\mathbf w_m\in\mathbb C^{N_t}$ 为第 $m$ 个单位范数预编码，
$m=0,\ldots,K_b-1$。第 $m$ 个波束域等效信道为

$$
\widetilde{\mathbf h}_m[s,k]=\mathbf H_s[k]\mathbf w_m\in\mathbb C^{N_r}.
$$

所有波束设计只使用 CDL profile、天面和长期统计；不得读取待评估 trial 的瞬时信道。

### 13.3 时延分支预编码方法

平台提供两种独立方法。两者输出相同 schema：$\mathbf W_b=[\mathbf w_0,\ldots,
\mathbf w_{K_b-1}]$、每列的设计区域/中心、归一化信息、方向图数据、逐波束 ray 功率、
逐波束 PDP 和完整 manifest。

#### 13.3.1 正规 DFT 码本法（默认）

“正规 DFT 码本法”是确定各时延分支预编码的一种方法，不是 32T2R 场景的固定 16 波束
定义。若配置未指定 `beam_design.method`，必须使用本方法。

从 resolved 天面和 TXRU 映射识别每极化的有效规则网格
$N_v\times N_h$。对 DFT index $(q_v,q_h)$ 构造

$$
c_{q_v,q_h}[m,n]
=\frac{1}{\sqrt{N_vN_h}}
\exp\left[-j2\pi\left(
\frac{q_vm}{O_vN_v}+\frac{q_hn}{O_hN_h}
\right)\right],
$$

其中 $O_v,O_h$ 为对应维的过采样倍数。无该空间维时令相应长度和过采样倍数为 1。
多极化端口默认复用同一空间向量，并按有效极化数再次归一化；端口 flatten 顺序、相位符号
和极化相位必须由配置与单射线测试确认。若天面或 TXRU 映射不是规则可分离网格，必须明确
报错 `REGULAR_DFT_UNSUPPORTED_GEOMETRY`，不得静默套用 32T 特例。

码本大小由 $(N_v,N_h,O_v,O_h)$ 自适应得到，不固定为 16。候选波束的长期功率定义为

$$
P_b=\mathbf b_b^H\mathbf R_t\mathbf b_b,
\qquad
\mathbf R_t=\mathbb E[\mathbf H^H\mathbf H],
$$

其中 $\mathbf R_t$ 必须由给定 CDL 的解析 ray 统计计算，并用小规模独立 realization 平均做
smoke 交叉验证。分支数为 $K_b$ 时，默认按 $(-P_b,b)$ 稳定排序选择前 $K_b$ 个码字；
也允许配置显式 beam index。选择只由长期统计或显式 index 决定，不使用 BLER/NMSE。
若候选数小于 $K_b$，配置校验失败，不循环复用波束。

#### 13.3.2 宽波束法

宽波束法采用固定的两级角域码本，不再根据 CDL ray 能量自适应搜索最小覆盖矩形。第一级
SSB 码本覆盖水平 AoD $[-60^\circ,60^\circ]$（等价方位表示为 $[300^\circ,360^\circ)\cup
[0^\circ,60^\circ]$）和垂直 ZoD $[90^\circ,110^\circ]$。角域按垂直 2 × 水平 4 均匀划分，
共生成 8 个 SSB 宽波束。每个 SSB 波束以其中心指向为中心，水平覆盖 $30^\circ$、垂直覆盖
$10^\circ$；边界采用闭区间，波束顺序为垂直优先、每个垂直行内水平角递增。

每个 SSB 宽波束在对应角域内以目标响应 1、角域外目标响应 0 进行确定性加权最小二乘合成，
施加 `weight_constraint` 后单位范数化。使用完整 CDL profile 的解析长期发射协方差
$\mathbf R_t$ 计算各候选的长期接收功率

$$
P_{{\rm SSB},b}=\mathbf w_{{\rm SSB},b}^H\mathbf R_t\mathbf w_{{\rm SSB},b},
$$

并以 $(-P_{{\rm SSB},b},b)$ 稳定排序选择接收功率最强的 SSB 波束。选择只使用长期统计，
不得读取待评估 trial 的瞬时信道，也不得用 BLER、CE NMSE 或后验 SNR 选择波束。

第二级只在最强 SSB 波束的 $30^\circ\times10^\circ$ 角域内生成 8 个 DFT 窄波束。该角域再次按
垂直 2 × 水平 4 均匀划分，每个窄波束指向对应子区域的中心；因此每个子区域宽
$7.5^\circ$、高 $5^\circ$。窄波束权值取该中心 AoD/ZoD 的 TXRU 阵列响应共轭并单位范数化，
即为离散空间频率上的 DFT/steering 波束。8 个窄波束必须完整、无缝覆盖最强 SSB 角域，
顺序同样为垂直优先、水平角递增；CDD 时延按该稳定顺序映射，不按窄波束后验功率重新排序。
配置、8 个 SSB 区域、最强波束 index 和功率、8 个窄波束区域及中心指向均写入 manifest。

#### 13.3.3 方向图与权值缓存

全部 8 个 SSB 宽波束、选中的最强 SSB 波束和其内部 8 个 DFT 窄波束必须画在方向图中：
至少包含归一化 dB 极坐标水平/垂直切面和 AoD--ZoD 二维图；最强 SSB 波束使用粗实线，
窄波束使用可区分颜色/线型，并标出两级区域边界、主瓣及同级栅瓣。水平切面必须注明固定 ZoD，垂直切面必须注明
固定 AoD；cut 的扫描角度、固定角度和精确方向图数组必须随绘图输入落盘。图必须由本地脚本
从已保存权值和网格数据生成。

平台按内容寻址缓存权值。cache key 至少由以下规范化输入的 SHA-256 组成：展开后的 CDL
profile/ray 表及版本、delay spread、载频、链路方向、天面几何、orientation、阵元方向图、
极化、AE-to-TXRU 映射、设计方法、$K_b$、两级角域与网格、过采样、权值约束和算法版本。
缓存目录保存 `weights.npz`、`manifest.json`、逐 ray/逐波束 CSV、方向图网格数据和 SHA-256。
命中时逐项校验 schema、shape、单位范数和 hash 后直接加载；任一字段或 hash 不同即重算，
不得只按“CDL-E + 32T”文件名复用。缓存只对应长期 CDL 定义；若未来以某个瞬时 realization
设计权值，realization seed/hash 必须额外进入 key，且不能与本 plan 的统计设计混用。

### 13.4 波束域/时延域功率可视化能力

新增统一分析入口，可对任意支持的 CDL、天面和已缓存波束集生成以下数据与图，不依赖链路
BLER 流程：

1. **正规 DFT 码本波束域功率分布：**输出全部自适应 DFT 码字的解析长期功率、归一化 dB、
   index 与空间频率；选中的 $K_b$ 个分支明确标记。
2. **两级波束功率分布：**输出 8 个 SSB 宽波束的解析长期接收功率并标出最强波束；同时
   输出最强 SSB 波束及其内部 8 个 DFT 窄波束的绝对捕获功率、占协方差总功率比例和方向图指标。
3. **原始 ray 角度—功率分布：**同一份原始 ray 表分别生成极坐标图和二维坐标图；二维图
   默认 AoD--ZoD，并可切换到 $(u_h,u_v)$。LoS specular 与 diffuse ray 使用不同标记，
   颜色或点面积表示线性功率；不得只画 cluster 中心。
4. **原始 PDP：**按公共物理时延合并 ray 功率后绘制，不含发射波束增益。
5. **逐波束 PDP：**正规 DFT 码本与最强 SSB/窄波束码本分图绘制，每个波束占一个子图；同一码本
   的子图共享以标准 CDL 最早物理路径为零点的时延轴，只画离散抽头，不连接相邻抽头。
   在选中正规 DFT、最强 SSB 波束和窄波束的全部绘制波束中，按 PDP 总绝对功率选出最强波束，
   并以该波束的总功率作为两张逐波束图唯一的公共归一化分母；不得对每个波束分别归一化。两张图
   的所有子图使用同一个纵轴范围，以便直接比较波束间时延域功率。数据表同时输出绝对功率
   和各波束自身归一化形状。

令原始 ray 的公共物理时延为 $\tau_\ell$，包含阵元方向图、极化、XPR、接收响应和
TXRU 映射的路径矩阵为 $\mathbf A_\ell$。第 $m$ 个波束在该路径上的长期功率为

$$
q_{m,\ell}
=\mathbb E\left[\|\mathbf A_\ell\mathbf w_m\|_2^2\right].
$$

原始 PDP 与第 $m$ 个波束的 PDP 分别为

$$
S_{\rm CDL}(\tau)=\sum_\ell p_\ell\delta(\tau-\tau_\ell),
\qquad
S_m^{\rm abs}(\tau)=\sum_\ell q_{m,\ell}\delta(\tau-\tau_\ell).
$$

因此不同波束 PDP 与原始 PDP **共享同一组物理时延取值** $\{\tau_\ell\}$；波束只改变相应
抽头功率 $p_\ell\rightarrow q_{m,\ell}$，不产生新的物理时延。只有后续叠加 CDD 人工时延
$d_m$ 时，该分支整体变为 $S_m^{\rm abs}(\tau-d_m)$。实现中必须先在公共时延轴聚合，禁止
逐波束把首径重新置零。所有图同时保存其绘图输入 CSV/NPZ、单位、归一化口径和复现命令。

### 13.5 原始 PDP、逐波束 PDP 与参考 PDP

原始 CDL PDP 从本轮实际 Sionna profile 读取 cluster delay、average power、20-ray 展开、
LoS specular 分量和 K-factor 规则，并按 RMS delay spread 换算为秒。必须保存标准原表、
Sionna 实例数组、最终展开 ray 表和相同 delay 合并表，检查 LoS 首径不重复计数、线性功率
归一化和单位一致。

对每个输出权值保存 $q_{m,\ell}$、未归一化总功率
$P_m=\sum_\ell q_{m,\ell}$ 及归一化形状

$$
S_m(\tau)=\frac{1}{P_m}
\sum_\ell q_{m,\ell}\delta(\tau-\tau_\ell).
$$

解析 ray-sum 是正式数值；独立随机初相平均只用于 smoke。必须验证解析总功率等于
$\mathbf w_m^H\mathbf R_t\mathbf w_m$。平台允许显式指定参考波束；未指定时，正规 DFT 法取
所选分支中长期功率最大者，宽波束法默认取第一级选中的最强 SSB 波束。参考选择规则和 hash 必须
进入 manifest，后续协方差构造不得按 SNR、trial、BLER 或 NMSE 改选。

### 13.6 坐标系、阵列和 profile 审计

任何权值缓存或可视化产物落盘前必须自动通过：

1. 明确 AoD/AoA 的零点、正方向、wrap 和 ZoD/ZoA 的天顶角语义；
2. 在 boresight、水平/垂直正负空间频率上做单射线校准，核对方向图与
   $\mathbf H\mathbf w$ 的最强响应；
3. 核对 AE/TXRU/极化 flatten 顺序、共轭、DFT 相位符号和下行 departure/arrival 方向；
4. 核对 $\mathbf T$ 的 shape、列范数、端口支持集和双极化隔离；
5. 用直接 AE 求和与可分离阵列因子两条路径核对主瓣、零点和栅瓣；不可分离天面只运行
   直接 AE 路径并明确记录；
6. 对 CDL 原表和展开 ray 数组检查 SHA-256、shape、单位、功率、K-factor、ray coupling、
   角度变换和 profile 版本；
7. 检查 $\mathbf R_t$ Hermitian、半正定、trace，以及 ray-sum、协方差二次型和小规模
   Monte Carlo 三者一致；
8. 检查全部输出权值有限、单位范数、缓存 round-trip 完全一致。

方向图存在栅瓣时必须全部显示，不得只保留与 CDL 强 ray 对齐的一个方向。

### 13.7 32T 到 $K_b$ 个等效分支及波束域 CDD

对任一方法输出的 $\mathbf W_b=[\mathbf w_0,\ldots,\mathbf w_{K_b-1}]$，第 $m$ 个分支为
$\widetilde{\mathbf h}_m[k]=\mathbf H[k]\mathbf w_m$。给定人工 delay index
$\mathbf j=(j_0,\ldots,j_{K_b-1})$ 和活动子载波数 $K_{\rm sc}$，最终预编码为

$$
\mathbf w_{\mathbf j}[k]
=\alpha[k]\sum_{m=0}^{K_b-1}
e^{-j2\pi k j_m/K_{\rm sc}}\mathbf w_m,
$$

其中

$$
\alpha[k]
=\left\|\sum_m e^{-j2\pi k j_m/K_{\rm sc}}\mathbf w_m\right\|_2^{-1}
$$

保证每个活动子载波单位总发射功率。只有当 $\mathbf W_b^H\mathbf W_b=\mathbf I$ 时，才可
把 $\alpha[k]$ 简化为 $1/\sqrt{K_b}$；宽波束分裂得到的子波束一般不预设正交，必须保留并
记录逐子载波归一化。接收端等效信道为

$$
\mathbf g_{\mathbf j}[k]
=\mathbf H[k]\mathbf w_{\mathbf j}[k]
=\alpha[k]\sum_m e^{-j2\pi k j_m/K_{\rm sc}}
\widetilde{\mathbf h}_m[k].
$$

实现必须用直接 $\mathbf H\mathbf w_{\mathbf j}$ 与先投影到波束域再叠加 CDD 两条路径逐元素
核对。delay 数组由调用方 plan 提供；本次只验证通用 B0/Sidon 小例子和相位分母，不选择
“最优” delay，也不运行 BLER。

### 13.8 信道估计协方差矩阵

实现参考仓库根目录 `CDL模式下PDP计算方法.md`，在公共接口中提供三种层级，并把所选方法
写入 receiver manifest：

1. `common_reference_pdp`：全部 CDD 分支复制同一参考 PDP；
2. `beam_specific_pdp_independent`：每个分支使用自身波束 PDP，忽略波束间相关项；
3. `beam_joint_covariance`：保留同一物理路径上不同波束的联合协方差，直接构造理想频域
   协方差；该方法不能伪称为普通非负标量 PDP。

对于方法二，若各分支等功率且人工时延为 $d_m$，频域协方差为

$$
R_f[k,k']
=\frac{1}{K_b}\sum_{m=0}^{K_b-1}\sum_\ell q_{m,\ell}
e^{-j2\pi(f_k-f_{k'})(\tau_\ell+d_m)}.
$$

对于方法三，令 $\mathbf C_\ell=\mathbb E[\mathbf g_\ell\mathbf g_\ell^H]$，并把实际
逐子载波归一化 $\alpha[k]$ 包含在 CDD 系数向量 $\mathbf c[k]$ 中，则

$$
R_f[k,k']
=\sum_\ell e^{-j2\pi(f_k-f_{k'})\tau_\ell}
\mathbf c^T[k]\mathbf C_\ell\mathbf c^*[k'].
$$

三种方法必须共享同一下行物理时延参考；不得逐波束重新置零。实现应接受任意子载波坐标，
输出 Hermitian 协方差及输入 PDP/联合统计 hash，并验证 Hermitian、半正定、对角功率、
从移位 PDP 直接求和与矩阵公式一致。时间协方差和二维 LMMSE 本轮只保留兼容接口与合成
小例子测试，不运行 CE 性能扫描。

### 13.9 实现位置、测试与 smoke

不得复制 CDL/PDSCH 主循环。可复用算法放入 `cdd_lls/`，专题配置和只读分析入口放入
`configs/` 与 `tools/`。实施前先检查现有模式 C、阵列构造、profile 提取和绘图代码，优先
扩展公共 schema。至少交付：

1. 天面自适应正规 DFT 码本生成与长期功率选择；
2. 固定 2×4 SSB 宽波束合成、长期接收功率最强波束选择，以及其角域内固定 2×4 DFT 窄波束生成；
3. 内容寻址的权值/统计缓存和严格失效规则；
4. 第 13.4 节统一波束域/时延域分析与本地绘图入口；
5. 公共 PDP、逐波束 PDP、参考 PDP 和三类频域协方差构造器；
6. 通用波束域 CDD 合成器及正交/非正交情况下的功率归一化；
7. 32T2R CDL-E 基准 YAML，以及不依赖 32T/CDL-E 的合成测试配置。

定向测试至少覆盖：

- 多种 $N_v\times N_h\times P$ 天面下 DFT 码本 shape、范数、相位、极化复制和过采样；
- 不规则天面拒绝路径，不得退化为固定 DFT16；
- SSB 总角域、2×4 区域边界和稳定顺序；最强 SSB 长期功率选择；其 $30^\circ\times10^\circ$
  角域内 2×4 窄波束的边界、中心指向、单位范数和主瓣指向；
- SSB/窄波束生成可重放，方向图输入数据、区域边界、中心指向和最强 SSB 标记齐全；
- cache 命中、任一物理/算法字段变化后的 miss、hash/shape 损坏拒绝；
- 原始/逐波束 PDP 共享 delay、只改变功率，LoS/K-factor 无重复计数；
- 三类协方差的手算小例子、Hermitian/半正定/对角功率和时延参考；
- 正交 DFT 分支使用 $1/\sqrt{K_b}$，非正交宽波束分支使用 $\alpha[k]$，全部子载波功率为 1；
- 直接 TXRU 合成与波束域 CDD 合成一致；
- 现有模式 C 相关测试和配置 validate 不回归。

smoke 仅生成少量 CDL-E realization，核对解析波束功率与样本均值的数量级、执行全部图表入口、
保存并重新加载 cache、构造三类协方差并完成一次合成 LS/LMMSE shape 测试。不得生成 BLER、
crossing、gain 或 formal trial 输出。

### 13.10 产物、验收与本轮停止条件

平台产物写入 `outputs/experiment037_cdl_platform/<run_id>/`，至少保存 resolved config、版本和
SHA-256、展开 CDL ray 表、阵列/TXRU manifest、DFT/宽波束/子波束权值、cache manifest、
逐 ray/逐波束功率、原始/逐波束 PDP、三类协方差审计和所有绘图输入。图至少包括：

- 全部 8 个 SSB 宽波束，以及最强 SSB 与其内部 8 个 DFT 窄波束的同图方向图；
- 正规 DFT 与宽/子波束的波束域功率分布；
- 原始 ray 角度—功率极坐标图和二维图；
- 原始 PDP 与全部选中波束 PDP 对照图。

本轮完成条件为：新增测试通过；现有相关回归通过；32T2R CDL-E validate 和无 BLER smoke
通过；cache round-trip 通过；上述数值表、图和复现命令齐全。不得因为平台 smoke 通过而声称
CDD 获得性能增益，也不得更新 `KNOWLEDGE.md` 的性能结论。

完成记录（2026-09-24）：CDL-E 与附加 CDL-C 场景的 validate、无 BLER smoke、全部图表入口、
三类频域协方差审计、合成 LS/LMMSE shape smoke 和 cache 二次命中均通过；定向测试及相关
`fixed_cdl_statistics`/Plan-037 回归共 22 项通过。复现配置为
`configs/plan037_cdl_c_e_beam_platform.yaml`，入口为 `tools/run_plan037_cdl_platform.py`，产物位于
`outputs/experiment037_cdl_platform/beam_patterns_e_c/`。本状态只表示平台开发验收完成，不表示
任何 BLER、CE NMSE 或 CDD 性能增益已经完成正式实验验证。

以下情况停止受影响能力并保留诊断：profile/ray 参数不能无损展开；坐标、共轭、极化或
AE-to-TXRU 审计失败；两级角域不能按配置完整无缝划分；最强 SSB 选择不能由解析长期功率唯一
重放；PDP delay 轴不一致；协方差非有限、
非 Hermitian 或显著非半正定；cache key 不完整；逐子载波单位总功率失败。研究者后续若要运行
BLER，必须另行确认场景、候选 delay、接收机方法、SNR 网格、trial 预算和正式判据。

## 附录 A：已废止的 CDL BLER 草案（仅供追溯，不执行）

> 本附录是 2026-09-24 重新规划前的历史文本，与当前第 13 节冲突时以第 13 节为准。
> 其中 top-8 固定 DFT16、BLER/CE prescan、formal 预算、Gate C 和 result 追加均已取消，
> 不属于本轮实现或验收范围。

### A.1 目的、状态与解释边界

> 历史状态（已废止）：原草案曾计划在研究者确认后实现和执行。本节使用
> `SIMULATION_PLATFORM_IMPLEMENTATION.md` 第 2.2 节的模式 C，同时运行 ideal CSI、
> estimated CSI 和 data-RE CE NMSE。现有模式 C 只有 ideal-CSI 链路；estimated-CSI
> 接收机、参考 PDP、二维时频 LMMSE 和 CE 落盘均属于本增补必须实现并通过 Gate C 的内容。

本增补回答三个问题：

1. 对固定 32 TXRU 天面和标准 CDL-E profile，根据独立长期统计样本得到的
   16 个 DFT 波束 RSRP 排序，是否与从 CDL-E 固定 cluster/ray 参数解析计算的
   波束化角度功率排序一致；
2. 若不一致，能否通过 AoD/ZoD 坐标、阵列索引、相位正负号、双极化端口顺序、
   AE-to-TXRU 映射和阵元方向图的定向审计解释；排除后仍不一致时，将其作为
   CDL profile 参数或实现异常，不继续正式 BLER；
3. 固结长期 RSRP 最大的 8 个波束后，在相同 8 维波束子空间内，原 8Tx
   `B0_QC` 和 `S0_SIDON` 时延是否比使用同一组波束的 `precoder cycling`
   获得 ideal-CSI 和 estimated-CSI 10%/1% BLER 目标 SNR 改善，以及改善是否伴随
   CE NMSE 或 ideal/estimated gap 的变化；
4. 标准 CDL-E 原始 PDP、选中8个波束各自的波束化 PDP 和最终接收机参考 PDP
   有何差异；分别使用8条真实波束 PDP、以及复制最强波束参考 PDP 的两类
   CDD-aware 接收机，其等效 PDP、CE NMSE 和 estimated-CSI BLER 有何差异。

本实验是 `GOALS.md` 第 6 节中“massive-MIMO 预编码降维后再叠加 CDD”的
研究者定向扩展，也直接检查 `KNOWLEDGE.md` K19 留下的几何 LoS/CDL 开放问题。
它不测试 32 元严格 Sidon；时延只作用于 8 个降维后的等效波束分支，
因此不与第 4.5 节的“32Tx 严格 Sidon 不可行”矛盾。

### A.2 冻结场景和局部符号

| 参数 | 冻结取值 |
|---|---|
| 场景 ID | `E100_NT32_NR2_V60_BEAM8` |
| 信道 | 模式 C `fixed_cdl_statistics`，Sionna 1.0.2 CDL-E，RMS delay spread 100 ns，下行 |
| 频率与速度 | 4 GHz，60 km/h；每个 realization 内按固定 ray Doppler 连续演化 |
| 发射天面 | $(M,N,P)=(8,8,2)$ AE；每极化 $(M_p,N_p)=(2,8)$ TXRU，共 32 TXRU |
| AE-to-TXRU | 每个 TXRU 用等相、单位范数权重驱动同一水平列的连续 4 个垂直 AE |
| BS 间距/方向图 | 水平 $0.5\lambda$，垂直 $0.8\lambda$，3GPP TR 38.901 阵元方向图，cross polarization |
| UE 阵列 | $(1,1,2)$，即单个空间位置的双极化 2R，omni，cross polarization |
| 角度变换 | 禁用 `mean_aod_deg` 平移及 AoD/AoA/ZoD/ZoA scale；使用 profile-native 固定角度和 ray offsets |
| 阵列姿态/拓扑 | 使用模式 C 当前固定的 BS/UT 位置、阵列 orientation 和速度方向；必须全量写入 manifest |
| 资源与链路 | 48 PRB，576 active SC，30 kHz SCS，FFT 4096，CP 288 samples，10 PDSCH symbols |
| DMRS/MCS | symbols `[2,7]`，comb-6；16QAM，NR 256QAM table MCS 8，目标码率 553/1024 |
| 接收机 | ideal CSI 与 estimated CSI，均为 2Rx MRC；cycling 使用透明逐 PRG 接收机，CDD 同时使用 `cdd_beamwise_pdp` 与 `cdd_reference_pdp_replicated` 两类全带接收机；estimated 对每根 Rx 独立做二维时频 LMMSE并报告 data-RE CE NMSE |
| 发射知识 | 只使用独立统计阶段固结的 8 波束集；不使用当前或过时瞬时 CSI |

本节的局部符号为：$\mathbf H_{d,s}[k]\in\mathbb C^{2\times32}$ 是 TXRU 域信道；
$\mathbf T\in\mathbb C^{128\times32}$ 是双极化 AE-to-TXRU 映射；
$\mathbf B=[\mathbf b_0,\ldots,\mathbf b_{15}]\in\mathbb C^{32\times16}$ 是本轮 DFT 码本；
$\mathcal I_8$ 是选中的 8 个波束 index 集；
$\mathbf B_8\in\mathbb C^{32\times8}$ 是按 index 升序排列的选中列。

### A.3 16 波束 DFT 码本的冻结定义

每极化的 TXRU 网格是 $2\times8$。对 $q_v\in\{0,1\}$、
$q_h\in\{0,\ldots,7\}$ 定义

$$
c_{q_v,q_h}[m,n]
=\frac{1}{4}
\exp\left[-j2\pi\left(\frac{q_vm}{2}+\frac{q_hn}{8}\right)\right],
\quad m\in\{0,1\},\ n\in\{0,\ldots,7\}.
$$

两个极化使用相同空间波束和相同相位：

$$
\mathbf b_{q_v,q_h}
=\frac{1}{\sqrt2}
\begin{bmatrix}
\mathbf c_{q_v,q_h}\\
\mathbf c_{q_v,q_h}
\end{bmatrix},
\qquad
b=q_v8+q_h.
$$

因此 $\mathbf B^H\mathbf B=\mathbf I_{16}$、$\|\mathbf b_b\|_2=1$。这是“水平 DFT8
$\times$ 垂直 DFT2，双极化共用空间波束”的明确定义。它不是当前模式 C 的
“16 个水平二倍过采样、垂直同相 secondary beams”；实现时必须新增独立
codebook type，不得复用现有 `secondary_horizontal_beams: 16` 的含义。

相位符号、$(m,n,极化)$ 的 flatten 顺序以 $\mathbf H\mathbf b$ 的实际发射语义为准，
并由第 13.6 节的单射线测试校准。由于每个垂直 TXRU 覆盖 4 个 $0.8\lambda$
间距 AE，垂直 TXRU 中心间距为 $3.2\lambda$，可能出现栅瓣；不得仅用
$\arcsin$ 把每个 DFT index 声明为唯一物理方向。正式 manifest 必须同时保存
TXRU 权重、$\mathbf T\mathbf b_b$ 的 AE 权重、主瓣与所有同级栅瓣方向。

### A.4 长期 RSRP 选择和解析角度功率交叉验证

统计阶段使用与正式链路不同的 `statistics_seed`。对第 $d$ 个独立
realization，先在该 realization 内的 $S$ 个 OFDM-symbol 时刻和 $K=576$ 个子载波上求

$$
\widehat{\mathbf R}_{t,d}
=\frac{1}{SK}
\sum_{s=1}^{S}\sum_{k=1}^{K}
\mathbf H_{d,s}^{H}[k]\mathbf H_{d,s}[k],
\qquad
P_{d,b}=\mathbf b_b^H\widehat{\mathbf R}_{t,d}\mathbf b_b.
$$

再以 $\overline P_b=D^{-1}\sum_dP_{d,b}$ 作为第 $b$ 个波束的长期 RSRP。
频域和时域 sample 因相关而不当作 $DSK$ 个独立样本；所有区间、bootstrap 和
选择稳定性都以独立 realization 向量 $(P_{d,0},\ldots,P_{d,15})$ 为重采样单位。

统计预算冻结为：初始 $D=2000$，每批追加 500，最大 $D=10000$；使用
2000 次 paired-realization bootstrap。只有同时满足以下条件才冻结 $\mathcal I_8$：

1. 按 $( -\overline P_b, b )$ 排序得到的 top-8 集合在连续 3 个累积检查点不变；
2. bootstrap 中每个选中波束进入 top-8 的概率不小于 0.99，每个未选波束的概率
   不大于 0.01；
3. 第 8 名与第 9 名之间 $\overline P_{(8)}-\overline P_{(9)}$ 的 paired-bootstrap
   99% 区间下界大于 0。

到 $D=10000$ 仍不满足时，标记 `TOP8_NOT_IDENTIFIABLE`，不用 beam index 人为打破
物理并列，不进入正式 BLER。最终选中集、排序、全部 16 个 RSRP、区间、
选中概率、统计样本数、seed 和 manifest SHA-256 必须在链路 prescan 前冻结。

解析参考不使用链路 Monte Carlo 结果。从 Sionna 1.0.2 的 CDL-E profile 读取并保存
cluster power/delay/AoD/AoA/ZoD/ZoA、20-ray offsets、XPR、ray coupling、LoS 分量和 K-factor，
使用实际阵元方向图、极化响应和 $\mathbf T$ 计算标准统计下的
$\mathbf R_t^{\rm ana}=\mathbb E[\mathbf H^H\mathbf H]$，并得到

$$
P_b^{\rm ana}=\mathbf b_b^H\mathbf R_t^{\rm ana}\mathbf b_b.
$$

等价地，可将 $P_b^{\rm ana}$ 实现为对 CDL-E 所有 cluster/ray 的功率加权和，
其中每一项使用完整发射波束 $\mathbf T\mathbf b_b$ 与 2R 接收响应，而不是只使用
cluster 中心角。该公式是主验证量。

本 plan 中的“角度功率”统一指：含阵元方向图、极化、XPR、子阵映射和
2R 接收响应的二维角度功率分布，经完整发射波束方向图加权积分后得到的
$P_b^{\rm ana}$。波束选择、排序和验收不使用波束中心角的点值。图上同时给出：

1. AoD--ZoD 二维角度功率图，标出 CDL-E 的离散 ray 与 16 个波束主瓣/栅瓣；
2. 16 个波束的 $P_b^{\rm ana}$ 与 $\overline P_b$ 并列图，Monte Carlo 点带 99% 区间；
3. top-8 选中概率及其随 $D$ 的稳定性图。

解析与 Monte Carlo 比较前都除以各自最强波束功率，去掉公共尺度。验收要求为：

- top-8 集合一致；若解析的第 8/9 名差小于 0.1 dB，则按该 0.1 dB 并列等价类比较；
- 16 点的 Spearman rank correlation 不小于 0.99；
- 每个波束的归一化功率差绝对值不超过 0.15 dB，或解析值落入该点的
  simultaneous 99% Monte Carlo 区间，两者满足其一即可。

不满足时先执行第 13.6 节坐标审计，不允许通过修改角度平移、缩放、波束集
或门槛追认一致。坐标审计通过但仍不一致时，标记 `CDL_ANGLE_POWER_MISMATCH`，
保存证据并停止 BLER。

### A.5 原始 PDP、波束化 PDP 与接收机参考 PDP

本节区分三个不能混用的对象。

1. **标准原始 PDP。**从本轮实际使用的 Sionna 1.0.2 CDL-E profile 读取标准化
   cluster delay、cluster average power 和 K-factor，并乘以 100 ns RMS delay spread。
   对 LoS profile，必须按 Sionna 实例化后的 K-factor 规则把 specular 分量与 diffuse
   分量拆开，再把数值相同的 delay 合并；不得把 LoS 首径功率重复计数。记最终离散测度为

   $$
   S_{m CDL}(\tau)=\sum_c p_c^{\rm CDL} \delta(\tau-\tau_c) ,
   \qquad  \sum_c p_c^{\rm CDL}=1.
   $$

   原始数值来自标准 CDL-E profile；实现仍须保存“标准原表、Sionna 原始数组、K-factor
   展开后的最终 delay/power 表”及三者的核对结果，不能只凭 profile 名称重写一份常数表。
   原始 PDP 不含阵列方向图、极化、TXRU 映射或发射波束增益。

2. **波束化 PDP。**把第 13.4 节解析角度功率的 ray-sum 保留到 delay 维。令 $\ell$
   表示 CDL-E 的 diffuse ray 或 LoS specular 分量，$\tau_\ell$ 为其物理 delay，
   $\mathbf A_\ell(\boldsymbol\phi)\in\mathbb C^{2\times32}$ 为已经包含标准 ray power、
   2R 响应、阵元方向图、极化/XPR、ray coupling 和 $\mathbf T$ 的 TXRU 域路径矩阵；
   $\boldsymbol\phi$ 表示标准规定的随机初相。第 $b$ 个波束在该分量上的长期功率定义为

   $$
   q_{b,\ell}:=
   \mathbb E_{\boldsymbol\phi}
   \left[\left\|\mathbf A_\ell(\boldsymbol\phi)\mathbf b_b\right\|_2^2\right].
   $$

   将相同 delay 的分量相加，得到未归一化波束化 PDP 及其归一化形状

   $$
   S_b^{\rm abs}(\tau)=\sum_\ell q_{b,\ell}\delta(\tau-\tau_\ell),
   \quad
   P_b^{\rm PDP}=\sum_\ell q_{b,\ell},
   \quad
   S_b(\tau)=\frac{S_b^{\rm abs}(\tau)}{P_b^{\rm PDP}}.
   $$

   必须验证 $P_b^{\rm PDP}=P_b^{\rm ana}$；数值实现用显式随机初相平均交叉验证，但
   正式 PDP 使用解析期望，不从某次瞬时 realization 的 $|H|^2$ 反推。对选中
   $\mathcal I_8$ 的 8 个波束分别保存 $q_{b,\ell}$、合并后的 delay/power 表、总功率和
   归一化 PDP；不得先把 8 个波束相加后只保留一条平均 PDP。

3. **接收机参考 PDP。**在冻结的 8 个波束中按第 13.4 节长期 RSRP 选择

   $$
   b_\star:=\arg\max_{b\in\mathcal I_8}\overline P_b,
   \qquad S_{\rm ref}(\tau):=S_{b_\star}(\tau).
   $$

   并列时只允许使用 Gate B 已冻结排序的第一个 beam ID；不得根据 BLER、NMSE、某次
   realization 或不同 SNR 重新选参考波束。manifest 保存 $b_\star$、$S_{\rm ref}$ 的
   delay/power 数组、$P_{b_\star}^{\rm PDP}$、$\overline P_{b_\star}$ 及各自 hash。
   `BEAM8_PRECODER_CYCLING` 的8个 PRG，以及两条 CDD 的
   `cdd_reference_pdp_replicated` 接收机所假设的8个分支，都共用这一条
   $S_{\rm ref}$；`cdd_beamwise_pdp` 则明确使用全部8条真实波束化 PDP。

   UE 测量参考波束 PDP 时，以该参考 PDP 的首个可检测径 $\tau_{\rm UE}$ 作为自身
   0 时延。正式 receiver manifest 因此保存原始物理 delay 和 UE-aligned delay
   $\tau_q^{\rm ref,UE}=\tau_q^{\rm ref}-\tau_{\rm UE}$ 两套数组，LMMSE 使用后者。
   本轮已冻结的最强 beam 8 在标准 0 ns 路径上有主导功率，因此
   $\tau_{\rm UE}=0$；若实现得到其他值，必须记录检测规则和阈值，不能静默改用数学上
   极小但接收机不可检测的抽头。

由参考 PDP 在 576 个 active subcarrier 上构造归一化频率协方差

$$
[\mathbf R_{f,\rm ref}]_{k\ell}
=\sum_q p_q^{\rm ref}
\exp[-j2\pi(k-\ell)\Delta f\,\tau_q^{\rm ref}],
\qquad \Delta f=30\ \mathrm{kHz}.
$$

为沿用此前移动 TDL 的二维 LMMSE 结构，本节还从同一个参考波束的解析 ray 权重和
每条 ray 的固定 Doppler $\nu_\ell$ 构造

$$
[\mathbf R_{t,\rm ref}]_{ss'}
=\frac{1}{P_{b_\star}^{\rm PDP}}
\sum_\ell q_{b_\star,\ell}
\exp[j2\pi\nu_\ell(t_s-t_{s'})].
$$

两根 Rx 各自独立估计，但共用这一个按 2R 功率聚合的 $\mathbf R_{t,\rm ref}$；同时保存
per-Rx 参考作为失配诊断。正式估计器使用可审计的可分离假设
$\mathbf R_{t,\rm ref}\otimes\mathbf R_f$。真实 CDL 的 delay--Doppler--angle 统计一般
非可分离，8 个波束也可能相关，因此本节的 estimated 接收机是预先规定的参考统计接收机，
不是“知道完整真实 CDL 协方差”的 fully matched receiver。

本地绘图必须在 Gate B 冻结后、链路 prescan 前生成并由数值表支撑：

- 标准原始 CDL-E PDP：delay 用 ns、功率相对该 PDP 总功率归一化；
- 选中 8 个波束的 8 条归一化 $S_b(\tau)$ 同图或分面图，并突出 $b_\star$；
- 选中 8 个波束的 $S_b^{\rm abs}$ 统一除以最强波束总功率后的图，用于同时观察
  shape 与波束总增益；原始 PDP 与波束化 PDP 不共享未经说明的绝对纵轴；
- $S_{\rm ref}$ 与 B0/Sidon 两条假设 CDD 等效 PDP 的对照图。人工 delay 使用 ns，
  并明确它是循环相位对应的等效 delay，不是新增传播路径或 CP 占用。

PDP 图同时保留两个明确标注的坐标口径：原始/波束化诊断图使用标准 CDL-E 最早物理
路径的 global 0 ns；接收机参考/等效 PDP 图使用 UE 测得参考波束首径的
$\tau_{\rm UE}=0$。8 个波束只能整体减去同一个 $\tau_{\rm UE}$，不得把每个波束各自
首径分别归零后再合成 CDD；否则会人为删除真实的波束间相对物理 delay。每条 PDP
另行计算并保存总功率、相对对应坐标零点的平均 delay
$\bar\tau=\sum_qp_q\tau_q$ 和 RMS delay spread

$$
\tau_{\rm rms}
=\sqrt{\sum_qp_q(\tau_q-\bar\tau)^2},
\qquad \sum_qp_q=1.
$$

RMS 围绕该 PDP 自身的功率加权均值计算，因此对上述 global/UE 公共平移不变；坐标
零点和各自均值分别服务于时间对齐/平均 delay 与 RMS，二者不冲突。CDD 等效 PDP 还必须按
$T_0=1/\Delta f$ 周期和圆周距离
$d_{T_0}(x,y)=\min_{n\in\mathbb Z}|x-y+nT_0|$ 定义
$\bar\tau_{\rm circ}=\arg\min_{\mu\in[0,T_0)}\sum_qp_qd_{T_0}(\tau_q,\mu)^2$ 及
$\tau_{\rm rms,circ}=\sqrt{\min_\mu\sum_qp_qd_{T_0}(\tau_q,\mu)^2}$；若本轮全部复合
delay 位于同一无 wrap 区间，则同时报告普通 RMS，并验证两种展开口径一致。
PDP 形状、平均 delay 或
RMS 的方案间差异只作接收机假设与机理诊断，不设相似性门槛，也不因 PDP 不同而停止
prescan/formal；只有 profile 读取失败、功率非正/未归一、解析与直接计算不一致或 hash
失配等实现正确性问题才阻断后续阶段。

### A.6 坐标系、阵列和 profile 审计

在任何性能仿真前必须自动完成以下测试：

1. **角度语义：**明确 AoD/AoA 的零点、正方向和 wrap 区间，明确 ZoD/ZoA 是天顶角
   而不是 elevation，并核对 downlink 下 BS 使用 departure angle、UT 使用 arrival angle；
2. **单射线校准：**分别在 boresight、水平正/负空间频率、垂直正/负空间频率
   注入单一已知射线，检查理论最强 beam index 与 $\mathbf H\mathbf b$ 直接计算一致；
3. **索引和共轭：**检查 vertical-major/horizontal-minor flatten 顺序、两极化端口顺序、
   $\mathbf H\mathbf b$ 与 $\mathbf b^H\mathbf R_t\mathbf b$ 的共轭关系，以及 DFT 正负号；
4. **TXRU 映射：**核对 $\mathbf T^H\mathbf T=\mathbf I_{32}$，每个 TXRU 只驱动预定 4 个垂直
   AE，每列权重为 $1/2$，双极化之间不串接；
5. **全阵列图：**用直接 AE 求和与“TXRU 阵列因子 $\times$ 4-AE 子阵因子”两种路径
   生成 16 个波束图，并验证主瓣、零点和栅瓣一致；
6. **profile 原样：**对保存的 CDL-E 原始表与 Sionna 实例内展开数组做 SHA-256、shape、
   单位、cluster/ray 数、功率归一化、LoS/K-factor 和 ray-coupling 对照，并确认未应用
   `mean_aod_deg` 或任何 angle scale；
7. **解析协方差：**对一个极小的合成 CDL 例子，用显式随机相位平均验证
   $\mathbf R_t^{\rm ana}$；再对正式 CDL-E 检查 Hermitian、半正定、trace 与 16 个
   $P_b^{\rm ana}$ 的直接 ray-sum 一致性。

图和数值中只使用通过上述测试的坐标标签。若单个 DFT beam 存在多个同级栅瓣，
则标出全部方向，不为了匹配 CDL cluster 而只保留其中一个。

### A.7 32T 到 8 个等效信道的严格表述

冻结 $\mathcal I_8$ 后，对其 index 升序排列得到
$\mathbf B_8=[\mathbf b_{i_0},\ldots,\mathbf b_{i_7}]$；不按带 Monte Carlo 噪声的 RSRP 名次
为 delay 重新排序。第 $m$ 个波束域等效信道为

$$
\widetilde{\mathbf h}_{m}[k]
=\mathbf H[k]\mathbf b_{i_m}\in\mathbb C^{2},
\qquad m=0,\ldots,7.
$$

对 delay index 向量 $\mathbf j=(j_0,\ldots,j_7)$，定义

$$
\mathbf v_{\mathbf j}[k]
=\frac{1}{\sqrt8}
\begin{bmatrix}
e^{-j2\pi k j_0/576}&\cdots&e^{-j2\pi k j_7/576}
\end{bmatrix}^{T},
\qquad
\mathbf w_{\mathbf j}[k]=\mathbf B_8\mathbf v_{\mathbf j}[k].
$$

因为 $\mathbf B_8^H\mathbf B_8=\mathbf I_8$，所以对每个活动子载波都有
$\|\mathbf w_{\mathbf j}[k]\|_2^2=1$，不应再使用依赖子载波或信道的功率归一化。
第 $r$ 根接收天线的最终等效信道为

$$
g_{r,\mathbf j}[k]
=\mathbf e_r^T\mathbf H[k]\mathbf w_{\mathbf j}[k]
=\frac{1}{\sqrt8}\sum_{m=0}^{7}
e^{-j2\pi k j_m/576}\widetilde h_{r,m}[k].
$$

这一表述将“先用 8 个预编码将 32T 降为 8 个等效分支，再对每个分支
施加一条人工时延”与单个 32 维、随频率变化的 rank-1 预编码严格等同。

### A.8 候选、接收机和公平性

三个方案只使用同一个冻结 $\mathbf B_8$：

| 方案 ID | 定义 |
|---|---|
| `BEAM8_B0_QC` | $\mathbf j=[0,9,18,27,36,45,54,63]$ |
| `BEAM8_S0_SIDON` | $\mathbf j=[0,1,3,7,12,20,30,65]$ |
| `BEAM8_PRECODER_CYCLING` | 48 PRB 按 6 RB 分为 8 个 PRG，第 $m$ 个 PRG 使用 $\mathbf b_{i_m}$ |

CDD 的 delay 与波束映射、cycling 顺序均按选中 beam index 升序冻结。此顺序
只为消除运行时任意性，不声称是非等功率波束的最优 delay permutation。如果后续要研究
permutation 敏感性，必须在查看正式 BLER 前另行冻结全部 permutation 或对称子集。

三条 ideal-CSI 曲线直接使用每个 data RE 上的真实 $g_{r}[s,k]$。estimated-CSI
接收机冻结为：

- **透明 precoder cycling：**接收机知道 6-RB PRG 边界，不知道当前 beam ID，也不使用
  该 beam 的真实 PDP。每个 PRG 内使用同一个 $\mathbf R_{f,\rm ref}$ 和
  $\mathbf R_{t,\rm ref}$，只用该 PRG 内两个 DMRS symbol 的 pilot 做二维时频 LMMSE，
  并只估计该 PRG 内的 data RE；不跨 PRG 插值或联合估计。60 km/h 下两个 DMRS
  不平均，均作为独立时频坐标的 noisy LS 观测输入。
- **CDD 非透明接收机：**B0 和 Sidon 各自同时运行两类知道本方案 delay vector 的
  CDD-aware 接收机：`cdd_beamwise_pdp` 使用8条真实波束化 PDP及其长期总功率；
  `cdd_reference_pdp_replicated` 只使用最强波束参考 PDP，并假设8个分支均具有该 PDP
  和参考波束总功率。两类均不使用瞬时 CSI 构造先验，也均不运行未知 CDD delay 的
  transparent-CDD estimated 曲线。

对 CDD delay $j_m$，人工 delay 为

$$
\tau_m^{\rm CDD}=\frac{j_m}{576\Delta f},
\qquad \Delta f=30\ \mathrm{kHz}.
$$

在离散频率采样意义下，delay 以 $1/\Delta f$ 为周期；本轮数值若未跨周期，仍须保留
同一 circular-delay 定义。两种非透明接收机的差别冻结如下。

1. **8波束真实 PDP 接收机 `cdd_beamwise_pdp`。**第 $m$ 个 CDD 分支使用所映射
   beam $i_m$ 的真实波束化 PDP。8个波束共同减去同一个 UE 定时原点
   $\tau_{\rm UE}$，并保留不同波束的长期总功率。其未归一化等效 PDP 为

$$
S_{{\rm beam},\mathbf j}^{\rm abs}(\tau)
:=\frac{1}{8}\sum_{m=0}^{7}P_{i_m}^{\rm PDP}
\sum_qp_{i_m,q}\delta_{1/\Delta f}
\left(\tau-[(\tau_{i_m,q}-\tau_{\rm UE})+\tau_m^{\rm CDD}]\right).
$$

其总先验功率为 $P_{\rm beam}=8^{-1}\sum_mP_{i_m}^{\rm PDP}$；画 PDP shape 时才除以
$P_{\rm beam}$ 归一化，构造 LMMSE 协方差时保留绝对功率。等价全带频率协方差为

$$
[\mathbf R_{f,\mathbf j}^{\rm beam}]_{k\ell}
=\frac{1}{8}\sum_{m=0}^{7}P_{i_m}^{\rm PDP}\sum_qp_{i_m,q}
\exp\{-j2\pi(k-\ell)\Delta f
[(\tau_{i_m,q}-\tau_{\rm UE})+\tau_m^{\rm CDD}]\}.
$$

2. **最强波束参考 PDP 复制接收机 `cdd_reference_pdp_replicated`。**接收机只测得
   $b_\star$ 的 $S_{\rm ref}$，并假设8个分支均具有该归一化 PDP、相同总功率
   $P_{\rm ref}$ 和相同 $\mathbf R_{t,\rm ref}$。其归一化等效 PDP 为

$$
S_{{\rm refrep},\mathbf j}(\tau)
=\frac{1}{8}\sum_{m=0}^{7}\sum_q p_q^{\rm ref}
\delta_{1/\Delta f}
\left(\tau-\tau_q^{\rm ref,UE}-\tau_m^{\rm CDD}\right),
$$

相同 circular delay 的功率必须合并。其带绝对功率的等价频率协方差为

$$
\mathbf R_{f,\mathbf j}^{\rm refrep}
=P_{\rm ref}\mathbf R_{f,\rm ref}\odot\mathbf A_{\mathbf j},
\qquad
[\mathbf A_{\mathbf j}]_{k\ell}
=\frac{1}{8}\sum_{m=0}^{7}
\exp\left[-j2\pi(k-\ell)\frac{j_m}{576}\right].
$$

上述 beamwise 等效 PDP 的归一化 shape 为

$$
S_{{\rm beam},\mathbf j}(\tau)
=\frac{1}{\sum_mP_{i_m}^{\rm PDP}}
\sum_{m=0}^{7}P_{i_m}^{\rm PDP}
\sum_qp_{i_m,q}
\delta_{1/\Delta f}
\left(\tau-[(\tau_{i_m,q}-\tau_{\rm UE})+\tau_m^{\rm CDD}]\right).
$$

这里 8 个波束共同减去同一个参考波束定时原点 $\tau_{\rm UE}$，因此保留波束间相对
物理 delay；不能分别把每个 $S_{i_m}$ 的首径归零。两类接收机都假设8个波束域分支
相互独立，并共用第 13.5 节的 $\mathbf R_{t,\rm ref}$。`cdd_beamwise_pdp` 的“真实”只表示
使用各波束真实长期 PDP/功率和相对 delay；真实 CDL CDD 的最终二阶统计还含
$\mathbb E[\widetilde h_m\widetilde h_{m'}^*]$ 交叉项，所以
$\mathbf R_{f,\mathbf j}^{\rm beam}$ 仍不是 fully matched covariance。两类接收机唯一允许的
先验差别是“8条 beam PDP/功率”与“最强 beam PDP/功率复制”；时间协方差、CDD delay、
pilot、噪声、numerical loading 和算法实现必须相同。必须报告两类等效 PDP 的
mean/RMS、频率协方差相对 Frobenius 误差、frequency-only analytic trace-NMSE，以及最终
CE NMSE/BLER 差；PDP 差异不作为停止条件。

实现必须分别验证两种“从移位 PDP 直接求协方差”与上述矩阵公式逐元素一致。B0 和
Sidon 的两类 CDD-aware 接收机都在全部576个 active subcarrier 上做二维时频 LMMSE。
对某个 receiver 的全带 CDD 或某个 cycling PRG 的 pilot/data 坐标集 $P,D$，令
$\mathbf R_a=\mathbf R_{t,\rm ref}\otimes\mathbf R_f^{(u)}$ 并按实际坐标截取，其中
$u\in\{\mathrm{beam},\mathrm{refrep}\}$；cycling 使用
$\mathbf R_f^{(u)}=P_{\rm ref}\mathbf R_{f,\rm ref}$。因此 beamwise receiver 的协方差
对角线为 $P_{\rm beam}$，reference-replicated receiver 和 cycling 为 $P_{\rm ref}$，
不得为方便比较把 beamwise 协方差重新缩放到 $P_{\rm ref}$。LMMSE 为

$$
\widehat{\mathbf g}_D
=\mathbf R_{a,DP}
\left(\mathbf R_{a,PP}+N_0\mathbf I\right)^{-1}\mathbf z_P.
$$

cycling 的 $P,D$ 只取当前 PRG，CDD 的 $P,D$ 取全带；两个 DMRS 保持不同 symbol 坐标，
不先平均。实现使用 Cholesky solve，不显式求逆，并保存 condition number、最小奇异值和
numerical loading。

对 absolute trial $t$，CE NMSE 在两根 Rx 和全部 data RE 上先求 trial 级线性能量比

$$
r_t=
\frac{\sum_{r=1}^{2}\sum_{(s,k)\in D}
|\widehat g_{t,r,s,k}-g_{t,r,s,k}|^2}
{\sum_{r=1}^{2}\sum_{(s,k)\in D}|g_{t,r,s,k}|^2}.
$$

每个 `scheme + receiver + SNR` 点在线性域报告 $\overline r=T^{-1}\sum_t r_t$、95% Monte Carlo
区间和 $10\log_{10}\overline r$。逐区间保存 error energy、true energy、$r_t$、和及平方和；
不得平均分段 dB 值，也不得用 $\sum E_t/\sum S_t$ 代替 $T^{-1}\sum r_t$。

每个 `SNR + absolute trial` 内三方案共享同一 CDL realization、payload、编码比特和
DMRS/data AWGN。三方案使用相同的 DMRS/data RE、MCS、译码器、2Rx MRC、
单位总发射功率和 trial 区间；同一方案的 ideal/estimated 使用相同真实信道、payload
和 data noise，同一 CDD 方案的两个非透明 receiver 还必须复用完全相同的 DMRS LS
观测和 absolute trials，只改变 LMMSE filter。统计选择样本与正式链路 realization seed
必须不同。

公共 SNR 参考功率定义为 $P_{\rm ref}=\max_{0\le b<16}\overline P_b$，三方案共用

$$
N_0=P_{\rm ref}10^{-\gamma_{\rm ref,dB}/10}.
$$

此定义取代当前模式 C 针对 8 个 parent SSB 的 $P_{\rm ref}$ 定义，只在本增补的
codebook type 下生效。正式输出必须同时保存 $P_{\rm ref}$、16 个 $\overline P_b$和实际
`noise_var`；禁止按方案或 trial 重新定标噪声。

### A.9 实现、测试和阶段门

不复制新的 CDL 或 PDSCH 主循环。在模式 C 的通用实现中扩展：

1. 新增显式 codebook type `dft_2x8_same_pol`、profile-native angle 模式、任意冻结 8 列的
   `BEAM8_PRECODER_CYCLING` 以及八分支 `BEAM8_CDD`；原有 parent/secondary 行为和 smoke 保持不变；
2. 新增解析 CDL ray-power/协方差审计入口，它只读 profile、阵列和 codebook，不读 BLER；
3. 新增第 13.5 节标准/波束化/参考/等效 PDP 的提取、RMS 汇总和本地绘图入口；
4. 扩展公共二维时频 LMMSE，使其支持参考统计下的逐 PRG transparent cycling、
   全带 `cdd_beamwise_pdp` 和全带 `cdd_reference_pdp_replicated` 三条路径，并复用
   现有 CE trial statistics、2Rx MRC 和 interval 落盘逻辑；
5. 新增长期 RSRP 分批统计、realization-level bootstrap、top-8 冻结 manifest 及本地绘图入口；
6. 新增 `configs/plan037_cdl_e_32tx_2rx_beam8_*.yaml` 和一个专题执行入口；配置必须
   通过冻结 manifest hash 加载波束，运行时不重新选择；
7. 扩展 result-037 分析入口，生成 RSRP/角度/PDP 核验、ideal/estimated BLER crossing、
   paired gain、ideal/estimated gap、CE NMSE 和数据完整性审计。

除第 13.6 节外，定向测试至少覆盖：

- $\mathbf B^H\mathbf B=\mathbf I_{16}$、$\mathbf B_8^H\mathbf B_8=\mathbf I_8$、双极化权重相同，
  以及 16 个 beam ID 与 $(q_v,q_h)$ 一一对应；
- 对全部 576 个子载波，B0、Sidon 和 cycling 的 $\|\mathbf w[k]\|_2^2=1$；
- 通过直接 $\mathbf H\mathbf w$ 与先算 8 个 $\widetilde{\mathbf h}_m$ 再加 delay 的结果相等；
- B0/Sidon delay 数组、`/576` 相位分母、index 升序映射、8 个 6-RB PRG 和无跨界 cycling；
- 标准 CDL-E PDP 的 K-factor/LoS 展开无重复计数，8 条 beam PDP 的总功率等于解析
  beam power，统一零点、各自均值 RMS 和 circular RMS 的手算小例子正确；
- strongest-beam $S_{\rm ref}$ 冻结不依赖 SNR/trial，cycling 的8个 PRG 使用同一个参考 PDP；
  B0/Sidon 的 beamwise 与 reference-replicated 两套移位 PDP 分别与各自频率协方差公式一致；
- UE timing 只从参考波束确定一次，8条 beam PDP 共同减去同一 $\tau_{\rm UE}$；
  beamwise receiver 保留波束间相对 delay 和长期功率，禁止逐 beam 重置首径或等功率化；
- cycling 滤波器不跨 PRG，CDD 滤波器覆盖全带；两个 60 km/h DMRS 不平均；两根 Rx
  独立估计后再 MRC；ideal 不生成 CE，estimated 的 trial NMSE 与手算结果一致；
- CDD 配置必须且只允许 `cdd_beamwise_pdp`、`cdd_reference_pdp_replicated` 两类
  receiver，拒绝 transparent-CDD；前者 manifest 必须含8条 beam PDP/power hash，后者
  只能含 reference-PDP hash，不得泄漏其他7条 PDP；两者均不得使用完整 CDL 交叉协方差；
- statistics/formal seed 分离、realization-level 统计、manifest hash 拒绝路径和三方案 paired identity；
- 原有两份 `fixed_cdl_statistics` smoke 及 `tests/test_fixed_cdl_statistics.py` 回归通过。

执行分为不可跨越的三个门：

1. **Gate A：几何与解析审计。**第 13.3、13.6 节全部测试通过；
2. **Gate B：长期统计。**第 13.4 节解析/Monte Carlo 验收通过且 top-8 可识别，
   完成第 13.5 节 PDP/RMS 数值产物并冻结8条 beam PDP/power、$b_\star$、$S_{\rm ref}$、
   两种 CDD 等效 PDP 构造规则和全部 hash，然后冻结 manifest；
   8 条 PDP 彼此不同不构成失败条件；
3. **Gate C：链路。**三条 ideal 曲线、cycling 的透明 estimated 曲线及两条 CDD 各自的
   两类 estimated 曲线完成 validate、滤波器、功率、配对、CE 统计和每条曲线2 trials
   的 smoke，才允许 prescan 和 formal。

### A.10 BLER/CE 预扫、正式预算和判定

prescan 使用公共网格 `[-4,-2,0,2,4,6,8,10,12,14] dB`。3条 ideal 加5条 estimated
共8条 BLER 曲线每点运行400个共同 paired trials；同一 CDD transmitter 的两个 receiver
不重新生成信道或噪声。若任一 BLER 曲线的 10% 或 1% crossing 没有双侧 bracket，只向缺失方向
按 2 dB 扩展，硬边界为 `[-10,22] dB`；到边界仍缺失则停止，不外推。

prescan 后为上述8条 BLER 曲线冻结共用 formal 网格：10%/1% crossing 附近间隔不大于
0.25 dB，其余过渡区不大于 0.5 dB。每点先运行 1000 trials，以 1000 为批次
共同追加，单点上限 50000。10% bracket 端点在 5%--20% BLER 时至少 200 errors；
1% 端点在 0.5%--2% 时至少 200 errors。使用真实相邻点之间的 log-BLER 线性插值，
不外推、平滑或单调修正。单点报 Wilson 95% 区间，crossing 和增益使用至少
2000 次 absolute-trial paired bootstrap。

主比较定义为

$$
G_{X\leftarrow\mathrm{cycling}}(p)
=\mathrm{SNR}_{\mathrm{BEAM8\_PRECODER\_CYCLING}}(p)
-\mathrm{SNR}_{X}(p),
\quad X\in\{\mathrm{BEAM8\_B0\_QC},\mathrm{BEAM8\_S0\_SIDON}\},
\quad p\in\{10\%,1\%\}.
$$

$G$ 对 ideal CSI 计算一次；对 estimated CSI 则分别以
$u\in\{\mathrm{beamwise},\mathrm{reference\ replicated}\}$ 计算
$G_{X,u\leftarrow\mathrm{cycling}}$，两者使用同一条 cycling 透明接收机曲线作为基线。
其95%区间完全大于0/跨0/完全小于0，分别表述为支持正增益/排序不确定/支持劣化；
两类 estimated-CSI $G$ 均为正式性能结论，ideal-CSI $G$ 是分集诊断。另报同一 CDD
transmitter 下 `beamwise - reference_replicated` 的 paired crossing 差和 CE NMSE 差，
B0 与 Sidon 之间的 paired 差、各 receiver 的 estimated-minus-ideal 目标 SNR penalty、
目标 crossing 邻域 CE NMSE 及相对 cycling 的 NMSE 差。CE NMSE 不设优劣硬门槛；
不用 RSRP、PDP/RMS、角度图或 ideal-CSI 均值接收功率代替 BLER 判定。

formal trial 1 前必须回填并由研究者确认：统计/formal seed、实际 $D$、
$\mathcal I_8$、排序、$b_\star$、8条 beam-PDP hash、reference-PDP hash、两类 CDD filter
manifest hash、Gate A/B/C 报告路径、
formal SNR 网格、batch size、
峰值 RSS、配置 SHA-256 和总预算上限。

### A.11 输出、result 和停止条件

本增补的所有产物写入
`outputs/experiment037_cdl_e_32tx_beam8/<run_id>/<stage>/E100_NT32_NR2_V60_BEAM8/`，
不覆盖第 10 节的 TDL 输出。至少保存：

- 原 YAML、resolved config、代码/工作区标识、Sionna 版本和全部 SHA-256；
- CDL-E profile 原表与展开 ray 参数、坐标系审计、$\mathbf T$、$\mathbf B$、
  $\mathbf R_t^{\rm ana}$、16 个解析波束功率和本地绘图数据；
- realization-level $P_{d,b}$、$\widehat{\mathbf R}_{t,d}$ 或其可无损重建汇总、bootstrap 结果、
  top-8 冻结 manifest 和选择稳定性记录；
- 原始 CDL-E PDP、8条波束化 PDP、$S_{\rm ref}$，以及 B0/Sidon 在 beamwise、
  reference-replicated 两种算法下的等效 PDP 精确 delay/power CSV；同时保存 global/UE
  两套时间原点、mean delay、普通/circular RMS、总功率、协方差误差、frequency-only
  trace-NMSE 和全部输入/输出 hash；
- 三类 receiver 的协方差/filter manifest、condition number、numerical loading 和解析
  posterior trace-NMSE 或可用的直接矩阵自检；
- 三方案完整 precoder manifest、逐子载波功率审计、$P_{\rm ref}$、`noise_var`、
  interval/trial metrics、ideal/estimated 错块计数、estimated CE trial arrays、配对审计、
  BLER/NMSE/crossing/gain CSV 和日志。

本地脚本生成第 13.4 节的三类图、第 13.5 节 global/UE-aligned PDP、RMS 与两种
CDD-aware 等效 PDP 对照图、三方案 ideal-CSI BLER、5条 estimated-CSI BLER、5条
CE NMSE 和10%/1% gain 图；同一 CDD 的两个 receiver 使用相同颜色、不同线型成对展示。
仿真数据和图不上传到云端服务。结果追加到现有成对
`research/result-037-PDSCH-Tx数影响.md` 和 `research/result-037-PDSCH-Tx数影响-text.md`，
必须将本 CDL-E 增补与原 TDL 4T1R/16T/32T 任务分开标注，不合并成跨场景 paired 结论。

出现以下任一条件立即停止受影响阶段：

- profile-native 参数不能完整取得，或 profile/angle transform 审计失败；
- Gate A 的坐标、共轭、极化、AE-to-TXRU 或波束图测试失败；
- 解析波束功率与 Monte Carlo 在排除坐标问题后仍不满足第 13.4 节门槛；
- $D=10000$ 时 top-8 仍不可识别；
- 正式运行的 codebook/选择/delay/PRG manifest hash 与冻结值不同；
- 标准/波束化/等效 PDP 的功率、LoS 展开、协方差等价或 reference hash 审计失败；
- estimated 滤波器非有限、PRG 越界、误把两个移动 DMRS 平均或 CE 统计口径错误；
- 任一方案逐子载波单位总功率失败，三方案的 `noise_var` 不同，或 trial 配对失败；
- prescan 在硬边界内缺 bracket，或最小 batch size 1 仍资源不可行。

未通过角度/RSRP 核验的负结果仍是本实验的有效诊断产物，必须在 result 中报告；
但不得在该情况下继续选波束或宣称 B0/Sidon 的 CDL-E 性能结论。
