# plan-028：comb-6 三类 CSI 曲线统一补充与 TDL-A 300 ns 方案对比

> 状态：已执行，新增结果待研究者确认。原计划已确认并执行；研究者于 2026-08-05 确认第 10 节“透明 PRG precoder cycling 基线”，其正式仿真、稳定性审计和 result-028 第 3.6 节三张同图增补现已完成。研究者于 2026-08-10 要求按第 11 节增加“CDD 发射不变、UE 使用底层物理协方差且不知道 CDD delay”的透明 CDD 基线；五组 delay set 的 estimated-CSI BLER、CE NMSE、逐 trial 审计及第 3.6 节同图增补现已完成。研究者于 2026-08-12 要求按第 16 节为固定小时延 CDD 补充 ideal CSI，以及协方差匹配的 estimated CSI/CE；正式运行、审计和第 3.6 节同图增补均已完成。研究者于 2026-08-03 明确要求：基于 result-027 的 A30 comb-6 与 A100 comb-6 正式数据，按统一大字号样式重绘 estimated-CSI BLER，并在相同 SNR 采样点补充真实 ideal-CSI 链路 BLER 与 CE NMSE；新增 TDL-A 300 ns、comb-6 的各方案对比，交付缩写图例和 `candidate NN` 图例版本。未单独给出的系统、统计预算、搜索与闭合规则按 A100 comb-6 增补保持不变。

## 1. 目的、研究问题与可证伪假设

本轮对应 `GOALS.md` 的物理展宽信道结构化 CDD 设计与 estimated-CSI 验证，并遵循 `DESIGN.md` 第 3.7、3.8 和第 6 节：ideal-CSI、matched CE NMSE 与 estimated-CSI BLER 是不同判据；同一场景内所有方案必须固定 DMRS、功率、MCS、数据 RE、接收机知识、trial 和随机样本。

研究问题：

1. result-027 的 A30 comb-6 与 A100 comb-6 各物理方案，在同一 SNR 采样网格下，ideal-CSI BLER、matched CE NMSE 与 estimated-CSI BLER 的曲线关系是什么？
2. 把 TDL-A RMS delay spread 增至 300 ns 后，按当前 PDP 与 comb-6 重新生成或搜索的方案，哪些能够形成 10%/1% estimated-CSI BLER 双侧 bracket，哪些受 CE 或 pilot rank 限制？
3. A300 中各搜索 family 相对同场景六基线最优者是否保留正的 estimated-CSI BLER 目标 SNR 增益？

预定可证伪假设：

- H1：真实 ideal-CSI 链路曲线通常不差于相同物理方案的 estimated-CSI 曲线；若有限样本下出现反向，必须报告误差计数和区间，不以该假设覆盖事实。
- H2：A300 的 `AP_RMS_T1`、`AP_TEPS_T1` 以及三个搜索 family 的精确时延会随 PDP 改变；直接复用 A100 的这些时延会造成不公平比较。
- H3：comb-6 的高导频密度不保证 A300 搜索候选在 estimated-CSI 下闭合或优于六基线最优者；任一目标未形成正式双侧 bracket 时不得外推。

## 2. 固定系统、物理定义与公平性

三个场景均固定：48 PRB、$K=576$、子载波间隔 30 kHz、8 Tx / 1 Rx、单层、static Sionna TDL-A、零速度、DMRS comb-6、两个 DMRS symbol、每 symbol 96 个 pilot、总计 192 pilot RE、5568 data RE、16QAM、MCS 8、码率 553/1024（取自 NR 256QAM MCS table）、Sionna LDPC 最多 8 次迭代。每个 pilot RE 功率不变。estimated-CSI 接收机为 two-DMRS averaged、known-$\mathbf V$/known-PDP、V-aware matched 全带频域 LMMSE；LLR 噪声方差不额外加入 CE 误差项。

CDD 相位固定为

$$
V_{k,n}=\exp\left(-j2\pi k\frac{j_n}{576}\right),
$$

相位分母为 576，参考点为首个 active subcarrier。$j_n$ 为允许连续值的 DFT 网格坐标，一个坐标单位为 $1/(576\Delta f)=57.870370\ldots$ ns。连续时延不得量化；`AP_TU_NTM1` 的末端 576 不得模回 0。

三个接收机/指标定义为：

1. **estimated-CSI BLER**：使用上述 matched LMMSE 估计值做均衡和译码；
2. **ideal-CSI BLER**：仍实际生成 TDL-A 信道、payload、编码、调制、AWGN、均衡、软解调和 LDPC 译码，但均衡器直接使用每个 data RE 的真实等效信道；禁止用 outage、互信息或 estimated-CSI 数据替代；
3. **CE NMSE**：data RE 上 matched LMMSE 估计误差功率与真实等效信道功率之比，先按 trial 在线性域求比值并平均，再换算为 dB。

同一场景、SNR、trial index 下，各方案共享 payload、TDL 信道样本和噪声。ideal 与 estimated 两种接收机也使用相同 seed 派生规则。跨 A30/A100/A300 不声称 paired。

## 3. A30/A100 已有 estimated-CSI 数据与新增采样

### 3.1 内容事实源

- A30：`outputs/experiment027_meff_sidon/20260726_main/e5_dense_dmrs/comb6/`；
- A100：`outputs/experiment027_meff_sidon/20260726_main/e6_a100_dense_dmrs/comb6/`。

必须先校验原 manifest SHA-256、逐点 CSV、逐 trial error flag、场景、comb、候选时延和 trial 数。原 estimated-CSI 链路不重跑；新结果通过只读引用 result-027 原始数据形成派生数据集。

### 3.2 统一曲线采样规则

对每个物理候选，将 result-027 的 `prescan`、`refine_10pct` 和 `refine_1pct` 逐点合并；同一 `candidate_id + snr_db` 重复时保留 trial 数最大的正式点，trial 相同则优先级为 `refine_1pct > refine_10pct > prescan`。所得候选专属 SNR 序列定义为本轮冻结曲线网格，并写入机器可读 schedule。

- estimated-CSI BLER 与 CE NMSE 直接取被选中的同一原始行；
- ideal-CSI BLER 对每个 `candidate_id + snr_db` 使用与该原始行相同的 trial 数运行真实链路；
- 三种曲线必须具有完全相同的候选集合与逐候选 SNR 集合；校验失败即停止出图。

## 4. A300 comb-6 候选、搜索与冻结

### 4.1 基线与参考

A300 只把 RMS delay spread 改为 300 ns，其余系统与 A100 comb-6 相同。链路 manifest 包含：

1. 六基线 `B0_QC`、`AP_RMS_T1`、`AP_TEPS_T1`、`AP_TU_NT`、`AP_TU_NTM1`、`AP_TALIAS_NT`；
2. 固定历史参考 `S0_SIDON`，不计入六基线包络；
3. `AP_T2_CTRL`、`GEO_T1_CTRL`、`MEFF_T2_CAND` 在 10%/1% ideal-CSI outage 上各自胜出代表的去重并集。

时延规则在看 BLER 前冻结：

- `B0_QC=[0,9,18,27,36,45,54,63]q`、`S0_SIDON=[0,1,3,7,12,20,30,65]q` 和 `AP_TU_NT=[0,72,144,216,288,360,432,504]q` 保持固定；
- `AP_TU_NTM1=[0,576/7,2\times576/7,\ldots,576]q` 保持连续坐标；
- `AP_TALIAS_NT=[0,12,24,36,48,60,72,84]q` 由 comb-6 的当前 $\tau_{\rm alias}/N_t$ 给出；
- `AP_RMS_T1` 使用连续 300 ns 公差；
- `AP_TEPS_T1` 从本轮实际展开 A300 PDP 重算覆盖 $1-\epsilon$ 能量的 $T_\epsilon$，再使用连续 $T_\epsilon$ 公差；$\epsilon$、PDP 离散化和取整口径与 plan-027 保持一致并在展开配置中显式保存；
- `AP_T2_CTRL`、`GEO_T1_CTRL`、`MEFF_T2_CAND` 使用 A300 完整物理协方差、comb-6 的 $N_p=96$，按 plan-027 E5/E6 相同可行性、起点、种子派生、$M_2^{\rm eff}/M_4^{\rm eff}$、pair/fold gap 和选择规则重新搜索。

A300 的硬厚 Sidon 可行性必须由当前支撑与装填上界重新判定；软几何候选不得标为硬厚 Sidon。

### 4.2 搜索、outage、manifest 与审批留痕

- `AP_T2_CTRL` 与 `MEFF_T2_CAND` 搜索使用 200,000 个唯一 T2 状态；
- `GEO_T1_CTRL` 使用 200,000 个确定性几何状态；
- ideal-CSI outage 使用 200,000 个共同 A300 信道样本和 1,000 次成对 bootstrap；
- 输出全部冻结候选的网格坐标、ns、residue、lift、pilot rank、condition number、pair/fold gap、选择理由和 ideal-outage 目标；
- 生成独立的 plan-028 A300 comb-6 manifest、SHA-256 和研究者指令回执。manifest 只由预定搜索/outage 规则决定，禁止按 BLER 后验改动。

## 5. 链路统计预算、SNR 网格与停止条件

A300 estimated-CSI 沿用 A100 comb-6：

1. 每候选 20 trials 集成 smoke；
2. 以候选 ideal-outage 目标加 5 dB 初始化，0.5 dB 步长、每点 400 paired trials 做 prescan；
3. 10%和1%目标分别以0.25 dB步长、每点3,000 paired trials细扫；
4. 每目标最多向所需方向扩展2 dB；正式目标未形成双侧 bracket 时不外推；目标附近累计错误块少于30时标记样本不足。

A300 的 estimated-CSI 正式运行结束后，按第 3.2 节相同“最大 trial、refine 优先”规则冻结逐候选曲线 SNR 网格。随后在完全相同的逐点 trial 数上运行 ideal-CSI 链路；CE NMSE 使用被选中的 estimated-CSI 同行。若 ideal-CSI 目标在冻结网格内形成双侧 bracket，可用与 estimated-CSI 相同的 logistic 方法报告10%/1%目标和95%区间；未 bracket 不外推。

以下任一情况停止相关阶段并记录，不静默修正：manifest/hash 不一致；候选少于7或多于13条物理曲线；SNR 网格不一致；时延非有限或分支数不是8；comb不是6；相位分母不是576；smoke出现非有限 CE/BLER；正式数据缺失逐 trial flag；A300 PDP 或 $T_\epsilon$ 未按展开配置重算。

## 6. 编号、共享样式与图表

### 6.1 candidate 编号

编号在每个场景内独立，从 `candidate 01` 开始；编号不跨场景表示同一物理时延。固定排序为：六基线的上述顺序、`S0_SIDON`、`AP_T2_CTRL`、`GEO_T1_CTRL`、`MEFF_T2_CAND`；同一搜索 family 有多个代表时按原 `candidate_id` 的数值后缀升序。物理去重后的 represented family 只占一个编号。

每个场景必须交付编号表，至少包含：candidate 编号、原缩写 `candidate_id`、family、represented families、角色、`delay_grid_coordinates`、精确 `delay_ns`、pilot rank、condition number。A30/A100 的表引用原 manifest 物理事实，A300 的表引用新 manifest。

### 6.2 共享样式

所有图调用一个随结果保存的共享 `curve_styles.json`，键为场景内物理 `candidate_id`，值至少包含 color、linestyle 和 marker。编号图与缩写图只改变 label，不改变样式；estimated-CSI BLER、ideal-CSI BLER、CE NMSE 和二维散点也不得改变同一方案的样式。所有图满足 `docs/agent/RESULT_SPEC.md`：标签、标题、图例至少16 pt，刻度至少14 pt，线宽至少2.5 pt，marker至少8 pt，并按约13 cm宽预览。

### 6.3 必须交付的图与对应数据

A30 comb-6、A100 comb-6 各交付编号图例版本：

1. estimated-CSI BLER–SNR；
2. 真实 ideal-CSI BLER–SNR；
3. matched CE NMSE(dB)–SNR。

A300 comb-6 交付两套版本：

- 缩写图例版：上述三张曲线，以及与 A100 comb-6 增补一致的全候选 prescan、正式10%/1% estimated-CSI 曲线、10%和1%“matched CE NMSE at target SNR–target SNR”二维散点；
- `candidate NN` 图例版：上述 estimated-CSI BLER、ideal-CSI BLER、CE NMSE 三张曲线。

每张图都有对应 CSV。BLER 为对数纵轴，零错误点以 `0.5/trials` 的绘图下界显示但 CSV 保留真实0和错误计数；CE NMSE 纵轴单位为 dB。图标题明确场景、comb和CSI类型。图例不得遮挡数据，必要时放到图外并使用多列。

## 7. 验收判据与解释边界

必须逐项报告：

1. 三个场景的候选数、逐候选 SNR 点数、trial 数与三类曲线网格一致性；
2. 每个点的 BLER、错误块数、Wilson 95%区间；目标拟合报告方法、点数、总 trial、目标附近错误数和95%区间；
3. A300 每个 family 的10%/1% estimated-CSI闭合状态，以及相对同目标六基线最优者的目标 SNR差和95%区间；
4. ideal 与 estimated 的差异只作为实际链路事实报告，不用 outage 代替 ideal-CSI BLER；
5. A300 的 exact delay、$T_\epsilon$、rank、condition number与未闭合/CE地板的对应关系；
6. 图例编号表、共享样式一致性、PPT缩放可读性检查；
7. 所有失败、中止、样本不足、反常排序和后验发现的偏差。

本轮仍只适用于 48 PRB、static TDL-A、300 ns及两个复用场景、comb-6、零速度、当前 MCS、matched receiver。不能推广到其他 TDL/CDL、移动性、空间相关、协方差失配、定时误差或相同 DMRS 开销比较。A30/A100 的 estimated-CSI 是原数据派生重绘，不应误写为本轮重跑。

## 8. 代码、测试、输出与复现

预计新增 `tools/run_plan028_csi_curves.py` 作为本轮搜索/manifest/estimated/ideal/analyze/plot 专题入口；可复用的链路和搜索逻辑优先调用 plan-027 现有函数，若需扩展只做向后兼容修改。新增测试至少覆盖：

- A300 场景参数、实际 PDP 支撑与动态时延；
- manifest/hash/comb/相位分母/8分支校验；
- 原始 estimated 数据合并优先级与逐候选 SNR 网格完全一致；
- ideal 接收机确实使用真实 data-RE 等效信道；
- candidate 编号双射、represented-family 去重；
- 同一物理方案跨所有图的 color/linestyle/marker 完全一致；
- 最小字体、线宽和 marker 约束。

输出写入 `outputs/experiment028_csi_curves/20260803_main/`，至少包含 `a30_comb6/`、`a100_comb6/`、`a300_comb6/`、`final/`、展开配置、manifest/hash/回执、schedule、逐点 CSV、逐 trial flags、日志、共享样式、编号表和图对应数据。Markdown 引用图复制到 `docs/figures/result-028/`。

建议复现入口：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage import-existing --scenario A30 --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage import-existing --scenario A100 --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage search --scenario A300 --run-id 20260803_main --search-states 200000 --geometry-states 200000
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage outage --scenario A300 --run-id 20260803_main --outage-samples 200000 --bootstrap-repeats 1000
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage manifest --scenario A300 --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage estimated-<smoke|prescan|refine-10pct|refine-1pct> --scenario A300 --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage ideal --scenario <A30|A100|A300> --run-id 20260803_main
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan028_csi_curves.py --stage analyze --scenario <A30|A100|A300> --run-id 20260803_main
```

## 9. result-028 必须回答

`research/result-028-comb6三类CSI-TDL-A300ns.md` 与 `research/result-028-comb6三类CSI-TDL-A300ns-text.md` 必须在同一个新增主章节中集中给出 A30 comb-6、A100 comb-6、A300 comb-6 的编号与时延表、三类曲线、原始采样表、目标和不确定性；图文版嵌入图，无图版不得嵌图但必须引用全部 CSV 和关键数值。两版结论、异常、验收、适用范围和证据路径必须一致，并更新 `research/README.md` 索引。未经研究者确认，不更新 `KNOWLEDGE.md` 已验证结论或 `GOALS.md` 阶段验收。

## 10. A100 透明 PRG precoder cycling 基线

### 10.1 目的、范围

本节为 result-028 第 3.6 节的补充任务，只在现有 A100、48 PRB、comb-6、static TDL-A 100 ns、8 Tx / 1 Rx / 1 layer、两个 DMRS symbol 和相同 MCS 下增加两条透明 precoder cycling 基线，不改变原五组 delay set 的数据、定义或编号。新增基线稳定标识为：

1. `A100_PRG_DFT8_4RB`：一个 precoding resource block group（PRG）为 4 RB，共 12 个 PRG；
2. `A100_PRG_DFT8_6RB`：一个 PRG 为 6 RB，共 8 个 PRG。

这里“透明”特指：在当前各发射天线独立、各天线具有相同 PDP、同一 PRG 内使用固定空域预编码向量的条件下，接收端不需要分离 8 根天线的底层信道，也不需要根据 DFT 向量修改 PDP 形状；接收端只需按 PRG 边界估计标量等效信道，并使用底层物理 PDP 形状及已知功率因子构造 LMMSE。透明不表示跨 PRG 的等效信道连续，也不表示该结论可推广到天线空间相关、各天线 PDP 不同、PRG 内切换预编码或移动信道。



### 10.2 发射端 DFT 码本、功率归一化与 PRG 映射

令天线索引和 DFT 向量索引均为 $n,m\in\{0,1,\ldots,7\}$。先定义酉归一化 8 点空域 DFT 矩阵的一列：

$$
V^{\rm unitary}_{n,m}
=
\frac{1}{\sqrt 8}
\exp\left(-j2\pi\frac{nm}{8}\right).
$$

当前 plan-028 的 CDD 链路使用每根天线预编码系数模长为 1 的约定，等效信道平均功率为单天线功率的 8 倍，并使用 `noise_variance = 8 / snr_linear`。为保持相同总发射功率、每个 DMRS RE 功率和 SNR 定义，本节实际接入链路的第 $m$ 个 DFT 向量固定为

$$
V^{\rm tx}_{n,m}
=
\sqrt 8 V^{\rm unitary}_{n,m}
=
\exp\left(-j2\pi\frac{nm}{8}\right),
\qquad
\left\|\mathbf V^{\rm tx}_{:,m}\right\|_2^2=8.
$$

因此不得把单位范数 DFT 向量直接接入后仍沿用 `noise_variance = 8 / snr_linear`。等价的另一种实现是使用单位范数向量并同时令 `noise_variance = 1 / snr_linear`；本轮为与既有 CDD 数据逐点比较，冻结采用前一种实现。所有数据和 DMRS RE 均通过同一个 PRG 向量发送，PRG 内不得只对数据而不对 DMRS 做预编码。

以 active-band 局部 RB 索引 $r=0,1,\ldots,47$ 和局部子载波索引 $k=0,1,\ldots,575$ 定义 PRG；PRG 边界从第一个 active RB 开始，不按物理 FFT bin 重新对齐。若 PRG 大小为 $B$ RB，则

$$
b=\left\lfloor\frac{r}{B}\right\rfloor
=
\left\lfloor\frac{k}{12B}\right\rfloor.
$$

两条基线的 DFT 向量映射规则冻结为：

- `A100_PRG_DFT8_6RB`：$B=6$，8 个 PRG 依次使用向量索引 `[0,1,2,3,4,5,6,7]`，每个 PRG 72 个连续 active subcarrier；
- `A100_PRG_DFT8_4RB`：$B=4$，前 8 个 PRG 在所有 trial 中依次使用 `[0,1,2,3,4,5,6,7]`；后 4 个 PRG 在每个 absolute trial 开始时从 8 个向量中重新做一次有序、无放回抽取。该 trial 抽到的 4 个索引在该 trial 的全部 OFDM symbol、两个 DMRS、全部 data RE、ideal/estimated 接收机中保持不变；下一个 absolute trial 重新抽取。不得逐 PRG 独立有放回抽样，也不得逐 OFDM symbol 或在 ideal/estimated 分支分别重抽。抽样 RNG 固定为 `numpy.random.Generator(PCG64)`，seed 由 `stable_seed(20260727, "A100_PRG_DFT8_4RB", snr_db, absolute_trial_index, "tail_prg_order")` 唯一派生。manifest 保存 RNG、base seed、派生字段和抽样算法，逐 trial 输出保存实际 4 元组及其哈希。

因此 `[4,7,1,0]` 只可作为 `numpy.random.default_rng(20260727).choice(...)` 的说明性单次样例，不再是正式仿真中冻结的尾部映射。正式比较评价的是“每 trial 尾部随机、trial 内固定”的随机 precoder cycling 策略对 BLER 的平均效果，而不是某一个尾部映射 realization。

对 PRG $b$ 内的子载波 $k$，标量等效信道为

$$
g_k
=
\sum_{n=0}^{7} h_{k,n}V^{\rm tx}_{n,m_b},
$$

其中 $m_b$ 是该 trial 内该 PRG 使用的 DFT 向量索引。发射端不根据 A100 的瞬时信道、PDP 或 SNR 自适应选择 $m_b$；6-RB 基线使用固定 cycling，4-RB 基线使用预先规定的“前 8 个固定 cycling、尾部 4 个逐 trial 随机”策略。

### 10.3 PRG 内等效协方差与 PDP 形状

本轮底层信道假设为 8 根发射天线相互独立且具有相同频域物理协方差：

$$
\mathbb E\!\left[h_{k,n}h_{\ell,n'}^*\right]
=
\begin{cases}
R_{{\rm phy},k\ell}, & n=n',\\
0, & n\ne n'.
\end{cases}
$$

因为一个 PRG 内 $\mathbf V^{\rm tx}_{:,m_b}$ 不随子载波改变，所以对同一 PRG 内任意 $k,\ell$，

$$
\begin{aligned}
\mathbb E[g_kg_\ell^*]
&=
\sum_{n,n'}V^{\rm tx}_{n,m_b}
\left(V^{\rm tx}_{n',m_b}\right)^*
\mathbb E[h_{k,n}h_{\ell,n'}^*]\\
&=
\left\|\mathbf V^{\rm tx}_{:,m_b}\right\|_2^2
R_{{\rm phy},k\ell}\\
&=8R_{{\rm phy},k\ell}.
\end{aligned}
$$

等价地，如果底层 PDP 的抽头功率为 $p_q$，则 PRG 内等效 PDP 为 $8p_q$：抽头时延、相对功率和归一化 PDP 形状均不变，只增加功率因子 8。接收端可以使用底层物理 PDP 构造 $\mathbf R_{\rm phy}$，再整体乘 8；也可以保持归一化 $\mathbf R_{\rm phy}$，同时把 LMMSE 中的观测噪声方差除以 8，两种实现必须数值等价。

该推导只对同一个 PRG 内成立。不同 PRG 使用不同 DFT 向量时，跨 PRG 协方差还含两个向量的内积；本轮接收机按研究者要求在每个 PRG 内独立估计，不利用跨 PRG 协方差，也不跨越 PRG 边界插值。4-RB 基线中重复使用同一向量的非相邻 PRG 仍分别估计，不做非连续 PRG 联合估计。

### 10.4 两个 DMRS 的 PRG 内 LMMSE 实现

两个 DMRS symbol 索引仍为 `[2,7]`，DMRS comb 间隔仍为 6 个子载波。4-RB PRG 每个 DMRS symbol 含 8 个 pilot，6-RB PRG 每个 DMRS symbol 含 12 个 pilot。对 PRG $b$，记其 pilot 子载波集合为 $P_b$，需要估计的本 PRG data 子载波集合为 $D_b$。对单位模已知 DMRS 解调并除去 DMRS 符号后，第 $i$ 个 DMRS symbol 的 LS 观测为

$$
\mathbf z_b^{(i)}
=
\mathbf g_{P_b}+\mathbf n_b^{(i)},
\qquad i\in\{1,2\}.
$$

A100 为 static、零速度信道，两个 DMRS 上的真实频域信道相同；两个独立噪声观测先做算术平均：

$$
\bar{\mathbf z}_b
=
\frac{\mathbf z_b^{(1)}+\mathbf z_b^{(2)}}{2}
=
\mathbf g_{P_b}+\bar{\mathbf n}_b.
$$

现有链路中单个 DMRS LS 观测的复噪声方差为

$$
\sigma_{\rm LS,one}^2
=
\frac{8}{10^{\mathrm{SNR}_{\rm dB}/10}},
$$

所以平均后进入 LMMSE 的噪声方差必须为

$$
\sigma_{\rm LS,avg}^2
=
\frac{\sigma_{\rm LS,one}^2}{2}
=
\frac{4}{10^{\mathrm{SNR}_{\rm dB}/10}}.
$$

从底层 A100 PDP 计算本 PRG 精确子载波坐标上的 $\mathbf R_{{\rm phy},P_bP_b}$ 与 $\mathbf R_{{\rm phy},D_bP_b}$，则 matched PRG-LMMSE 为

$$
\widehat{\mathbf g}_{D_b}
=
8\mathbf R_{{\rm phy},D_bP_b}
\left(
8\mathbf R_{{\rm phy},P_bP_b}
+
\sigma_{\rm LS,avg}^2\mathbf I
\right)^{-1}
\bar{\mathbf z}_b.
$$

实现使用 Cholesky 线性求解，不得显式求逆；数值 jitter 或 diagonal loading 必须写入输出。若对 $8\mathbf R_{\rm phy}$ 形式加入 loading，则改写为归一化 $\mathbf R_{\rm phy}$ 形式时 loading 也必须除以 8；以下两式是不含数值 loading 的物理定义。等价实现

$$
\widehat{\mathbf g}_{D_b}
=
\mathbf R_{{\rm phy},D_bP_b}
\left(
\mathbf R_{{\rm phy},P_bP_b}
+
\frac{\sigma_{\rm LS,avg}^2}{8}\mathbf I
\right)^{-1}
\bar{\mathbf z}_b
$$

必须在单元测试中与上式一致。每个 PRG 独立求得 $\widehat{\mathbf g}_{D_b}$ 后，按 active-band 原位置拼接；不得把某个 PRG 的 pilot 用于另一个 PRG。由于信道 static，同一个频率估计供所有 data-bearing OFDM symbol 使用。estimated-CSI 均衡、LLR 与 LDPC 口径保持第 2 节不变，LLR 噪声方差不额外加入 CE 误差项；ideal-CSI 分支直接使用相同 trial 的真实 $\mathbf g_D$。

实现还必须用 LMMSE 后验误差协方差做解析校验。对 PRG $b$，

$$
\mathbf C_{e,b}
=
8\mathbf R_{{\rm phy},D_bD_b}
-
8\mathbf R_{{\rm phy},D_bP_b}
\left(
8\mathbf R_{{\rm phy},P_bP_b}
+
\sigma_{\rm LS,avg}^2\mathbf I
\right)^{-1}
8\mathbf R_{{\rm phy},P_bD_b}.
$$

解析 trace NMSE 与大样本 Monte Carlo NMSE 应在预先固定的统计容差内一致；容差在 smoke 前写入测试，不得看到正式结果后调整。

### 10.5 信道估计误差与 CE NMSE 统计口径

CE 误差只评价 estimated-CSI 分支，评价集合是当前 grid 的全部 5568 个 data RE，不包括 192 个 DMRS RE。对 trial $t$，记 data RE 集合为 $D$；static 信道的频率估计按实际 data RE 位置展开后，先在线性域计算该 trial 的误差能量、真实信道能量和比值：

$$
E_t
=
\sum_{(s,k)\in D}
\left|
\widehat g_{t,s,k}-g_{t,s,k}
\right|^2,
$$

$$
S_t
=
\sum_{(s,k)\in D}
\left|g_{t,s,k}\right|^2,
$$

$$
\mathrm{NMSE}_t
=
\frac{E_t}{\max(S_t,10^{-30})}.
$$

这里 $s$ 是 OFDM symbol 索引，$k$ 是 active subcarrier 索引；因此各 PRG 自动按其实际 data RE 数参与同一个 trial 的总误差和总信号能量。主曲线不得先对每个 PRG 的 NMSE 做等权平均，也不得用 pilot RE 误差替代 data RE 误差。可额外输出每 PRG 的 $E_{t,b}/S_{t,b}$ 作为边界诊断，但它不是 result-028 主 CE 曲线的统计量。

对一个 `baseline_id + SNR` 点的 $N$ 个 trial，先对 trial 级线性 NMSE 做算术平均，再转换为 dB：

$$
\overline{\mathrm{NMSE}}
=
\frac{1}{N}\sum_{t=1}^{N}\mathrm{NMSE}_t,
\qquad
\mathrm{CE\ NMSE}_{\rm dB}
=
10\log_{10}\overline{\mathrm{NMSE}}.
$$

禁止采用以下两种不同口径替代：$10\log_{10}(\mathrm{NMSE}_t)$ 的 trial 平均，以及 $\sum_tE_t/\sum_tS_t$。输出至少保存 `ce_error_energy`、`ce_true_energy`、`ce_nmse_trial` 的逐 trial 可核验值或等价无损统计、`ce_nmse_sum`、`ce_nmse_mean` 和 `ce_nmse_mean_db`。追加 trial 时只在线性域累加 `ce_nmse_sum` 后除以总 trial 数，最后一次性转 dB。

该口径与 result-028 第 3.6 节现有五组 delay set 完全一致：现有实现也是先在每个 trial 的全部 data RE 上计算 $E_t/S_t$，再在线性域跨 trial 求平均并转换为 dB。新增透明基线不得改用 PRG 等权平均、总能量比或 dB 域平均，否则不能与现有 CE NMSE 曲线同图比较。

### 10.6 SNR 网格、trial、配对和统计判据

result-028 第 3.6 节当前五条曲线在 16 dB 以内所用 SNR 点的去重并集冻结为

$$
[14,14.25,14.5,14.75,15,15.5,15.75,16]\ \mathrm{dB}.
$$

两条新增透明基线的 ideal-CSI BLER、estimated-CSI BLER 和 CE NMSE 都必须覆盖以上全部 8 个点；不得因为某点零错误而删除该点，也不得新增 15.25 dB 或对单条基线后验加点。BLER 为零时 CSV 保留 0 和错误计数，图上仍以 `0.5/trials` 显示。

正式预算为：batch size 20；每个 SNR 至少运行 10,000 个共同 paired trials；之后两条基线继续使用相同 absolute trial 区间追加，直到两条基线各自的 ideal 与 estimated 四个接收机点均累计至少 200 个错误块，或共同达到 50,000 trials。达到 50,000 后任一点仍少于 30 个错误块，则该点标记“样本不足”，保留 Wilson 95%区间和零错误绘图下界，不外推。2026-08-05 在没有任何正式 trial 区间落盘前发现 batch size 100 会造成内存/解码吞吐退化；因此在正式数据生成前把纯计算分块改为 20。运行环境还对单条命令施加约 9 分钟外层时限；4,000/2,000-trial 区间在负载波动下仍过于接近该限制，因此在没有正式区间落盘前，把每次可恢复落盘区间最终冻结为最多 1,000 trials，并通过重复调用继续绝对 trial 区间。正式后台运行已稳定落盘 6 个区间后，将每个 Python 进程可连续处理的 SNR task 数由 1 调为 8；两个已排队的单任务进程结束后，该调度从第 9 个区间起生效，以复用同一进程中的 TensorFlow、信道和译码器初始化。每个 task 仍单独落盘不超过 1,000 trials。以上变更都只影响计算分块和调度，不改变 absolute trial 集合、物理配置、随机策略或总统计预算。

初始运行结束后必须执行曲线稳定性审计，而不是只检查总 trial 数：

1. 报告每点错误块数和 Wilson 95%区间；少于 200 个错误块的点显式标记；
2. 检查相邻 SNR 点是否出现高 SNR 经验 BLER 高于低 SNR 经验 BLER 的反向；禁止用单调拟合、删点或平滑隐藏反向；
3. 对视觉上控制曲线走向的点，报告 Wilson 区间宽度和相邻点区间重叠；若反向点任一侧少于 200 个错误块，优先追加这两个点；
4. CE NMSE 同时报告 trial 级线性 NMSE 的均值、标准误和 95% Monte Carlo 区间，用于区分真实趋势与有限 trial 波动。

如果研究者审图后仍认为曲线波动过大，允许在不重跑旧 trial 的前提下继续追加：新配置保持物理参数、SNR 网格、seed 派生和 baseline ID 不变，只提高 `target_errors`、`max_total_trials` 或指定需要补样的 SNR 点；例如下一轮可提高到 400 个错误块、最多 100,000 paired trials。被指定的 SNR 点仍对两条基线和 ideal/estimated 接收机追加相同 absolute trial 区间，以保留成对比较。每次追加必须以已有 CSV 为 `base_csv`，从 `previous_trial_end + 1` 开始连续运行；旧、新 trial 的错误数直接相加，CE NMSE 通过 `ce_nmse_sum` 在线性域合并。任一 absolute trial 区间重叠、缺口、seed 规则变化或 PRG 尾部映射不可重放时停止合并。

同一 SNR、absolute trial index 下，两条新增基线共享 payload、底层 A100 TDL 信道和原始噪声样本；同一基线的 ideal 与 estimated 接收机也使用同一 trial。6-RB 映射始终固定；4-RB 尾部映射按第 10.2 节在每个 absolute trial 重抽一次，并在该 trial 内固定。seed 仍为 `20260727`：共享 payload、信道和噪声的 seed 只包含场景、SNR、absolute trial index 与用途标签，不包含 baseline ID；只有 baseline 专属随机量（本轮即 4-RB 尾部映射）的 seed 才包含 baseline ID。增量运行只能继续未使用的绝对 trial 区间。

每个 BLER 点保存 trial 数、错误块数、BLER 和 Wilson 95%区间。第 10 节只要求同图曲线，不做 10%/1% logistic 目标拟合，也不据此改变 plan-028 原候选的排序或验收结论。

### 10.7 实现、测试、输出与 result-028 修改

优先扩展固定曲线入口，不复制链路主循环：

- 在 `cdd_lls/phy/precoding.py` 增加可复用的 8 Tx 单层 DFT codebook 和显式 PRG-to-vector 映射，保持旧 2/4 Tx 行为兼容；
- 在 `cdd_lls/phy/estimators.py` 把现有 PRG RMMSE 扩展为由 `prg_size_rb` 驱动的 4-RB/6-RB 通用实现，显式处理等效协方差功率因子；
- 扩展 `tools/run_bler_curves.py` 的配置 schema，使其可声明透明 PRG baseline、冻结 SNR schedule 和 paired ideal/estimated 运行；新增稳定 YAML 配置，不新建另一份链路主循环；
- 扩展 `tools/plot_result028_a100_subset.py`，把两条透明基线加入现有三张图，并为其分配与原五条曲线可区分且跨三图一致的 color、linestyle 和 marker；原五条曲线样式不变。

测试至少覆盖：8 点 DFT 列正交和每列功率 8；4-RB/6-RB PRG 边界；6-RB 固定映射；4-RB 尾部逐 absolute trial 重抽、trial 内固定、四元素无重复、同 seed 可重放且不同 trial 实际出现多种映射；DMRS 与 data 使用相同向量；PRG 内 Monte Carlo 协方差等于 $8\mathbf R_{\rm phy}$；scaled-covariance 与 scaled-noise 两个 LMMSE 公式等价；每个 PRG 只消费自己的 8 或 12 个 averaged pilot；两个 DMRS 平均后的噪声方差减半；解析与 Monte Carlo NMSE 一致；trial 级线性 NMSE 的合并顺序；8 点 SNR 网格和 paired absolute trial 区间无重叠；从 `base_csv` 追加后 BLER、Wilson 区间与 CE NMSE 等于一次性运行相同 absolute trial 并集的结果。

新增原始产物写入 `outputs/experiment028_csi_curves/20260803_main/transparent_prg_baselines/a100/`，至少包含展开配置、PRG 随机策略 manifest/hash、逐 trial 4-RB 尾部映射或可无损重放的等价记录、逐点 CSV、逐 trial error flags、CE 误差统计、每轮绝对 trial 区间、日志和样式表。最终把两条基线并入 result-028 第 3.6 节现有三张同图：estimated-CSI BLER、ideal-CSI BLER 和 CE NMSE；同步更新 `research/result-028-comb6三类CSI-TDL-A300ns.md`、`research/result-028-comb6三类CSI-TDL-A300ns-text.md`、对应图 CSV 和第 3.6 节文字说明。不更新 `KNOWLEDGE.md` 或 `GOALS.md`，直至研究者确认新增 result。

## 11. A100 透明 CDD：UE 使用底层物理信道协方差

### 11.1 目的、范围与稳定标识

本节在第 3.6 节已有五组 CDD delay set 上增加接收机失配基线。发射端、底层信道、DMRS、data、噪声归一化、MCS 和 SNR 采样均不改变；唯一变化是 UE 不知道 CDD delay set，信道估计器不使用 CDD 造成的等效频域协方差形状，而只使用底层 static TDL-A 100 ns 物理信道协方差及已知的 8 路总接收功率。这里“透明 CDD”是接收机知识口径，不是新的发射预编码方案。

五条新增曲线与原五组 delay set 一一对应，稳定标识为：

| result-028 标签 | 新 candidate ID | 对应现有发射 candidate |
|---|---|---|
| delay set 1, transparent CDD | `A100_B0_QC_TRANSPARENT_CDD` | `A100_B0_QC` |
| delay set 2, transparent CDD | `A100_AP_RMS_T1_TRANSPARENT_CDD` | `A100_AP_RMS_T1` |
| delay set 3, transparent CDD | `A100_AP_TU_NT_TRANSPARENT_CDD` | `A100_AP_TU_NT` |
| delay set 4, transparent CDD | `A100_S0_SIDON_TRANSPARENT_CDD` | `A100_S0_SIDON` |
| delay set 5, transparent CDD | `A100_MEFF_T2_06_TRANSPARENT_CDD` | `A100_MEFF_T2_06` |

本节只新增 estimated-CSI BLER 和 CE NMSE，不运行、不输出透明 CDD 的 ideal-CSI BLER。原五组 matched receiver 数据只读复用，不重跑旧 trial，不改变原 candidate ID、样式或结论。

### 11.2 发射端与真实等效信道

对每个现有 delay set，继续使用原 CDD 发射矩阵

$$
V_{k,n}=\exp\!\left(-j2\pi k\frac{j_n}{576}\right),\qquad n=0,\ldots,7,
$$

其中 $j_n$ 为该 candidate 在 source manifest 中冻结的 `delay_grid_coordinates`。8 根天线各系数幅度均为 1，不做 $1/\sqrt 8$ 归一化，所以单层总发射功率和等效信道平均功率因子均为 8；噪声方差继续使用

$$
\sigma_w^2=\frac{8}{10^{SNR_{dB}/10}}.
$$

对底层第 $n$ 根天线频域信道 $h_{k,n}$，实际用于 DMRS、data、信道估计误差真值和均衡的等效信道均为

$$
g_k=\sum_{n=0}^{7}V_{k,n}h_{k,n}.
$$

因此每条曲线的真实协方差 $\mathbf R_t$ 必须由该条 CDD delay set 与底层 TDL-A PDP 共同构造，等价于现有 matched receiver 使用的 candidate-specific 等效协方差。不得为了“透明”而改变发射矩阵或用底层物理信道直接充当误差真值。

### 11.3 UE 假设协方差与失配 LMMSE

令 $\mathbf R_{\rm phy}$ 为一根发射天线在全部 576 个 active subcarrier 上、由 static TDL-A 100 ns 真实 PDP 构造的频域协方差，且其对角线功率归一为 1。UE 不知道任何 $j_n$，因此不得使用 CDD 相位矩阵修正协方差形状；但 UE 知道 8 根独立同分布发射分支的总接收功率，所以估计器采用

$$
\mathbf R_a=8\mathbf R_{\rm phy}.
$$

采用 $8\mathbf R_{\rm phy}$ 而不是未缩放的 $\mathbf R_{\rm phy}$，是为了只研究“未知 CDD delay 导致的协方差形状失配”，不额外混入约 9 dB 的总功率标定错误。该总功率知识不泄露任何 delay set 信息。

两个 static DMRS symbol 的 LS 观测先做复数域平均。由于两个 DMRS 的等效信道相同，平均后 pilot 噪声方差为

$$
\sigma_{\rm LS}^2=\frac{\sigma_w^2}{2}=\frac{4}{10^{SNR_{dB}/10}}.
$$

令 $P$、$D$ 分别表示 pilot 与 data 子载波局部索引，$\mathbf z_P$ 为 averaged LS pilot。透明 CDD 接收机对五条曲线均使用同一个、与 delay set 无关的全带失配 LMMSE 滤波器

$$
\mathbf W_a
=\mathbf R_{a,DP}
\left(\mathbf R_{a,PP}+\sigma_{\rm LS}^2\mathbf I\right)^{-1},
\qquad
\widehat{\mathbf g}_D=\mathbf W_a\mathbf z_P.
$$

与第 10 节的 PRG 内估计不同，CDD 的相位随 active subcarrier 连续变化，本节在全部 576 个 active subcarrier 上一次性做频域 LMMSE，不切 PRG。估计后的单层 data RE 继续使用现有 estimated-CSI MRC/标量均衡、Max-Log 解调与 LDPC 解码链路。

对于实现自检，给定某条发射 candidate 的真实协方差 $\mathbf R_t$，失配估计器在 data 子载波上的理论误差协方差为

$$
\begin{aligned}
\mathbf C_e={}&\mathbf R_{t,DD}
-\mathbf W_a\mathbf R_{t,PD}
-\mathbf R_{t,DP}\mathbf W_a^H\\
&+\mathbf W_a\left(\mathbf R_{t,PP}+\sigma_{\rm LS}^2\mathbf I\right)\mathbf W_a^H,
\end{aligned}
$$

相应的解析参考值为

$$
NMSE_{\rm analytic}=\frac{\operatorname{tr}(\mathbf C_e)}
{\operatorname{tr}(\mathbf R_{t,DD})}.
$$

解析值只用于验证滤波器、索引、噪声方差和失配协方差实现，不替代 Monte Carlo 曲线。

### 11.4 CE NMSE 的逐 trial 统计口径

对 absolute trial $t$，只在全部 data RE 上计算

$$
E_t=\sum_{k\in\mathcal D}|\widehat g_k^{(t)}-g_k^{(t)}|^2,
\qquad
S_t=\sum_{k\in\mathcal D}|g_k^{(t)}|^2,
\qquad
r_t=E_t/S_t.
$$

每个 SNR、每条 transparent-CDD candidate 的正式 CE 点定义为

$$
\overline r=\frac{1}{T}\sum_{t=1}^{T}r_t,
\qquad
NMSE_{dB}=10\log_{10}(\overline r).
$$

必须逐区间保存 $E_t$、$S_t$、$r_t$，并累计 `sum(r_t)`、`sum(r_t^2)`、样本标准差、标准误和线性域 95%区间。追加 trial 时按 absolute-trial 区间合并这些充分统计量；不得平均分段 dB 值，不得改成 $\sum E_t/\sum S_t$。该定义与第 3.6 节现有 matched CDD 和透明 PRG 曲线完全一致。

### 11.5 SNR、配对、停止规则与可追加性

estimated-CSI BLER 和 CE NMSE 固定使用第 3.6 节相同的 8 点 SNR 网格：

$$
[14,14.25,14.5,14.75,15,15.5,15.75,16]\ {\rm dB}.
$$

同一 SNR 和 absolute trial 下，五条 transparent-CDD candidate 必须共享 payload、底层 8 天线物理信道 realization、averaged-LS noise 和 data noise；仅 CDD 发射矩阵及其导致的真实等效信道不同。五条曲线使用共同的连续 absolute-trial 区间，避免因独立停止造成比较噪声。

正式停止规则为：每个 SNR 至少 10,000 个共同 trials；之后只有当五条 estimated-CSI BLER 点均累计至少 200 个误块时才停止；任何一点不足 200 个误块则继续共同追加，单点上限 50,000 trials。可恢复区间长度不超过 1,000 trials，内部 batch size 为 20。每个区间完成即原子落盘 error flags、CE trial arrays 和 supplemental CSV；中断后从已验证的最大连续 absolute trial 继续，不重跑 trial 1，不与旧区间相加两次。研究者认为曲线仍抖时，可提高目标误块数或最小 trial 数，从现有 `base_csv`/supplemental 区间继续追加，而无需重新全跑。

每个 BLER 点保存 trial 数、误块数、BLER 和 Wilson 95%区间；最终执行相邻 SNR 单调性审计。若出现肉眼可见的反向波动，先追加共同 trials 并复核误块区间，而不是平滑或手工修改曲线。

### 11.6 实现、测试、输出与 result-028 修改

继续扩展 `tools/run_bler_curves.py`，增加独立的 `transparent_cdd_physical_covariance` scene mode；复用 source manifest、CDD precoder、信道生成、解调解码与 CSV 合并逻辑，不复制链路主循环。新增正式及 smoke YAML 配置，并提供可重复调用的控制脚本，直至所有 SNR 满足共同停止规则。

测试至少覆盖：五个新 ID 与 source candidate 一一映射；发射 precoder 与原 candidate 完全相同；真实协方差使用 candidate-specific $\mathbf R_t$；五条估计器滤波器均严格使用同一个 $8\mathbf R_{\rm phy}$；UE 配置和落盘 manifest 不包含 CDD delay-aware filter；两个 DMRS 平均后的噪声方差减半；理论失配误差公式与直接矩阵/Monte Carlo 结果一致；CE 逐 trial 统计和追加合并正确；五条曲线的 absolute-trial 区间相同、连续且无重叠；8 点 SNR 网格冻结；该 mode 不创建 ideal receiver 数据。

新增正式原始产物写入 `outputs/experiment028_csi_curves/20260803_main/transparent_cdd_physical_covariance/a100/`，至少包含展开配置与 hash、透明 CDD manifest/hash、source candidate 映射、滤波器及理论失配 NMSE 诊断、逐点 CSV、逐 trial error flags、逐 trial CE arrays、区间清单、稳定性审计、日志和样式表。

最终将五条 transparent-CDD 曲线加入 result-028 第 3.6 节现有 estimated-CSI BLER 和 CE NMSE 同图，使用与各自 source delay set 相同的颜色和统一的透明 CDD 线型/标记形成配对；ideal-CSI 图和数据保持不变。同步更新 `research/result-028-comb6三类CSI-TDL-A300ns.md`、`research/result-028-comb6三类CSI-TDL-A300ns-text.md`、精确图数据 CSV、证据 hash 和复现命令。不更新 `KNOWLEDGE.md` 或 `GOALS.md`，直至研究者确认新增 result。

## 12. A100 小时延透明 CDD：CE NMSE 低于 -15 dB 且 estimated BLER 闭合

### 12.1 目标与范围

本节在第 11 节五组强失配透明 CDD 之外增加一条新的小时延 CDD 发射 candidate。目标不是最大化人工频率分集，而是在 UE 仍不知道 CDD delay、仍只使用 $\mathbf R_a=8\mathbf R_{\rm phy}$ 的条件下，同时满足：

1. 正式采样网格上每个 CE NMSE 点均低于 `-15 dB`；
2. estimated-CSI BLER 在正式采样范围内同时形成 10%和 1%目标的双侧 bracket；
3. 发射 delay 不全为零，保留可复现但较弱的人工频率分集；
4. 正式结果继续支持 absolute-trial 追加，不因曲线波动而从 trial 1 重跑。

本节只新增 estimated-CSI BLER 和 CE NMSE，不要求 ideal-CSI BLER。系统、物理 TDL-A 100 ns、48 PRB、comb-6、8 Tx / 1 Rx、MCS、两个 static DMRS、功率和噪声归一化均与第 11 节一致。

### 12.2 候选选择与冻结发射定义

先在不运行译码链路的解析阶段扫描等间隔 delay family

$$
\mathbf j(\Delta q)=[0,\Delta q,2\Delta q,\ldots,7\Delta q].
$$

对每个 $\Delta q$，真实等效协方差仍由该 CDD 发射矩阵和 TDL-A 100 ns PDP 构造，接收机固定使用 $8\mathbf R_{\rm phy}$，按第 11.3 节失配误差协方差计算 14/16 dB 的解析 trace-NMSE。选择规则为：在 14 dB 解析 NMSE 不高于 `-16.5 dB` 的候选中，选择最大的 $\Delta q$，为 Monte Carlo 波动和实现差异相对 `-15 dB` 留至少 1.5 dB 设计余量，同时尽可能保留人工时延跨度。

2026-08-10 的冻结解析扫描中，$\Delta q=0.5$ 在 14/16 dB 分别为 `-16.925/-17.462 dB`，而 $\Delta q=0.75$ 已为 `-14.385/-14.892 dB`，不满足门限；因此先把 $\Delta q=0.5$ 送入链路 prescan。该候选在 13--19 dB 的 200-trial prescan 中 BLER 仅从 `0.49` 降至 `0.075`，未闭合 1%目标，故在正式 trial 1 前明确淘汰，不进入正式结果。

随后比较 $\Delta q=0.25$ 与 $\Delta q=0.1$。两者均满足 CE 余量和双目标 bracket；为在满足估计性能的前提下保留尽可能大的非零人工时延跨度，正式冻结较大的 $\Delta q=0.25$。其 15--22 dB、每点 200 trials 的 prescan 中，Monte Carlo CE NMSE 为 `-19.61` 至 `-22.92 dB`；BLER 在 16.0/16.5 dB 分别为 `0.125/0.080`，在 19.5/20.0 dB 分别为 `0.015/0.005`，因而同时观察到 10%和 1%的双侧 bracket。新增稳定 candidate 定义为：

- candidate ID：`A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD`；
- 图例：`small-delay transparent CDD`；
- `delay_grid_coordinates=[0,0.25,0.5,0.75,1,1.25,1.5,1.75]`；
- 物理人工 delay：`[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273] ns`；
- 发射系数：$V_{k,n}=\exp(-j2\pi k j_n/576)$，幅度 1，不做 $1/\sqrt8$ 归一化；
- UE 假设协方差：$8\mathbf R_{\rm phy}$，`knows_cdd_delays=false`。

### 12.3 BLER prescan 与正式 SNR 冻结规则

正式运行前先做不进入正式累计的短 prescan。prescan 使用与正式链路相同的 candidate、接收机、MCS、seed derivation 和 SNR 定义，初始扫描 13–19 dB、0.5 dB 间隔，每点 200 trials；若未同时观察到 BLER 高于 10%、低于 10%、高于 1%和低于 1%的点，则只向缺失方向按 0.5 dB 扩展。prescan 只定位门限，不与正式误块数或 CE 统计相加。

根据上述 prescan，在正式 trial 1 前冻结正式 SNR 网格为

$$
[15.5,15.75,16,16.25,16.5,17,17.5,18,18.5,19,19.25,19.5,19.75,20,20.25,20.5,21]\ {\rm dB}.
$$

其中 16 dB 附近和 19.5--20 dB 附近采用 0.25 dB 间隔，分别保证 10%和 1%交叉附近的采样分辨率；21 dB 作为 1%以下的额外闭合保护点。正式 YAML、plan 本节和落盘 manifest 必须使用完全相同的点集；正式运行后不得按结果删除或移动 SNR 点。prescan 只定位候选和门限，不与正式误块数或 CE 统计相加。

### 12.4 CE、BLER 统计与正式停止规则

CE NMSE 继续按第 11.4 节逐 trial 计算 $r_t=E_t/S_t$，在线性域平均后转 dB，保存 $E_t$、$S_t$、$r_t$、和、平方和、标准误及 95%区间。正式验收要求每个采样点的 Monte Carlo `ce_nmse_mean_db < -15 dB`；解析 trace-NMSE只作实现自检，不能替代该判据。

每个正式 SNR 至少 10,000 trials；达到最小预算后，只有当该 candidate 累计至少 200 个误块时才停止，否则继续追加到最多 50,000 trials。单个可恢复区间不超过 1,000 trials，内部 batch size 为 20。每个 BLER 点保存误块数、BLER 和 Wilson 95%区间。10%和 1%闭合均以实际采样点形成双侧 bracket 为准，不做样本范围外外推；曲线若出现肉眼可见反向，则追加相邻点的 absolute trials，不平滑数据。

### 12.5 实现、测试与交付

扩展 `transparent_cdd_physical_covariance` mode，使其除 approved source candidate 映射外，还能接收在配置和策略 manifest 中显式冻结的 8 维 `delay_grid_coordinates`；approved A100 source manifest 仍只提供物理系统定义和 hash，不被修改。新增 prescan/正式 YAML、可恢复控制脚本和候选审计。测试至少覆盖：显式 delay 的维度/有限值/单调性；实际 precoder delay 与 `[0,0.25,...,1.75]` 完全一致；五组旧透明 CDD 映射保持兼容；接收机滤波器不消费新增 delay；无 ideal 数据；SNR bracket、CE 门限、区间连续性和追加合并。

正式产物写入 `outputs/experiment028_csi_curves/20260803_main/transparent_cdd_small_delay/a100/`。最终把该 candidate 加入 result-028 第 3.6 节 estimated-CSI BLER 和 CE NMSE 图，ideal 图保持不变；同步更新 paired result、精确图 CSV、解析/Monte Carlo CE 对照、BLER bracket 审计、hash 和复现命令。不更新 `KNOWLEDGE.md`、`GOALS.md` 或 Git checkpoint，直至研究者确认。

## 13. 小时延 CDD 的 ideal 分集筛选与透明接收机联合验收

### 13.1 目标、事实边界与参考曲线

本节补充多组小时延 CDD 的 ideal-CSI BLER，目标是找到一条同时满足以下条件的新 candidate：

1. ideal-CSI BLER 相对第 10 节 transparent PRG precoder cycling 不构成可见劣化；
2. UE 不知道 CDD delay、只用 $8\mathbf R_{\rm phy}$ 时，CE NMSE 在正式 estimated 网格上全部低于 `-15 dB`；
3. estimated-CSI BLER 在实际采样点上同时闭合 10%和 1%，且曲线不因误块不足出现可见反向；
4. ideal 与 estimated 的原始 trial 均可按 absolute-trial 区间继续追加。

“分集性能”只用 ideal-CSI BLER 衡量，不用 delay span、协方差秩或 CE NMSE 代替。参考数据固定为第 10 节 `A100_PRG_DFT8_4RB` 和 `A100_PRG_DFT8_6RB` 的正式 ideal 曲线；共同比较网格为 `[14,14.25,14.5,14.75,15,15.5,15.75,16] dB`。本节不重跑或修改已有 PRG 点。

### 13.2 冻结候选池与解析 CE 预筛

候选均显式冻结 8 维 `delay_grid_coordinates`，发射系数仍为 $V_{k,n}=\exp(-j2\pi k j_n/576)$，幅度 1、总功率 8。首轮候选池为：

| ID 后缀 | delay-grid coordinates | 设计目的 |
|---|---|---|
| `QSTEP0P25` | `[0,0.25,0.5,0.75,1,1.25,1.5,1.75]` | 第 12 节已完成的 CE/estimated 参考 |
| `QSTEP0P30` | `[0,0.3,0.6,0.9,1.2,1.5,1.8,2.1]` | 均匀增大人工时延 |
| `QSTEP0P35` | `[0,0.35,0.7,1.05,1.4,1.75,2.1,2.45]` | 均匀 trade-off 中点 |
| `QSTEP0P40` | `[0,0.4,0.8,1.2,1.6,2,2.4,2.8]` | CE 门限内较大的均匀 span |
| `G34MAX2P1` | `[0,0.061765,0.247059,0.617647,1.111765,1.420588,1.544118,2.1]` | 以 34 刻度 Golomb 型差集提高 pairwise-delay 多样性 |
| `SIDON65MAX2P8` | `[0,0.043077,0.129231,0.301538,0.516923,0.861538,1.292308,2.8]` | 非均匀长尾，在限制 RMS delay 时保留较大最大 span |

首轮 6 条候选的 300-trial ideal 公共网格 prescan 全部未满足第 13.3 节的 PRG 非劣点估计门限，因此在任何正式 trial 1 前冻结第二轮零均值循环移位候选。CDD ideal 接收机只依赖各分支 delay 的相对差；给全部分支减去同一个常数只给每个子载波的等效信道乘公共相位，不改变 ideal-CSI BLER。透明接收机却会受到该公共线性相位的协方差失配影响。因此第二轮去掉不贡献分集的公共 delay 偏置，使用以下关于零对称的循环 delay/advance；负坐标表示相对于共同 OFDM 定时参考的循环 advance，不解释为负物理传播时延：

| ID 后缀 | delay-grid coordinates | 与非中心化等价的 ideal delay 差 |
|---|---|---|
| `CSTEP0P50` | `[-1.75,-1.25,-0.75,-0.25,0.25,0.75,1.25,1.75]` | `[0,0.5,...,3.5]` |
| `CSTEP0P75` | `[-2.625,-1.875,-1.125,-0.375,0.375,1.125,1.875,2.625]` | `[0,0.75,...,5.25]` |
| `CSTEP1P00` | `[-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5]` | `[0,1,...,7]` |
| `CSTEP1P25` | `[-4.375,-3.125,-1.875,-0.625,0.625,1.875,3.125,4.375]` | `[0,1.25,...,8.75]` |
| `CSTEP1P50` | `[-5.25,-3.75,-2.25,-0.75,0.75,2.25,3.75,5.25]` | `[0,1.5,...,10.5]` |

第二轮仍逐条执行 14/16/20 dB 的 `-16.5 dB` 解析 CE 预筛和相同的 ideal 公共网格 prescan；不得因中心化而跳过透明 estimated/CE 验收。

第二轮解析结果表明，零均值并不能使当前固定 DMRS/LMMSE 的预测相位参考保持不变：即使 `CSTEP0P50` 在 14 dB 也只有 `-12.214 dB`，五条均未通过 `-16.5 dB`，故不进入译码链路。第三轮在任何新增 link trial 前冻结为二维解析设计族

$$
\mathbf j(a,q)=a+q[0,1,\ldots,7],
$$

其中相对步长 $q$ 决定 ideal 分集，公共 offset $a$ 不改变 ideal-CSI 的每 RE 信道幅度或 BLER，但会改变固定物理协方差 LMMSE 的预测相位。先对 $q\in\{0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2\}$、$a\in[-4,2]$（0.25 栅格）计算 14/16/20 dB 解析 NMSE；每个 $q$ 只保留最小化三点中最差 NMSE 的 offset，并仅把三点全部不高于 `-16.5 dB` 的 `(a,q)` 冻结为新增 link candidate。该 offset 搜索不运行信道 realization、payload 或译码，不计为 link prescan。

冻结解析搜索结果为：`q=0.5,a=-1`、`q=0.6,a=-0.75`、`q=0.7,a=-0.75` 通过，14 dB 解析 NMSE 分别为 `-18.545/-17.426/-16.703 dB`；`q=0.8` 的最优 offset 在 14 dB 只有 `-15.893 dB`，其余更大步长更差，均淘汰。进入 ideal 公共网格 prescan 的三条新增 ID 和 delay 数组固定为：

- `A100_SMALL_CDD_OQ0P50_AM1P00_TRANSPARENT_CDD`：`[-1,-0.5,0,0.5,1,1.5,2,2.5]`；
- `A100_SMALL_CDD_OQ0P60_AM0P75_TRANSPARENT_CDD`：`[-0.75,-0.15,0.45,1.05,1.65,2.25,2.85,3.45]`；
- `A100_SMALL_CDD_OQ0P70_AM0P75_TRANSPARENT_CDD`：`[-0.75,-0.05,0.65,1.35,2.05,2.75,3.45,4.15]`。

三条各 300 trials 的公共网格 prescan 后，`q=0.5/0.6` 整体偏弱，`q=0.7` 与 PRG 4-RB 最接近但尾点仅有 7--11 个错误块，尚不能区分真实差异与抽样波动。因此只对 `q=0.7` 使用上述 300-trial CSV 为 base，从 absolute trial 301 追加到 3,000；不得重跑前 300。该 refine 仍属 prescan，只有 pointwise 差值和 crossing 审计后才决定是否进入 transparent estimated prescan。

`q=0.7 @ 14 dB` 追加到 3,000 trials 后为 `404/3000=0.13467`，相对 PRG 4-RB 的 `0.1157` 已确认偏弱，故均匀 family 不进入 estimated prescan。下一轮冻结为 deterministic 非均匀搜索：seed `20260727`，生成 2,000 组严格递增的 8 维相对 delay，首坐标为 0、总 span 在 `[4.9,10.5]`，7 个增量先独立采样 `Uniform[0.2,2.0]` 再按目标 span 等比缩放；每组在 `[-3,2]`、0.25 栅格解析选择最优公共 offset。只保留 14/16/20 dB 解析 NMSE 全部不高于 `-16.5 dB` 者，并按 A100 物理协方差加权的 $M_4^{eff}$ 从小到大排序。$M_4^{eff}$ 只作筛选，排名前 3 的精确 delay 数组必须在 link trial 1 前回写本节，最终非劣仍只由 ideal BLER 判定。

搜索共 75 组通过 CE 门限，link trial 1 前冻结 $M_4^{eff}$ 最小的三组为：

- `A100_SMALL_CDD_NU419_TRANSPARENT_CDD`：`[-1,-0.6021138567581448,0.0895996186153063,0.6239161297150366,1.5984488322422608,2.6801431206827058,3.826494734287368,4.322142395582422]`；
- `A100_SMALL_CDD_NU1100_TRANSPARENT_CDD`：`[-0.75,-0.29809293777231693,-0.05505748890107687,0.9387998633766634,1.873540817584035,3.0169058861647016,4.1646556275916815,4.391081388131776]`；
- `A100_SMALL_CDD_NU126_TRANSPARENT_CDD`：`[-1,-0.7414683241142626,-0.17562137895638907,0.7033701780596064,1.2344736058932755,2.3922171517916433,3.405211310731164,4.454550929516178]`。

三组在 14 dB 的解析 NMSE 分别为 `-16.617/-16.661/-16.697 dB`，对应 $M_4^{eff}$ 为 `25720.47/25885.49/25886.78`；只有完整 ideal BLER 公共网格通过者才进入 estimated prescan。

300-trial ideal 公共网格中，`NU1100` 的 10% crossing 约为 `14.20 dB`，相对 PRG 4-RB 约 `14.12 dB` 晚 `0.08 dB`，通过 crossing 预筛；高 SNR 点只有 10--11 个误块，pointwise 非劣仍待追加确认。故 `NU1100` 是唯一 provisional finalist，允许先运行 14--24 dB、0.5 dB 间隔的 300-trial transparent estimated prescan以确认 CE 和 10%/1%闭合；这不等同于正式选中，ideal 公共点仍须从 absolute trial 301 追加并通过正式非劣审计。

`NU1100` 的 transparent estimated prescan 在 21--24 dB 仍为约 `2.7%--4.7%`，虽 Monte Carlo CE 为约 `-16.3` 至 `-17.7 dB`，但未闭合 1%，因此淘汰且不再追加 ideal。为直接针对该失败机理，最后一轮 deterministic 非均匀搜索把解析 CE 门限收紧为 14/16/20 dB 全部不高于 `-19.5 dB`，seed 仍为 `20260727`，样本数 10,000，总 span 改为 `[1.75,5.5]`，其余生成、offset 优化和 $M_4^{eff}$ 排序规则不变。该门限来自已正式闭合的 q=0.25 候选解析/Monte Carlo CE 量级，只作更保守的闭合预筛；排名前 3 仍必须跑完整 ideal 与 estimated 链路。

严格搜索共 3,428 组通过，link trial 1 前冻结 $M_4^{eff}$ 最小的三组：

- `A100_SMALL_CDD_STRICT512_TRANSPARENT_CDD`：`[-0.75,-0.5881588162481517,-0.3934827742815943,-0.12846160007940122,0.5645272524222211,1.0108750753299287,1.8030895119760455,2.7408216662551523]`；
- `A100_SMALL_CDD_STRICT2717_TRANSPARENT_CDD`：`[-0.75,-0.6286838856150396,-0.03690256624247046,0.0756489126928106,0.45856291237065805,0.9847526186824276,1.7239788978808308,2.827984167958107]`；
- `A100_SMALL_CDD_STRICT8263_TRANSPARENT_CDD`：`[-0.75,-0.6441578172090183,-0.4232914387586246,-0.20176604413630228,0.5936518750446993,1.1776850367200484,1.9863260227898762,2.4551325054894164]`。

三组 14 dB 解析 NMSE 为 `-19.536/-19.595/-19.578 dB`，$M_4^{eff}$ 为 `34780.96/35176.49/35195.21`。它们进入最后一轮 300-trial ideal 公共网格 prescan。

严格 CE 组在已完成的 14/14.25/14.5 dB 配对点上，最好候选 `STRICT512` 的 ideal BLER 为 `0.1267/0.1267/0.1000`，对应 PRG 4-RB 为 `0.1157/0.0861/0.0725`；后两点劣化远超 `0.005`，因此提前淘汰并停止剩余 prescan。至此本节冻结并执行的均匀、offset 优化、宽 span 非均匀和严格 CE 非均匀候选中，没有一条同时通过 ideal PRG 非劣与 transparent estimated 1%闭合：

- q=0.25 的正式 transparent estimated 曲线闭合 1%，但其 300-trial ideal 公共点明显弱于 PRG；
- q=0.7 在 14 dB refine 至 3,000 trials 后仍为 `0.13467`，弱于 PRG 4-RB 的 `0.1157`；
- `NU1100` 的 ideal 10% crossing 接近 PRG，但 transparent estimated 在 21--24 dB 保持约 `2.7%--4.7%`，形成不可接受的尾部地板；
- 把解析 CE 收紧到 `-19.5 dB` 后，最优非均匀候选的 ideal 公共点再次明显弱于 PRG。

因此本轮不冻结正式 candidate，不启动 10,000-trial 正式 ideal/estimated 运行，也不向 result-028 第 3.6 节加入一条不满足用户联合要求的曲线。该结论是当前已搜索候选空间的负结果，不证明任意透明 CDD 设计在数学上都不可能满足两项要求；若继续研究，必须新增不同发射结构或放宽“PRG pointwise 非劣”“UE 完全不知道 CDD”中的至少一项，并先形成新的 plan 节。

候选进入链路 prescan 前，先在 14、16、20 dB 计算第 11.3 节的解析失配 trace-NMSE；任一点高于 `-16.5 dB` 的候选淘汰，为正式 Monte Carlo `-15 dB` 门限保留至少 1.5 dB 设计余量。解析量只作预筛，不能替代正式 trial-ratio NMSE。

### 13.3 ideal/estimated prescan 与选择判据

对通过解析预筛的候选先运行 ideal-only 公共网格 prescan。所有候选使用相同 seed derivation、payload、物理信道和 data noise，每点 300 trials，固定覆盖参考 PRG 的 `[14,14.25,14.5,14.75,15,15.5,15.75,16] dB`。在该公共网格上已经明显违反下述点估计非劣门限的候选立即淘汰，不继续浪费高 SNR 译码预算。只有通过公共网格的候选才把 ideal prescan 向高 SNR 延伸到观察到 1%以下，并运行 14–22 dB、0.5 dB 间隔的 transparent estimated-only prescan；两种 receiver 在重叠点继续使用相同 absolute seed。estimated prescan 必须同时观察到 10%与 1%的上下侧采样点。prescan 不计入正式累计。

候选的 ideal 非劣判据冻结为：在八个 PRG 共同 SNR 点上，点估计相对两条 PRG 中较小 BLER 的最大绝对劣化不超过 `0.005`，且按独立 Wilson 区间构造的差值单侧 95%上界不超过 `0.01`；同时 log-BLER 线性插值得到的 10% crossing 不晚于较优 PRG 超过 `0.1 dB`。正式选择时先排除不满足 CE/estimated 闭合者，再在满足 ideal 非劣者中选择解析及 Monte Carlo CE NMSE 最低者；若首轮池没有合格者，只能在本节追加并冻结新 delay 数组后继续 prescan，不得事后放宽判据。

### 13.4 正式运行、停止规则与交付

在正式 trial 1 前，把 prescan 选中的唯一 candidate ID、delay 数组、ideal SNR 网格和 estimated SNR 网格回写本节及正式 YAML。ideal 正式网格包含八个 PRG 共同点，并扩展到 1%以下；estimated 正式网格在 10%与 1%附近不大于 0.25 dB，其余不大于 0.5 dB。ideal 与 estimated 分开落盘但在重叠 SNR 使用相同 absolute seed。

每个正式点至少 10,000 trials、目标至少 200 errors、最多 50,000 trials；可恢复区间不超过 1,000，batch size 20。ideal 执行 pointwise PRG 非劣审计；estimated 执行 `CE NMSE < -15 dB`、10%/1%双 bracket 和相邻反向审计。出现可见反向时只追加相关点，不平滑。逐 trial error flags、CE 的 $E_t/S_t$ 数组、Wilson 区间、解析/Monte Carlo CE 对照、candidate manifest/hash 和续跑日志全部保存。

代码继续扩展 `transparent_cdd_physical_covariance` mode，使显式 CDD candidate 可按配置单独运行 `ideal` 或 `estimated` receiver；旧五组 estimated-only 配置和第 12 节结果必须保持兼容。正式产物分别写入 `outputs/experiment028_csi_curves/20260803_main/transparent_cdd_small_delay_diversity_<ideal|estimated>/a100/`。最终把选中 candidate 的 ideal 曲线加入第 3.6 节 ideal 同图，并把 estimated 与 CE 加入对应同图；同步更新 paired result 和精确 CSV。不更新 `KNOWLEDGE.md`、`GOALS.md` 或 Git checkpoint，直至研究者确认。

## 14. A100 8Tx/4Rx 首批 ideal/estimated-CSI BLER 增补

> 状态：研究者于 2026-08-11 确认执行。首批只运行 delay set 2、3、4、transparent PRG 6 RB 和 small-delay transparent CDD；生成首批图并写入 result-028 新节后暂停，其他方案须等待研究者后续要求。

### 14.1 目标与首批候选

本增补把第 3.6 节的 A100、8Tx/1Rx 比较扩展到 8Tx/4Rx，回答接收 MRC 分集下五条指定发射方案的 ideal-CSI BLER、estimated-CSI BLER 和 data-RE CE NMSE 曲线如何变化。首批候选固定为：

| 图例 | candidate ID | 发射与 estimated 接收机 |
|---|---|---|
| delay set 2 | `A100_AP_RMS_T1` | 原 candidate 02 CDD；known-delay matched 全带 LMMSE |
| delay set 3 | `A100_AP_TU_NT` | 原 candidate 04 CDD；known-delay matched 全带 LMMSE |
| delay set 4 | `A100_S0_SIDON` | 原 candidate 07 CDD；known-delay matched 全带 LMMSE |
| transparent PRG 6 RB | `A100_PRG_DFT8_6RB` | 8 个 6-RB PRG 固定循环 DFT 向量；逐 PRG matched physical-covariance LMMSE |
| small-delay CDD | `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | `j=[0,0.25,...,1.75]`；UE 不知道 CDD delay，只用 $8\mathbf R_{\rm phy}$ 的全带失配 LMMSE |

五条方案都运行真实 ideal-CSI 链路和各自上述 estimated-CSI 链路。small-delay CDD 的 ideal 接收机只使用相同发射矩阵和真实 data-RE 等效信道，不使用 transparent 接收机假设。不得复用 1Rx BLER trial，也不得把第 13 节其他 small-delay 搜索候选加入本轮。

### 14.2 系统、公平性与 4Rx 定义

固定 48 PRB、576 个有效子载波、30 kHz SCS、8 Tx / 4 Rx / 1 layer、static Sionna TDL-A 100 ns、零速度、各 Tx/Rx 分支独立且同 PDP、DMRS comb-6、两个 static DMRS symbol、每 symbol 96 个 pilot、5568 data RE、16QAM、NR 256QAM MCS table index 8、码率 553/1024、Sionna LDPC 最多 8 次迭代。各方案的发射矩阵、总功率 8、DMRS/data 映射和协方差知识沿用第 3.6、10、12 节。

横轴 SNR 继续定义为每根 Rx 分支的平均接收 SNR。噪声方差固定为

$$
\sigma_n^2=\frac{8}{10^{\mathrm{SNR}_{dB}/10}},
$$

不得因 $N_r=4$ 再乘或除 4。因此 4Rx MRC 同时包含约 6.02 dB 的平均合并功率增益和接收分集，不能把横轴解释为四分支合计 SNR。ideal 与 estimated 均按 Rx 维执行 MRC；estimated 的 LLR 噪声方差仍不额外加入 CE 误差项。

同一 SNR 和 absolute trial index 下，五条方案尽可能共享 payload、底层 TDL 信道和原始噪声；同一方案的 ideal/estimated 使用相同 seed 派生和重叠 absolute trial 区间。base seed 固定为 `20260727`。所有新产物必须记录 `n_rx=4`，且输出目录不得读取或合并任何 1Rx `base_csv`。

### 14.3 CE NMSE 与解析自检

estimated trial $t$ 的主 CE 统计在四根 Rx 和全部 data RE 上先求能量比：

$$
r_t=
\frac{\sum_{r=1}^{4}\sum_{(s,k)\in D}|\widehat g_{t,r,s,k}-g_{t,r,s,k}|^2}
{\sum_{r=1}^{4}\sum_{(s,k)\in D}|g_{t,r,s,k}|^2}.
$$

一个 `candidate + SNR` 点先在线性域计算 $\overline r=N^{-1}\sum_t r_t$，再报告 $10\log_{10}\overline r$、标准误和 95% Monte Carlo 区间。必须保存可核验的 `ce_nmse_sum`、`ce_nmse_sumsq`；透明 PRG 和 transparent CDD 继续保存逐 trial error energy、true energy 与 ratio 数组。

独立同分布 Rx 下，LMMSE 解析 trace-NMSE 的分子和分母都随 $N_r$ 等比例增加，所以固定 SNR 的解析值应与 1Rx 相同；该性质只作实现校验。主图仍使用本轮 4Rx trial-ratio Monte Carlo 统计，不复制 1Rx CE 曲线。

### 14.4 prescan、正式网格与预算

代码和测试通过后，先在五条方案、两个 receiver 上使用公共粗网格 `[4,6,8,10,12,14,16,18,20] dB`，每点固定 300 trials。prescan 只用于寻找每条曲线的 10%和 1%双侧 bracket、发现 CE/BLER 地板以及确定正式采样范围，不写入正式累计。

prescan 后、任何正式 trial 1 前，冻结以下网格。三条 manifest CDD 的 ideal 网格为 `[4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25] dB`；estimated 网格为上述十点再加 `[8,12,16,20] dB`，其中高 SNR 四点用于确认 delay set 3 的地板，三条 candidate 仍共享同一网格以保持配对。transparent PRG 6 RB 的 ideal/estimated 共用上述十点。small-delay CDD 的 ideal 使用上述十点，estimated 使用 `[4,4.5,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,8] dB`，用于闭合其较晚的 1% crossing。不得只对结果有利的 candidate 删除困难点。

正式点至少 10,000 trials、目标至少 200 个错误块、最多 50,000 trials；单次可恢复 interval 不超过 1,000。4Rx throughput smoke 表明 batch size 100 会因单个约 44.9 GiB complex128 张量而 OOM，batch size 25 完成 300 trials 且比 20 更快，因此在正式 trial 1 前把执行 batch size 冻结为 25；该值只改变执行分批和 absolute seed 的批次起点，不改变物理定义、trial 数或配对规则。零错误点 CSV 保留真实 0，图上使用 `0.5/trials` 下界。出现相邻反向时检查 Wilson 区间并只追加相关点，不平滑原始数据。10%/1% crossing 只在形成双侧 bracket 时报告，不外推。

### 14.5 代码、测试、输出与停止条件

扩展 `tools/run_bler_curves.py` 的 `bler-curve-runner-v1` scene schema，增加正整数 `n_rx`，默认值保持 1 以兼容既有配置；普通 manifest CDD、transparent PRG 和 transparent CDD 三条路径都必须由该字段驱动信道与噪声张量。测试至少覆盖：默认仍为 1Rx；非法 `n_rx` 被拒绝；4Rx 信道/噪声/估计形状；ideal 与 estimated MRC 在 Rx 维求和；固定样本下 4Rx 手算均衡一致；4Rx CE trial ratio 跨全部 Rx 求和；展开配置和策略 manifest 明确保存 `n_rx=4`；旧 1Rx 配置 validate 不变。

prescan 与正式产物写入 `outputs/experiment028_csi_curves/20260811_4rx_subset/` 下互不混合的目录，保存原 YAML、展开配置及 hash、源/策略 manifest 及 hash、逐点 CSV、error flags、CE 统计、absolute-trial 区间、日志、共享样式和图对应精确 CSV。交付三张无标题首批同图：4Rx ideal-CSI BLER、4Rx estimated-CSI BLER、4Rx CE NMSE；图写入 `docs/figures/result-028/`，字体、线宽、marker 和约 13 cm 预览遵循 `docs/agent/RESULT_SPEC.md`。

最终在 `research/result-028-comb6三类CSI-TDL-A300ns.md` 与 `research/result-028-comb6三类CSI-TDL-A300ns-text.md` 新增独立的 4Rx 首批结果节，逐点或按预先说明的覆盖规则给出原始数值、trial/error、Wilson 区间、CE 区间、SNR 定义、hash、复现命令、异常和适用范围。完成五条首批曲线、证据审计和文档后停止，不运行 delay set 1/5、transparent PRG 4 RB 或其他 transparent CDD。结果未经研究者确认前不更新 `KNOWLEDGE.md`、`GOALS.md`，不创建 Git checkpoint。

### 14.6 2026-08-12 快速 estimated-CSI 三曲线预览

研究者要求暂停耗时的五方案完整正式队列，优先查看 delay set 4、transparent PRG 6 RB 和 small-delay CDD 的 4Rx estimated-CSI BLER 对比。已落盘的完整 1,000-trial 区间均保留，后台停止最多只丢弃尚未完成、未落盘的当前区间，后续仍可按 absolute trial 继续。

快速图复用 delay set 4 已完成的正式数据；PRG 6 RB 和 small-delay CDD 在 `[4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7] dB` 各运行固定 1,000 trials、batch size 25、seed `20260727`。PRG 当前入口保持 ideal/estimated 成对执行，但快速交付只画 estimated。该图是用于确定趋势和后续预算的快速诊断，不作为第 14.4 节 10,000-trial 正式验收；必须报告 Wilson 95%区间和不同曲线 trial 数不相同，不得据此声称小于约 1%的细微差异已确认。

### 14.7 服务器执行与紧凑证据回传

`run_bler_curves.py` 不设置 `tf.device`、`CUDA_VISIBLE_DEVICES` 或多 GPU distribution strategy；Sionna 信道生成和 LDPC 译码中的 TensorFlow 算子沿用 TensorFlow 自动设备放置。兼容 GPU 可见时，支持的算子可落到单 GPU；没有兼容 GPU 时自动使用 CPU。NumPy、Python 调度及 TensorFlow 到 NumPy 的数据回传仍在 CPU，故本节不把当前实现描述为端到端 GPU 并行，也不承诺未实测的加速比。服务器正式执行前须先用相同配置做小规模 smoke，记录 TensorFlow 可见设备、实际 placement probe、显存是否满足当前 batch size 和每 trial 用时；若因显存降低 batch size，trial 区间、seed、物理配置和停止规则不变。

服务器运行结束并完成 `--stage merge` 后，使用 `tools/export_bler_result_bundle.py` 把多个配置导出为一个 JSON 和一个纯文本摘要。JSON 必须包含展开配置、source/transparent manifest hash、代码 hash、Python/TensorFlow/Sionna 版本、GPU 可见性与 placement probe、逐点 trials/errors/BLER/Wilson 95%区间、CE NMSE 及其区间、absolute-trial 连续性审计、逐 SNR 耗时、关键产物 hash 和原始 NPY 的聚合 hash；TXT 提供可直接复制的逐点表。原始 NPY 不通过文本复制，但服务器端保留，出现审计问题时再按聚合 hash 定位并补传。

## 15. A100 8Tx/1Rx 第 3.6 节两条曲线的 SNR 延伸

> 状态：已完成正式运行、合并、审计和第 3.6 节更新，待研究者确认。只扩展第 3.6 节已有的小时延透明 CDD estimated-CSI/CE 和 transparent PRG 6 RB estimated-CSI/CE；不改变第 14 节 4Rx 工作，不重跑已有正式点。

### 15.1 目标、范围与事实源

小时延透明 CDD `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` 在原正式网格之外新增 `14/14.5/15 dB` 三个 estimated-CSI BLER 与 CE NMSE 点。发射 delay、UE 假设协方差、MCS、信道、seed 和统计口径与第 12 节完全相同；第 12 节已有 17 个正式点只读保留，新点从各自 absolute trial 1 开始，不把同 SNR prescan trial 并入正式结果，不生成 ideal-CSI 数据。

transparent PRG 6 RB `A100_PRG_DFT8_6RB` 复用第 10 节的固定 8 点 DFT codebook、6-RB PRG 边界、`[0,1,2,3,4,5,6,7]` 映射和 PRG 内 matched physical-covariance LMMSE，只扩展 estimated-CSI BLER 与 CE NMSE，目标是在原 `16 dB` 点以下继续形成 1% BLER 双侧 bracket。原 4-RB PRG 和 6-RB ideal-CSI 曲线均不扩展；原目录和原八点数据只读保留。

### 15.2 SNR 网格、预算与可追加性

小时延透明 CDD 的完整正式 estimated/CE 网格更新为

$$
[14,14.5,15,15.5,15.75,16,16.25,16.5,17,17.5,18,18.5,19,19.25,19.5,19.75,20,20.25,20.5,21]\ \mathrm{dB}.
$$

只计算新增前三点；每点至少 10,000 trials、目标至少 200 errors、最多 50,000 trials，单个可恢复区间不超过 1,000，batch size 20。由于新增点预计处于高 BLER 区，满足最小 trial 后即可停止，但仍必须保存 Wilson 95%区间、逐 trial CE 三数组和解析失配 trace-NMSE。

6-RB PRG 的新点先用不进入正式累计的 `[16.5,17,17.5,18] dB`、每点 200 trials 短 prescan 核对 1% crossing 位置。2026-08-12 的 prescan BLER 分别为 `0.03/0.015/0.005/0`，没有生成 ideal 数据；据此在任何正式 trial 1 前把正式网格冻结为

$$
[16.25,16.5,16.75,17,17.25,17.5,17.75]\ \mathrm{dB}.
$$

它与原 16 dB 点共同用于 1% bracket，17.75 dB 是低于目标的保护点；18 dB 只保留 prescan，不进入正式累计。正式点至少 10,000 trials、目标至少 200 errors、最多 50,000 trials，单个可恢复区间不超过 1,000，batch size 20；若 `17.75 dB` 仍未低于 1%，停止并记录未闭合，不在正式结果后临时改网格。新增输出使用相同 base seed `20260727` 和 absolute-trial seed derivation，独立目录中的每个新 `candidate+SNR` 从 trial 1 开始；后续只可从该点 `previous_trial_end+1` 追加。

### 15.3 实现、测试、输出与验收

固定曲线入口增加范围受限的第 15 节 transparent PRG 运行：只允许 1Rx、`A100_PRG_DFT8_6RB`、estimated receiver 和严格递增的新 SNR 网格；旧第 10/14 节配置行为不变。estimated-only 路径不得构造、译码或落盘 ideal receiver 数据。测试覆盖第 15 节候选/receiver约束、单 receiver paired-task 预算、无 ideal 输出和旧配置 validate 兼容。

小时延新增点继续写入 `outputs/experiment028_csi_curves/20260803_main/transparent_cdd_small_delay/a100/`；6-RB 新点写入独立的 `outputs/experiment028_csi_curves/20260803_main/transparent_prg_6rb_1pct_extension/a100/`。两处均保存展开配置、manifest/hash、逐点/逐 trial 数据、interval、日志和稳定性审计。最终将两处新增 estimated/CE 点与原第 3.6 节图数据只读合并，更新 estimated-CSI BLER、CE NMSE 两图及对应精确 CSV；ideal-CSI 图保持不变。两版 result 必须给出新增逐点 trials/errors/BLER/Wilson 95%区间、CE 95%区间、1% bracket、异常和复现命令。未经研究者确认，不更新 `KNOWLEDGE.md`、`GOALS.md`，不创建 Git checkpoint。

## 16. A100 小时延 CDD 的 ideal CSI 与 matched-covariance estimated CSI

> 状态：已完成正式运行、合并、审计和第 3.6 节三张同图更新，待研究者确认。发射端固定为第 12/15 节的 `j=[0,0.25,...,1.75]`，新增 ideal-CSI BLER，以及 UE 知道该 CDD 并使用真实等效信道协方差的 estimated-CSI BLER/CE NMSE；原透明接收机结果只读保留。

### 16.1 研究对象与接收机定义

发射端候选仍为

$$
\mathbf j=[0,0.25,0.5,0.75,1,1.25,1.5,1.75],
\qquad
V_{k,n}=\exp(-j2\pi k j_n/576).
$$

8 个系数幅度均为 1，总发射功率为 8；DMRS、data、48 PRB、comb-6、static TDL-A 100 ns、1Rx、16QAM MCS 8、噪声方差 `8/SNR_linear`、两个 static DMRS 平均、base seed `20260727` 和 absolute-trial seed derivation 均不改变。本节只改变接收机知识：

- ideal CSI：在 data RE 直接使用真实等效标量信道 $g_k=\sum_nV_{k,n}h_{k,n}$ 做均衡，不运行 CE，不生成 NMSE；
- matched estimated CSI：UE 知道 CDD delay，并对该候选使用真实等效频率协方差

$$
\mathbf R_g=\mathbf R_{\rm phy}\odot(\mathbf V\mathbf V^H)
$$

构造全带 frequency-LMMSE。令 $P,D$ 分别为 comb-6 导频和 data 子载波索引，两个 DMRS 平均后的 LS 噪声方差仍为 $4/SNR_{\rm linear}$，则

$$
\widehat{\mathbf g}_D
=\mathbf R_{g,DP}
\left(\mathbf R_{g,PP}+\frac{4}{SNR_{\rm linear}}\mathbf I\right)^{-1}
\bar{\mathbf z}_P.
$$

CE NMSE 保持第 10–15 节口径：每 trial 先在全部 data RE 上计算 $r_t=\sum|\widehat g-g|^2/\sum|g|^2$，在线性域跨 trial 求均值和 95% Monte Carlo 区间，最后转 dB；同时保存误差能量、真实能量和 $r_t$ 三类逐 trial 数组。matched 曲线与原透明曲线使用不同 candidate ID，避免把不同 receiver 假设误合并。

### 16.2 SNR 网格、预算与冻结规则

已有、不进入正式累计的 300-trial ideal prescan 对本候选在 `14/14.5/15/15.5/16/16.5/17/17.5/18/18.5 dB` 给出 BLER `0.1733/0.1433/0.1367/0.0933/0.0833/0.0367/0.03/0.0167/0.0133/0.00667`。据此在正式 trial 1 前冻结 ideal 网格为

$$
[14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19]\ \mathrm{dB}.
$$

matched estimated 在独立目录运行 `[14,15,16,17,18,19,20] dB`、每点 200 trials 的 prescan，BLER 为 `0.19/0.175/0.09/0.045/0.015/0.01/0`，CE NMSE 为 `-24.28/-24.76/-26.20/-26.82/-27.45/-28.54/-29.60 dB`；这些 trial 只用于定位 crossing，不进入正式累计。据此在任何正式 trial 1 前冻结正式网格为

$$
[14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20]\ \mathrm{dB}.
$$

其中 19/19.5/20 dB 覆盖 1% crossing 两侧和一个低于 1%的保护点；正式运行后不再后验改变网格。

两条正式曲线均采用每点至少 10,000 trials、目标至少 200 errors、最多 50,000 trials、单个可恢复区间不超过 1,000、batch size 20。达到最小 trial 且错误数达到 200 后停止；否则按 absolute trial 区间继续，达到 50,000 仍不足 200 errors 时保留该点并报告样本限制。正式数据不得混入 prescan trial。

### 16.3 实现、输出与验收

`tools/run_bler_curves.py` 的 explicit-delay CDD 路径增加显式 `covariance_mode`：原透明 estimated 配置保持 `physical` 默认行为；本节 matched estimated 必须使用 `matched_effective`，逐候选从真实 $\mathbf R_g$ 构造滤波器。manifest 和 filter diagnostics 必须记录 UE 是否知道 CDD、使用的协方差类型、候选特定滤波器、condition number、解析 matched trace-NMSE 和实际配置 hash。测试至少覆盖：旧透明配置行为不变、matched filter 权重确由真实 $\mathbf R_g$ 构造、ideal 不生成 CE 数据、第 16 节候选/网格/receiver 策略约束。

prescan、ideal 正式运行和 matched estimated 正式运行使用互相隔离的输出目录。正式产物保存展开配置、manifest/hash、逐点 CSV、逐 trial arrays、absolute-trial 区间、日志、稳定性与 BLER bracket 审计。最终把 matched estimated BLER/CE 加入第 3.6 节现有 estimated/CE 图，把小时延 CDD ideal BLER 加入 ideal 图；透明小时延曲线继续保留，ideal 图横轴按新网格扩展。两版 result 必须给出全部新增正式点的 trials/errors/BLER/Wilson 95%区间、estimated CE NMSE/95%区间、区间连续性、1% bracket、异常和精确证据路径。未经研究者确认，不更新 `KNOWLEDGE.md`、`GOALS.md`，不创建 Git checkpoint。
