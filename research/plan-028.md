# plan-028：comb-6 三类 CSI 曲线统一补充与 TDL-A 300 ns 方案对比

> 状态：已确认，可执行。研究者于 2026-08-03 明确要求：基于 result-027 的 A30 comb-6 与 A100 comb-6 正式数据，按统一大字号样式重绘 estimated-CSI BLER，并在相同 SNR 采样点补充真实 ideal-CSI 链路 BLER 与 CE NMSE；新增 TDL-A 300 ns、comb-6 的各方案对比，交付缩写图例和 `candidate NN` 图例版本。未单独给出的系统、统计预算、搜索与闭合规则按 A100 comb-6 增补保持不变。

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

三个场景均固定：48 PRB、$K=576$、子载波间隔 30 kHz、8 Tx / 1 Rx、单层、static Sionna TDL-A、零速度、DMRS comb-6、两个 DMRS symbol、每 symbol 96 个 pilot、总计 192 pilot RE、5568 data RE、256QAM MCS 8、Sionna LDPC 最多 8 次迭代。每个 pilot RE 功率不变。estimated-CSI 接收机为 two-DMRS averaged、known-$\mathbf V$/known-PDP、V-aware matched 全带频域 LMMSE；LLR 噪声方差不额外加入 CE 误差项。

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

`research/result-028.md` 与 `research/result-028-text.md` 必须在同一个新增主章节中集中给出 A30 comb-6、A100 comb-6、A300 comb-6 的编号与时延表、三类曲线、原始采样表、目标和不确定性；图文版嵌入图，无图版不得嵌图但必须引用全部 CSV 和关键数值。两版结论、异常、验收、适用范围和证据路径必须一致，并更新 `research/README.md` 索引。未经研究者确认，不更新 `KNOWLEDGE.md` 已验证结论或 `GOALS.md` 阶段验收。
