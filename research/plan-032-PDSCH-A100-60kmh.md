# plan-032：A100 TDL-A 60 km/h 下的 CDD、透明 PRG 与过时 CSI 闭环预编码

状态：研究者于 2026-09-06 明确要求执行。本计划先完成实现、smoke 和预扫描；正式 trial 1 开始前，必须把预扫描选出的正式 SNR 网格追加并冻结在本文件中。

## 1. 目的与研究问题

本轮回答：result-028 中 A100 的多种 CDD 时延设计准则在 60 km/h 的 TDL-A 移动性信道下，使用已知速度的二维时频 LMMSE 信道估计时性能如何；相对不要求 UE 知道预编码向量的 transparent 6-RB precoder-cycling 基线，CDD 是否仍保持 estimated-CSI BLER 增益；再增加一条由过时下行 CSI 驱动的闭环预编码基线，量化 CSI 老化后的性能。

本轮只研究 TDL-A，不研究 CDL。结论只适用于本计划冻结的 8Tx/1Rx、单层、A100、60 km/h、30 kHz、comb-6 和无 ICI 链路，不外推至其他速度、反馈周期、天线相关性或实际码本 PMI。

## 2. 曲线、基线与公平性

### 2.1 CDD 候选

从 result-028 的 A100 source manifest 原样复用以下 10 个 CDD transmitter，不重新搜索、不改变 delay、相位参考或功率：

1. `A100_B0_QC`；
2. `A100_AP_RMS_T1`；
3. `A100_AP_TEPS_T1`；
4. `A100_AP_TU_NT`；
5. `A100_AP_TU_NTM1`；
6. `A100_AP_TALIAS_NT`；
7. `A100_S0_SIDON`；
8. `A100_AP_T2_01`；
9. `A100_MEFF_T2_04`；
10. `A100_MEFF_T2_06`。

CDD 仍按

$$
V_{k,n}=\exp\!\left(-j2\pi k j_n/576\right)
$$

构造，首个 active subcarrier 为相位参考，每个频点的 8 维向量平方范数为 8；数字循环移位不计入传播时延或 CP。UE 知道实际 CDD 矩阵和 60 km/h 速度，使用匹配的等效频域协方差及匹配的时间协方差。

### 2.2 transparent 基线

透明参考固定为 result-028 的 `A100_PRG_DFT8_6RB`：48 PRB 分成 8 个 6-RB PRG，每个 PRG 72 个连续子载波，依次使用 8 个未归一化 DFT 向量 `[0,1,2,3,4,5,6,7]`，每个向量平方范数为 8。UE 不需要知道 DFT 向量索引；信道估计在每个 PRG 内独立执行二维时频 LMMSE，使用 $8\mathbf R_{\rm phy}$，不跨 PRG 插值。

### 2.3 过时 CSI 闭环基线

新增 `A100_AGED_CSI_MRT_PRG6_SLOTS10`。30 kHz SCS 对应每 slot 0.5 ms，CSI report 周期取 `slots10 = 5 ms`。`slots10` 是 3GPP TS 38.331 的合法 `CSI-ReportPeriodicityAndOffset`；3GPP TS 38.508-1 的 FR1 测试配置也使用 `slots10`。本轮把它称为“标准允许且具有代表性的 FR1 配置”，不声称它是所有现网的统一常用值。

本基线使用一整个周期以前的物理信道快照。当前固定 CP 链路的 OFDM symbol 时长为 35.677083333 us；Sionna TDL 以相同 symbol 间隔连续生成，第 0 个样本作为反馈 CSI，第 140–149 个样本作为当前 10-symbol PDSCH，因此实际 CSI age 为 4.994791667 ms，与名义 5 ms 的相对差为 0.1042%。该实现误差必须写入展开配置和 result。

为同时贴近 PMI 的 PRG 粒度和研究者要求的“精确匹配滤波”，每个 6-RB PRG 用过时的 1Rx×8Tx 频域信道构造

$$
\mathbf G_b=\sum_{k\in b}\mathbf h_{\rm old,k}^{H}\mathbf h_{\rm old,k},\qquad
\mathbf w_b=\sqrt 8\,\mathbf u_{\max}(\mathbf G_b),
$$

并在当前 PDSCH 的该 PRG、全部 10 个 symbol 上保持 $\mathbf w_b$ 不变。这里使用未量化的精确主特征向量，属于“过时、未量化 subband CSI 的 MRT 上界型基线”，不是严格的 3GPP Type-I 码本 PMI；不得把二者混称。接收端与 transparent PRG 使用相同的 PRG 内二维时频 LMMSE 假设 $8\mathbf R_{\rm phy}$，因此它同时包含发射端 CSI 老化和当前 DMRS 信道估计的影响。

### 2.4 共同样本与归一化

同一 SNR/absolute trial 的 12 条曲线共享当前 TDL realization、payload、编码比特、data AWGN 和 DMRS AWGN；闭环基线额外从同一个连续 TDL realization 读取旧 CSI。所有预编码向量平方范数为 8，噪声方差统一为 $8/\mathrm{SNR}_{\rm linear}$，SNR 定义为每 Rx branch 平均接收 SNR。不同 SNR 使用由全局 seed、SNR 和 absolute trial 唯一派生的可重放随机流。

## 3. 移动信道与接收机

- Sionna 1.0.2 TDL-A，RMS delay spread 100 ns，carrier frequency 3.5 GHz，UE speed 60 km/h，20 sinusoids，8Tx/1Rx，各 Tx branch 独立且同统计；
- 只采用 OFDM symbol 采样的时变频域乘法信道。假设定时和载波同步，不模拟 CFO、ICI 或 ISI；多普勒扩展只表现为 symbol 间复信道变化和信道老化；
- 48 PRB、576 active subcarriers、30 kHz SCS、4096 FFT、CP 288 samples、10 PDSCH symbols；DMRS symbols `[2,7]`，comb-6，每个 DMRS symbol 96 pilots，共 192 个独立 LS 观测，5568 data RE；
- 两个 DMRS 不再先平均。二维时频 LMMSE 直接使用 192 个 `(symbol, subcarrier)` 观测估计全部 data RE；单个 LS 观测噪声方差仍为 $8/\mathrm{SNR}_{\rm linear}$；
- 时间协方差由 UE 已知的 60 km/h、3.5 GHz 和实际 OFDM symbol duration 通过 Sionna `tdl_time_cov_mat` 构造；频域协方差由 TDL-A 100 ns PDP 构造；
- CDD 的匹配频域协方差为 $\mathbf R_g=\mathbf R_{\rm phy}\odot(\mathbf V\mathbf V^H)$。这里 $\mathbf V$ 的行范数平方为 8，因而无需另乘 $1/8$；transparent PRG 和闭环 PRG 在各 PRG 内使用 $8\mathbf R_{\rm phy}$；
- 16QAM、NR 256QAM MCS table 的 MCS 8、目标码率 553/1024、Sionna LDPC、最多 8 次迭代，与 result-028 保持一致；
- CE NMSE 在 5568 个 data RE 上逐 trial 计算线性误差能量/真实能量，再对 trial 求均值并转 dB。

## 4. 实现范围与测试

1. 在 `cdd_lls/phy/channel_tdl.py` 增加可复用的 Sionna TDL active-subcarrier 历史/当前抽样接口：一次 SoS realization 内生成至最大时间索引，只把请求的 CIR 时间样本转换为 active-subcarrier 响应；旧接口和零速行为不变。
2. 在 `cdd_lls/phy/estimators.py` 增加任意 pilot/data coordinate 子集的二维 RMMSE 构造，并据此实现 PRG-local 二维滤波；全带现有接口保持兼容。
3. 新增范围固定的 `tools/run_plan032_tdl_mobility.py` 专题入口；复用 `tools/run_bler_curves.py` 的 LDPC decode、CE 统计和 MRC helpers，以及 plan-027 的资源、manifest 与 seed helpers，不复制这些已验证算法。移动性所需的历史信道和二维滤波放在 `cdd_lls/phy/` 可复用层；旧固定静态入口行为不变。
4. 新增 `configs/bler_curves_result032_a100_60kmh_<smoke|prescan|formal>.yaml` 和 result-032 专用分析入口。
5. 单元测试至少覆盖：历史/当前样本索引与 shape；60 km/h 下经验 symbol 相关和 Sionna 理论时间相关同量级；0 km/h 的二维滤波与“两 DMRS 平均后频域 LMMSE”等价；CDD 协方差功率和相位；PRG 不跨边界；过时 CSI 与当前 H 来自同一 realization；MRT PRG 权值范数为 8 且确实最大化旧 CSI 的 PRG 平均功率；旧 result-028 配置 validate 回归。
6. smoke 使用 20 trials，检查无噪声或高 SNR 译码、CE 数组有限、12 条曲线均落盘、trial pairing、seed 重放、维度、功率和内存峰值。

## 5. 预扫描、正式预算与统计方法

预扫描固定 SNR 网格为 `[12,14,16,18,20,22,24] dB`，每点 400 个共同 paired trials；预扫描只用于选择正式网格，不并入正式结果。

正式网格必须在预扫描后追加至本计划，覆盖各曲线 10% 和 1% BLER 的双侧 bracket；目标附近间隔不大于 0.25 dB，过渡段可用 0.5 dB。正式公共网格每点先运行至少 1,000 个共同 trials；所有实际 10%/1% bracket 的端点追加至至少 3,000 trials。沿用 result-028 的目标带最低样本判据：BLER 为 0.5%–2% 的端点至少需要 30 个错误块；不足时按不重叠 absolute-trial 区间继续追加，最多 10,000 trials。若达到上限仍不足，保留点并标注样本边界。曲线在出现首个 BLER 低于 0.3% 且 Wilson 95%上界低于 1%后可停止向更高 SNR 扩展；冻结范围内未双侧 bracket 的目标不外推。执行中已先完成的 14、14.25、16.5、18.75 dB 为 3,000 trials，属于超过最低预算的合法正式点。

逐点 BLER 报告 Wilson 95%区间；10%/1%目标 SNR 只在线性 SNR 邻点双侧 bracket 内做对数 BLER 插值，并用 paired bootstrap 给出差值区间。每个 CDD 候选分别相对 transparent PRG6 和 aged-CSI MRT 报告目标 SNR 差；同时报告其相对 result-028 static 同名曲线的目标 SNR 变化。由于 0 km/h 与 60 km/h 的 trial 不配对，跨速度差值不得宣称 paired 方差缩减。

## 6. 输出、复现与审计

正式产物写入 `outputs/experiment032_tdl_mobility/20260906_main/`，按 `smoke/`、`prescan/`、`formal/`、`final/` 隔离。至少保存：原配置及 SHA-256、展开配置、源 result-028 manifest/hash、12 条曲线定义、反馈周期及实际 CSI age、Sionna/软件版本、时间索引、seed 规则、任务 schedule、逐区间 CSV、逐 trial error flags、逐 trial CE 误差/真实能量数组、闭环旧/当前信道相关诊断、滤波器条件数、运行日志、稳定性审计、目标 bracket 与 bootstrap 数据、共享样式和图对应精确 CSV。

核心入口预期为：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan032_tdl_mobility.py --config configs\bler_curves_result032_a100_60kmh_smoke.yaml --stage validate
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan032_tdl_mobility.py --config configs\bler_curves_result032_a100_60kmh_smoke.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan032_tdl_mobility.py --config configs\bler_curves_result032_a100_60kmh_prescan.yaml --stage validate
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan032_tdl_mobility.py --config configs\bler_curves_result032_a100_60kmh_prescan.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan032_tdl_mobility.py --config configs\bler_curves_result032_a100_60kmh_formal.yaml --stage validate
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan032_tdl_mobility.py --config configs\bler_curves_result032_a100_60kmh_formal.yaml --stage run
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_result032_tdl_mobility.py
```

## 7. 验收、停止和 result 要求

验收要求：

1. 0 km/h 数学回归证明二维实现与 result-028 的 two-DMRS averaged frequency-LMMSE 口径等价；60 km/h 时不再平均两个 DMRS；
2. 经验时间相关、闭环实际 CSI age、旧/当前信道相关和最大 Doppler 诊断与配置一致；
3. 12 条正式曲线完成冻结网格和预算，数据可由 absolute trials 无损续跑且无重叠/缺口；
4. 每条 CDD 对 transparent PRG6 和 aged MRT 的 10%/1%差值均按“已 bracket / 未 bracket”明确报告，不能因目标未闭合而外推；
5. 结果明确区分：移动信道事实、二维估计影响、透明/非透明协方差知识、闭环 CSI 老化、以及“未量化 MRT 并非标准码本 PMI”的适用边界；
6. 图与数据满足 `docs/agent/RESULT_SPEC.md`，并按约 13 cm PPT 半页宽完成可读性检查。

若 smoke 出现功率、维度、seed pairing、历史/当前连续性或静态等价失败，立即停止正式运行并修复；若预扫描表明 24 dB 仍不能形成必要上界，只能在正式 trial 1 前扩展并冻结网格；正式开始后不得根据结果回选候选、改变反馈周期或改变接收机。

最终同时生成 `research/result-032-PDSCH-A100-60kmh.md` 与 `research/result-032-PDSCH-A100-60kmh-text.md`；图放 `docs/figures/result-032/`，无图版必须含足以独立核验的逐点关键表和精确原始路径。结果经研究者确认前，不更新 `KNOWLEDGE.md`、`GOALS.md`，不创建 Git checkpoint。

## 8. 标准依据

- 3GPP TS 38.331 / ETSI TS 138 331 V18.8.0，`CSI-ReportPeriodicityAndOffset` 以 slot 配置，包含 `slots10`：https://www.etsi.org/deliver/etsi_ts/138300_138399/138331/18.08.00_60/ts_138331v180800p.pdf
- 3GPP TS 38.508-1 的 FR1 测试参数示例包含 `CSI-ReportPeriodicityAndOffset = slots10`；该证据只支持“标准测试中使用”，不支持“所有现网普遍采用”。

## 9. 正式 SNR 网格冻结区

2026-09-06 在任何 formal trial 开始前冻结。预扫描证据为 `outputs/experiment032_tdl_mobility/20260906_main/prescan/final/estimated_csi_bler_points.csv`：可闭合曲线的 10%/1% 过渡集中在 14–20 dB；`AP_TU_NT`、`AP_T2_01`、`MEFF_T2_04`、`MEFF_T2_06` 至 24 dB 仍呈地板。正式网格固定为 `[14,14.25,14.5,14.75,15,15.25,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,19.75,20,22,24] dB`。其中 14–20 dB 的 0.25 dB 网格用于目标 bracket，22/24 dB 只审计高 SNR 地板；不增加 24 dB 以上点。冻结配置为 `configs/bler_curves_result032_a100_60kmh_formal.yaml`。

执行记录补充：14 dB 首点完成并落盘后，为降低墙钟时间，剩余冻结 SNR 按不重叠集合拆为 `formal_shard1/2/3` 三个并行目录；分片配置分别为 `configs/bler_curves_result032_a100_60kmh_formal_shard1.yaml`、`...shard2.yaml`、`...shard3.yaml`。分片只改变执行分组和输出目录，`seed=20260727`、batch 20、每点 3,000 trials、候选、接收机和全部物理参数不变；不同 SNR 本来就使用独立 seed，因此不改变任一点样本。最终必须核验四个目录的 SNR 并集恰好等于上述冻结网格且交集为空。

统计执行补充（2026-09-06，剩余分片第二个 SNR 尚未落盘时）：三分片并行实测每个 3,000-trial 点约 17 分钟，若所有远离目标的地板点也一律运行 3,000，会产生大量不影响目标 bracket 的计算。故公共冻结网格的初始最低预算改为 1,000 trials；所有 10%/1% bracket 的两个端点必须再追加至至少 3,000 trials，若任一相关曲线在端点处 BLER 为 0.5%–20% 且错误块少于 100，则继续追加，最多 10,000 trials。14、14.25、16.5、18.75 dB 已完成的 3,000 trials 原样保留。该变更只调整预先声明的自适应预算，不删除正式样本、不改变 SNR 网格/候选/物理口径，也不依据候选优劣回选；被中断且未形成完整区间的部分计算不落盘、不计入结果。初始 1,000-trial 分片使用 `...formal_shard[1-3]_initial1000.yaml`，原 3,000-trial 配置和已完成点的展开配置保留用于审计。

统计执行补充（2026-09-07，公共网格全部完成、refine 尚未开始）：为与 result-028 的 1% 目标样本判据保持一致，上段“0.5%–20% 均追到 100 errors”修订为“0.5%–2% 至少 30 errors”；10%端点在 3,000 trials 下自然远超 30 errors。该修订统一作用于所有候选和两条基线，不按优劣选择。公共网格识别出的 refine SNR 为 shard1 `[14.5,14.75,15,15.25,15.5,15.75,16.25]`、shard2 `[16.5,16.75,17,17.25,17.5,17.75]`、shard3 `[19.25,19.5]`；16.5 已有 3,000 trials，其余端点从 absolute trial 1001 追加。配置为 `configs/bler_curves_result032_a100_60kmh_formal_refine_shard[1-3].yaml`。

低 BLER 错误数补充（2026-09-07，全部 bracket 端点达到 3,000 trials 后）：审计发现 16.5、17、17.75、19.5 dB 的一个或多个 0.5%–2% 端点只有 19–29 个错误块，未满足上述 30-error 判据。统一追加 16.5/17/17.75 dB 至 5,000 trials，19.5 dB 至 6,000 trials；配置为 `configs/bler_curves_result032_a100_60kmh_formal_refine_lowbler_shard2.yaml` 和 `...shard3.yaml`。选择总量只由最低错误数和已测 BLER 决定，全部 12 条曲线仍共享相同 absolute-trial 区间。

第一次低 BLER 追加后，`AP_TU_NTM1` 的 1% 实测 bracket 因 17 dB 合并 BLER 从 0.933% 更新为 1.06%，由 `[16.75,17]` 移至 `[17,17.25]`；新高侧 17.25 dB 在 3,000 trials 下为 22 errors/0.733%。因此按同一判据用 `configs/bler_curves_result032_a100_60kmh_formal_refine_lowbler_17p25.yaml` 将 17.25 dB 统一追加至 5,000 trials，再重新冻结最终统计。

1%可视化降波动补充（2026-09-07，研究者审阅首版图后）：首版图的剩余抖动来自 1%邻域内仍只有 1,000 trials 的非 bracket 中间点，而不是信道或译码异常。按“只平滑至 1%目标、不追逐 0.1%”的要求，将 16、18、18.25、18.5、19 dB 从 absolute trial 1,001 追加至 3,000 trials；配置为 `configs/bler_curves_result032_a100_60kmh_formal_smooth_to_1pct_shard[1-3].yaml`。19.75/20 dB 时所有可闭合曲线均已低于 1%，仍高于 1%的只有四条未闭合地板曲线，故不追加。分析时使用 `--skip-ce-plot`，保持已交付 NMSE 图文件不重画。BLER 图保留未经修改的原始 marker，同时叠加 trials 加权的 log-BLER 非增 PAVA 与保单调 PCHIP 展示线；拟合不参与目标 SNR、差值或置信区间计算。

result-028 样本量级补充（2026-09-07，研究者认为上述曲线仍有明显有限样本波动后）：result-028 第 3.7 节中用于 1%邻域比较的 transparent PRG6、small-delay transparent CDD 和 matched CDD 均以每点至少 10,000 trials、目标 200 errors 自适应追加；1%邻域的实际点为 12,040–38,000 trials，matched 尾点上限为 50,000 trials。相比之下，本轮此前 1%邻域只有 1,000–6,000 trials 和 6–83 errors，统计量级不足。故对 14–20 dB 的全部点先提高到至少 10,000 trials；再根据追加前原始 BLER，对 0.5%–2% 内的任一曲线按 $N_{\rm target}=\lceil 200N/E\rceil_{20}$ 计算该 SNR 的共同 paired-trial 目标，上限 50,000。首轮冻结目标合计为 491,880 个 SNR-level trials，即在现有 82,000 个基础上追加 409,880 个；12 条曲线等价追加 4,918,560 candidate-trials。完成后重新检查 0.5%–2% 点，若仍少于 200 errors，则沿相同绝对 trial 序列继续追加至满足判据或 50,000 上限。22/24 dB 只用于地板审计且不参与 1%图示，因此不追加。对应配置为 `configs/bler_curves_result032_a100_60kmh_result028_scale_<formal|shard1|shard2|shard3>.yaml`；专题入口新增逐 SNR 总 trial 目标字段，只改变统计预算，不改变物理链路、seed、候选或接收机。

手动编排与恢复补充（2026-09-08）：新增 `tools/run_result032_result028_scale.py` 作为上述四份配置的唯一手动入口。入口启动时读取每个分片的 `intervals.csv`，验证每个 candidate/SNR 的 absolute-trial 区间从 1 连续且 12 条曲线保持 paired，打印当前/目标/剩余预算，并以最多三个子进程调度分片。内层 `run_plan032_tdl_mobility.py` 每完成一个 SNR 点即原子更新 `intervals.csv` 及分片合并 CSV，同时保存该点 error flags、CE NMSE 数组和滤波诊断；因此异常中断最多丢失正在计算但尚未落盘的单个 SNR 点，已完成点不会重算。重复执行同一命令会从首个未达到目标的 SNR 续跑。编排器使用由操作系统释放的单实例文件锁防止两个窗口并发写同一输出，并把状态摘要和逐分片控制台日志保存在 `outputs/experiment032_tdl_mobility/20260906_main/result028_scale_orchestration/`。这项修改仅涉及外围预算、调度、日志与 checkpoint 检查，不改变任何信道、估计、预编码、译码或随机数定义。

200-error 二次验收补充（2026-09-09）：首轮 491,880 个 SNR-level trials 完成后，31 个最终 BLER 位于 0.5%–2% 的点中仍有 9 个少于 200 errors，分别位于 16.25、16.5、17、17.5、18.5、19、19.25、19.75 dB。按上一补充中已经冻结的公式和 20-trial 对齐重新计算共同目标，将这些 SNR 分别提高至 19,300、32,920、38,140、14,760、14,040、16,900、21,080、33,340 trials；其中同一 SNR 有多条不足曲线时取所需目标的最大值。二次追加共 39,180 个 SNR-level trials，等价 470,160 candidate-trials；仍低于每点 50,000 上限。该追加是预先声明的样本充分性闭环，不改变候选、SNR 网格、目标定义或物理口径。二次完成后必须再次以最终合并 BLER 检查 0.5%–2% 内所有点；未达到 200 errors 的点继续按同一规则追加或在 50,000 上限处标注。

200-error 最终边界补充（2026-09-09）：二次追加后，上述 9 个点中的 8 个达到 200 errors；唯一例外为 aged-MRT 19.75 dB，最终为 179/33,340，即 0.5369%，仍处于 0.5%–2% 审计带。若再次只按 $200N/E$ 设置目标，新样本下错误数仍有约一半概率低于 200，因此为避免重复追加并遵守既定上限，直接把该 SNR 的共同目标提高至 50,000 trials。新增 16,660 个 SNR-level trials，等价 199,920 candidate-trials；完成后无论该曲线是否达到 200 errors，均按“达到 50,000 上限”的预定停止条件结束并报告实际错误数。

最终执行与作图补充（2026-09-09）：aged-MRT 19.75 dB 在 50,000 trials 下为 270 errors/0.540%，因此 31 个最终 BLER 位于 0.5%–2% 的点全部达到至少 200 errors。14–20 dB 共完成 547,720 个共同 trials；连同不追加的 22/24 dB 各 1,000 trials，最终分析覆盖 549,720 个共同 trials、6,596,640 candidate-trials。研究者要求增加样本后直接显示原始数据，故最终 BLER 图移除先前的 PAVA/PCHIP 展示拟合，marker 与相邻线段均直接使用原始 Monte Carlo BLER；图窗仍为 14–20 dB、BLER 不低于 0.5%，NMSE 图不重画。最终统计、图和成对 result 由 `tools/analyze_result032_tdl_mobility.py --skip-ce-plot` 生成。
