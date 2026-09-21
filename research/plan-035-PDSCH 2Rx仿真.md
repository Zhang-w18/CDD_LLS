# plan-035：PDSCH 2Rx、60 km/h 下 4/8Tx 的 estimated CSI、ideal CSI 与 CE NMSE 曲线

> 状态：**草案，待研究者确认后执行**。研究者于 2026-09-17 提出：时延曲线和系统设置与 plan-033 相同，只把接收天线改为 2Rx；当前阶段只运行 60 km/h，并输出 estimated-CSI BLER、ideal-CSI BLER 和 NMSE 曲线。本 plan 将“plan-033 的设置”展开为下列显式参数；第 8 节统计预算和停止条件仍待研究者确认，确认前不启动正式 trial。

## 1. 目的与交付目标

本轮在 plan-033 的 PDSCH 单层移动 TDL-A 链路上，把接收天线数从 4 改为 2，回答：

1. 在 `8Tx/2Rx/1 layer` 与 `4Tx/2Rx/1 layer`、60 km/h 下，六个冻结方案的 estimated-CSI BLER、ideal-CSI BLER 和 data-RE CE NMSE 随 SNR 如何变化；
2. ideal CSI 与 estimated CSI 的 BLER 间隔如何随方案变化，哪些排序主要来自理想分集，哪些排序同时受到 DMRS 信道估计影响；
3. 2Rx 下 `S0_SIDON` 相对 `B0_QC`、small-delay matched 相对 transparent、aged MRT 相对 transparent PRG6 的 10%/1% BLER 目标 SNR 差异如何；
4. 在与 plan-033 相同的每根 Rx 分支 SNR 定义下，2Rx 结果与对应 4Rx/60 km/h 结果有何变化。该跨 plan 比较不是 paired comparison，只报告独立不确定性。

本轮不搜索新 delay、不改变发射端候选、不运行 3 km/h、不引入空间相关、CDL、ICI/CFO、反馈量化或多层检测。estimated-CSI BLER 是正式性能判据；ideal-CSI BLER 和 CE NMSE 是机理诊断，不能单独作为方案胜出依据。

本轮对应 `GOALS.md` 的多 Rx 适用范围检查，并遵循 `DESIGN.md` 第 3.7--3.9、5.4--5.5 和第 6 节关于单位总发射功率、接收机知识、ideal/estimated CSI 区分和公平比较的要求。

## 2. 冻结系统条件

| 参数 | 冻结取值 |
|---|---|
| 场景 | `A100_NT8_NR2_V60`、`A100_NT4_NR2_V60` |
| 天线与层数 | `8Tx/2Rx/1 layer` 与 `4Tx/2Rx/1 layer` |
| 波形 | 48 PRB，576 个连续有效子载波，30 kHz SCS，FFT 4096，CP 288 samples，10 个 PDSCH symbols |
| 信道 | Sionna 1.0.2 TDL-A，RMS delay spread 100 ns，3.5 GHz，60 km/h，20 sinusoids；Tx--Rx 分支独立、同 PDP、无空间相关 |
| 时变模型 | OFDM-symbol 采样的时变频域乘法；假设定时与载波同步，不模拟 CFO、ICI 或 ISI |
| DMRS | symbols `[2,7]`，comb-6；每 symbol 96 pilot RE，共 192 个独立 LS 观测；5568 data RE |
| estimated-CSI 接收机 | 两个 DMRS 不平均；使用真实 60 km/h 时间协方差及候选指定的频域协方差做二维时频 LMMSE；两根 Rx 独立估计，data RE 上做 2Rx MRC |
| ideal-CSI 接收机 | 在每个 data RE 直接使用当前真实等效信道做 2Rx MRC；不使用 DMRS 估计值，也不统计 CE NMSE |
| 调制编码 | 16QAM；NR 256QAM MCS table 的 MCS 8；目标码率 553/1024；Sionna LDPC 最多 8 次迭代 |
| CSI age | 10 slots，名义 5 ms；旧 CSI 为第 0 个 symbol sample，当前 PDSCH 为第 140--149 个 sample，实际 age 4.994791667 ms |
| 随机性 | base seed `20260727`；absolute-trial 可重放 |
| 执行设备 | CPU-only；运行时保存 CPU placement、逐 batch RSS 峰值和墙钟耗时 |

`ideal CSI` 只改变接收端。对 aged-CSI MRT，发射权值仍由 5 ms 旧 CSI 产生；不得用当前 CSI 重新生成发射权值。这样 ideal/estimated 差值只反映接收端 CSI，而不同时改变发射策略。

## 3. 归一化、SNR、MRC 与 NMSE 定义

发射端使用单位总功率预编码。若恒模相位矩阵为 $\mathbf V$，实际预编码为

$$
\mathbf W=\frac{\mathbf V}{\sqrt{N_t}},
\qquad
\lVert\mathbf W_k\rVert_2^2=1.
$$

每根 Rx 分支的噪声方差为

$$
\sigma_n^2=\frac{1}{10^{\mathrm{SNR}_{\rm dB}/10}}.
$$

横轴 SNR 是每根 Rx 分支的平均接收 SNR，不因 $N_r=2$ 再乘或除 2。2Rx MRC 包含约 3.01 dB 的平均合并功率增益和接收分集。matched CDD 使用 $\mathbf R_g=\mathbf R_{\rm phy}\odot(\mathbf W\mathbf W^H)$；transparent PRG、transparent CDD 和 aged MRT 使用 $\mathbf R_{\rm phy}$。data AWGN 和每个 DMRS LS 观测的噪声方差都采用上述 $\sigma_n^2$。

estimated-CSI MRC 在每个 data RE 使用

$$
\widehat x=
\frac{\sum_{r=0}^{1}\widehat g_r^*y_r}
{\sum_{r=0}^{1}|\widehat g_r|^2},
\qquad
\widehat\sigma_{\rm eff}^2=
\frac{\sigma_n^2}{\sum_{r=0}^{1}|\widehat g_r|^2}.
$$

ideal-CSI MRC 用 $g_r$ 替换 $\widehat g_r$。LLR 噪声方差不额外加入 CE-error-aware 项，以保持 plan-033 的链路口径。

CE NMSE 只对 estimated CSI 定义。每个 trial 在两根 Rx 和全部 5568 个 data RE 上先计算

$$
\mathrm{NMSE}_t=
\frac{\sum_{r=0}^{1}\sum_{(s,k)\in D}
|\widehat g_{t,r,s,k}-g_{t,r,s,k}|^2}
{\sum_{r=0}^{1}\sum_{(s,k)\in D}|g_{t,r,s,k}|^2},
$$

再在线性域跨 trial 求均值，最后转换为 dB。保存 trial 级线性 NMSE、和、平方和、均值及 95% Monte Carlo 区间；禁止先对 Rx、PRG 或 RE 分别转 dB 后平均。

## 4. 冻结的六个方案与时延曲线

CDD 相位统一为

$$
V_{k,n}=\exp\left(-j2\pi k j_n/576\right),
\qquad
\tau_n=\frac{j_n}{576\times30\,\mathrm{kHz}},
$$

其中 $k$ 是从首个 active subcarrier 起算的局部索引，一个 delay-grid 单位为 57.870370 ns。人工 delay 是数字相位/循环移位，不计入真实传播时延或 CP。

| $N_t$ | 方案 | delay grid $\mathbf j$ / PRG 规则 | 接收端频域知识 |
|---:|---|---|---|
| 8 | `B0_QC` | `[0,9,18,27,36,45,54,63]` | matched CDD covariance |
| 8 | `S0_SIDON` | `[0,1,3,7,12,20,30,65]` | matched CDD covariance |
| 8 | `transparent PRG6` | 8 个单位范数 DFT8 向量，PRG index `[0,1,2,3,4,5,6,7]` | 每 PRG physical covariance |
| 8 | `aged-CSI MRT PRG6` | 每 trial、每 PRG 用旧信道生成单位范数主特征向量 | 每 PRG physical covariance |
| 8 | `small-delay CDD, transparent` | `[0,0.25,0.5,0.75,1,1.25,1.5,1.75]` | physical covariance，UE 不知道 delay |
| 8 | `small-delay CDD, matched` | 与上一行相同发射向量 | matched CDD covariance |
| 4 | `B0_QC` | `[0,24,48,72]` | matched CDD covariance |
| 4 | `S0_SIDON` | `[0,1,3,7]` | matched CDD covariance |
| 4 | `transparent PRG6` | 4 个单位范数 DFT4 向量，PRG index `[0,1,2,3,0,1,2,3]` | 每 PRG physical covariance |
| 4 | `aged-CSI MRT PRG6` | 每 trial、每 PRG 用旧信道生成单位范数主特征向量 | 每 PRG physical covariance |
| 4 | `small-delay CDD, transparent` | `[0,0.25,0.5,0.75]` | physical covariance，UE 不知道 delay |
| 4 | `small-delay CDD, matched` | 与上一行相同发射向量 | matched CDD covariance |

物理时延分别为 grid 坐标乘以 57.870370 ns；正式 manifest 必须保存精确 grid 坐标、ns、FFT-sample 等效值、comb-6 residue、无序二元和、pair/fold gap、pilot rank、condition number 和选择规则。4Tx 数组是独立冻结配置，不允许从 8Tx 数组在运行时隐式截断。

两根 Rx 的 aged MRT 对每个 6-RB PRG 构造

$$
\mathbf G_b=\sum_{k\in b}\mathbf H_{{\rm old},k}^{H}\mathbf H_{{\rm old},k},
\qquad
\mathbf H_{{\rm old},k}\in\mathbb C^{2\times N_t},
$$

并取 $\mathbf G_b$ 的单位范数主特征向量。它是过时、未量化 subband CSI 的 MRT 上界型基线，不是 3GPP Type-I 码本 PMI。

## 5. 公平性与 paired trial

每个场景包含 6 个发射方案 × 2 个接收模式，共 12 条 BLER 曲线；另有 6 条 estimated-CSI CE NMSE 曲线。每个 `scenario + SNR + absolute trial` 必须共享：

- 当前时变 TDL realization 和 aged-MRT 所需的旧信道样本；
- payload、编码比特和原始 data AWGN；
- 六个 estimated-CSI 方案使用的原始 DMRS AWGN。

同一发射方案的 ideal/estimated 接收端必须使用相同发射波形、当前信道、payload 和 data AWGN。ideal 分支不得生成新的信道或 data noise；estimated 分支独有的 DMRS noise 不进入 ideal 分支。六个方案之间以及同一方案的 ideal/estimated 差值可使用 paired bootstrap。4Tx 与 8Tx 是不同维度场景，跨 $N_t$ 比较不声明 paired 方差缩减；035 的 2Rx 与 033 的 4Rx 比较也按独立场景处理。

## 6. 实现范围与兼容性

现有 `tools/run_plan033_tdl_mobility_mimo.py` 把 `n_rx=4`、`receiver=estimated` 以及部分输出 identity 写死。035 不能只复制 YAML。实现应复用 033 已验证的候选、移动 TDL、二维时频 LMMSE、PRG、aged MRT、MRC、LDPC、interval 合并和 absolute-trial seed 逻辑，不复制第二套链路主循环。

预计修改：

1. 把 033 runner 中与 plan 编号无关的单层移动 MIMO trial 核心抽为可复用函数，显式接收 `n_rx` 和 `receiver_modes`；保留 033 schema、candidate ID、输出字段、seed 和恢复行为；
2. 新增 `tools/run_plan035_pdsch_2rx.py`，只接受 `n_rx=2`、`speed_kmh=60`、`n_tx in {4,8}` 和 `receiver_modes=[estimated,ideal]`；所有物理参数显式写入展开配置；
3. 输出 interval identity 增加 `receiver`，estimated 与 ideal 的 error flags 分开保存；CE 数组只允许出现在 estimated record，ideal record 的 CE 字段不得伪填零值；
4. 新增 `tools/analyze_result035_pdsch_2rx.py`，读取本地 CSV/NPY，完成 bracket、Wilson 区间、paired bootstrap、estimated/ideal gap、CE NMSE 汇总和本地图生成；
5. 新增 smoke、prescan 配置；formal YAML 只能在 prescan 完成并冻结 SNR 网格后生成。

现有 033 正式实验仍在进行中。重构必须先用现有 033 配置做 `--stage validate`，并通过定向回归，禁止改变已存在的 033 resolved-config hash 语义、candidate 定义、trial seed 或输出合并规则。若无法在不改变 033 行为的条件下安全抽取公共核心，则先在 `cdd_lls/` 增加公共函数，再让 033/035 两个薄入口分别调用；不得复制完整 trial 主循环。

## 7. 测试与 smoke

新增测试只覆盖 035 增量：

1. 035 配置只接受 `(n_tx,n_rx,speed)=(4,2,60)` 或 `(8,2,60)`，拒绝隐式默认 Rx 数；
2. 4Tx/8Tx 的 delay、DFT-PRG 映射和 candidate manifest 与第 4 节完全一致；
3. 2Rx aged-MRT Gram 对两根 Rx 和 PRG 子载波求和正确，且权值平方范数为 1；
4. ideal MRC 使用真实 data-RE 等效信道，estimated MRC 使用 LMMSE 输出；零噪声且注入精确信道时二者译码输入一致；
5. 一个小型 paired trial 验证 12 条 BLER 记录共享 trial key、channel/payload/data-noise seed，只有 estimated 记录含 CE 数组；
6. interval 合并按 `candidate_id + receiver + SNR` 分组，absolute-trial 连续且不会把 ideal/estimated 混合；
7. 033 的现有定向测试和一份 033 V60 配置 validate 通过，确认兼容性。

smoke 分别运行 `A100_NT4_NR2_V60` 和 `A100_NT8_NR2_V60`，每场景取 `[4,18] dB`、每点 20 个 paired trials。检查 12 条 BLER 记录、6 条有限 NMSE 记录、维度、功率、old/current 同 realization、MRC denominator、CPU placement、峰值 RSS、墙钟耗时和所有数组引用。`batch_size` 从 25 起步；若 RSS 超出 plan-033 smoke 的可用范围，只允许降低 batch，不改变 seed、trial 或物理定义，并把理由和实测值写入展开配置。

## 8. 预扫描、正式预算与停止条件（待确认）

### 8.1 建议的预扫描

每个场景先运行 400 个 paired trials/点，初始公共网格建议为

`[4,6,8,10,12,14,16,18] dB`。

该范围依据 plan-033 的 4Rx/60 km/h 初始 crossing 约位于 5--7.6 dB，并考虑从 4Rx 降为 2Rx 后约 3 dB 的平均 MRC 功率差；这只是预扫范围推断，不是性能结论。若任一 estimated/ideal 曲线未同时观察到 10%和1% BLER 的上下侧点，只向缺失方向按 2 dB 扩展，边界暂定 `[-2,22] dB`。达到边界仍未 bracket 时暂停并报告，不外推正式 crossing。

### 8.2 建议的正式网格与预算

prescan 后分别为两个场景冻结一套 12 条 BLER 曲线共同使用的 SNR 网格：10%和1% crossing 邻域最大间隔 0.25 dB，过渡区最大 0.5 dB，只保留 bracket 所需范围及两侧最多一个保护点。正式 trial 1 前把网格、配置 SHA-256、预算和预计总 trial 数写入本 plan 的冻结区并由研究者确认。

正式预算采用预声明的分阶段自适应规则：

1. 每个冻结 SNR 点先运行 1,000 个共同 paired trials；检查点间隔固定为 1,000 trials，单点上限 50,000 trials；
2. runner 在所有点达到 1,000 trials 后，从保存的真实 BLER 自动识别每条 estimated/ideal 曲线的 10%和1%相邻双侧 bracket；缺少任一 bracket 时停止并报告，不外推；
3. 10% bracket 端点的观测 BLER 位于 5%--20%时要求至少 200 个错误块；1% bracket 端点位于 0.5%--2%时要求至少 200 个错误块。区间外的保护端点不追 200 errors；
4. 同一 SNR 只要任一曲线的端点条件未满足，该 SNR 的 12 条 BLER 曲线就共同追加 1,000 trials；非 bracket 点在完成初始 1,000 trials 后停止；
5. 达到 50,000 trials 仍不足时保留原始点并标注样本上限，不平滑、不替换伪计数；
6. 10%/1%目标只在相邻真实采样点的双侧 bracket 内按 log-BLER 线性插值，不使用 PAVA、PCHIP、单调修正或外推；
7. 每条曲线均取得 10%和1%真实双侧 bracket，且所有需追加的端点达到 200 errors 或 50,000 trials，即完成该场景；
8. CE NMSE 与 estimated 曲线使用完全相同的 trial，不单独追加 NMSE-only trials。

Wilson 95%区间用于单点 BLER；同场景方案差值和同方案 ideal/estimated 差值使用 absolute-trial paired bootstrap。跨 4Tx/8Tx、2Rx/4Rx 的差值使用独立场景不确定性。

自适应检查只发生在上述固定 1,000-trial 检查点，停止原因和每轮目标写入 `adaptive_status.json`。Wilson 95%区间作为最终点估计的不确定性描述，不作为逐轮停止判据；result 必须披露分阶段停止规则。上述 1,000/200/50,000 预算、BLER band 以及“12 条曲线共同追加”在 formal trial 1 前冻结，之后不得依据结果更改并宣称正式验收。

## 9. 输出、图和 result 要求

正式产物写入

`outputs/experiment035_pdsch_2rx/<run_id>/<stage>/<scenario_id>/`。

至少保存：原 YAML、展开配置和 SHA-256、代码版本/工作区变更标识、candidate/receiver manifest 和 hash、delay audit、interval CSV、estimated/ideal error flags、estimated CE trial arrays、合并 BLER/NMSE CSV、absolute-trial/seed/pairing audit、滤波器诊断、MRC/功率诊断、old/current 时间相关、运行日志、耗时、CPU placement 与 RSS 峰值。

每个场景至少生成三张主图：

1. 六方案 estimated-CSI BLER--SNR；
2. 六方案 ideal-CSI BLER--SNR；
3. 六方案 estimated-CSI data-RE CE NMSE--SNR。

另生成一张同方案 ideal/estimated 目标 SNR 间隔汇总表；只有表格过密时才生成汇总图。同一物理方案在三张图中复用 color/marker，CSI 模式用固定 linestyle 区分。BLER 图使用对数纵轴，零误块点的 CSV 保留真实零值，仅绘图时使用 `0.5/trials` 下界。所有图由本地脚本读取已保存 CSV 生成，图数据另存 CSV，并生成约 13 cm 宽预览供研究者人工核对；Agent 不加载结果图片。

正式完成后成对生成：

- `research/result-035-PDSCH 2Rx仿真.md`；
- `research/result-035-PDSCH 2Rx仿真-text.md`。

两版必须逐项回答第 1 节问题，报告 10%/1% crossing、真实 bracket、错误数、trials、Wilson 区间、bootstrap 区间、CE NMSE 及 95%区间，并明确 ideal CSI、estimated CSI 和 NMSE 的证据边界。结果未经研究者确认前不更新 `KNOWLEDGE.md`/`GOALS.md`，不创建 Git checkpoint。

## 10. 执行顺序与验收条件

1. 研究者确认第 8 节统计预算、停止条件以及“4Tx 和 8Tx 均运行”的范围；
2. 实现公共 trial 核心、035 入口、分析入口和定向测试；
3. 运行新增测试、033 回归 validate 和两个 smoke；smoke 审计失败时暂停，不进入 prescan；
4. 运行两个 2Rx/60 km/h prescan，生成 bracket audit；
5. 把正式 SNR 网格、配置 hash、总预算写入本 plan 冻结区，经研究者确认后开始 formal trial 1；
6. 完成初始预算和 1%端点追加，运行数据一致性与配对审计；
7. 本地生成图、成对 result 和复现命令，交研究者确认。

代码验收要求：测试通过；033 兼容性不变；两场景 smoke 的维度、功率、seed/pairing、数组引用和数值有限性全部通过。实验验收要求：12 条 BLER 曲线均有 10%/1%真实 bracket，需追加端点达到预定错误数或样本上限，6 条 estimated NMSE 曲线与对应 BLER trial 完全对齐，result 的原始数字可从保存的 CSV/NPY 独立复算。

## 11. 正式冻结区

prescan 已于 2026-09-17 完成。首次 `20260917_main` 运行因 estimated error-flag
保存键错误而作废；修复后使用新目录 `20260917_fix1` 重跑 smoke 和 prescan，旧目录仅保留作
审计证据。以下冻结提案已通过 runner `--stage validate`，**尚待研究者在 formal trial 1 前确认**：

- 两个场景共同使用 25 点 SNR 数组
  `[8,8.25,8.5,8.75,9,9.25,9.5,9.75,10,10.25,10.5,10.75,11,11.25,11.5,11.75,12,12.25,12.5,12.75,13,13.25,13.5,13.75,14] dB`。
  修正版 prescan 中所有 estimated/ideal 曲线的 10% 与 1% crossing 均在 8--14 dB 内取得双侧点；
  0.25 dB 全区间网格满足 crossing 邻域最大间隔要求。
- `A100_NT4_NR2_V60`：`configs/plan035_nt4_nr2_v60_formal.yaml`，SHA-256
  `95e30cee2eb1079478a7ad9526505dc686a72d0305b852097ba96ca1c0580b42`。
- `A100_NT8_NR2_V60`：`configs/plan035_nt8_nr2_v60_formal.yaml`，SHA-256
  `c6e4d5fd8a9342deac504eb30e081b832c6de6c12bad030a572adb19db60c0d4`。
- 初始预算为每点 1,000 个共同 paired trials；每个自适应检查/可恢复 interval 为 1,000 trials，
  `batch_size=25`。每场景初始 25,000 trials，两场景合计 50,000 trials。
- 10%/1% bracket 端点的自适应追加规则按第 8.2 节执行：相关 SNR 的 12 条 BLER 曲线共同
  追加至至少 200 个错误块或 50,000 trials 上限。按所有 25 点均达到上限计算的严格上界为每场景
  1,250,000 trials、两场景合计 2,500,000 trials；实际只追加命中的 bracket 端点。
- 4Tx smoke audit：
  `outputs/experiment035_pdsch_2rx/20260917_fix1/smoke/A100_NT4_NR2_V60/smoke_audit.json`，
  SHA-256 `dc33e538887d122c64056f0908c2bdf3a4f7c5592ef6dc0339dac6c611982a3a`。
- 8Tx smoke audit：
  `outputs/experiment035_pdsch_2rx/20260917_fix1/smoke/A100_NT8_NR2_V60/smoke_audit.json`，
  SHA-256 `9fd3bf3be426a639d3c4857b70372fadefb62f66198ac866ae8767acc66b1920`。
- 4Tx prescan interval evidence：
  `outputs/experiment035_pdsch_2rx/20260917_fix1/prescan/A100_NT4_NR2_V60/intervals.csv`，
  SHA-256 `b19f467128134eff7157ff56c67f9f600dd0f646a0a8231761fb897d60eef067`。
- 8Tx prescan interval evidence：
  `outputs/experiment035_pdsch_2rx/20260917_fix1/prescan/A100_NT8_NR2_V60/intervals.csv`，
  SHA-256 `7a1506fbde8041a7b7c4434be58aa790ed4f492a49c43f629a4f5b67facfba79`。
- 研究者确认日期：待确认。

## 12. 增补：CDL-D 对照运行

在 035 链路上增加一组 CDL-D（偏 LoS）对照，使用 `tools/run_plan035_cdl_pdsch_2rx.py`。
除信道由 TDL-A 改为 Sionna CDL-D 外，第 2--10 节的天线数、60 km/h、100 ns RMS delay
spread、六个候选及人工时延、CSI age、DMRS、estimated/ideal 接收机、SNR/trial 规则和交付
产物均不变。CDL 阵列固定为单极化 V、omni、半波长间距的 `1×Nt` Tx ULA 和 `1×Nr` Rx
ULA；真实信道保留 CDL 阵列几何空间相关性，035 接收机仍逐 Rx 使用仅含时频统计的协方差，不利用
空间协方差。新场景 ID 使用 `D100_NT{4|8}_NR2_V60`，输出写入
`outputs/experiment035_cdl_d_pdsch_2rx/`，不得覆盖原 TDL-A 结果。先完成两个场景的 validate 和
smoke；随后按第 8--10 节执行 prescan、冻结正式网格、formal 与同规格 result 交付。
