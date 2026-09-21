# plan-036：PDCCH 4T4R、1-symbol CORESET 的 CDD 时延与 precoder-cycling 基线

> 状态：**已确认，待执行**。研究者于 2026-09-18 确认 AL2 不使用原 8 元 Sidon 的前四项，而从中选择 `[1,12,20,65]`，消除公共循环移位后冻结为 `[0,11,19,64]`。本 plan 已在第 4 节记录选择规则、精确时延和可辨识性审计。

## 1. 目的与交付目标

本轮比较 1-symbol PDCCH 中三类发射方案相对 4Tx DFT precoder cycling 基线的 DCI BLER：严格 Sidon CDD、指定小时延 QC CDD，以及 transparent DFT cycling。每个物理发射波形同时运行 estimated CSI 和 ideal CSI；estimated 分支另报告 data-RE 信道估计 NMSE。只运行 AL1 和 AL2，并分别覆盖最大 Doppler 5 Hz 与 1100 Hz。

每个 `(AL, Doppler)` 的 estimated-CSI 交付四条 BLER 和四条 CE NMSE 曲线：Sidon non-transparent、QC-delay non-transparent、同一 QC-delay transparent、DFT4 cycling transparent。ideal-CSI 交付三条 BLER 曲线：Sidon、QC-delay 和 DFT4 cycling。QC transparent/non-transparent 的逐 RE 发射矩阵完全相同，ideal CSI 下不存在接收机协方差知识差异，因此只运行和绘制一条 `QC-delay ideal-CSI`，不得把同一波形重复计为两个独立结果。

报告完整 `-10:1:5 dB` 原始 BLER、estimated-CSI data-RE CE NMSE、相对同 CSI 模式 DFT baseline 的 10%/1% BLER 目标 SNR 差，以及同一波形的 estimated-to-ideal 接收代价和统计区间。正增益统一定义为

$$
G_{p,c}=\mathrm{SNR}_{p,\mathrm{DFT\ baseline},c}-\mathrm{SNR}_{p,\mathrm{candidate},c},
\qquad p\in\{0.1,0.01\}.
$$

其中 $c\in\{\mathrm{estimated},\mathrm{ideal}\}$。同一物理波形的信道估计代价定义为

$$
\Delta_{\rm CE,p}
=\mathrm{SNR}_{p,\mathrm{estimated}}
-\mathrm{SNR}_{p,\mathrm{ideal}},
$$

正值表示 estimated CSI 需要更高 SNR。

本轮不搜索新时延，不运行 AL4/8，不引入空间相关、CFO、ICI、定时误差或多层检测。结论只适用于表列信道、资源、接收机知识和独立 4Rx 模型。

## 2. 冻结系统条件

| 参数 | 冻结取值 |
|---|---|
| SCS / FFT / CP | 30 kHz / 4096 / 288 samples |
| CORESET | 48 RB，duration 1 symbol，first CCE 0 |
| CCE-to-REG mapping | non-interleaved |
| REG bundle | `L=6 REG`，在本场景等于 `6 RB × 1 symbol` |
| Aggregation level | AL1、AL2 |
| DCI / CRC | 40 information bits / CRC24C 24 bits；RNTI `0xFFFF` |
| 调制与编码 | QPSK；NR PDCCH Polar，list size 8 |
| DMRS | PDCCH comb-4，频域密度 25%；data/DMRS 使用相同预编码 |
| 信道 | Sionna 1.0.2 TDL-C，RMS delay spread 300 ns，4 GHz，20 sinusoids |
| Doppler | 最大 Doppler 5 Hz、1100 Hz；配置中保存对应速度及回算误差 |
| 天线与层数 | 4Tx / 4Rx / 1 layer |
| Rx 空间模型 | 四个 Rx 分支独立同分布、同 PDP、噪声独立；无空间相关 |
| 接收模式 | estimated CSI：REG-bundle-based RMMSE；ideal CSI：直接使用真实 data-RE 等效信道 |
| SNR | `[-10,-9,...,5] dB`，固定 1 dB 间隔，不按结果增删点 |
| 功率 | 每个 RE 的 4Tx 总发射功率为 1 |

配置接口当前以速度而非 Doppler 为输入。应使用同一光速常数把 $f_D$ 换算为 $v=f_Dc/f_c$，约为 1.35 km/h 与 296.79 km/h；展开配置必须同时保存请求的 Doppler、实际速度、回算 Doppler和所用光速常数，不得只保存四舍五入速度。

4Rx 的横轴 SNR 定义为单位总发射功率相对单个 Rx 分支噪声功率，故每个分支的复高斯噪声方差为

$$
\sigma_n^2=10^{-\mathrm{SNR}_{\rm dB}/10},
$$

不因 `n_rx=4` 再除以 4。四个分支分别估计等效信道后，在 data RE 上做 coherent MRC。

## 3. 资源计数与相位定义

一个 PDCCH candidate 含 $6\,\mathrm{AL}$ 个 REG。duration 为 1 symbol 时，每个 REG 对应 1 RB，因此

$$
N_{\rm RB,occ}=6\,\mathrm{AL},\qquad
K=72\,\mathrm{AL},\qquad
N_p=K/4.
$$

| AL | 占用 RB | $K$ | $q_K=1/(K\Delta f)$ | DMRS RE | data RE | coded bits $E$ | 6-REG bundle 数 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 6 | 72 | 462.962963 ns | 18 | 54 | 108 | 1 |
| 2 | 12 | 144 | 231.481481 ns | 36 | 108 | 216 | 2 |

CDD 使用 candidate 内从 0 开始的局部子载波索引 $k$：

$$
V_{k,n}=\frac12\exp\left(-j2\pi k j_n/K\right),
\qquad
j_n=\tau_nK\Delta f.
$$

人工时延是数字线性相位，不并入 TDL-C 的真实传播时延，也不改变 CP/ISI 模型。每行满足 $\sum_n|V_{k,n}|^2=1$。

## 4. 候选与时延表

### 4.1 AL1

| family / candidate | 接收口径 | $\mathbf j$ | 物理人工时延 / ns |
|---|---|---|---|
| `SIDON` / `SIDON_AL1_NT` | non-transparent | `[0,1,3,7]` | `[0,462.963,1388.889,3240.741]` |
| `QC-delay-112` / `QC_DELAY_112_NT` | non-transparent | `[0,0.24192,0.48384,0.72576]` | `[0,112,224,336]` |
| `QC-delay-112` / `QC_DELAY_112_TRANSPARENT` | transparent | 同上 | 同上 |
| `DFTcodebook` / `DFT4_CYCLING_BASELINE` | transparent，逐 bundle | 不适用 | DFT4 index `[0]` |

AL1 的 Sidon 值等于 plan031 的 2-symbol AL2、4Tx `S0_SIDON` 物理时延，因为两者均有 $K=72$。AL1 只有一个 6-REG bundle，所以 DFT baseline 实际固定使用单位范数 DFT4 第 0 列 `[1,1,1,1]/2`；名称保留 precoder cycling，但本 AL 内没有发生 bundle 间切换。

### 4.2 AL2

研究者给出的完整 Sidon 输入为

`[0,231.481,694.444,1620.370,2777.778,4629.630,6944.444,15046.296] ns`，

即 8 元坐标 `[0,1,3,7,12,20,30,65]`。该集合的 36 个无序二元和在模 $K=144$ 下全部唯一，且模 $N_p=36$ 的 8 个 fold residue `[0,1,3,7,12,20,30,29]` 全部互异，因此任意四元素子集都继承严格 Sidon 和基本导频可辨识性，但不同子集的人工时延 RMS 与 fold 几何不同。

研究者确认从原集合选择 `[1,12,20,65]`。按公共循环移位等价性减去 1，正式配置只使用规范化坐标

$$
\mathbf j_{\rm SIDON,AL2}=[0,11,19,64].
$$

这里的公共移位消除不改变人工时延 RMS、严格 Sidon、pair-sum 圆周间距或 fold 圆周间距。对应精确物理人工时延为

$$
\boldsymbol\tau
=\frac{[0,11,19,64]}{144\times30\,\mathrm{kHz}}
=[0,2546.296296,4398.148148,14814.814815]\ \mathrm{ns}.
$$

正式表格和图例显示到 3 位小数，但 YAML、manifest 和审计回执必须保存上述整数坐标，并由坐标反算双精度 ns，不得把 3 位小数重新作为相位输入。该候选的 equal-power 人工时延 RMS 为 `5633.377 ns`；fold residues 为 `[0,11,19,28]`，最小圆周间距为 8；10 个无序二元和为 `[0,11,19,22,30,38,64,75,83,128]`，模 144 全部唯一，最小圆周间距为 3。完整 comb-4 导频矩阵应为 rank 4、condition number 1。

| family / candidate | 接收口径 | $\mathbf j$ | 物理人工时延 / ns |
|---|---|---|---|
| `SIDON` / `SIDON_AL2_NT` | non-transparent | `[0,11,19,64]` | `[0,2546.296,4398.148,14814.815]` |
| `QC-delay-96` / `QC_DELAY_96_NT` | non-transparent | `[0,0.41472,0.82944,1.24416]` | `[0,96,192,288]` |
| `QC-delay-96` / `QC_DELAY_96_TRANSPARENT` | transparent | 同上 | 同上 |
| `DFTcodebook` / `DFT4_CYCLING_BASELINE` | transparent，逐 bundle | 不适用 | DFT4 index `[0,1]` |

`QC-delay-112` 只用于 AL1，`QC-delay-96` 只用于 AL2；图例、CSV 和 manifest 必须保留这两个精确 family 名。QC 的 transparent/non-transparent 两条曲线必须使用完全相同的逐 RE 发射矩阵，唯一有意差异是接收机假设的协方差。

## 5. 接收机口径

- Sidon non-transparent：`receiver_covariance_mode: matched_effective`。接收机知道真实 CDD $V$，用 $R_g=R_{phy}\odot(VV^H)$ 在整个 candidate 占用带宽构造频域 RMMSE。
- QC non-transparent：同样使用 `matched_effective` 和真实等效 PDP/协方差。
- QC transparent：`receiver_covariance_mode: physical_fullband`。接收机不知道人工 delay 或 $V$，只用原始 TDL-C 物理 PDP 推导的 $R_{phy}$，但仍在整个 candidate 带宽估计。
- DFT cycling baseline：`receiver_covariance_mode: physical_prg`。接收机只用原始 $R_{phy}$，并严格限制在每个 6-REG bundle，即 72 个连续子载波内构造 RMMSE；不得跨 bundle 插值或联合求逆。
- ideal CSI：在每个 data RE、每个 Rx 分支直接令 $\widehat g_r=g_r$，再用真实等效信道做 4Rx coherent MRC；DMRS RE 仍保留且不承载数据，但不参与信道估计，不构造 RMMSE filter。ideal CSI 只改变接收端，不改变发射预编码、物理信道、DMRS 开销或噪声。

CE NMSE 只对 estimated CSI 具有统计意义。每个 trial 在四个 Rx 和全部 data RE 上先在线性域计算

$$
\mathrm{NMSE}_t=
\frac{\sum_r\sum_{k\in D}|\widehat g_{t,r,k}-g_{t,r,k}|^2}
{\sum_r\sum_{k\in D}|g_{t,r,k}|^2},
$$

再跨 trial 在线性域平均并转换为 dB。ideal CSI 的逐元素误差必须精确为 0，保存为实现审计，但不把恒零值包装成独立 Monte Carlo NMSE 曲线，也不与 estimated NMSE 做 dB 差。

因为 CORESET 只有一个 symbol，RMMSE 只在频域使用同一 symbol 的 DMRS/data 坐标，不存在跨 symbol 时间插值。现有链路是 OFDM-symbol 采样的频域乘法模型，不模拟 symbol 内信道演化和 ICI；因此 5 Hz 与 1100 Hz 在单 symbol 上具有相同边缘分布，理论上不应产生系统性 BLER 差异。两档 Doppler仍按要求运行，但只能作为模型不变性检查，不能据此声称 1100 Hz 移动性无损。若研究目标是观察高 Doppler 的时间选择性损失，需要另建多 symbol 或含 symbol 内时变/ICI 的计划。

## 6. 公平性、随机性与统计预算

每个 `(AL, Doppler, SNR, absolute trial)` 的所有 candidate 和 CSI 模式共享底层 4Tx×4Rx TDL realization、payload、编码比特和 data AWGN；四条 estimated 曲线还共享 DMRS AWGN。QC transparent/non-transparent 必须逐 trial 共享完全相同的发射波形和样本，使差值只反映接收协方差知识。同一物理方案的 estimated/ideal 分支也必须共享发射波形、信道、payload 和 data AWGN；ideal 分支不得重新生成信道或 data noise。两个 Doppler 场景分别统计，不声明 paired 方差缩减。

SNR 网格固定为 16 点 `-10:1:5 dB`。每个 `(AL,Doppler,SNR)` 以 1,000 个共同 trial 为检查步长，最少 2,000、最多 50,000 个共同 trial。达到最少样本后，仅当四条 estimated BLER 和三条不重复的 ideal BLER 各自满足下列任一条件时停止该点：

1. 累计至少 200 个 DCI block errors；
2. Wilson 95% 上限已经低于 0.01，可将该点可靠归类为低于 1% BLER；
3. 达到 50,000 trials。

这样在真实 BLER 约 0.01 时通常需要约 20,000 trials，而高 BLER点和明显低于 0.01 的点可以较早停止。停止判断只能在预定 1,000-trial 边界执行，并保存每轮 `trials/errors/Wilson interval/stop_reason`；不得依据方案排序修改阈值。全部 estimated/ideal 曲线在同一场景/SNR 使用相同 trial 数，以保留 paired comparison。estimated CE NMSE 使用完全相同的 trial，不单独追加 NMSE-only trial。

单点报告原始 BLER 与 Wilson 95% 区间。10%/1%目标只允许在相邻的 1 dB 原始点形成双侧 bracket 时，按 `log10(BLER)` 线性插值；不做单调修正或外推。相对 baseline、QC transparent/non-transparent 差值及同波形 estimated/ideal 差值使用 absolute-trial paired bootstrap 4,000 次；若自举样本中有效双侧 bracket 少于 95%，只报告原始 bracket，不给目标差值区间。

“看到 0.01 BLER”的验收定义为每条曲线在 `[-10,5] dB` 内至少有一个原始点满足 `BLER<=0.01`，且 Wilson 区间、错误数和 trials 均已保存。若 5 dB 仍未达到，不扩大 SNR 范围、不外推，明确标记该曲线未满足本轮范围内的 1%验收。

## 7. 实现、测试与 smoke

现有 `cdd_lls/sim/pdcch_cdd.py` 已支持任意正整数 `n_rx` 的独立分支 MRC，但场景校验仍限制为 `8Tx/1-symbol/AL2,4,8` 或 `4Tx/2-symbol/AL1,2,4`。本轮应扩展同一公共实现以接受 `4Tx/4Rx/1-symbol/AL1,2/TDL-C300`，不得复制链路主循环。

预计新增 `configs/pdcch_plan036_{al1,al2}_{fd5,fd1100}_formal.yaml` 和薄编排/分析入口；若现有 plan031 runner 的 adaptive stop 能直接复用，只扩展 schema、配置与分析，不新增专用 trial 核心。定向测试至少覆盖：

1. AL1/2 的 REG、RB、$K$、DMRS、data 和 $E$ 计数；
2. 4Tx Sidon、QC ns↔坐标换算、fold residue、严格 Sidon pair-sum 和逐 RE 单位总功率；
3. AL1 `[0]`、AL2 `[0,1]` 的 DFT4 bundle 映射，data/DMRS 在 bundle 内使用同一向量；
4. 三种 covariance mode 的滤波器分别等于直接构造的 matched-effective、physical-fullband 和 72-subcarrier physical-PRG RMMSE；ideal 分支不得构造或调用 RMMSE filter；
5. 4Rx 分支独立、每分支噪声方差不缩放、MRC 分母和 LLR 有效噪声正确；
6. 同一物理波形的 estimated/ideal 分支共享 channel/payload/data-noise key，ideal 分支逐元素 $\widehat g=g$、CE NMSE 精确为 0；
7. QC transparent/non-transparent 在 ideal CSI 下的 LLR、判决和 error flag 逐 trial 完全相同，并只输出一条正式 ideal 曲线；
8. A=40/CRC24C Polar 的无噪声编解码、所有 candidate 与 CSI 模式在高 SNR 无 NaN/Inf；
9. adaptive stop 的最少/错误数/Wilson 上限/上限样本三类路径及断点续跑；
10. 原 plan031 的 8Tx/1-symbol 和 4Tx/2-symbol validate 回归不变。

smoke 对四个 `(AL,Doppler)` 场景各运行 `[-10,5] dB × 20 shared trials`，同时覆盖 estimated/ideal 分支，只检查维度、资源、时延、协方差模式、4Rx MRC、CSI 共享键、ideal 精确 CSI、seed 对齐、输出字段和数值有限性，不用于性能判断。另以相同底层随机键比较 5/1100 Hz 的单-symbol 边缘统计；若出现超出数值/Monte Carlo 波动的系统差异，应先审计信道生成器的时间采样语义，不进入正式运行。

## 8. 输出、结果与执行顺序

输出写入 `outputs/experiment036_pdcch_4t4r_1symbol/<run_id>/`，至少保存原 YAML、展开配置和 SHA-256、代码版本/工作区标识、candidate/CSI manifest、精确 delay/坐标/DFT 映射、Doppler↔速度回执、estimated/ideal 逐 trial error flags、estimated CE trial NMSE、ideal 精确零误差审计、逐点 CSV、adaptive status、seed/pairing audit、pilot rank/condition、解析 zero-noise CE floor、日志、耗时和环境信息。

本地分析脚本为每个 `(AL,Doppler)` 生成三张主图：四方案 estimated-CSI BLER、三个不重复物理波形的 ideal-CSI BLER、四方案 estimated-CSI CE NMSE。ideal CSI 的 NMSE 恒为 0，只进入审计表，不单独绘图。同一物理 family 跨图复用颜色，transparent/non-transparent 和 CSI 模式用固定线型区分。图只读取已保存的本地 CSV，零错误点保留原始零值，绘图时才使用 `0.5/trials` 下界。按规范生成约 13 cm 预览，但 Agent 不加载结果图片。

正式完成后成对生成 `research/result-036-PDCCH-4T4R-1symbol.md` 与 `research/result-036-PDCCH-4T4R-1symbol-text.md`。两版必须回答：estimated/ideal CSI 下 Sidon 和 QC 相对各自 DFT baseline 的 10%/1%差异；同一 QC 波形的透明性代价；各物理波形的 estimated-to-ideal 目标 SNR 代价；5/1100 Hz 是否仅表现为单-symbol 模型下的统计不变；每条 BLER 曲线是否在规定 SNR 范围内取得原始 `BLER<=0.01` 点；estimated CE NMSE、pilot rank/condition 和 zero-noise floor 是否解释了异常排序。结果确认前不更新 `KNOWLEDGE.md`/`GOALS.md`，不创建 Git checkpoint。

执行顺序：

1. 按已确认的 AL2 Sidon `[0,11,19,64]` 扩展 estimated/ideal 校验、配置、编排与分析入口，运行定向测试和 plan031 回归；
2. 运行四场景 smoke 并审计 Doppler 语义；
3. 冻结 estimated/ideal 正式配置、hash、seed namespace、跨 CSI 共享规则和预算后执行 `-10:1:5 dB` 自适应 trial；
4. 完成 trial/seed/停止条件/数值一致性审计，本地生成图和成对 result；
5. 交研究者确认，未达 1% 的曲线作为范围内负结果保留，不后验扩网格。
