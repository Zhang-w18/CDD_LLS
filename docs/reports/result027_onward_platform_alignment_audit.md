# result-027 起仿真平台与 target-BLER SNR 对齐审计

## 1. 范围与结论

本审计覆盖 `research/result-027` 至 `research/result-036` 及其 plan、展开配置、候选 manifest、目标门限 CSV 和相关实现。`result-033` 仍是部分正式结果，`result-036` 只完成 5 Hz 缩减范围；二者不能作为最终完整实验使用。

结论如下：

1. 未发现同一入口内信道生成、CDD 相位、发射功率、噪声注入、MRC、解调或译码的明显实现错误。定向回归测试在工作区内以独立 pytest 临时目录运行，结果为 `98 passed, 14 warnings`；warning 均为 Matplotlib/PyParsing deprecated API。
2. 各实验并非使用完全相同的仿真平台口径。除用户预期的信道估计器和 SNR 数值归一化外，target-SNR 提取方法、PDCCH payload、资源映射和随机流配对方式也发生过变化。CSV 已把这些字段显式列出。
3. PDSCH 的两种功率数值归一化在物理 SNR 上等价：`||W||^2=8, noise=8/SNR` 与 `||W||^2=1, noise=1/SNR` 都定义为“每 RE 总发射功率 / 单个 Rx 分支噪声方差”。多 Rx 不把噪声除以 Rx 数。
4. 相同名称不保证相同时延集合。必须同时比较 `n_tx`、`delay_coordinates`、`delay_ns` 和相位分母。特别是 `S0_SIDON`、`B0_QC`、small-delay/QC 和不同 AL 下的严格 Sidon，名称相同但集合可能不同。
5. 没有发现“完整仿真设置相同、时延集合相同、但由两次独立正式运行得到互相矛盾 target SNR”的实例。现有完全相同的数据主要是只读复用或断点续跑，不是独立重复实验；因此可以确认复用一致性，但不能据此估计跨独立重复运行的再现误差。

机器可读汇总为 `docs/reports/result027_onward_target_bler_snr_gain.csv`。该表共 552 行、40 个比较组；每行是一条方案在一个 target BLER 下的记录。未闭合目标保留空值，不外推。

## 2. 平台分支与实现对齐

| 范围 | 信道与发射端 | 接收端 | 解调/译码 | SNR 口径 | 对齐判定 |
|---|---|---|---|---|---|
| result-027/028 PDSCH static | Sionna TDL-A；8Tx/1Rx；CDD `exp(-j2πkj/576)`；预编码平方范数 8 | 两个静态 DMRS 观测平均后做频域 LMMSE；另有 ideal、transparent PRG、transparent/matched small-delay 分支 | 16QAM；NR MCS 8；Sionna LDPC，最多 8 次迭代 | `noise=8/SNR` | 主链路实现一致；028 的 A30/A100 主曲线直接复用 027 原始数据，A300 和补充分支是新增场景/接收机 |
| result-032 PDSCH mobility | 与 result-028 相同的 8Tx CDD 集合和相位定义；TDL-A 100 ns、60 km/h | 两个 DMRS 不平均，做二维时频 RMMSE；CDD matched，PRG/MRT transparent | 与 027/028 相同 | `noise=8/SNR` | 发射、调制、译码和 SNR 物理口径对齐；估计器因移动性有意改变 |
| result-033/035 PDSCH multi-Rx | TDL-A；单位范数预编码；4/8Tx、2/4Rx；相同 8Tx 名义候选沿用相同坐标 | 各 Rx 独立二维时频 RMMSE 或 ideal CSI，随后 coherent MRC | 与 027/028/032 相同 | `noise=1/SNR`，每 Rx 分支相同 | 与 result-032 物理 SNR 等价；多 Rx 和 4Tx 子集是不同场景。033 尚未完成，035 已完成 |
| result-029/030 PDCCH | TDL-A 100 ns、8Tx/1Rx、2-symbol；单位范数 DFT8 bundle cycling | bundle 内跨两个 symbol 的 LS，再做频域线性插值 | QPSK；A=41；CRC24C；Polar list 8；030 为 AL8 repetition | `noise=1/SNR` | 029/030 内部对齐；与 031 的 LMMSE 接收机不对齐 |
| result-031 A100 PDCCH | TDL-A 100 ns、8Tx/1Rx、1-symbol；CDD/PRG | CDD matched frequency LMMSE；PRG physical-PRG；small-delay 有 matched/transparent | QPSK；A=41；Polar list 8 | `noise=1/SNR` | 发射/编译码实现沿用 029/030；资源、符号数和 CE 方法不同，不能把跨实验差值归因于 CDD 单因素 |
| result-031 C300 PDCCH | TDL-C 300 ns、4Tx、2-symbol、3 km/h；1Rx/2Rx | 每 Rx 二维时频 LMMSE 或 ideal CSI，随后 MRC | QPSK；A=41；Polar list 8 | `noise=1/SNR`，不随 Rx 缩放 | 同一 C300 分支内对齐；1Rx 数据在后续补充中复用。AL、Rx 和 receiver mode 必须分组比较 |
| result-034 PDCCH | TDL-C 300 ns、4Tx/4Rx、2-symbol、3 km/h；固定 DFT0 | 每 Rx 二维时频 LMMSE，随后 4Rx MRC | QPSK；A=40；Polar list 8 | `noise=1/SNR` | 与 031 C300 的信道/MRC 框架一致，但 payload 从 41 改为 40，且只跑固定 DFT0，不能直接并入 031 target 比较 |
| result-036 PDCCH | TDL-C 300 ns、4Tx、4Rx/2Rx、1-symbol、5 Hz；CDD/QC/DFT4 | 每 Rx 单-symbol frequency LMMSE 或 ideal CSI，随后 MRC | QPSK；A=40；Polar list 8 | `noise=1/SNR` | 与 034 的 payload 和 SNR 对齐，但 symbol 数、Doppler、DMRS、候选集合不同；1100 Hz 未运行 |

## 3. 需要特别注意的不对齐项

### 3.1 信道估计器

信道估计器不是只增加了多 Rx：

- result-027/028 static PDSCH：两个 DMRS 观测先平均，再做一维频域 LMMSE；
- result-032/033/035 mobility PDSCH：保留两个 DMRS 的独立观测，使用二维时频 RMMSE；
- result-029/030 PDCCH：bundle LS 与频域线性插值；
- result-031 A100 1-symbol PDCCH：频域 LMMSE；
- result-031 C300 与 result-034：二维时频 LMMSE；
- result-036：单 symbol 频域 LMMSE；
- transparent、matched、ideal 也是不同接收机知识口径，不能只按发射时延集合合并。

### 3.2 SNR 与功率归一化

固定 8Tx PDSCH 入口使用 `||W||²=8`、`noise=8/SNR`；PDCCH、plan-033/035 和通用多 Rx PDSCH 使用 `||W||²=1`、`noise=1/SNR`。二者的物理定义一致。当前多 Rx 入口保持每个 Rx 分支噪声方差不变，MRC 后产生阵列与接收分集增益，因此跨 Rx 数的 target SNR 不能解释为发射端 CDD 单独增益。

### 3.3 target SNR 提取

- result-027/028 主 `target_summary.csv` 使用局部 logistic fit；
- result-031 使用 Jeffreys 平滑、按 trial 加权的单调递减 isotonic 曲线，再在局部按 `log10(BLER)` 插值；
- result-032/033/035/036 使用相邻真实点的 `log10(BLER)` 插值，不平滑、不外推；
- result-029/030 和 result-034 原 result 没有发布统一的 crossing 表。汇总 CSV 对 029/030 及 034 使用相邻原始点补算，其中 034 明确标为 `derived_in_audit`。

因此，小于约一个采样间隔的跨 result 差异可能同时包含门限提取方法差异。严格复现实验应固定同一个 crossing 算法后重新分析原始 BLER 点。

### 3.4 时延方案身份

CSV 同时保存 `delay_coordinates` 和 `delay_ns`。判断“同一个 CDD 方案”时采用如下身份键：

`(n_tx, occupied bandwidth/SCS, phase denominator, ordered delay_coordinates or delay_ns, receiver_knowledge)`。

只比较 candidate 名称会产生以下误合并：

- 8Tx `S0_SIDON=[0,1,3,7,12,20,30,65]` 与 4Tx `S0_SIDON=[0,1,3,7]`；
- result-028 的 8Tx small-delay `[0,14.468,...,101.273] ns` 与 4Tx C300 small-delay `[0,14.468,28.935,43.403] ns`；
- result-036 AL1/AL2 的 Sidon 和 QC 集合彼此不同；
- result-031 不同 AL 下，固定整数坐标与固定物理时延是两种不同冻结规则。

## 4. CSV 口径

基线按用户指定规则选择：同一 `comparison_group + target_bler` 内优先选择 precoder cycling；若不存在则选择 B0QC；两者都不存在则标记 `NA`。增益定义为

`gain_vs_baseline_db = baseline target SNR - candidate target SNR`，

正值表示 candidate 达到同一 BLER 所需 SNR 更低。若基线曲线存在但 target 未闭合，保留基线 ID，增益为空，不使用其他候选替代。

重要状态：

- `bracketed`：正式目标已闭合；
- `preliminary_bracketed`：result-033 初始结果，尚未完成端点追加和 paired bootstrap；
- `unbracketed` / `unqualified`：不满足原 result 的 bracket 判据，target 与增益为空；
- `derived_in_audit`：原 result 未发布 crossing，本审计只用相邻原始点补算；
- result-036 的全部行仍缺少原 plan 要求的 paired bootstrap，只有点估计。

生成脚本为 `tools/build_target_bler_gain_table.py`。脚本只读取已保存的本地 CSV、JSON、YAML，不读取结果图片。
