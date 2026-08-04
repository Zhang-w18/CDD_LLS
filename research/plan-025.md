# plan-025：信道模型与信道估计器重构

## 目标

把当前静态、自写指数 PDP 信道替换为 Sionna 1.0.2 的 3GPP TR 38.901 TDL 信道，并把当前一维频域 RMMSE 扩展为使用匹配时频协方差的二维 RMMSE。

重构后需要：

- 支持 TDL profile、RMS delay spread、载频和 UE 速度配置；
- 支持 8 RB、48 RB，10 个 OFDM symbols，2 个 DMRS symbols；
- 不发送真实 DMRS，直接构造除掉已知 DMRS 符号后的 LS 信道观测；
- known-delay 分支使用 CDD-shifted 频率协方差；
- unknown-delay 分支只使用未加 CDD 的基准 TDL 频率协方差；
- 两个分支使用相同且匹配的时间协方差和噪声方差；
- 使用零速度、极小 delay spread 的 Sionna TDL 重新运行 plan-024 的 QC/Sidon 对比，观察原结论在更实际信道模型下是否保持。

本轮包含 plan-024 条件向近似平坦 TDL 的迁移验证，但不预设必须复现原来的 0.33/0.85 dB 数值。3/60 km/h 下的正式移动 TDL BLER 扫描另行讨论。

> 2026-07-23 修订说明：第二部分第 5 节已改写为“保持 plan-024 物理 CDD 时延与链路定义、仅替换底层信道生成器”的补做实验草案。该补做尚未执行，须经研究者审核确认后才能修改代码或启动仿真。2026-07-22 已完成的 `/4096` 整数采样延迟结果保留为首次执行记录，不作为单独归因于 TDL-A 的证据。

## 一、代码修改计划

### 1. 建立本地运行环境

创建独立虚拟环境，初始固定：

- NumPy 1.26.4
- SciPy 1.15.3
- Matplotlib 3.10.3
- PyYAML 6.0
- Sionna 1.0.2
- TensorFlow 2.15.1

Sionna 1.0.2 的 PHY 后端是 TensorFlow，PyTorch 不是本轮必需依赖。环境建立后保存 Python、TensorFlow、Sionna 和 GPU/CPU 信息。

### 2. 扩展配置

在 `cdd_lls/core/config.py` 增加：

```yaml
channel:
  backend: sionna_tdl
  tdl_profile: A
  delay_spread_ns: 5
  carrier_frequency_hz: 3.5e9
  ue_speed_kmh: 0
  num_sinusoids: 20

resource:
  n_prbs: 8
  pdsch_n_symbols: 10
  dmrs_symbol_indices: [2, 7]
  dmrs_spacing_sc: 6
  dmrs_offset_sc: 0
  cyclic_prefix_length: 288
```

要求：

- `n_prbs` 至少支持 8 和 48；
- `ue_speed_kmh` 支持 0、3 和 60；
- 通用入口支持 8 Tx；
- TDL profile、载频、delay spread、DMRS 图样均保持可配置；
- 保留旧配置字段，避免旧脚本立即失效。

### 3. 把资源栅格改成显式二维时频坐标

修改 `cdd_lls/phy/resource_grid.py`：

- 显式保存 pilot RE 的 `(symbol, subcarrier)` 坐标；
- 显式保存 data RE 的 `(symbol, subcarrier)` 坐标；
- 保留当前 DMRS symbol 上非 pilot RE 继续承载数据的语义；
- 增加 OFDM symbol duration；
- 明确项目有效子载波与 Sionna full-FFT 子载波的映射。

后续信道、估计和数据抽取都必须同时使用 symbol 与 subcarrier 索引，不能再忽略时间维。

### 4. 使用 Sionna 生成 TDL 时频信道

重构 `cdd_lls/phy/channel_tdl.py`：

- 使用 `sionna.phy.channel.tr38901.TDL`；
- 使用 `GenerateOFDMChannel` 生成每个 OFDM symbol 的频域信道；
- 设置 `min_speed=max_speed=ue_speed_kmh/3.6`；
- 输出统一形状：

```text
[batch, n_rx, n_tx, n_symbols, n_active_subcarriers]
```

- 不做每个 realization 的额外归一化；
- 按项目有效子载波位置从 full FFT 中抽取信道；
- Sionna 对象按场景缓存，不在每个 trial 重建。

### 5. 扩展 CDD 等效信道

修改 `cdd_lls/phy/precoding.py`，使 CDD 支持时间维：

```text
H[n_rx, n_tx, n_symbol, n_sc]
    ->
g[n_rx, n_symbol, n_sc]
```

新路径继续使用单位范数 CDD：

$$
C_{k,n}=\frac{1}{\sqrt{N_t}}
\exp\left(-j\frac{2\pi k d_n}{N_{FFT}}\right).
$$

CDD delay 继续以采样点表示，同时在输出中记录换算后的秒数。

### 6. 直接构造二维 LS 观测

修改 `cdd_lls/sim/orchestrator.py`：

在每个配置的 pilot RE 上直接构造

$$
z_{s,k}=g_{s,k}+w_{s,k},
\qquad w_{s,k}\sim\mathcal{CN}(0,N_0).
$$

这表示单位功率已知 DMRS 已被除掉。

两个 DMRS symbols 的观测分别保留，不能再提前平均，也不能再设置 `noise_var_ls=N_0/2`。两个 DMRS 之间的信息由时间协方差和 RMMSE 联合利用。

### 7. 构造匹配的 TDL 时频协方差

重构 `cdd_lls/phy/estimators.py`：

- 使用 Sionna `tdl_time_cov_mat` 构造时间协方差；
- 使用 Sionna `tdl_freq_cov_mat` 构造基准频率协方差；
- 从 full-FFT 协方差中抽取有效子载波部分；
- 使用与仿真器相同的 profile、delay spread、速度、载频和 OFDM symbol duration。

CDD known-delay 频率协方差为

$$
R_f^{known}(k,l)=R_f^{TDL}(k,l)
\frac{1}{N_t}\sum_n
\exp\left[-j\frac{2\pi(k-l)d_n}{N_{FFT}}\right].
$$

unknown-delay 分支使用

$$
R_f^{unknown}=R_f^{TDL}.
$$

unknown-delay 协方差构造接口不接收真实 CDD delay，防止无意泄漏。

### 8. 实现二维 RMMSE

把 pilot 和目标 data RE 按固定 `(symbol, subcarrier)` 顺序展开，计算

$$
\hat{g}_D=R_{DP}(R_{PP}+N_0I)^{-1}z_P.
$$

实现要求：

- known/unknown 使用相同的 `z_P`、`R_t`、`N_0` 和 pilot 坐标；
- 使用 Cholesky 或线性求解，不显式求逆；
- 缓存每个场景和 SNR 的矩阵分解；
- 主 NMSE 在 data RE 上计算；
- 输出 condition number、数值 jitter 和协方差类型；
- 保留旧估计器接口所需的兼容层，但 plan-024 的新验证必须使用重构后的二维估计链路。

### 9. 调整仿真编排和输出

修改 `cdd_lls/sim/orchestrator.py`：

- 数据信道按二维 data RE 坐标抽取；
- known/unknown 分支共享信道和噪声；
- 同时设置 NumPy、TensorFlow 和 Sionna seed；
- 输出实际 profile、速度、载频、delay spread、OFDM symbol duration、DMRS 坐标、CDD delay 和噪声方差；
- 新旧配置使用明确的信道和估计器名称，避免混用归一化方式。

## 二、验证计划

### 1. 基础单元测试

增加以下测试：

- 8 RB 对应 96 个有效子载波，48 RB 对应 576 个；
- 两种带宽均为 10 symbols、两个 DMRS symbols；
- pilot/data 二维坐标不重复且数量正确；
- 3/60 km/h 正确转换为 m/s；
- Sionna TDL 输出维度、能量和 seed 重放正确；
- full-FFT 到有效子载波的索引与相位符号正确；
- CDD 时间维扩展不改变静态单 symbol 结果；
- 两个 DMRS 的 LS 噪声方差均为 `N_0`，没有提前除以 2。

### 2. 协方差验证

使用 TDL-A、100 ns、3.5 GHz，分别测试：

```text
RB:       8, 48
速度:     3, 60 km/h
CDD:      全零、QC、Sidon
```

检查：

- 理论协方差 Hermitian、半正定、对角为 1；
- 全零 CDD 时 known 与基准 TDL 频率协方差一致；
- 60 km/h 的长时间 lag 相关性低于 3 km/h；
- Sionna 信道样本的经验时间/频率相关性与理论值一致；
- QC/Sidon 的经验有效信道协方差与 CDD-shifted 理论协方差一致。

经验相关性的目标误差暂定：代表性相关元素 RMSE 不超过 0.03，最大绝对误差不超过 0.06。若有限样本不足，增加 realization，不调整公式或事后放宽阈值。

### 3. RMMSE 验证

在 8/48 RB、3/60 km/h、多个 SNR 下，用 Monte Carlo NMSE 对照闭式 MSE：

- matched known-delay 的经验 NMSE 与闭式值差不超过 0.2 dB；
- mismatched unknown-delay 的经验 NMSE 与一般线性 mismatch MSE 公式一致；
- known-delay 的期望 MSE 不高于 unknown-delay；
- 全零 CDD 时 known/unknown 的矩阵和输出一致；
- 两分支的时间协方差、噪声方差和 pilot 坐标完全相同。

该阶段只验证信道估计，不运行新的 TDL LDPC BLER。

### 4. 软件兼容性回归

软件回归与后面的物理模型扩展分开执行。软件回归只回答重构是否破坏旧接口和旧结果。

在修改功能代码前，使用新建的同一本地环境保存一个小型回归基准：

- 运行现有全部单元测试；
- 运行 plan-024 E1 的闭式 segment 计算；
- 运行 plan-024 E2 的固定小样本：SNR `14.25,14.5,14.75,16.25,17.25` dB，每点 100 trials，`batch_size=1`；
- 保存 CSV、JSON、测试日志、包版本和命令。

重构完成后，以相同环境、seed、参数和 batch size 回放。要求：

- 旧 YAML 配置仍能加载；
- 旧专题脚本仍能运行；
- 现有单元测试继续通过；
- E1 候选、pilot 数和判定完全一致，NMSE 数值差不超过 `1e-12`；
- E2 每个 SNR、每个候选的 TB 错误数完全一致；
- BLER、CE NMSE 和条件数差不超过 `1e-12`。

旧路径允许通过兼容层继续使用原来的一维静态信道和 DMRS 平均语义，但不能影响新路径。软件回归失败时先修复兼容性，不开始新的 TDL 对比。

### 5. 补做：保持 plan-024 物理 CDD 时延，仅用 Sionna TDL-A 替换底层信道

#### 5.1 研究问题、假设和变量隔离

研究问题：在保持 plan-024 E2 的 QC/Sidon 预编码矩阵、物理 CDD 时延、资源分配、DMRS、功率与噪声定义、编码调制、接收机知识和统计方法不变时，仅把底层平坦分支信道替换为 Sionna 1.0.2 的 3GPP TR 38.901 TDL-A、5 ns、零速度信道，Sidon 相对 QC 的 estimated-CSI BLER 优势是否保持？

可证伪假设：

- `H1`：Sidon 的 10% BLER 目标 SNR 相对 QC 改善不低于 0.15 dB，且保守 95% 区间不支持反向排序；
- `H2`：样本充分时，Sidon 的 1% BLER 目标 SNR 不高于 QC；若目标附近任一候选累计错误块少于 30，只报告先导结果；
- `H3`：补做实验使用的 QC/Sidon 频域预编码矩阵与 plan-024 对应矩阵逐元素一致，最大绝对误差不超过 `1e-12`。若 `H3` 不通过，不得启动正式 BLER 扫描。

本补做只允许一个主动变化：底层分支信道从 plan-024 的频率平坦 `h_n~CN(0,1)` 改为 Sionna TDL-A 5 ns 频率响应 `H_{s,k,n}`。matched 接收机的真实频率协方差必须随物理信道改为 TDL-A 与 CDD 的复合协方差，这是保持“双方均为 known-delay matched/oracle 接收机”的必要变化，不作为第二个独立实验变量。不得同时改变 CDD 相位斜率、归一化、导频处理、数据映射、LLR 噪声定义或 MCS。

#### 5.2 物理时延与相位矩阵的精确定义

固定有效子载波数 `K=576`、FFT 长度 `N_FFT=4096`、子载波间隔 `Delta_f=30 kHz`。plan-024 的 DFT 栅格索引 `j_n` 必须先换算成物理时延和允许为非整数的 FFT 采样延迟：

$$
\tau_n=\frac{j_n}{K\Delta f},
\qquad
d_n=\tau_n N_{FFT}\Delta f
=j_n\frac{N_{FFT}}{K}
=j_n\frac{64}{9}.
$$

频域预编码必须以有效带宽第一个子载波为相位参考，使用 plan-024 的局部有效子载波索引 `m=0,...,575`：

$$
V_{m,n}=\exp(-j2\pi m\Delta f\tau_n)
=\exp\left(-j\frac{2\pi m j_n}{K}\right).
$$

不得把原始 `j_n` 直接代入 `/4096` 公式。若底层接口使用中心化子载波索引，必须补偿每端口的常相位，使最终 `V` 与上述局部索引公式逐元素一致；不能只依赖“统计分布等价”。本补做保留 plan-024 的未归一化矩阵 `|V_{m,n}|=1`，不使用 `1/sqrt(8)` 归一化。


| 候选 | DFT 栅格索引 `j_n` | FFT 采样延迟 `d_n=j_n*64/9`，samples | 物理时延 `tau_n=j_n/(576*30 kHz)`，ns |
|---|---|---|---|
| QC | `[0,9,18,27,36,45,54,63]` | `[0,64,128,192,256,320,384,448]` | `[0,520.833333,1041.666667,1562.5,2083.333333,2604.166667,3125,3645.833333]` |
| Sidon | `[0,1,3,7,12,20,30,65]` | `[0,7.111111,21.333333,49.777778,85.333333,142.222222,213.333333,462.222222]` | `[0,57.870370,173.611111,405.092593,694.444444,1157.407407,1736.111111,3761.574074]` |



输出必须同时保存 `j_n`、`d_n`、`tau_n_s`、相位参考、`K`、`N_FFT` 和 `Delta_f`，不得只保存名为 `delay` 的无单位数组。

#### 5.3 完整仿真条件

| 类别 | 参数 | 取值 | 与 plan-024 的关系 |
|---|---|---|---|
| 信道 | 生成器 | Sionna 1.0.2 `TDL` + `GenerateOFDMChannel` | 唯一主动变化 |
| 信道 | profile / RMS delay spread | TDL-A / 5 ns | 新物理信道条件 |
| 信道 | UE 速度 / 载频 | 0 km/h / 3.5 GHz | 零多普勒；载频仅为 Sionna 必需配置 |
| 信道 | realization 归一化 | 不做 per-realization normalization | 保留 TDL 自然功率波动 |
| 天线 | 发射 / 接收 / 层数 | 8 / 1 / 1 | 相同 |
| 资源 | SCS / FFT / CP | 30 kHz / 4096 / 288 samples | 相同 SCS 和 FFT；CP 供 Sionna OFDM 使用 |
| 资源 | 分配 | 48 PRB，576 active SC，中心连续映射 | 相同 |
| 资源 | PDSCH symbols | 10 | 相同 |
| DMRS | symbols / comb / offset | `[2,7]` / 24 / 0 | 相同 |
| DMRS | pilot RE / data RE | 48 / 5712 | 必须与 plan-024 坐标逐项相同 |
| 候选 | QC / Sidon | 第 5.2 节两组 `j_n` 及其物理时延 | 相同物理时延与相位矩阵 |
| 预编码 | 幅度 | 未归一化，`|V|=1` | 与 plan-024 相同 |
| 接收机 | 知识 | 双方 known-delay matched/oracle | 相同公平性等级 |
| 接收机 | DMRS 处理 | 两个零速度 DMRS 等效平均，平均后 LS 噪声方差 `N0/2` | 与 plan-024 相同，不使用首次执行的 2D 联合接口 |
| 接收机 | 频率估计 | 全带 matched LMMSE；使用各候选真实 TDL-A + CDD 复合协方差 | 与信道匹配所需变化 |
| 链路 | 调制 / MCS / 码率 | 16QAM / MCS 8 / 553/1024 | 相同 |
| 链路 | LDPC / LLR clip | 最多 8 次迭代 / 50 | 相同 |
| 噪声 | 数据噪声方差 | `N0=8/SNR` | 与未归一化 `V` 的 plan-024 定义相同 |
| 噪声 | LLR 有效噪声 | 只使用 `N0`，不加入 CE-error-aware 项 | 相同 |
| 随机性 | seed / batch size | 20260716 / 20 | 复用 plan-024 seed 与正式扫描 batch |

对每个 SNR 和 trial，QC 与 Sidon 必须共享同一 TDL realization、payload、平均后 LS noise 和 data noise。

零速度下两个 DMRS symbol 的底层 TDL tap realization 应相同。实现必须先验证这一点，再按 plan-024 语义生成等效平均 LS 观测：

$$
\bar z_P=g_P+\bar w_P,
\qquad
\bar w_P\sim\mathcal{CN}(0,N_0/2).
$$

matched 频率协方差使用实际 TDL-A 基准协方差与未归一化 CDD 的复合形式：

$$
R_g(k,l)=R_{TDL}(k,l)
\sum_{n=0}^{7}V_{k,n}V_{l,n}^{*}.
$$

#### 5.4 正式扫描前的最小验证

第二部分第 2–4 节以及首次执行已经验证 Sionna TDL 生成器、通用 TDL 协方差、RMMSE 和旧路径回归。本补做不重复 4000-realization 功率检查、经验协方差 Monte Carlo 或多 SNR smoke，只验证本次新增的分数采样延迟及 024 兼容链路：

1. 确定性 preflight：检查两组 `j_n` 到 `d_n/tau_n` 的换算；QC/Sidon `V` 与 `tools/run_experiment024_segment_sidon_qc.py::cdd_V` 的逐元素最大绝对误差不超过 `1e-12`；active、pilot、data 坐标逐项相同且 pilot/data RE 为 48/5712；未归一化 `|V|=1`、数据噪声 `N0=8/SNR`、平均 LS 噪声 `N0/2`；`R_g` Hermitian、半正定、对角为 8，matched LMMSE 矩阵数值有限；
2. 自动测试：运行现有单元测试，并新增分数采样延迟、相位参考和 `/576` 等价测试；不重复已经通过且本次代码未涉及的独立大样本验证；
3. 单点端到端 smoke：SNR `14.5 dB`、20 trials、batch size 20，只验证 TDL、CDD、平均 DMRS、matched LMMSE、LDPC、共同随机数、输出字段和相同 seed 重放一致，不据此下性能结论。

preflight、自动测试或单点 smoke 任一失败即停止，修复后重新执行这三项；全部通过后直接进入粗扫。

#### 5.5 正式预算、统计方法与停止条件

粗扫和精扫沿用 plan-024 的两阶段方法，不设置精扫前的人工确认点：

1. 粗扫：SNR `13.0:0.5:18.0 dB`，每点 400 trials，QC/Sidon 成对运行；
2. 对 10% 和 1% 目标分别使用粗扫的单调化 BLER 序列，为 QC 和 Sidon 找到跨越目标的相邻 0.5 dB 区间；取两个候选区间的并集，以 0.25 dB 间隔生成该目标的精扫网格；
3. 精扫网格在运行任何精扫 trial 前写入 `refinement_grid.json`，随后直接执行，每点 3000 trials，不再等待研究者确认；
4. 若任一目标未在 `13–18 dB` 内形成跨越区间，粗扫按相应方向以 0.5 dB 步长自动扩展，每点仍为 400 trials，最多扩展 2 dB；达到扩展上限仍未跨越时，报告目标超出扫描范围，不对该目标拟合或继续增加预算；
5. 正式点不因结果有利或不利而删除；BLER 非单调点保留并作为有限样本现象报告；
6. 1% 目标附近任一候选累计错误块少于 30 时，只报告先导 1% 结果，不自动增加 trials；10% 主判据不受此条影响。

每个 SNR 点报告 TB 错误数、trial 数、BLER 和 Wilson 95% 区间。10%/1% 目标 SNR 使用预定精扫区间全部点的二项 logit 拟合；目标 SNR 区间使用 delta method。Sidon 改善定义为：

$$
G_{Sidon}=SNR_{QC}-SNR_{Sidon},
$$

正值表示 Sidon 更好。主区间继续采用不利用正配对协方差的保守 95% 近似，同时保存逐点配对四格计数并报告 McNemar 精确检验作为方向性佐证。主结论以目标 SNR 及其预定区间为准，不以 McNemar 汇总替代。

判定分为：

- **保持 plan-024 主结论**：10% 改善至少 0.15 dB，且保守 95% 区间下限不小于 0；
- **方向保持但未达到原判据**：10% 点估计为正，但小于 0.15 dB或区间跨 0；
- **持平/不可判定**：点估计接近 0 且区间跨 0，或目标区间数据不足；
- **显著反向**：10% 改善的保守 95% 区间上限小于 0；
- 1% 结果单独报告，不替代 10% 主判据。

不要求补做结果数值接近 result-024 的 `+0.33/+0.85 dB`；必须报告相对这两个历史值的差，但不得把差异事后归因于未单独验证的机制。

#### 5.6 代码范围、输出和复现入口

审核通过后预计增量修改：

- `cdd_lls/core/config.py`：CDD delay 支持浮点采样延迟并明确单位；旧整数采样配置语义保持兼容；
- `cdd_lls/phy/precoding.py`：支持分数采样延迟和显式相位参考，输出 `j_n/d_n/tau_n_s` 元数据；
- `cdd_lls/phy/estimators.py`：确认 matched TDL + 未归一化 CDD 协方差的幅度标度；
- `tools/run_plan025_delay_matched_tdl.py`：补做实验入口，不覆盖首次执行输出；
- `tools/analyze_plan025_sionna_e2.py`：允许读取补做目录并生成同口径目标拟合；
- `tests/test_precoding.py`、`tests/test_rmmse_time_frequency.py`：增加第 5.4 节对应测试。

固定输出根目录：

```text
outputs/experiment025_sionna_tdl_rmmse/20260723_delay_matched/
```

至少包含：

- `validation/preflight.json`
- `validation/tests.log`
- `smoke/`
- `prescan/sidon_qc_bler.csv`
- `prescan/paired_error_counts.csv`
- `refinement_grid.json`
- `refine_10pct/` 和 `refine_1pct/`
- `final/final_summary.json`
- 每个运行目录下的 `resolved_experiment.json`、`environment.json`、`commands.json` 和 UTF-8 日志。

`sidon_qc_bler.csv` 至少包含 `design`、`snr_db`、`trials`、`tb_errors`、`bler`、`bler_wilson95_lo`、`bler_wilson95_hi`、`ce_nmse_mean`、`ce_nmse_mean_db`、`estimator_condition_number` 和 `estimator_min_singular_value`。`paired_error_counts.csv` 至少包含每个 SNR 的双方都正确、仅 QC 错误、仅 Sidon 错误、双方都错误四格计数。所有 delay 字段必须在 `resolved_experiment.json` 中带单位保存。

计划复现命令；具体参数名可在实现时调整，但最终命令必须保存在 `commands.json` 并同步回写 result：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan025_delay_matched_tdl.py --stage validate --run-id 20260723_delay_matched
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan025_delay_matched_tdl.py --stage smoke --snrs 14.5 --trials 20 --batch-size 20 --run-id 20260723_delay_matched
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan025_delay_matched_tdl.py --stage link --snrs 13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18 --trials 400 --batch-size 20 --run-id 20260723_delay_matched\prescan
```

粗扫完成后由分析入口按第 5.5 节自动生成并保存 `refinement_grid.json`；精扫入口读取该文件直接运行，不接受人工临时修改网格。最终实际命令必须保存到 `commands.json`。

#### 5.7 result 必须回答的问题

补做完成后，成对更新 `research/result-025.md` 和 `research/result-025-text.md`，两版结论和数字必须一致，并至少回答：

1. `V` 是否与 plan-024 逐元素一致，物理时延与分数采样转换是否通过；
2. 除 Sionna TDL-A 5 ns 信道及其 matched 协方差外，是否存在任何未预定的配置变化；
3. QC/Sidon 的 10% 和样本充分时的 1% 目标 SNR、错误计数、置信区间与 Sidon 改善；
4. 导频矩阵秩、条件数、最小奇异值，以及 data-RE CE NMSE 差；
5. 相对 result-024 `+0.33/+0.85 dB` 和首次 `/4096` 执行 `-0.161/-0.582 dB` 的差别；
6. `H1/H2/H3` 和第 5.4 节验证项逐条是否通过；
7. 异常、样本不足、适用范围、完整证据路径和复现命令。

## 三、预计修改文件

```text
requirements-sionna102.txt
cdd_lls/core/config.py
cdd_lls/phy/channel_tdl.py
cdd_lls/phy/resource_grid.py
cdd_lls/phy/precoding.py
cdd_lls/phy/estimators.py
cdd_lls/sim/orchestrator.py
configs/smoke_sionna_tdl_rmmse.yaml
tools/validate_tdl_rmmse.py
tools/replay_plan024_regression.py
tools/run_plan024_sionna_tdl.py
tests/test_channel_tdl.py
tests/test_resource_grid_time_frequency.py
tests/test_rmmse_time_frequency.py
```

## 四、执行顺序

1. 建立并验证本地 Sionna 1.0.2 环境。
2. 运行现有测试并建立旧路径小型回归基准。
3. 修改配置和二维资源栅格。
4. 接入 Sionna TDL。
5. 修改 CDD、LS 观测和数据 RE 抽取。
6. 实现 TDL 协方差和二维 known/unknown RMMSE。
7. 运行单元测试、协方差验证和 RMMSE 验证。
8. 回放旧路径回归基准，确认软件兼容性通过。
9. 研究者审核并确认第二部分第 5 节补做规格、预算和停止条件；确认前不得执行后续步骤。
10. 实现分数采样延迟、相位参考和补做入口，不改变旧整数采样配置的既有语义。
11. 完成第 5.4 节全部验证、单元测试和 smoke；任一失败则停止并修复。
12. 运行预定 400-trial 粗扫，保存完整配置、计数和日志。
13. 按粗扫结果自动生成并固化 0.25 dB 精扫网格，随后运行 3000-trial 精扫，不设置人工确认点。
14. 成对更新 result-025 两版，逐项回答第 5.7 节问题；研究者确认结果后再决定是否扩展到更大 delay spread。

## 五、2026-07-22 首次执行范围确认（已完成，保留为历史记录）

研究者确认本轮正式链路扫描只执行 plan-024 E2 的 QC/Sidon delay 对比：

- 只运行 known-delay matched/oracle 二维 RMMSE；不运行 unknown-delay 正式 BLER；
- 粗扫复用 plan-024 的 `13.0:0.5:18.0 dB`，每点 400 trials；
- 粗扫先确认 5 ns TDL-A 下 10% 和 1% BLER 交叉区域；
- 若交叉仍落在原区间，10% 区间使用 `14.25,14.5,14.75 dB`，每点 3000 trials；
- 若交叉仍落在原区间，1% 区间使用 `16,16.25,16.5,16.75,17,17.25,17.5 dB`，每点 3000 trials；
- 若 TDL-A 使交叉移出上述精化区间，先报告粗扫结果，再由研究者确认新的精化点；
- 1% 目标附近任一候选累计错误块少于 30 时，只报告先导结果，不做确定性结论。

首次执行把 QC/Sidon 的整数数组解释为 FFT 采样点并使用 `/4096` 相位分母，因此没有保持 plan-024 的物理时延和预编码矩阵。对应结果只说明该联合配置下的性能，不能回答“仅替换为 TDL-A 5 ns 后 plan-024 排序是否保持”。第二部分第 5 节的补做草案用于消除该混杂。

## 六、2026-07-23 补做审核状态

- 状态：研究者已在 2026-07-23 明确要求按第二部分第 5 节开始实验，视为审核通过；补做已执行完成；
- 已确认方向：正式扫描前采用最小验证；粗扫和精扫沿用 plan-024 两阶段方法；精扫网格由粗扫自动生成，不再单独请求研究者确认；
- 执行回执：preflight、自动测试和 smoke 通过后，已完成 400-trial 粗扫、自动精扫网格固化、每点 3000-trial 精扫和结果分析；
- 结果位置：`research/result-025.md`、`research/result-025-text.md` 和 `outputs/experiment025_sionna_tdl_rmmse/20260723_delay_matched/`；
- 后续状态：等待研究者审核补做结果；确认前不更新 `KNOWLEDGE.md` 和 `GOALS.md`。
