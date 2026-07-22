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

### 5. 用近似平坦 Sionna TDL 重新运行 plan-024 对比

这里的“复现 plan-024”是复现其比较问题和仿真条件，不是复现旧平坦信道代码或要求逐比特一致。

使用重构后的统一信道和估计链路，配置：

```text
信道：          Sionna TDL-A
RMS delay spread：5 ns
UE 速度：       0 km/h
载频：          3.5 GHz
发射/接收：     8 Tx / 1 Rx
资源：          48 RB × 10 symbols
DMRS symbols：  [2, 7]
DMRS comb：     24
QC delay：      [0,9,18,27,36,45,54,63]
Sidon delay：   [0,1,3,7,12,20,30,65]
接收机：        known-delay matched/oracle 二维 RMMSE
链路：          16QAM、MCS 8、码率 553/1024、8 次 LDPC
```

两个候选使用相同的 TDL realization、payload、LS noise 和 data noise。两个 DMRS symbols 不提前平均，由零速度匹配时间协方差完成联合估计。

先进行小规模 SNR 扫描定位 10% 和 1% BLER 区域，再在交叉区域增加 trial。报告：

- QC 和 Sidon 的 10% BLER 目标 SNR；
- 样本足够时的 1% BLER 目标 SNR；
- Sidon 相对 QC 的 SNR 差和置信区间；
- 两者的 data-RE CE NMSE；
- 与 result-024 平坦模型中 0.33 dB 和 0.85 dB 改善的差别。

判读方式：

- 若 Sidon 仍显著优于 QC，说明 plan-024 的排序至少保持到 5 ns、零速度 TDL-A；
- 若优势缩小、消失或反向，按实际结果记录，说明原平坦模型结论不能直接推广；
- 不把“接近 0.33/0.85 dB”设为代码验收条件；代码正确性由前面的协方差和 RMMSE 验证保证。

unknown-delay 分支也使用相同 realization 运行并单独报告，但不与 known-delay 结果混合拟合。它用于验证接收机不知道 CDD delay 时的性能变化。

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
9. 使用 5 ns、零速度 TDL-A 小规模重跑 plan-024 的 QC/Sidon 对比。
10. 定位 BLER 交叉区域后增加 trial，完成近似平坦 TDL 对比。
11. 根据结果决定下一步是否扩展到 3/60 km/h 和更大 delay spread。

## 五、2026-07-22 执行范围确认

研究者确认本轮正式链路扫描只执行 plan-024 E2 的 QC/Sidon delay 对比：

- 只运行 known-delay matched/oracle 二维 RMMSE；不运行 unknown-delay 正式 BLER；
- 粗扫复用 plan-024 的 `13.0:0.5:18.0 dB`，每点 400 trials；
- 粗扫先确认 5 ns TDL-A 下 10% 和 1% BLER 交叉区域；
- 若交叉仍落在原区间，10% 区间使用 `14.25,14.5,14.75 dB`，每点 3000 trials；
- 若交叉仍落在原区间，1% 区间使用 `16,16.25,16.5,16.75,17,17.25,17.5 dB`，每点 3000 trials；
- 若 TDL-A 使交叉移出上述精化区间，先报告粗扫结果，再由研究者确认新的精化点；
- 1% 目标附近任一候选累计错误块少于 30 时，只报告先导结果，不做确定性结论。
