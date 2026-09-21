# 仿真平台信道、SNR 与分集统计实现说明

## 1. 检查信息与适用范围

- 最后更新：2026-09-17 11:34:36 +08:00（Asia/Singapore）
- 代码基线：Git `83b25bc` 加当前工作区 CDL 实现；本轮尚未创建 checkpoint
- 本机依赖：Sionna 1.0.2、NumPy 1.26.4
- 检查范围：PDCCH、PDSCH、多个接收天线、Sionna TDL/CDL 信道、预编码功率归一化、AWGN、最大比合并（maximum-ratio combining, MRC）以及用于观察分集增益的接收功率统计。
- 主要实现位置：`cdd_lls/phy/channel_tdl.py`、`cdd_lls/phy/precoding.py`、`cdd_lls/sim/orchestrator.py`、`cdd_lls/sim/pdcch.py`、`cdd_lls/sim/pdcch_cdd.py`、`tools/run_bler_curves.py` 和 `tools/run_plan033_tdl_mobility_mimo.py`。

本文记录检查当日的实现事实。后续修改信道后端、预编码归一化、噪声注入或接收合并方式后，应更新检查日期和代码基线，并重新核对本文结论。

## 2. 结论摘要

1. Sionna TDL 的 PDP 是固定的统计参数：抽头时延和抽头平均功率由所选 profile 决定。每次信道实现中的复抽头系数是随机的，瞬时抽头功率不是固定值。
2. 非视距 TDL 抽头的复包络由有限个正弦波叠加产生。当前 `num_sinusoids=20`，复包络近似复高斯，幅度近似瑞利分布，瞬时功率近似指数分布；每个抽头的平均功率由 PDP 固定。
3. 当前 C300 场景实际配置为 `tdl_profile: C` 和 `delay_spread_ns: 300.0`，使用 24 抽头 TDL-C profile，并将归一化时延乘以 300 ns。它不是 Sionna 名为 `C300` 的固定 12 抽头 profile。
4. $N_t$ Tx、$N_r$ Rx 时，代码只构造一个 TDL 对象并一次输出形状为 $[B,N_r,N_t,S,K]$ 的信道张量，不是循环创建 $N_tN_r$ 个 TDL 对象。未传入空间相关矩阵；所有天线对共享同一个 PDP、路径时延集合和 Doppler sinusoid 频率结构，但每个天线对、每条路径、每个 sinusoid 都有独立随机初相位。该实现不是“所有天线对只有一个不同总相位、衰落幅度完全相同”的退化模型。
5. Sionna 路径先生成连续时延的 CIR 路径系数与时延，再直接按傅里叶和计算各子载波频响。这里不是先把抽头量化到离散采样点后再执行数值 FFT。
6. SNR 的物理定义是：单个 Rx 分支上，单个 RE 的总发射功率与该 Rx 分支单个复噪声样本方差之比。噪声方差不随 $N_r$ 缩放。
7. 平台有两种数值归一化。PDCCH、通用 PDSCH 和 plan-033 使用每 RE 总发射功率 1、噪声方差 $1/\mathrm{SNR}$；固定 8Tx PDSCH BLER 入口保留每 RE 总发射功率 8、噪声方差 $8/\mathrm{SNR}$。两者的 SNR 定义相同。
8. 通用 `run.py` 现可通过 `channel.backend` 选择 `sionna_tdl` 或 `sionna_cdl`。CDL-A～E 使用 Sionna `PanelArray` 的实际阵列几何生成空间相关 MIMO 信道；现阶段只支持单极化阵列。
9. CDL 支持 `IDEAL`、`TF_RMMSE_KNOWN` 和 `TF_RMMSE_UNKNOWN` 接收机。两种 TF-RMMSE 的时间/频率协方差由 CDL cluster PDP 和 Doppler 构造，但没有纳入 CDL 的空间相关，因此是 spatial-unaware 近似，不是完全 matched CDL LMMSE。

### 2.1 当前平台支持范围与使用方法

通用 PDSCH 入口为：

```powershell
python run.py --config <config.yaml>
```

信道选择由 `channel.backend` 决定：

| `channel.backend` | 信道 | profile 字段 | 适用入口 |
|---|---|---|---|
| `legacy_exponential` | 旧的离散指数 PDP | 不使用 3GPP profile | 通用 `run.py` |
| `sionna_tdl` | 3GPP TR 38.901 TDL | `tdl_profile: A`～`E` | 通用 `run.py`；现有固定 TDL 实验入口 |
| `sionna_cdl` | 3GPP TR 38.901 CDL | `cdl_profile: A`～`E` | 通用 `run.py` |

现有 TDL YAML 不需要修改。一个最小 CDL 配置如下，完整可运行示例见 `configs/smoke_sionna_cdl.yaml`：

```yaml
antenna:
  n_tx: 2
  n_rx: 1

channel:
  backend: sionna_cdl
  model: 3gpp_tr38901_cdl
  cdl_profile: A
  delay_spread_ns: 100.0
  carrier_frequency_hz: 3.5e9
  ue_speed_kmh: 3.0
  cdl_direction: downlink
  cdl_tx_array_rows: 1
  cdl_tx_array_cols: 2
  cdl_rx_array_rows: 1
  cdl_rx_array_cols: 1
  cdl_polarization: single
  cdl_polarization_type: V
  cdl_antenna_pattern: omni
  cdl_element_vertical_spacing: 0.5
  cdl_element_horizontal_spacing: 0.5

channel_estimation:
  ce_method: IDEAL
```

阵列行列数必须满足

$$
N_{{\rm row,tx}}N_{{\rm col,tx}}=N_t,
\qquad
N_{{\rm row,rx}}N_{{\rm col,rx}}=N_r.
$$

`cdl_tx_array_cols: 0` 或 `cdl_rx_array_cols: 0` 表示按对应天线数自动建立单行均匀线阵；显式配置非零列数时，行列乘积必须与 `antenna.n_tx` 或 `antenna.n_rx` 一致。元素间距的单位是波长。`cdl_direction` 决定移动端（user terminal, UT）与基站（base station, BS）阵列在 Sionna CDL 中的角色：下行时 Tx 是 BS、Rx 是 UT；上行时相反。

`ue_speed_kmh` 接受任意非负速度，平台换算为 m/s 后同时设置 Sionna 的 `min_speed` 与 `max_speed`，因此每个配置使用固定速度、随机三维运动方向。

当前 CDL 接收机选项为：

- `IDEAL`：数据 RE 直接使用真实等效信道，适合验证信道、预编码和理想 CSI 性能；
- `TF_RMMSE_KNOWN`：接收机知道人工 CDD delay，并使用 CDL PDP/Doppler 的 spatial-unaware 协方差；
- `TF_RMMSE_UNKNOWN`：接收机不使用人工 CDD delay，使用未叠加人工 delay 的 CDL PDP/Doppler spatial-unaware 协方差。

输出的 `resolved_config.yaml` 保存展开配置；`summary.csv`/`summary.json` 中保存 `channel_backend`、`channel_model`、`tdl_profile`、`cdl_profile`、`covariance_type` 等字段。CDL estimated-CSI 行的 `covariance_type` 会包含 `cdl_pdp_doppler_spatial_unaware`，用于防止把该接收机误认为空间统计完全匹配。

本次扩展没有改变冻结的专题实验入口范围。`tools/run_bler_curves.py`、plan-031/032/033/034 与 PDCCH 专用入口仍按各自已有 schema 使用 TDL；需要在这些入口中开展 CDL 正式实验时，应分别扩展其 schema、固定物理范围和测试，不能只修改 YAML 字段。

## 3. TDL 抽头与 PDP

### 3.1 固定量与随机量

对第 $r$ 根 Rx、第 $n$ 根 Tx、第 $\ell$ 条路径，当前 Sionna 非视距 TDL 实现可写为

$$
a_{r,n,\ell}(t)
=
\sqrt{p_\ell}\frac{1}{\sqrt{N_{\rm sos}}}
\sum_{q=1}^{N_{\rm sos}}
\exp\!\left(j\left[\omega_{\ell,q}t+\phi_{r,n,\ell,q}\right]\right),
$$

其中：

- $\tau_\ell$ 和 $p_\ell$ 是 profile 给出的固定路径时延和归一化平均功率；
- $N_{\rm sos}=20$ 是当前配置的正弦波数量；
- $\phi_{r,n,\ell,q}$ 是随机初相位；
- 移动场景中的 $\omega_{\ell,q}$ 由载频、速度及随机到达角产生。

因此，PDP 描述的是集合平均功率

$$
\mathbb E\!\left[|a_{r,n,\ell}(t)|^2\right]=p_\ell,
\qquad
\sum_\ell p_\ell=1,
$$

而不是要求每次实现都满足 $|a_{r,n,\ell}(t)|^2=p_\ell$。准确地说，近似瑞利的是幅度 $|a|$，对应的瞬时功率 $|a|^2$ 近似指数分布。由于当前使用有限的 20 个正弦波，这只是对复高斯/瑞利过程的有限和近似。

Sionna 在载入 profile 时先把 dB 功率转换到线性域，再归一化到总和为 1。平台调用 `GenerateOFDMChannel(..., normalize_channel=False)` 或 `cir_to_ofdm_channel(..., normalize=False)`，不会把每次信道实现再次强制归一化到单位功率。因此单次实现的全带平均功率可以随机波动。YAML 中的 `channel.normalize: true` 当前不控制 Sionna TDL 的这一行为；该字段对 Sionna 路径实际上没有生效，判断归一化时应以 `normalize_channel=False` 和 profile 的集合归一化为准。

旧的 `legacy_exponential` 后端不是 3GPP TDL-A～E。它先在 FFT 采样间隔上构造指数 PDP，然后对每个 Tx--Rx--tap 直接生成独立的复高斯变量。PDCCH 当前不使用该后端；通用 PDSCH 入口仍允许显式选择它。

### 3.2 当前使用的 TDL-C、300 ns 参数

当前 PDCCH C300 配置使用以下参数：

- `tdl_profile: C`
- `delay_spread_ns: 300.0`
- `carrier_frequency_hz: 4.0e9`
- `ue_speed_kmh: 3.0`
- `num_sinusoids: 20`

TDL-C 是非视距 profile。表中“相对功率”是 Sionna profile 文件中的原始 dB 值；“归一化线性功率”是 Sionna 转为线性域并除以所有路径功率之和后的实际 $p_\ell$。行的次序保持 profile 原顺序，时延与功率必须逐行配对，不应分别排序。

| $\ell$ | 归一化时延 $\tau_\ell/\mathrm{DS}$ | 300 ns 下时延（ns） | 相对功率（dB） | 归一化线性功率 $p_\ell$ |
|---:|---:|---:|---:|---:|
| 0 | 0.0000 | 0.00 | -4.4 | 0.061806 |
| 1 | 0.2099 | 62.97 | -1.2 | 0.129130 |
| 2 | 0.2219 | 66.57 | -3.5 | 0.076038 |
| 3 | 0.2329 | 69.87 | -5.2 | 0.051408 |
| 4 | 0.2176 | 65.28 | -2.5 | 0.095726 |
| 5 | 0.6366 | 190.98 | 0.0 | 0.170227 |
| 6 | 0.6448 | 193.44 | -2.2 | 0.102572 |
| 7 | 0.6560 | 196.80 | -3.9 | 0.069347 |
| 8 | 0.6584 | 197.52 | -7.4 | 0.030976 |
| 9 | 0.7935 | 238.05 | -7.1 | 0.033192 |
| 10 | 0.8213 | 246.39 | -10.7 | 0.014489 |
| 11 | 0.9336 | 280.08 | -11.1 | 0.013214 |
| 12 | 1.2285 | 368.55 | -5.1 | 0.052605 |
| 13 | 1.3083 | 392.49 | -6.8 | 0.035565 |
| 14 | 2.1704 | 651.12 | -8.7 | 0.022963 |
| 15 | 2.7105 | 813.15 | -13.2 | 0.008148 |
| 16 | 4.2589 | 1277.67 | -13.9 | 0.006935 |
| 17 | 4.6003 | 1380.09 | -13.9 | 0.006935 |
| 18 | 5.4902 | 1647.06 | -15.8 | 0.004477 |
| 19 | 5.6077 | 1682.31 | -17.1 | 0.003319 |
| 20 | 6.3065 | 1891.95 | -16.0 | 0.004276 |
| 21 | 6.6374 | 1991.22 | -15.7 | 0.004582 |
| 22 | 7.0427 | 2112.81 | -21.6 | 0.001178 |
| 23 | 8.6523 | 2595.69 | -22.8 | 0.000893 |

Sionna 1.0.2 还提供名为 `C300` 的另一份固定参数文件。它有 12 条路径，时延为 $[0,65,70,190,195,200,240,325,520,1045,1510,2595]$ ns，相对功率为 $[-6.9,0,-7.7,-2.5,-2.4,-9.9,-8.0,-6.6,-7.1,-13.0,-14.2,-16.0]$ dB。当前被检查的 PDCCH C300 YAML 均写为 profile `C` 加 300 ns，而不是 profile `C300`；两者不能互换名称或混用参数。

## 4. 多 Tx、多 Rx 信道的生成方式

调用接口向同一个 Sionna `TDL` 对象传入 `num_tx_ant=N_t` 和 `num_rx_ant=N_r`。仓库取回的频域信道形状为

$$
[B,N_r,N_t,S,K],
$$

其中 $B$ 是 batch 大小，$S$ 是 OFDM 符号数，$K$ 是子载波数。实现没有为每个天线对循环实例化一个独立 TDL 对象，也没有传入 `spatial_corr_mat`、`rx_corr_mat` 或 `tx_corr_mat`。

### 4.1 Doppler sinusoid 与初相位的共享关系

Sionna 1.0.2 的相关随机张量维度如下。表中省略了值为 1 的逻辑 `num_rx`、`num_tx` 外层链路维，保留其原始位置以便与源码核对。

| 随机量 | Sionna 中的维度 | 天线对之间的关系 |
|---|---|---|
| 最大 Doppler | $[B,1,1,1,1,1,1,1]$ | 同一 batch realization 内由所有 Tx--Rx 天线对共享 |
| sinusoid 到达角偏移 $\theta$ | $[B,1,1,1,1,L,1,N_{\rm sos}]$ | 每条路径、每个 sinusoid 可不同，但由所有天线对共享 |
| 初相位 $\phi$ | $[B,1,N_r,1,N_t,L,1,N_{\rm sos}]$ | 对每个天线对、每条路径、每个 sinusoid 分别随机生成；仅沿时间维保持不变 |

因此，第 $r$ 根 Rx、第 $n$ 根 Tx、第 $\ell$ 条路径的系数是

$$
a_{r,n,\ell}(s)
=
\sqrt{\frac{p_\ell}{N_{\rm sos}}}
\sum_{q=1}^{N_{\rm sos}}
\exp\!\left(
j\left[\omega_{\ell,q}sT+\phi_{r,n,\ell,q}\right]
\right),
$$

其中当前 $N_{\rm sos}=20$。不同天线对共享 $\omega_{\ell,q}$，但使用不同的整组

$$
\left\{\phi_{r,n,\ell,q}:q=1,\ldots,N_{\rm sos}\right\}.
$$

这里的“各自独立抽取正弦波初相位”不是给每个天线对只抽取一个总相位 $\phi_{r,n}$。当前实现一般不能写成

$$
a_{r,n,\ell}(s)=e^{j\phi_{r,n}}a_{0,0,\ell}(s).
$$

只有当两个天线对的所有 sinusoid 初相位之差都恰好相同，即

$$
\phi_{r,n,\ell,q}-\phi_{r',n',\ell,q}
=
\Delta\phi
\quad\text{对所有 }q\text{ 成立},
$$

才能提出一个公共总相位。当前各 $\phi_{r,n,\ell,q}$ 分别随机生成，连续分布下出现上述条件的概率为零。因此，不同天线对通常同时具有不同的合成幅度和合成相位，而不是满足

$$
|a_{r,n,\ell}(s)|=|a_{r',n',\ell}(s)|.
$$

每个 sinusoid 自身的幅度为 1，但路径系数的幅度由 20 个相量的矢量和产生：

$$
|a_{r,n,\ell}(s)|
=
\sqrt{\frac{p_\ell}{N_{\rm sos}}}
\left|
\sum_{q=1}^{N_{\rm sos}}
e^{j(\omega_{\ell,q}sT+\phi_{r,n,\ell,q})}
\right|.
$$

平台没有再独立抽取一个显式的“瑞利幅度变量”。近似瑞利的随机幅度来自独立初相位相量的叠加，并满足

$$
\mathbb E\!\left[|a_{r,n,\ell}(s)|^2\right]=p_\ell.
$$

速度为 0 时，$\omega_{\ell,q}=0$，所以同一次 realization 中每个天线对的路径幅度在 OFDM 符号间保持不变，但不同天线对的幅度通常仍不相同。速度非零时，共享的 Doppler 频率与各自独立的初相位共同产生不同的随时间幅度轨迹。

对于不同天线对 $(r,n)\ne(r',n')$，独立初相位使路径系数的空间互相关为零：

$$
\mathbb E\!\left[
a_{r,n,\ell}(s)a_{r',n',\ell}^*(u)
\right]=0.
$$

因此，对平台的二阶统计和接收机实现，正确表述是“$N_tN_r$ 个空间相关矩阵为单位阵、具有相同 PDP 的分支”，平台模型记为 `identity_no_spatial_correlation`。但若“完全独立的 MIMO fading”要求每个天线对连 Doppler sinusoid 频率集合也分别独立生成，则当前实现不满足该更强定义：不同天线对共享 Doppler 频率结构，完整时间过程可能共享高阶结构，只是空间二阶互相关为零，而且固定时刻的合成幅度通常不相同。

在多 Rx 接收端，每根 Rx 分支分别进行相同的信道估计，不进行跨 Rx 联合信道估计。数据检测阶段再沿 Rx 维执行相干 MRC。

## 5. 从时延路径到频域信道

Sionna 先生成每条连续时延路径的系数 $a_{r,n,\ell}(s)$ 和时延 $\tau_\ell$，随后直接在每个 OFDM 子载波频率 $f_k$ 上计算

$$
H_{r,n,s,k}
=
\sum_{\ell=0}^{L-1}
a_{r,n,\ell}(s)
\exp\!\left(-j2\pi f_k\tau_\ell\right).
$$

`generate_sionna_tdl_channel()` 通过 Sionna `GenerateOFDMChannel` 完成这一步；`generate_sionna_tdl_channel_active()` 为节省内存，仅在使用中的 active subcarriers 上调用 `cir_to_ofdm_channel()` 计算同一个傅里叶和。PDCCH CDD 入口使用 active-only 路径，通用 PDSCH 和固定静态 PDSCH BLER 入口主要使用 full-FFT 生成后再截取 active subcarriers 的路径。

因此可概括为“先生成时延域的物理路径，再转换为频域响应”，但转换是按连续时延直接求和，不是将路径先放到整数离散 tap 栅格后执行数值 FFT。当前 OFDM 链路直接在频域施加 $y=Hx+n$，不包含时域波形卷积、循环前缀截断、ISI 或 ICI；移动性只使 CIR 在所请求的 OFDM 符号时刻发生变化。

## 6. 预编码、SNR 与多 Rx 合并

### 6.1 信号模型与 SNR 定义

对单层传输，第 $q$ 个 RE 上第 $r$ 根 Rx 的等效信道和接收信号为

$$
g_{t,r,q}
=
\sum_{n=0}^{N_t-1}H_{t,r,n,q}W_{t,q,n},
\qquad
y_{t,r,q}=g_{t,r,q}x_{t,q}+z_{t,r,q},
$$

其中 $\mathbb E[|x|^2]=1$，$z_{t,r,q}\sim\mathcal{CN}(0,\sigma_n^2)$，且各 Rx 分支的噪声独立。定义每 RE 总发射功率

$$
P_{{\rm tx},q}=\sum_{n=0}^{N_t-1}|W_{t,q,n}|^2.
$$

平台的输入 SNR 定义为

$$
\mathrm{SNR}
=
\frac{P_{{\rm tx},q}}{\sigma_n^2}.
$$

它是单个 Rx 分支、单个 RE 上的口径，不把 $N_r$ 个分支的接收功率预先相加，也不把噪声方差除以 $N_r$。

### 6.2 不同入口的数值归一化

| 入口 | $P_{{\rm tx},q}$ | $\sigma_n^2$ | 说明 |
|---|---:|---:|---|
| PDCCH `pdcch-bler-v1` / `pdcch-cdd-bler-v2` | 1 | $10^{-\mathrm{SNR}_{\rm dB}/10}$ | 预编码向量逐 RE 单位范数 |
| 通用 PDSCH `run.py` / `cdd_lls.sim.orchestrator` | 1 | $10^{-\mathrm{SNR}_{\rm dB}/10}$ | CDD 和 PRG 预编码均按总功率 1 构造 |
| plan-033 多 Rx PDSCH | 1 | $10^{-\mathrm{SNR}_{\rm dB}/10}$ | 明确要求 `precoder_normalize: true` |
| 固定 8Tx PDSCH `tools/run_bler_curves.py` | 8 | $8\,10^{-\mathrm{SNR}_{\rm dB}/10}$ | 历史未归一化预编码，$\|W_q\|^2=8$ |

固定 8Tx PDSCH 入口中 DMRS 对两个符号做平均，所以该平均后的 LS 观测噪声方差为 $\sigma_n^2/2$。这不改变原始单个 DMRS RE 的噪声方差定义。

### 6.3 MRC

estimated-CSI MRC 为

$$
\widehat x_{t,q}
=
\frac{\sum_r\widehat g_{t,r,q}^*y_{t,r,q}}
{\sum_r|\widehat g_{t,r,q}|^2},
\qquad
\widehat\sigma_{{\rm eff},t,q}^2
=
\frac{\sigma_n^2}{\sum_r|\widehat g_{t,r,q}|^2}.
$$

ideal CSI 时用真实 $g$ 代替 $\widehat g$。对应的瞬时合并后 SNR 为

$$
\gamma_{t,q}^{\rm MRC}
=
\frac{\sum_r|g_{t,r,q}|^2}{\sigma_n^2}.
$$

增加独立 Rx 分支会提高 MRC 合并功率并提供接收分集，但不会改变配置的单 Rx SNR 或每个 Rx 分支上的噪声方差。

## 7. 用频域平均接收功率 CDF 观察分集

### 7.1 建议统计量

该思路是合理的，但建议把指标称为“归一化频域平均有效信道功率”，而不是直接称为 RSRP。只有当集合 $\mathcal Q$ 明确定义为参考信号 RE、并且接收合并口径也固定时，才适合使用 RSRP 名称。

令 $\mathcal Q$ 为同一传输块内用于统计的 RE 集合，可以只取频率维，也可以包含多个 OFDM 符号。对第 $t$ 次信道实现，定义

$$
\overline P_t(W)
=
\frac{1}{N_r|\mathcal Q|P_{\rm tx}}
\sum_{q\in\mathcal Q}\sum_{r=0}^{N_r-1}
\left|
\sum_{n=0}^{N_t-1}H_{t,r,n,q}W_{q,n}
\right|^2,
$$

其中候选的 $P_{\rm tx}=\sum_n|W_{q,n}|^2$ 应在所有 $q$ 上相同。PDCCH、通用 PDSCH 和 plan-033 取 $P_{\rm tx}=1$；固定 8Tx PDSCH BLER 入口取 $P_{\rm tx}=8$。再转为 dB：

$$
X_t(W)=10\log_{10}\overline P_t(W).
$$

使用 $T$ 次独立或按共同随机数配对的信道实现，经验 CDF 为

$$
\widehat F_W(x)
=
\frac{1}{T}\sum_{t=1}^{T}
\mathbf 1\!\left\{X_t(W)\le x\right\}.
$$

若只研究单 Rx 的发射分集，应固定 $N_r=1$。若研究 Tx 预编码与 Rx 分集的合成效果，可使用上式，但比较不同 $N_r$ 时必须保留 $1/N_r$ 归一化，否则平均值会因接收天线数线性增长。若目标就是观察合并后的绝对接收功率，则可去掉 $1/N_r$，但不能再把均值移动解释为纯粹的“稳定性改善”。

### 7.2 如何判定“更集中”和“低尾增益”

推荐同时报告以下量：

$$
Q_p(W)=\inf\{x:\widehat F_W(x)\ge p\},
$$

$$
\Delta_p(W;W_0)=Q_p(W)-Q_p(W_0),
\qquad p\in\{0.01,0.10\},
$$

以及

$$
D_{90-10}(W)=Q_{0.90}(W)-Q_{0.10}(W).
$$

在总发射功率、$N_t$、$N_r$、统计 RE 和物理信道样本完全相同的前提下：

- $\Delta_{0.01}>0$ 或 $\Delta_{0.10}>0$ 表示低尾接收功率提高；
- 更小的 $D_{90-10}$ 或更小的 $\operatorname{std}[X_t]$ 表示分布更集中；
- 应同时核对 $\mathbb E[\overline P_t]$，避免把功率归一化错误造成的均值移动误判为分集增益；
- 最好对所有候选使用相同的底层 $H_t$，即共同随机数配对，以便直接统计每次实现上的差值。

在独立同分布 Tx 分支、单位范数确定性预编码下，单个 RE 的 $g_{r,q}$ 边缘分布通常不因 $W_q$ 改变。频率变化的预编码主要改变不同 RE 之间的相关性，因此可能使块内平均功率 $\overline P_t$ 更集中，并抬高低尾分位数。这个结论不是对所有预编码自动成立，必须由 CDF 或分位数验证。

统计时应使用无噪声的真实 $H$ 和 $g=HW$。若直接从 $|y|^2$ 统计，会把 AWGN、调制符号幅度和信道估计误差混入指标。若确实要模拟测量型 RSRP，应固定参考信号幅度，并单独说明是否去噪、是否按 Rx 分支平均以及是否在 MRC 后统计。

### 7.3 与链路性能的关系

平均接收功率 CDF 是直观的机理诊断，但不是频率选择性编码链路的最终分集判据。相同的平均功率仍可能对应不同的深衰落位置和不同 BLER。建议同时统计每次实现的理想 CSI 块平均调制互信息

$$
I_t(W)
=
\frac{1}{|\mathcal Q|}
\sum_{q\in\mathcal Q}
I_{\rm QAM}\!\left(
\frac{\sum_r|g_{t,r,q}|^2}{\sigma_n^2}
\right),
$$

并比较 $I_t$ 的 1% 和 10% outage 分位数。接收功率 CDF 用于说明“总有效功率是否稳定”，互信息 outage 用于说明“这种稳定是否转化为编码链路可利用的频率分集”，最终仍以相同接收机条件下的 BLER--SNR 曲线为准。

## 8. 已识别的实现注意点

1. `channel.normalize` 在 Sionna TDL 路径中目前是无效配置字段，容易让使用者误以为每次 realization 被归一化。当前真实行为是 PDP 集合功率归一化、单次 realization 不归一化。
2. `tdl_profile: C` 加 `delay_spread_ns: 300` 与 `tdl_profile: C300` 是不同模型；当前 C300 场景采用前者。
3. PDSCH 不同入口的预编码数值尺度不同。比较噪声方差、信道估计 NMSE 或接收功率时，必须同时读取 $P_{\rm tx}$，不能只比较 `noise_variance`。
4. 当前多 Rx 模型没有空间相关，也没有跨 Rx 联合信道估计。若未来加入相关矩阵，本文关于独立分支、CDF 集中程度和 MRC 增益的解释需要重检。
5. 上一条只适用于 TDL。CDL 已通过 Tx/Rx `PanelArray` 几何引入空间相关，但接收端仍逐 Rx 估计后执行 MRC，且当前 TF-RMMSE 不使用空间协方差。
6. CDL 当前仅支持 `single` polarization，类型为 `V` 或 `H`；尚未开放双极化端口映射、阵列朝向配置、路径损耗、阴影衰落和多链路拓扑。

## 9. 复查记录

| 日期 | Git 基线 | 结论 |
|---|---|---|
| 2026-09-16 | `83b25bc` | 首次检查并记录当前 PDCCH/PDSCH、多 Rx、TDL、SNR 和分集统计实现；补充 Doppler sinusoid 与天线对初相位的精确共享关系。 |
| 2026-09-17 11:34:36 +08:00 | `83b25bc` + 工作区修改 | 通用 `run.py` 增加 Sionna CDL-A～E、单极化 ULA/URA、上下行方向、IDEAL 与 spatial-unaware TF-RMMSE；增加 CDL smoke 配置和回归测试。 |
