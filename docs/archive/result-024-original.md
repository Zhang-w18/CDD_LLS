# result-024：分段边界已知接收机诊断，以及 Sidon / QC 的 V-aware 链路终审

> 对应 `research/plan-024.md`。代码：`tools/run_experiment024_segment_sidon_qc.py`。原始数据：`outputs/experiment024_segment_sidon_qc/20260716_main/`。

## 1. TL;DR

**E1：分段对齐估计不是B3/B4的普遍补救。** 在接收机已知候选分段边界、每段独立使用均匀PDP MMSE的有利条件下，48PRB只有6/16个B3获得正改善（最大1.27 dB），24PRB只有7/16个B3获得正改善（最大1.63 dB）；B4在两个带宽各100个候选全部恶化。重新代入plan-023的H1门槛后，48PRB仍0个通过，24PRB仅`B3_cc_nseg8_T6_seq`诊断性通过。由于R4使用了候选边界，这个点属于半透明实现，不恢复原Track B“UE不知道V结构”的结论。

**E2：Sidon优势传递到真实有限码长链路，且强于outage预判。** 在相同V-aware matched LMMSE下，binomial-logit拟合得到Sidon的10% BLER SNR优势 **0.33 dB**（保守95% CI `[0.20,0.45]` dB），1% BLER SNR优势 **0.85 dB**（保守95% CI `[0.63,1.07]` dB）。1%区间每点3000 trials；Sidon在16.5 dB为31/3000错，QC在17.25 dB为30/3000错。两者matched CE NMSE在全部加密点最大差仅0.06 dB，优势来自分集/有限码长结构而非估计条件不公平。

## 2. 实验配置回执

### 2.1 E1：R4分段对齐MMSE

- 系统、候选、DMRS、SNR归一化与plan-023一致。
- 全部B3/B4候选；48PRB与24PRB。
- R4使用候选真实`N_seg`和名义边界，但不知道候选相位值；每个名义段只使用段内comb导频，目标为段内全部非导频RE。
- 假设协方差仍为均匀PDP `[0,300 ns]`，真实MSE使用`R_g=VV^H`闭式计算。
- B3过渡带位于边界两侧，过渡带RE未删除。

### 2.2 E2：V-aware matched estimated-CSI BLER

- 48PRB，`K=576`，`N_t=8`，单Rx，平坦`h~CN(0,I_8)`。
- `QC_arith_s9`：`j=[0,9,18,27,36,45,54,63]`。
- `Sidon`：`j=[0,1,3,7,12,20,30,65]`。
- DMRS comb `S_f=24`、两个DMRS符号平均；两者均用自己的真实`R_g=VV^H`构造全带matched LMMSE。
- 16QAM、MCS 8（码率553/1024）、8次LDPC迭代。
- 两候选逐trial共用`h`、payload、导频噪声、数据噪声；粗扫400 trials/点，1%区间3000 trials/点。
- 10%区间`{14.25,14.5,14.75}`和1%区间`{16,16.25,...,17.5}`均为3000 trials/点。
- 加密使用batch size 20：每批随机生成一个payload，批内20个独立信道/噪声realization复用该码字；不同批重新随机payload。两候选始终成对比较。

## 3. E1结果：分段估计是否降低NMSE

定义

```text
gain_R4_vs_R1 = NMSE_R1(dB) - NMSE_R4(dB),
```

正值表示R4更好。

| 带宽 | 家族 | 候选数 | R4改善数 | 改善中位数 | 最大改善 | 最差变化 |
|---|---:|---:|---:|---:|---:|---:|
| 48PRB | B3 | 16 | 6 | −0.53 dB | +1.27 dB | −3.43 dB |
| 48PRB | B4 | 100 | 0 | −2.46 dB | −2.02 dB | −2.93 dB |
| 24PRB | B3 | 16 | 7 | −0.43 dB | +1.63 dB | −3.24 dB |
| 24PRB | B4 | 100 | 0 | −4.41 dB | −3.06 dB | −5.97 dB |

最明显的正改善点：

| 带宽 | 候选 | R1 NMSE | R4 NMSE | 改善 | 每段导频数 | H1诊断复判 |
|---|---|---:|---:|---:|---:|---|
| 48PRB | `B3_cc_nseg16_T6_bitrev` | −6.05 | −7.31 | +1.27 | 1–2 | 不通过 |
| 48PRB | `B3_cc_nseg16_T6_seq` | −6.59 | −7.85 | +1.26 | 1–2 | 不通过 |
| 24PRB | `B3_cc_nseg8_T6_seq` | −6.50 | −8.13 | +1.63 | 1–2 | **通过** |
| 24PRB | `B3_cc_nseg8_T6_bitrev` | −5.97 | −7.57 | +1.59 | 1–2 | 不通过 |

图：`outputs/experiment024_segment_sidon_qc/20260716_main/figures/segment_nmse_gain.png`。

### 3.1 机理解释

1. R4避免把分段两侧的稳态区域放进同一个4RB窗口，所以对“段较短、过渡带较窄”的部分B3有效。
2. R4仍必须估计边界两侧的过渡带RE；相位连续不等于局部协方差平稳，因此改善上限有限。
3. 分段越短，每段不同频率导频越少。1–2个导频限制了降噪和局部变化辨识能力；避免跨界得到的偏差改善可能被观测数下降带来的方差损失抵消。
4. B4段内虽然斜率固定，但局部信道本来已较平滑；R1的重叠4RB窗可从相邻频率借用更多导频。R4把观测截断到每段后减少导频，且B4边界只有斜率变化、没有相位跳变，因此“禁止跨界”的收益小于导频损失，导致全部恶化。

### 3.2 对H1的含义

- 48PRB：即使授予候选边界信息，仍无候选越过原透明基线门槛。
- 24PRB：一个B3点诊断性通过，说明窗口/边界不对齐确实是该点的重要损失来源。
- 该点需要UE知道`N_seg`和边界；若这些信息未标准化或未信令，就不属于原Track B。若继续研究，应把问题改写为“分段结构与接收机/DMRS联合设计”，并加入同等边界信息下的PRG/CDD基线，而不是继续称为透明R1胜出。

## 4. E2结果：Sidon与QC的V-aware BLER

目标SNR使用各目标区间内所有3000-trial点做二项logit拟合：

```text
logit(BLER) = a + b*SNR_dB.
```

由拟合反解10%/1%目标SNR；单候选区间用Fisher信息矩阵和delta method计算。候选逐trial成对，但原始文件只保留聚合误块数，无法恢复误块联合表，因此候选差值的区间忽略正配对协方差，属于保守近似。

| 目标 | QC SNR (95% CI) | Sidon SNR (95% CI) | Sidon优势 | 优势保守95% CI |
|---|---:|---:|---:|---:|
| 10% BLER | 14.65 `[14.56,14.74]` | 14.32 `[14.23,14.41]` | **0.33 dB** | **`[0.20,0.45]`** |
| 1% BLER | 17.23 `[17.06,17.40]` | 16.38 `[16.24,16.52]` | **0.85 dB** | **`[0.63,1.07]`** |

1%附近的直接计数：

| SNR | QC错误/3000 | QC BLER | Sidon错误/3000 | Sidon BLER |
|---:|---:|---:|---:|---:|
| 16.00 | 119 | 3.967% | 48 | 1.600% |
| 16.25 | 77 | 2.567% | 30 | 1.000% |
| 16.50 | 64 | 2.133% | 31 | 1.033% |
| 16.75 | 58 | 1.933% | 19 | 0.633% |
| 17.00 | 41 | 1.367% | 12 | 0.400% |
| 17.25 | 30 | 1.000% | 13 | 0.433% |
| 17.50 | 19 | 0.633% | 8 | 0.267% |

16.25/16.5 dB处Sidon的轻微非单调（30→31错）在二项统计范围内；拟合使用全部7点而不是以单点累计最小值强制插值。

图：`outputs/experiment024_segment_sidon_qc/20260716_main/figures/sidon_qc_bler_refined.png`。

### 4.1 与H3 outage的关系

result-023的outage预测为：Sidon相对等差CDD在10%处好0.210 dB、1%处好0.368 dB。本轮真实estimated-CSI BLER分别为0.33和0.85 dB，方向一致，但深尾优势更大。

CE公平性检查：加密SNR点上`NMSE_Sidon(dB)-NMSE_QC(dB)`的最大绝对值为0.060 dB，且正负均有，不支持“Sidon因估计更准而赢”的解释。差异来自`V`导致的频域增益联合结构及其与有限码长译码的相互作用。一个合理但尚未单独证实的机理是：等差集合保留大量低阶加性共振，Sidon消除这些四阶共振，使坏信道实现中的信息量分布更均匀；有限码长译码在1%尾部放大了这一区别。

### 4.2 判定与边界

- 10%主判据`≥0.15 dB`通过，且保守区间下界0.20 dB仍为正。
- 1%具有≥30错误事件的直接锚点，拟合优势0.85 dB，保守区间不跨0，结论成立。
- 本轮“真实BLER”指真实16QAM/LDPC/软解调/matched CE链路，但物理分支信道仍采用H3一致的带内平坦模型。尚未验证5 ns TDL展宽、定时误差和V失配。Sidon的最小delay栅格间距小于QC，下一轮最重要的风险检查是5–100 ns物理PDP下的导频条件数、NMSE和BLER。

## 5. Sidon名称与本项目中的严格含义

Sidon集（也称`B_2`集）来自加性组合学，以数学家Simon Sidon命名。整数集合`A`若满足

```text
a+b=c+d
```

时只有`{a,b}={c,d}`这种平凡解，则称为Sidon集。对delay索引，非平凡四元关系

```text
j_a+j_c=j_b+j_d (mod K)
```

会产生四阶Weyl共振；等差集合因`j_{n-1}+j_{n+1}=2j_n`具有大量此类关系。Sidon型集合通过减少非平凡成对和重复，把最低阶残差压低。

本项目使用`[0,1,3,7,12,20,30,65]`。8个元素含自配对共有`8×9/2=36`个无序成对和；检查结果为：整数域36个和全部不同，模`K=576`后仍为36个不同余数。因此它在本实验的全带循环群内也没有非平凡四阶加性关系，可以严格称为Sidon集。模导频周期24后的成对和并不唯一，但8个单独的`j_n mod 24`仍互异，因此不破坏导频域列正交。

## 6. 原始数据与复现

```text
outputs/experiment024_segment_sidon_qc/20260716_main/
├── segment_nmse.csv
├── segment_summary.json
├── figures/segment_nmse_gain.png
├── sidon_qc_bler.csv                         # 400-trial粗扫
├── sidon_qc_summary.json
├── sidon_qc_bler_refined.csv                 # 10%/1%加密点合并
├── final_link_summary.json                    # logit目标SNR、CI、CE公平性
├── figures/sidon_qc_bler_refined.png
├── refine_10pct/                              # 14.25/14.5 dB，各3000
├── refine_10pct_14p75/                        # 14.75 dB，3000
└── refine_1pct/
    ├── sidon_qc_bler.csv                     # 3000-trial深尾加密
    ├── sidon_qc_summary.json
    └── figures/sidon_qc_bler.png
```

复现：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/cdd_lls_matplotlib \
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python \
tools/run_experiment024_segment_sidon_qc.py --stage segment

MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/cdd_lls_matplotlib \
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python \
tools/run_experiment024_segment_sidon_qc.py --stage link \
--snrs 16,16.25,16.5,16.75,17,17.25,17.5 --trials 3000 --batch-size 20 \
--out outputs/experiment024_segment_sidon_qc/20260716_main/refine_1pct

MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/cdd_lls_matplotlib \
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python \
tools/analyze_experiment024.py
```
