# CDD等效PDP构造方法总结

## 1. 统一模型

考虑 \(K\) 个 CDD 分支。第 \(m\) 个分支对应一个空间波束 \(\mathbf w_m\)，并施加人工时延 \(d_m\)。

设物理信道包含 \(L\) 条时延路径，其公共物理时延为

$$
\tau_1,\tau_2,\ldots,\tau_L.
$$

对于第 \(m\) 个波束，第 \(\ell\) 条路径经过空间投影后的复信道系数记为

$$
g_{m,\ell}.
$$

相应的波束特定路径功率为

$$
p_{m,\ell}
=
\mathbb E\left[|g_{m,\ell}|^2\right].
$$

第 \(m\) 个波束对应的 PDP 为

$$
P_m(\tau)
=
\sum_{\ell=1}^{L}
p_{m,\ell}\delta(\tau-\tau_\ell).
$$

CDD 人工时延会将第 \(m\) 个分支的全部物理路径整体平移 \(d_m\)，因此其等效路径时延为

$$
\tau_{\ell,m}^{\rm eff}
=
\tau_\ell+d_m.
$$

以下三种方法的区别主要在于：如何获得各个 CDD 分支对应的路径功率，以及是否考虑不同波束分支之间的相关性。

对于方法 2 和方法 3，默认各个波束测得的 PDP 使用统一的下行时延参考，因此不同波束之间的路径时延轴已经对齐。

---

## 2. 方法一：基于公共 PDP

### 2.1 基本思想

所有 CDD 分支均采用同一个参考 PDP，例如所关联 SSB 波束对应的 PDP：

$$
P_{\rm ref}(\tau)
=
\sum_{\ell=1}^{L}
p_{{\rm ref},\ell}
\delta(\tau-\tau_\ell).
$$

因此假设每一个 CDD 分支具有相同的路径功率分布：

$$
p_{m,\ell}
=
p_{{\rm ref},\ell},
\qquad
m=1,\ldots,K.
$$

不同分支之间仅通过人工 CDD 时延 \(d_m\) 区分。

### 2.2 等效 PDP

若各个 CDD 分支采用等功率分配，则等效 PDP 为

$$
\boxed{
P_{\rm eq}^{\rm common}(\tau)
=
\frac{1}{K}
\sum_{m=1}^{K}
P_{\rm ref}(\tau-d_m)
}
$$

展开后为

$$
\boxed{
P_{\rm eq}^{\rm common}(\tau)
=
\frac{1}{K}
\sum_{m=1}^{K}
\sum_{\ell=1}^{L}
p_{{\rm ref},\ell}
\delta\left(
\tau-\tau_\ell-d_m
\right)
}
$$

对于 \(K=2\)，且

$$
d_1=0,\qquad d_2=d,
$$

有

$$
\boxed{
P_{\rm eq}^{\rm common}(\tau)
=
\frac12 P_{\rm ref}(\tau)
+
\frac12 P_{\rm ref}(\tau-d)
}
$$

即复制同一个参考 PDP，并按照不同的 CDD 人工时延进行平移。

### 2.3 对应频域协方差

对应的频域协方差为

$$
R_f[k,k']
=
\frac{1}{K}
\sum_{m=1}^{K}
\sum_{\ell=1}^{L}
p_{{\rm ref},\ell}
e^{-j2\pi(f_k-f_{k'})(\tau_\ell+d_m)}.
$$

该方法不需要获得各个窄波束独立的 PDP，因此不存在不同波束 PDP 之间的参考时延对齐问题。

---

## 3. 方法二：基于波束特定 PDP，不考虑波束间相关性

### 3.1 基本思想

每一个 CDD 分支采用其自身空间波束对应的 PDP：

$$
P_m(\tau)
=
\sum_{\ell=1}^{L}
p_{m,\ell}
\delta(\tau-\tau_\ell).
$$

不同窄波束可能对不同物理路径产生不同的空间增益，因此一般有

$$
p_{1,\ell}
\neq
p_{2,\ell}.
$$

该方法保留各个波束自身的路径功率分布，但假设不同波束分支之间互不相关：

$$
\mathbb E
\left[
g_{m,\ell}g_{n,\ell}^{*}
\right]
=
0,
\qquad m\neq n.
$$

默认各个波束 PDP 使用统一的下行时延参考，因此相同的 \(\tau_\ell\) 表示相同的物理传播时延位置。

### 3.2 等效 PDP

加入 CDD 人工时延后，第 \(m\) 个波束的 PDP 整体平移 \(d_m\)：

$$
P_m^{\rm CDD}(\tau)
=
P_m(\tau-d_m).
$$

因此等效 PDP 为

$$
\boxed{
P_{\rm eq}^{\rm beam}(\tau)
=
\frac{1}{K}
\sum_{m=1}^{K}
P_m(\tau-d_m)
}
$$

展开为

$$
\boxed{
P_{\rm eq}^{\rm beam}(\tau)
=
\frac{1}{K}
\sum_{m=1}^{K}
\sum_{\ell=1}^{L}
p_{m,\ell}
\delta\left(
\tau-\tau_\ell-d_m
\right)
}
$$

对于 \(K=2\)：

$$
\boxed{
P_{\rm eq}^{\rm beam}(\tau)
=
\frac12P_1(\tau-d_1)
+
\frac12P_2(\tau-d_2)
}
$$

若

$$
d_1=0,
$$

则

$$
\boxed{
P_{\rm eq}^{\rm beam}(\tau)
=
\frac12P_1(\tau)
+
\frac12P_2(\tau-d_2)
}
$$

### 3.3 对应频域协方差

其频域协方差为

$$
\boxed{
R_f[k,k']
=
\frac{1}{K}
\sum_{m=1}^{K}
\sum_{\ell=1}^{L}
p_{m,\ell}
e^{-j2\pi(f_k-f_{k'})(\tau_\ell+d_m)}
}
$$

对于 \(K=2\)：

$$
R_f[k,k']
=
\frac12
\sum_\ell
p_{1,\ell}
e^{-j2\pi\Delta f_{kk'}(\tau_\ell+d_1)}
+
\frac12
\sum_\ell
p_{2,\ell}
e^{-j2\pi\Delta f_{kk'}(\tau_\ell+d_2)}.
$$

其中

$$
\Delta f_{kk'}=f_k-f_{k'}.
$$

该方法相比公共 PDP 方法能够描述不同窄波束对多径功率分布的不同加权，但忽略了两个波束来源于同一个空间信道而产生的交叉相关项。

---

## 4. 方法三：理想方法，考虑波束间相关性

### 4.1 基本思想

进一步保留不同 CDD 波束之间的联合统计特性。

对于第 \(\ell\) 条物理路径，定义波束域信道向量

$$
\mathbf g_\ell
=
\begin{bmatrix}
g_{1,\ell}\\
g_{2,\ell}\\
\vdots\\
g_{K,\ell}
\end{bmatrix}.
$$

其波束域路径协方差矩阵为

$$
\boxed{
\mathbf C_\ell
=
\mathbb E
\left[
\mathbf g_\ell
\mathbf g_\ell^{H}
\right]
}
$$

即

$$
\mathbf C_\ell
=
\begin{bmatrix}
p_{1,\ell} & c_{12,\ell} & \cdots\\
c_{21,\ell} & p_{2,\ell} & \cdots\\
\vdots & \vdots & \ddots
\end{bmatrix},
$$

其中

$$
c_{mn,\ell}
=
\mathbb E
\left[
g_{m,\ell}g_{n,\ell}^{*}
\right].
$$

方法二仅保留 \(\mathbf C_\ell\) 的对角元素。方法三同时保留非对角元素。

### 4.2 CDD 分支系数

第 \(m\) 个分支在频率 \(f_k\) 上的 CDD 系数为

$$
c_m[k]
=
\frac{1}{\sqrt K}
e^{-j2\pi f_k d_m}.
$$

定义

$$
\mathbf c[k]
=
\frac{1}{\sqrt K}
\begin{bmatrix}
e^{-j2\pi f_kd_1}\\
e^{-j2\pi f_kd_2}\\
\vdots\\
e^{-j2\pi f_kd_K}
\end{bmatrix}.
$$

若发射端还进行了已知的逐子载波功率归一化，可将对应归一化系数直接包含在 \(\mathbf c[k]\) 中。

### 4.3 等效频域协方差

等效信道为

$$
H_{\rm eq}(f_k)
=
\sum_{\ell=1}^{L}
e^{-j2\pi f_k\tau_\ell}
\mathbf c^{T}[k]
\mathbf g_\ell.
$$

因此其频域协方差为

$$
\boxed{
R_f[k,k']
=
\sum_{\ell=1}^{L}
e^{-j2\pi(f_k-f_{k'})\tau_\ell}
\mathbf c^{T}[k]
\mathbf C_\ell
\mathbf c^{*}[k']
}
$$

这是考虑所有波束间相关项后的完整表达式。

对于 \(K=2\)：

$$
\mathbf C_\ell
=
\begin{bmatrix}
p_{1,\ell} & c_{12,\ell}\\
c_{12,\ell}^{*} & p_{2,\ell}
\end{bmatrix}.
$$

展开后：

$$
\begin{aligned}
R_f[k,k']
=
\frac12\sum_\ell
e^{-j2\pi\Delta f_{kk'}\tau_\ell}
\Big[
& p_{1,\ell}
e^{-j2\pi(f_kd_1-f_{k'}d_1)}
\\
+&p_{2,\ell}
e^{-j2\pi(f_kd_2-f_{k'}d_2)}
\\
+&c_{12,\ell}
e^{-j2\pi(f_kd_1-f_{k'}d_2)}
\\
+&c_{12,\ell}^{*}
e^{-j2\pi(f_kd_2-f_{k'}d_1)}
\Big].
\end{aligned}
$$

前两项正是方法二保留的两个独立波束 PDP 项。

后两项为波束间交叉相关项。

### 4.4 关于理想 PDP 的表述

严格来说，当

$$
c_{mn,\ell}\neq0
$$

时，已经无法用一个普通的非负标量 PDP

$$
P_{\rm eq}(\tau)
$$

完整描述等效信道的二阶统计特性。

原因在于交叉项同时依赖 \(f_k\) 和 \(f_{k'}\)，不能简单表示成

$$
\sum_i p_i
e^{-j2\pi(f_k-f_{k'})\tau_i}.
$$

因此方法三更准确的名称应为：

$$
\boxed{\text{基于波束联合统计的理想等效频域协方差}}
$$

若继续使用理想 PDP 这一名称，需要明确其中实际保留了超出普通 PDP 的波束间联合协方差信息。

---

## 5. 三种方法的关系

三种方法可以统一写成以下层级关系。

### 方法一：公共 PDP

假设

$$
p_{m,\ell}
=
p_{{\rm ref},\ell}
$$

且忽略波束间相关性：

$$
c_{mn,\ell}=0.
$$

等效 PDP：

$$
\boxed{
P_{\rm eq}^{\rm common}(\tau)
=
\frac1K
\sum_m
P_{\rm ref}(\tau-d_m)
}
$$

### 方法二：波束特定 PDP

允许

$$
p_{m,\ell}
$$

随波束变化，但仍假设

$$
c_{mn,\ell}=0,\qquad m\neq n.
$$

等效 PDP：

$$
\boxed{
P_{\rm eq}^{\rm beam}(\tau)
=
\frac1K
\sum_m
P_m(\tau-d_m)
}
$$

### 方法三：理想联合统计

同时保留

$$
p_{m,\ell}
$$

和

$$
c_{mn,\ell}.
$$

完整统计量为

$$
\boxed{
R_f[k,k']
=
\sum_\ell
e^{-j2\pi(f_k-f_{k'})\tau_\ell}
\mathbf c^T[k]
\mathbf C_\ell
\mathbf c^*[k']
}
$$

其中方法二相当于令

$$
\mathbf C_\ell
\rightarrow
\operatorname{diag}
\left(
p_{1,\ell},
\ldots,
p_{K,\ell}
\right).
$$

方法一进一步令

$$
p_{1,\ell}
=
p_{2,\ell}
=
\cdots
=
p_{K,\ell}
=
p_{{\rm ref},\ell}.
$$

---

## 6. 时延参考假设

对于方法二和方法三，默认不同波束对应的 PDP 或路径统计量具有统一的时延参考。

即对于不同波束 \(m\)，同一个 \(\tau_\ell\) 均表示相对于同一下行 OFDM timing reference 的物理传播时延。

因此无需额外引入 beam-specific reference offset：

$$
\Delta\tau_m=0.
$$

在该假设下，CDD 人工时延可以直接叠加到各波束的物理路径时延上：

$$
\tau_\ell
\rightarrow
\tau_\ell+d_m.
$$

若实际测量过程中对每个波束的 PDP 独立重新归零，则必须额外估计不同波束之间的时延参考偏差，此时不能直接使用上述方法二和方法三公式。
