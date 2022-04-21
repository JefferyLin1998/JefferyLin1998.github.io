### **多种协方差计算、检验及调整**

此节我们计算多种协方差矩阵，同时对这些协方差矩阵的预测能力进行偏差检验，并且使用偏差调整的方法对协方差矩阵进行调整.

协方差类型：

- 样本协方差矩阵(sample covariance matrix)

- 协方差压缩矩阵(shrink covariance matrix)：利用样本协方差矩阵和一个目标正定矩阵做加权

- 单因子模型协方差矩阵(single covariance matrix)：

  选一个股指，将 index return 作为自变量去对个股return进行回归
  $$
  r_{i}(t)=\alpha_{i}(t)+\beta_{i}(t) x(t)+\epsilon_{i}(t),
  $$
  然后以 $\Sigma=\sigma_{0} \boldsymbol{\beta} \boldsymbol{\beta}^{\boldsymbol{T}}+\boldsymbol{D}$ 作为协方差矩阵的估计量，其中 $\sigma_{0}$ 是index return的方差， $\boldsymbol{\beta}$ 是 $\beta_{i}$ 组成的向量， $\boldsymbol{D}$ 是 $\epsilon_{i}$ 的方差组成的对角矩阵.

- 多因子模型-Barra协方差矩阵(Barra covariance matrix)

- 多因子模型-统计风险主成分模型(Stats covariance matrix)

偏差检验组合：

- 组合一：正态分布随机仓位，和为0
- 组合二：alpha信号多头

调整方法：

- 方法一：[特征值调整方法](https://jefferylin1998.github.io/676c5c7d353dd6f0232a19fe12a7b011c600a8bd/Risk_model/Eigenfactor_Adjusted_Covariance_Matrices_May2011.pdf)
- 方法二：滚动偏差检验调整

最终优化尝试：

- 使用统计风险主成分模型+滚动偏差检验调整进行组合优化

样本协方差矩阵与统计风险主成分模型的介绍见[Stats_model](https://jefferylin1998.github.io/676c5c7d353dd6f0232a19fe12a7b011c600a8bd/index.html#/Risk_model/Stats_model). 

#### 一、偏差检验

偏差检验的介绍见[Bias_statistic](https://jefferylin1998.github.io/676c5c7d353dd6f0232a19fe12a7b011c600a8bd/index.html#/Risk_model/Bias_statistic). 下图中，偏差值可以理解为组合实际标准差比上协方差矩阵预测标准差，我们希望为1，两条虚线代表95%置信区间。

组合一偏差检验

<img src="..\CovBiasTest\scm-sto.png" style="zoom: 50%;" /><img src="..\CovBiasTest\shrink-sto.png" style="zoom: 50%;" /><img src="..\CovBiasTest\single-sto.png" style="zoom: 50%;" /><img src="..\CovBiasTest\barra-sto.png" style="zoom: 50%;" /><img src="..\CovBiasTest\stats-sto.png" style="zoom: 50%;" />

组合二偏差检验

<img src="..\CovBiasTest\scm-alpha.png" style="zoom: 50%;" /><img src="..\CovBiasTest\shrink-alpha.png" style="zoom: 50%;" /><img src="..\CovBiasTest\single-alpha.png" style="zoom: 50%;" /><img src="..\CovBiasTest\barra-alpha.png" style="zoom: 50%;" /><img src="..\CovBiasTest\stats-alpha.png" style="zoom: 50%;" />

上述五个组合在alpha多头上的风险预测都明显高估或低估，在随机组合上只有Barra模型表现不佳。因而对协方差矩阵进行偏差调整很有必要。

#### 二、偏差调整(以Stats covariance matrix为例)

##### 1. 特征值偏差调整

首先, 我们对Stats covariance matrix(后续记为$V_0$)做对角化$$\mathbf{D}_{0}=\mathbf{U}_{0}^{\prime} \mathbf{V}_{0} \mathbf{U}_{0}$$,得到对角矩阵 $\mathbf{D}_{0}$ , 特征矩阵为$U_0$.然后，采样收益率序列，满足均值为 0 ，方差为 $D_{0}(k)$ ，得到 $\mathbf{b}_{m}$ ，那么可以得到仿真的个股收益率序列。
$$
\mathbf{f}_{m}=\mathbf{U}_{0} \mathbf{b}_{m}
$$
接下来我们就可以计算仿真的协方差矩阵 $\mathbf{V}_{m}$
$$
\mathbf{V}_{m}=\frac{\mathbf{f}_{m} \cdot \mathbf{f}_{m}^{\prime}}{T-1}
$$
接着我们求仿真的特征组合，对 $\mathbf{V}_{m}$ 进行对角化
$$
\mathbf{D}_{m}=\mathbf{U}_{m}^{\prime} \mathbf{V}_{m} \mathbf{U}_{m}
$$
其中 $\mathbf{U}_{m}$ 的每一列对应一个仿真特征组合。我们利用 $\mathbf{V}_{0}$ 求解仿真特征组合的真实方差
$$
\tilde{\mathbf{D}}_{m}=\mathbf{U}_{m}^{\prime} \mathbf{V}_{0} \mathbf{U}_{m}
$$
因为 $\mathbf{U}_{m}$ 和 $\mathbf{V}_{0}$ 是不匹配的，所以 $\tilde{\mathbf{D}}_{m}$ 不是对角矩阵。图(2)说明了仿真特征矩阵的预测方差 $\mathbf{D}_{m}$ 是有偏的，即 $E\left[D_{m}(k)\right] \neq \tilde{D}_{m}(k)$ 。我们可以计算仿真波动率偏差
$$
\lambda(k)=\frac{1}{M} \sum_{m} \sqrt{\frac{\tilde{D}_{m}(k)}{D_{m}(k)}}
$$
$\mathrm{M}$ 表示仿真的次数。接下来我们进行尺度变换，来解决正态性和平稳性的问题。
$$
\gamma(k)=a[\lambda(k)-1]+1
$$
我们将仿真的波动率偏差相对于1做了一个偏离， $a$ 是一个常数。再利用
$$
\tilde{\mathbf{V}}_{0}=\mathbf{U}_{0} (\gamma^2{\mathbf{D}}_{0}) \mathbf{U}_{0}^{\prime}
$$
调整特征值，得到新的协方差矩阵。

我们进行实验得到使用调整后的协方差矩阵的表现如下

<img src="..\CovBiasTest\xie.png" style="zoom:67%;" /><img src="..\CovBiasTest\ping.png" style="zoom: 67%;" />

通过特征调整仅将不同持仓数目的alpha多头展平，对于整体偏差的下降作用不大。

##### 2. 滚动偏差检验

由于上述图中的点对应的偏差是整个时间段的偏差，我们可以利用前63/126个工作日的滚动偏差调整当日的风险预测

<img src="..\CovBiasTest\ROLL_BIAS.png" style="zoom:67%;" />

可以发现滚动偏差在一段时间内比较稳定。于是我们使用滚动偏差调整预测风险
$$
b_{t}=\frac{R_{t}}{\sigma_{t}*rollbias_t}
$$
alpha多头偏差检验如下：

<img src="..\CovBiasTest\adj.png" style="zoom:67%;" />

可以看到，预测偏差落入置信区间内。同时我们观察alpha多头未调整预测风险与已调整预测风险曲线

<img src="..\CovBiasTest\adj_risk.png" style="zoom:67%;" />

可以发现2021年的风险已显著拔高.

#### 三、信号优化

这里我们比较三个组合：

- bench: 取信号值前5%
- opt: 按照下述优化

$$
\begin{gathered}
\max : \operatorname{sig}^{\prime} \mathrm{w} \\
\mathrm{st}: \operatorname{SUM}(\mathrm{w}) \leq 1 \\
\mathrm{w} \leq 0.015,(\text { 控制最大持仓 } 1.5 \%) \\
\mathrm{w}[: 200] \geq 0.001,(\text { 最少持仓 } 200 \text { 只 }) \\
\frac{(w)^{\top} \cdot V_{\mathrm{i}} \cdot V_{\mathrm{i}}^{\top}(\mathrm{w})}{\Psi_{\mathrm{p}}} \leq b_{i}(\text { 控制主成分风险) } \\
\frac{(\mathrm{w})^{\top} \cdot V \cdot(\mathrm{w})}{\Psi_{\mathrm{p}}} \leq B \text { (控制总风险) }
\end{gathered}
$$

- adj-opt: 对风险进行调整再进行优化
  $$
  \begin{gathered}\max : \operatorname{sig}^{\prime} \mathrm{w} \\\mathrm{st}: \operatorname{SUM}(\mathrm{w}) \leq 1 \\\mathrm{w} \leq 0.015,(\text { 控制最大持仓 } 1.5 \%) \\\mathrm{w}[: 200] \geq 0.001,(\text { 最少持仓 } 200 \text { 只 }) \\\frac{(w)^{\top} \cdot V_{\mathrm{i}} \cdot V_{\mathrm{i}}^{\top}(\mathrm{w})}{\Psi_{\mathrm{p}}} \leq b_{i}*rollbias_t(\text { 控制调整后主成分风险) } \\\frac{(\mathrm{w})^{\top} \cdot V \cdot(\mathrm{w})}{\Psi_{\mathrm{p}}} \leq B*rollbias_t \text { (控制调整后总风险) }\end{gathered}
  $$
  此处的$rollbias_t$使用预先bench组合的滚动偏差





结果：

| port    | 2020-RET% | 2020-MDD% | 2021-SHARPE | 2021-RET% | 2021-MDD% | 2021-SHARPE | 2020-STD% | 2021-STD% |
| ------- | --------- | :-------- | ----------- | --------- | --------- | ----------- | --------- | --------- |
| bench   | 93.64     | 2.68      | 8.33        | 74.48     | 8.67      | 5.26        | 12.45     | 14.62     |
| opt     | 100.64    | 2.55      | 8.86        | 82.15     | 7.04      | 6.10        | 12.58     | 13.93     |
| adj-opt | 99.32     | 2.31      | 8.82        | 76.96     | 5.76      | 6.44        | 12.46     | 12.36     |

累计收益图：

<img src="..\CovBiasTest\cumrtn.png" style="zoom:67%;" />
