## 风险模型

#### 风险
对于风险的定义，我们有以下几点：
- 风险是收益率的标准差
- 风险不具有可加性
- 许多机构投资者关注主动风险和残差风险甚于总风险
- 主动风险主要依赖于主动头寸的规模，而不是基准头寸的规模
- 风险造成的成本与方差成正比
- 风险模型识别重要的风险来源，并把风险分解为多个组成部分

#### 风险模型
风险模型主要常见的有两种，一种为结构化风险模型以Barra模型常见，另一种为统计风险模型。

Barra模型的表示如下：

$V_{n, m}=\sum_{k_{1}, i_{1}=1}^{K} X_{n_{, k_{1}}} \cdot F_{k_{1}, k_{1}} \cdot X_{m, k_{1}}+\Delta_{n, m}$

其中$X_{n_{, k_{1}}}$为因子载荷,$F_{k_{1}, k_{1}}$为因子协方差矩阵，$\Delta_{n, m}$为股票特质性方差。

统计风险模型依赖于历史数据的样本方差和样本协方差，在该模型下我们使用T个时期的样本来估计一个N*N的协方差矩阵。