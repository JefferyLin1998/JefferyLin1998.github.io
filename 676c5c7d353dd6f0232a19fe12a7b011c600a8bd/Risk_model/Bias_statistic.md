## **协方差矩阵偏差统计**
这里主要参考[《Eigen-Adjusted Covariance Matrices》](https://jefferylin1998.github.io/676c5c7d353dd6f0232a19fe12a7b011c600a8bd/Risk_model/Eigenfactor_Adjusted_Covariance_Matrices_May2011.pdf)以及[《协方差矩阵估计：特征调整》](https://zhuanlan.zhihu.com/p/82047494)中的解读。

### **定义**
先定义
$$
f_{m}=r_{m}-R_{t}^{M}
$$
表示第m个股票的主动收益率，$R_{t}^{M}$ 是t天的市场收益率。一般会定义为股票的市值加权收益率。根据性质， $f_{m}$的市值加权求和的结果等于0。

根据股票协方差矩阵的定义:
$$
\mathbf{V}_{0}(m n)=\frac{1}{T-1} \sum_{t=1}^{T}\left(f_{m t}-\bar{f}_{m}\right)\left(f_{n t}-\bar{f}_{n}\right)
$$
其中样本的时间窗口为T =200。

最后定义偏差统计量，这是一个用于投资组合的指标，用于判断协方差矩阵预测投资组合的风险能力强弱。定义 $R_{t}$ 为投资组合第 $\mathrm{t}$ 的收益率， $\sigma_{t-1}$ 为t-1天收盘时给出的组合标准差预测。利用 $\sigma_{t-1}$ 对 $R_{t}$ 做标准化，得到 $\mathrm{t}$ 天的z值。
$$
b_{t}=\frac{R_{t}}{\sigma_{t}}
$$
z值的标准差就是偏差统计量
$$
B=\sqrt{\frac{1}{\tau-1} \sum_{t=1}^{\tau}\left(b_{t}-\bar{b}\right)^{2}}
$$
其中 $\tau$ 是测试窗口的大小，可以设置成滚动窗口。
偏差统计量可以看成是实际的标准差和预测的标准差之比，我们希望 $B \approx 1$ 。但是这往往不现 实，因此我们会设置一个置信区间，希望它落在这个置信区间内，95\%置信区间对应的B值范围是 $1 \pm \sqrt{2 / \tau}$ 。当B值大于1时，表明预测的标准差低于实际的标准差，即低估了。反之则高估 了。

