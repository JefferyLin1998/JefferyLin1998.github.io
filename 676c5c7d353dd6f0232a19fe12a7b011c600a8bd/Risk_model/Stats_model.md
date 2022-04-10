## **Statistical Risk Models**

这里我们主要参考Z.Kakushadze和W.Yu的论文[Statistical Risk Model](https://jefferylin1998.github.io/676c5c7d353dd6f0232a19fe12a7b011c600a8bd/Risk_model/Kakushadze.pdf)中对统计风险的描述。

### **样本协方差矩阵**

记资产i在时间s的收益率为$R_{is}$,那么我们的样本协方差矩阵则为
$$
C_{i j}=\frac{1}{M} \sum_{s=1}^{M+1} X_{i s} X_{j s}
$$
其中$X_{i s}=R_{i s}-\bar{R}_{i}$ 同时$\bar{R}_{i}=\frac{1}{M+1} \sum_{s=1}^{M+1} R_{i s}$。

由于A股现有4000+支股票，因此我们对M<N时的情况感兴趣，而此时我们会发现$X_{is}$是奇异的。同时，我们观察到X的最后一列$X_{i,M+1} = −\sum_{s=1}^MX_{is}。$
因此我们仅取X的前M列对$C_{ij}$进行分解如下：
$$
C_{i j}=\sum_{s, s^{\prime}=1}^{M} X_{i s} \phi_{s s^{\prime}} X_{j s^{\prime}}
$$

### **主成分方法**
这里我们考虑相关系数矩阵$\Psi_{ij}$:
$$
\Psi_{i j}=\sum_{s, s^{\prime}=1}^{M} Y_{i s} \phi_{s s^{\prime}} Y_{j s^{\prime}}
$$
其中$Y_{i s}=X_{i s} / \sigma_{i}$，$\sigma_i$取自$C_{ij}$的对角线。
如果我们直接使用上述$\phi_{ss^\prime}$作为我们的相关系数矩阵，我们就没有特质方差。

因此我们选择使用主成分分析的方法进行建模，记$V_i^{(a)},a=1,...,N$是$\Psi_{i j}$的主成分：
$$
\begin{aligned}
&\sum_{j=1}^{N} \Psi_{i j} V_{j}^{(a)}=\lambda^{(a)} V_{i}^{(a)} \\
&\sum_{i=1}^{N} V_{i}^{(a)} V_{i}^{(b)}=\delta_{a b}
\end{aligned}
$$
并且满足$\lambda^{(1)}>\lambda^{(2)}>\cdots$,因而有
$$
\Psi_{i j}=\sum_{a=1}^{M} V_{i}^{(a)} \lambda^{(a)} V_{j}^{(a)}
$$
接下来我们选取K<M得到新的因子模型如下：
$$
\begin{aligned}
&\Gamma_{i j}=\xi_{i}^{2} \delta_{i j}+\sum_{A=1}^{K} \lambda^{(A)} V_{i}^{(A)} V_{j}^{(A)} \\
&\xi_{i}^{2}=1-\sum_{A=1}^{K} \lambda^{(A)}\left(V_{i}^{(A)}\right)^{2}
\end{aligned}
$$
那么我们写成因子模型的形式如下：
$$
\Gamma_{i j}=\xi_{i}^{2} \delta_{i j}+\sum_{A, B=1}^{K} \Omega_{i A} \Phi_{A B} \Omega_{j B}
\\
\begin{aligned}
&\Omega_{i A}=\sqrt{\lambda^{(A)}} V_{i}^{(A)}, \quad A=1, \ldots, K \\
&\Phi_{A B}=\delta_{A B}.
\end{aligned}
$$