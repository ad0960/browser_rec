# browser_rec
浏览器网址导航推荐功能

**主要内容**

1. 使用Pig执行一些数据ETL操作，得到训练集。Pig Latin类似SQL语句，比用MapReduce编写更加方便，更多资料可以看[Pig][1]
2. 使用LBFGS算法求解逻辑回归，相比SGD算法迭代次数从10000次降到14次，想要对LBFGS和LineSearch等算法有更深入了解请看[Numerical Optimization][2]
3. 使用AUC作为二分类算法的评估指标，相比accuracy、precision、recall等更完善，更多资料可以看[AUC][3]


**核心代码**

- LBFGS.scala  LBFGS算法  
- StrongWolfe.scala  LBFGS中使用的LineSearch方法  
- AreaUnderCurve.scala  计算AUC
- script.pig  提供了一些Pig操作的简单例子




[1]:http://www.codelast.com/?p=3621  "Pig"
[2]:www.bioinfo.org.cn/~wangchao/maa/Numerical_Optimization.pdf  "Numerical Optimization"
[3]:http://alexkong.net/2013/06/introduction-to-auc-and-roc/  "AUC"
