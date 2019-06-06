# somaticSignatures
## 1.原理

    根据六种碱基突变类型A>T,A>C,A>G,T>C,T>G,G>C及其上下游各一个碱基序列，构成6*4*4组三碱基组合，统计所有三碱基组合的frequency，构成了多样本三碱基frequency频率矩阵，然后根据NMF(非负矩阵分解法）或PCA解析出能够代表多个样本共同特征的三碱基频率分布作为体细胞突变特征。再将该特征与cosmic数据库中已知30个肿瘤的体细胞突变特征进行比较，从而解释候选突变特征表示的含义。

## 2.条件
   
    1）要求多个样本
    2）只使用点突变数据
    
## 3.NMF算法

    对于数据较大时该算法比较耗计算资源。
    The NMF decomposes M with the constraint of positive components in W and H [7]. The method was used [11] for the identification of mutational signatures, and can be computationally expensive for large data sets.
    https://www.cnblogs.com/pinard/p/6812011.html?utm_source=itdadao&utm_medium=referral

## 4.PCA算法

    The PCA employs the eigenvalue decomposition and is therefore suitable for large data sets [13]. While this is related to the NMF, no constraint on the sign of the elements of W and H exists.

## 5.评估最优突变特征数目

    使用assessNumberSignatures探索我们到底应该把ｓｐｅｃｔｒｕｍ分成多少个ｓｉｇｎａｔｕｒｅ
    n_sigs = 2:8
    gof_nmf = assessNumberSignatures(sca_mm, n_sigs, nReplicates = 5)
    gof_pca = assessNumberSignatures(sca_mm, n_sigs, pcaDecomposition)
    plotNumberSignatures(gof_nmf)　## 可视化展现

## 6.结果可视化

http://www.bioconductor.org/packages/release/bioc/vignettes/SomaticSignatures/inst/doc/SomaticSignatures-vignette.html
http://www.bio-info-trainee.com/1623.html
