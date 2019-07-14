# CNVseq方法介绍

## 1.WGS

### 1.1 低深度WGS：5x以下，分辨率低，对于大的CNV敏感度高，但成本低、样本周转周期短、对计算机硬性条件低
    分辨率在几十kb以上，因为了保证每个窗口有足够的reads，一般设置窗口大小为20kb，这样分辨率为60kb，对于小的CNV检测不到。

### 1.2 高深度WGS：5x以上，如30X或50X,分辨率高，敏感度高，问题是成本高、样本周转周期长（测序、分析等）、对计算机硬性条件和分析人员要求高。

### 1.3 Array VS low WGS vs high WGS（比较low WGS是否能代替Array）
    参考链接：https://mp.weixin.qq.com/s?__biz=MzA4MDQ0NzM1MQ==&mid=2650190427&idx=1&sn=2e6ba5f4bbce220fbcbdedd782c05f69&chksm=87a60043b0d1895511b18bebd8927db0990b055a49c4c5f7bb3d00f460642ee690bfbf8fa0d5&mpshare=1&scene=1&srcid=0618CdFLteORM9Uq0fm9KD6U&from=singlemessage&clicktime=1563013462&ascene=7&devicetype=android-24&version=27000536&nettype=WIFI&abtest_cookie=BgABAAgACgALABIAEwAVAAYAnYYeACOXHgBWmR4AzpkeAPaZHgAMmh4AAAA%3D&lang=zh_CN&pass_ticket=Qqq9pMkDQZov5Lp7nyDrhZrjSqrruJ2qB1B43ThVXkvNTKO%2Bz8DsZZsZ%2BuHs9om%2B&wx_header=1
    参考文献：Whole-genome sequencing analysis of genomic copy number variation (CNV) using low-coverage and paired-end strategies is highly efficient and outperforms array-based CNV analysis

    1)low WGS vs Array: 
    对多种不同检测CNV的低深度WGS方法进行综合分析，同时研究了350 bp 短插入文库构建法WGS（short-insert）,3kb配对文库构建法WGS（3kb-insert mate-pair）, 及5kb配对文库构建法WGS（ 5kb-insert mate-pair）分别在1X, 3X, 及 5X测序深度时检测CNV的有效性。以NA12878样本中经验证的CNV作为分析的金标准（Gold Standard ，GS)，用以确定这些CNV是否能在每种低深度WGS方法中检测到。为比较金标准确定的CNV（GS-CNV），我们采取了Haraksingh et al 23在研究高密度寡聚物阵列（high-density oligomer arrays）时的标准方法,此后我们将低测序深度的WGS与之前arrays的研究结果进行比较
    *CNV分析方法：
     一、read depth：此种方法通过CNVnator (version 0.2.7)来实现，算法过程中大小为5000bp的二进制文件用于生成柱状图及统计分析等，产出的CNV数据将与从UCSC下载的人类参考基因组hg19进行比对，接着会滤过大于300kb且与基因间隙区重叠超过50%的数据及X染色体上的主要组织相容性复合物区。二、disconcordant read pairs：此种方法通过 LUMPY (version 0.6.11)来实现。从UCSC Genome Browser下载的segmentation duplication和reference gaps从分析的基因组区域去掉，仅筛选出缺失及重复CNV用于后续分析。大于150kb的CNV也被滤除。
    *CNV分析标准：
    合并read depth和disconcordant read pair法检出的CNV。如果两种方法检测出来的CNV有至少50%的重叠区，则由disconcordant read pairs法检测出的数据代替。与Haraksingh 等描述的方法一致，合并后的CNV将与由NA12878基因组数据设定的金标准CNV进行比较，结果共分两类：（1）有大于或等于50%的相应重叠区；（2）有大于或等于10%却小于50%的相应重叠区。合并后的CNV依照银标准也可划分为两组：（1）有大于50%的相应重叠区；（2）有小于50%的相应重叠区。不同方法检出的CNV的敏感性通过计算得出：即检测到的与金标准CNV大于50%重叠的个数除以先前已鉴定的所有金标准的CNV个数。
    *Array数据的分析：
    所有样本的array测序数据获得后，仅信号>10的连续探针及最大Log BAF > 10 的CNV数据被用于分析。
    
    *总结：
    low WGS比array敏感度高，在权衡短插入文库构建与长插入文库构建的利弊后，控制每份样本的花费及测序周转时间后，low WGS性价比更高。
    ![](low_wgs_vs_array.png)
    
    
    
    

## 2.WES

### 2.1 问题

### 2.2 解决方法

## 3.Panel

### 3.1 问题

### 3.2 解决方法
