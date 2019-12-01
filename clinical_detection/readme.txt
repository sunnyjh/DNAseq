# NGS临床检测的流程

## 1.[数据拆分](https://www.jianshu.com/p/0eaa6bce82b2)

### 1.1 文库结构

![](library.png)
    
    说明如下：
    1)P5和P7是接头序列
    2）I1和I2是index序列，长度都为7nt，用于区分同一条lane的不同子文库。
    3）B1和B2是Barcode序列，又称UMI，长度为8nt，用来区分同一子文库中不同的样本。
    4）R1和R2包括barcode序列信息，bcl2fastq按照I1和I2拆分得到子文库，split_fastq按照B1和B2拆分得到样本fastq，此时fastq包含Barcode序列，利用fastp将Barcode序列从R1和R2中移除，并放到read ID上。如下图：

![](barcode_umi.png)
    
    
### 1.2 fastq数据的一级拆分
     BCL转换成fastq后进行，基于index将fastq数据进行一级拆分，得到子文库类型的fastq文件，该子文库包含多个样本的fastq数据，这么做的原因是为了保证GC平衡，另一份方面也是为了保证一条lane能上足够多的样本。具体原理如下：
    根据每条lane所有子文库index序列的比对结果来确定拆分是否容错：当同lane index间最小差异碱基数大于等于3时，允许1个碱基错配，其余情况不容错
    
### 1.3 fastq数据的二级拆分（针对子文库的fastq数据）
    根据barcode信息将子文库fastq数据被拆分成单个样本的fastq数据，该阶段不容错，拆分后barcode序列还是保留着。

### 1.4 fastp切除barcode
    fastp去掉低质量reads和接头，然后将barcode序列从read序列上切除，放在read ID上，如UMI_ATGCTAGG_GTCAGTAA

## 2.bwa比对、排序和合并
    1）Bwa mem比对，过滤低质量比对reads(-q 10), 过滤未比对reads（-F 4）。
    2）Samtools merge对比对结果进行合并。

## 3.去重和校正（consensus bam）
   一种barcode包括15种UMI(长度为8bp），双端UMI一共为15*15=225组合，利用[gencore](https://github.com/OpenGene/gencore)进行去重和校正：
   1）去重：比对位置起始、终止位置都相同，UMI相同的reads为duplication
   2）校正：利用高深度测序中的duplication对PCR和测序错误进行校正，具体算法参考gencore。
   3）对于肿瘤样本，使用gencore按照同barcode，同起始，同终止的原则对*.sort.bam进行去重与碱基校正。默认血浆使用>=2x raw reads支持，alt_ratio>=0.8；组织使用>=1x raw reads支持，alt_ratio >=0.9。对于白细胞对照样本，使用sambamba进行去重。
   
   gencore的算法：
   1）clusters the reads by their mapping positions and UMIs (if UMIs are applicable).
   2）for each cluster, compares its supporting reads number (the number of reads/pairs for this DNA fragment) with the threshold specified by supporting_reads. If it passes, start to generate a consensus read for it.
   3）if the reads are paired, finds the overlapped region of each pair, and scores the bases in the overlapped regions according their concordance and base quality.
   4）for each base position at this cluster, computes the total scores of each different nucleotide (A/T/C/G/N).
   5）if there exists a major nucleotide with good quality, use this nucleotide for this position; otherwise, check the reference nucleotide from reference genome (if reference is specified).
   6）when checking the reference, if there exists one or more reads are concordant with reference genome with high quality, or all reads at this positions are with low quality, use the reference nucleotide for this position.

  
## 4.变异检测
    1）samtools mpileup（-B -q 20 -Q 0）分别建立肿瘤及白细胞样本pileup文件。
    2）MutLoci检出原始突变。白细胞保留所有突变用作对照，肿瘤样本保留初筛阈值线以上突变（vaf >=0.001/0.01, alt>=2, uniq_depth >=100 ）
    3）链特异性过滤（fisher‘s exact p value<0.001,倍数>=10）。
    4）对于normal中存在alt支持的突变，要求alt<=10, tumor_vaf >= 10 * normal_vaf。
    5）位点合并，annovar注释，transvar校正，过滤人群频率>=0.05的点。
    6）alt碱基位置与一致性过滤，热点与非热点检出。对于候选alt位点，过滤碱基位置位于比对首尾及PE reads不一致的reads, 热点tumor_vaf>=0.001/0.01, alt>=3; 对于非热点tumor_vaf>=0.005/0.03, alt>=3。
    7）检出结果一般只保留热点及非同义突变位点。



## 5.阴性背景池构建
    过滤背景信号，457血浆阴性背景池包含27个正常人的血浆样本，457组织阴性背景池包含70个正常人的组织样本，31组织阴性背景池包含30个正常人的组织样本，构建方法如下：
### 5.1 得到consensus bam
    同上方式
### 5.2 检测单个样本的背景突变
    samtools mpileup -f hg19.fasta -l LC_CRC_456.cnv.bed -d 100000  -B -q 20 -Q 0 -a YF1918P.cons.bam | python getBackgroundMutation.v2.1.py YF1918P.bgm.cons.xls

### 5.3 合并多个样本的背景突变
    python combineMutation.py conf.xls cons.total.bgm.xls
    
    conf.xls格式为：样本名 样本名..bgm.cons.xls。如下：
        YF1918P YF1918P.bgm.cons.xls

### 5.3 背景模型拟合
    利用[fitdistrplus](https://www.cnblogs.com/ywliao/p/6297162.html)包判断背景突变的vaf符合哪种统计学模型，然后通过qqplot判断实际数据与该统计学模型数据的相关性大小和pvalue值。([判断数据是否服从某一分布（一）](https://www.cnblogs.com/ywliao/p/6265945.html))
    perl FitBGMModel.pl -i cons.total.bgm.xls -o cons.total.bgm.xls.fit

### 5.4 利用背景突变过滤掉背景噪音
     perl filterBGM.pl -i annovar.combined.xls.strand -n info.txt -o annovar.combined.xls.test -bp sort.total.bgm.xls.fit -bt sort.total.bgm.xls.fit
     当只有阴性背景池中一个样本有突变，采用的二项分布检验；当只有阴性背景池中2-4个样本有突变，采用的z检测；其他情况为weibull分布。
     
## 6.变异过滤

### 6.1 链偏好性过滤[getStrandInfo.py]
    过滤条件：1）plus_alt = 0 或者 minus_alt = 0
             2）plus_alf_vaf = plus_alt/(plus_ref+plus_alt),minus_alf_vaf = minus_alt/(minus_ref+minus_alt),plus_alt_vaf/minus_alt_vaf>10或者1/10.

### 6.2 

## panel设计
    优化panel设计算法，提高各种突变类型的检测性能
        Indel 区（关联cosmic数据库,掺入突变型探针）
        Fusion区（增加探针密度，并找到最佳间隔距离）
        低GC区（增加探针数量，并找到增加最佳值）
        重复区（mismatch对捕获效率的影响，去除重复区的原则 ） 



## 1.位点合并

## 2.SNP突变阴性背景池构建

## 3.


## 参考链接
[gencore](https://github.com/OpenGene/gencore)
[fitdistrplus](https://www.cnblogs.com/ywliao/p/6297162.html)
[weibull distribution](http://reliawiki.org/index.php/The_Weibull_Distribution)
