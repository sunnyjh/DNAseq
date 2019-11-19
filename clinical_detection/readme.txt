# NGS临床检测的流程

## 1.数据拆分（参考链接：https://www.jianshu.com/p/0eaa6bce82b2 ）

### 1.1 文库结构
    ![](library.png)
    
    说明如下：
    1)P5和P7是接头序列
    2）I1和I2是index序列，长度都为7nt，用于区分同一条lane的不同子文库。
    3）B1和B2是Barcode序列，又称UMI，长度为8nt，用来区分同一子文库中不同的样本。
    4）R1和R2包括barcode序列信息，bcl2fastq按照I1和I2拆分得到子文库，split_fastq按照B1和B2拆分得到样本fastq，此时fastq包含Barcode序列，利用fastp将Barcode序列从R1和R2中移除，并放到read ID上。如下图：
    ![](barcode_umi.png)
    
    
    
### 1.2 fastq数据的一级拆分
    BCL转换成fastq后进行，基于index将fastq数据进行一级拆分，得到子文库类型的fastq文件，该子文库包含多个样本的fastq数据，这么做的原因是为了保证GC%在正常范围内。具体原理如下：
    根据每条lane所有子文库index序列的比对结果来确定拆分是否容错：当同lane index间最小差异碱基数大于等于3时，允许1个碱基错配，其余情况不容错
    
### 1.3 fastq数据的二级拆分（针对子文库的fastq数据）
    二级拆分后子文库fastq数据被拆分成单个样本的fastq数据，该阶段不容错，但是还保留二级barcode的序列信息。

### 1.4 fastp切除barcode

## 2.


## 3.


## 4.

## 5.去重和校正（onsensus bam）
   1.








## 1.位点合并

## 2.SNP突变阴性背景池构建

## 3.
