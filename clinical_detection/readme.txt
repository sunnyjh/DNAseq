# NGS临床检测的流程

## 1.数据拆分（参考链接：https://www.jianshu.com/p/0eaa6bce82b2 ）

### 1.1 文库结构
    ![](文库结构图.png)
    
### 1.2 fastq数据的一级拆分
    BCL转换成fastq后进行，基于index将fastq数据进行一级拆分，得到子文库类型的fastq文件，该子文库包含多个样本的fastq数据，这么做的原因是为了保证GC%在正常范围内。具体原理如下：
    根据每条lane所有子文库index序列的比对结果来确定拆分是否容错：当同lane index间最小差异碱基数大于等于3时，允许1个碱基错配，其余情况不容错
    

### 1.3 fastq数据的二级拆分（针对子文库的fastq数据）
    二级拆分后子文库fastq数据被拆分成单个样本的fastq数据，该阶段不容错，但是还保留二级barcode的序列信息，需要利用fastp进一步切除barcode信息。

### 

## 2.


## 3.


## 4.

## 5.consensus bam
   1.








## 1.位点合并

## 2.SNP突变阴性背景池构建

## 3.
