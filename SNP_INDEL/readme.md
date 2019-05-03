# 突变检测软件原理研究

## 1.[unifiedgenotyper](https://gatkforums.broadinstitute.org/gatk/discussion/1237/how-unified-genotyper-works-retired)
   https://software.broadinstitute.org/gatk/download/archive
   https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef

### 1.1 简介

    UnifiedGenotyper是集合多种变异检测方法而成的一种Variants Caller，既可以用于单个样本的变异检测，也可以用于群体的变异检测。UnifiedGenotyper使用贝叶斯最大似然模型，同时估计基因型和基因频率，最后对每一个样本的每一个变异位点和基因型都会给出一个精确的后验概率。
    UnifiedGenotyper为了保证结果的敏感度，输出的raw vcf结果中具有很高的假阳性，需要用VQSR（variant quality scores reclibration)进行过滤。
    现在该软件已经被haplotypercaller取代，它检测结果比UnifiedGenotyper更准确，特别是indel，而且适用于更大样本量的数据。
    
### 1.2 步骤

    https://gatkforums.broadinstitute.org/gatk/discussion/1237/how-unified-genotyper-works-retired
    https://software.broadinstitute.org/gatk/media/docs/Multiallelic.pdf
    
