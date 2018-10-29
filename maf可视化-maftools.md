# [MAF可视化-maftools](http://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html)

MAF为常用的体细胞突变注释文件格式，可通过maftools对MAF文件进行可视化操作（如果希望可视化vcf文件，可将vcf转换成MAF，[参考链接](https://github.com/mskcc/vcf2maf))，如下：
 - 读入annovar文件转换为maf——annovarToMaf
 - 读入maf文件——read.maf
 - 绘制maf文件的摘要——plotmafSummary
 - 绘制瀑布图——oncoplots
 - 绘制箱线图—— titv
 - 绘制棒棒图——lollipopPlot

## 1.maftools 安装
    source("http://bioconductor.org/biocLite.R")
    biocLite("maftools")
    library(maftools)
   
## 2.maftools功能
### 2.1 读入annovar文件转换为maf——annovarToMaf
    #produce maf
    var.annovar.maf = annovarToMaf(annovar = "all_annovar3", Center = 'NA', refBuild = 'hg19', tsbCol = 'Tumor_Sample_Barcode', table = 'ensGene',sep = "\t")
    write.table(x=var.annovar.maf,file="var_annovar_maf",quote= F,sep="\t",row.names=F)
    
### 2.2 绘制maf文件的摘要——plotmafSummary
     #先去掉NA或者unkown的突变
     sed 's/^NA/unknown/' var_annovar_maf > var_annovar_maf2
     grep -v "^NA" var_annovar_maf | grep -v -P "\tUNKNOWN\t"> var_annovar_maf2

     var_maf = read.maf(maf ="var_annovar_maf2")
     plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median'))

![](https://github.com/sunnyjh/DNAseq/edit/master/images/maftools_plotsummary.png)
### 2.3 绘制瀑布图——oncoplots]
     oncoplot(maf = var_maf, top = 400, writeMatrix=T,removeNonMutated = F,showTumorSampleBarcodes=T)
    
### 2.4 绘制箱线图—— titv
### 2.5 绘制棒棒图——lollipopPlot
