# samtools学习计划

## 1.samtools 基本参数

## 2.samtools flags

### 2.1 samtools flagstat
![samtools flagstat](flagstat.png)

    total：bam文件所有alignment数（即bam文件行数）,包括没有比对到参考基因组的、secondary和supplementary alignment，所以总行数大于或等于fastq的read总数。(分为QC-passed reads和 QC-failed reads两类)，与samtools view -c得到结果相同
    secondary: secondary alignment数，由于reads的multiple mapping产生的，与samtools view -f 256 得到结果相同
    supplementary: supplementary alignment数，由于chimeric reads产生的，与samtools view -f 2048得到结果相同
    duplicates:PCR重复？？
    mapped：比对上参考基因组的alignment数(该数占bam文件所有alignment的比列），包含secondary和supplementary alignment，与比对到基因组的read数是不一样的，与samtools view -F 4得到的结果相同。
    paired in sequencing：bam文件中成对的reads总数，包括比对上和没有对上的read（注：不是aligment数）
    read1：bam文件中属于reads1的reads数量
    read2：bam文件中属于reads2的reads数量
    properly paired：比对到基因组上正确配对的reads数量（即read1和read2都比对到同一染色体上，同时insert size在正常范围内。弄清楚代码原理），与samtools view -f 2得到结果相同
    with itself and mate mapped：比对到基因组上配对的reads数，包括read1和read2比对到不同染色体上或者insert size较大的成对的read数。
    singletons：只有单条reads比对上的reads数
    以上计数均以reads条数计，一对reads计为两条。
    
    samtools view -q 与samtools view -F 256区别：
    前者是去除MAPQ低的alignment，不仅包括唯一alignment reads还包括多重alignment reads，但是一般情况下具有多重alignment的read，其大部分alignment的MAPQ比较低。而后者则只是过滤掉具有多重alignment的reads中除primary alignment外的所有secondary alignment。


### 2.2 samtools tview 
    1.DNA重测序比对结果说明
    “.” 比对到正链;
    “，” 表示比对到负链;
    “<”或“>” 表示reference skip   RNA-seq当中内含子剪切;
    "ATCGN"  表示正向mismatch;
    "atcgn"  表示反向mismatch;
    ‘+[0-9]+[ACGTNacgtn]+’ insertion;
    ‘-[0-9]+[ACGTNacgtn]+’ 表示deletion;
    “^”标记reads起始;
    “$”标记reads segment结尾;
    
    2.以下为RNA转录组比对结果展示：

    正链:
        CATCACTGGTTTAAAGACAAACTTGCATTGTGAGATTCCAAAATAACAACAACAAAAAACAATTTGCATTGAGAACATTTTGAAG

        .........A.......
        .....>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        .....>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
    负链：
        TTTCATTTGCAAGTAATCGATTTAGGTTTTTGATTTTAGGGTTTTTTTTTGTTTTGAACAGTCCAGTCAAAGTACAAATCGAGAG
        ...KK....KKK..KK.K.K...K........K....K..................KKKK.........K...K...........
        <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,,,,,,,,t,,,c,,,,,,,,,,
        
        

    参考链接：https://www.omicsclass.com/article/416

## 3.samtools cigar

## 4.samtools tags 

## 5.samtools source code 


