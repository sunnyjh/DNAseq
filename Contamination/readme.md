# 样本污染

## 1.样本污染的原因
         
    以HiSeq3000/4000、HiSeq X Ten和NovaSeq为代表的测序平台在运行单index文库时均存在一定程度的样本间数据污染问题。尽管造成污染的可能原因有多种，如接头制备时的交叉污染、样品间的交叉污染、建库实验的交叉污染，以及捕获实验多杂一等，但真正让数据污染问题变得普遍到让人闻风丧胆的“罪魁祸首”，却是一个叫做标签跳跃（Index Hopping）的“新手”。众所周知，为提高测序产出通量，上述测序平台均采用了规则流动槽（Patterned Flow Cell Technology, PFCT）芯片和排他性扩增（Exclusive Amplification, ExAmp）成簇两种新技术，然而利剑有双刃，也正是这两个新技术使得pooling在一起的文库更容易发生标签跳跃，导致标签错配（Index Misassignment），进而造成样品间数据污染。该现象及其原理早在Illumina于2017年4月发布的官方白皮书中就有过详细介绍，但时至今日，这一问题似乎并未引起足够的重视。
    
    参考链接：https://www.illumina.com/science/education/minimizing-index-hopping.html
    
## 2.样本污染类型

### 2.1 样本错配

### 2.2 样本交叉污染



## 参考链接：

    http://www.berrygenomics.com/ke-yan-dong-tai/udi/
    http://www.seqchina.cn/3276.html
    


