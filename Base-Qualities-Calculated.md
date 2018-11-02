# [Probability of Incorrect Base Call](https://www.ecseq.com/support/ngs/how-are-base-qualities-calculated-and-stored)

## 影响Base calling准确性的因素： 
 - the strength of the light signal 
 - the appearance of the light spot
 - other characteristic features are measured

## 计算Probability of Incorrect Base Call 
 - evaluates the detected light signal for each base call for every cluster, on every tile, for every sequencing cycle 
simultaneous with a sequencing run. 
 - measures various aspects correlating with the quality of the base call like single-to-noise ration and light intensity 
 - Based on these parameters a quality predictor value (QPV) is calculated
