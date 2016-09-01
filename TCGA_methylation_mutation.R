# read max-stddev data
#首先需要用Excel把下下来的数据打开，删除第二行，然后把数据存成.csv的格式，在R中读取csv
setwd('/Users/xiayukun/Desktop/')
methylation <- read.csv('GBM.meth.by_max_stddev.data.csv')
data=na.omit(methylation)

library('gplots')
heatmap.2(as.matrix(methylation[,-1]), col=redgreen)