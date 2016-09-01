# read max-stddev data
#首先需要用Excel把下下来的数据打开，删除第二行，然后把数据存成.csv的格式，在R中读取csv
methylation <- read.csv('/Users/xiayukun/Desktop/GBM.meth.by_max_stddev.data.csv')
#将数据框变形
library(reshape2)
methylation_m <- melt(methylation)
#根据两个变量的名称作图
library(ggplot2)
ggplot(data = methylation_m, aes(x=variable, y=Hybridization.REF, fill=value)) + 
    geom_tile()