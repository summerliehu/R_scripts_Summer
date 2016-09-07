# read max-stddev data
#首先需要用Excel把下下来的数据打开，删除第二行，然后把数据存成.csv的格式，在R中读取csv
setwd('/Users/xiayukun/Desktop/')
methylation <- read.csv('GBM.meth.by_max_stddev.data.csv')
data=na.omit(methylation)

library('gplots')
heatmap.2(as.matrix(methylation[,-1]), col=redgreen)

#-----------------

# from LP
# HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

samples <- read.csv("samples.csv", as.is = T)
samples.with.mutations <- samples[, 1]
samples.with.mutations <- gsub("\\.", "-", samples.with.mutations)
CEBPA.samples <- samples.with.mutations[samples[, 2] != ""]
DNMT3A.samples <- samples.with.mutations[samples[, 3] != "" & samples[, 6] =="" & samples[, 4] =="" & samples[, 5] ==""]
IDH.samples <- samples.with.mutations[(samples[, 4] != "" | samples[, 5] != "") & samples[, 3] == "" & samples[, 6] == ""]
TET2.samples <- samples.with.mutations[samples[, 6] != "" & samples[, 3] =="" & samples[, 4] == "" & samples[, 5] == ""]
IDHDNMT3A.samples<-samples.with.mutations[(samples[, 3] != "")&(samples[, 4] != ""|samples[, 5] != "")]
DNMT3ATET2.samples<-samples.with.mutations[samples[, 3] != "" & samples[, 6] != ""]

hm450 <- dir("/Users/peng/Documents/AML/450K/JHU_USC__HumanMethylation450/Level_3")
hm450 <- gsub("A-01D-074[1-3]-05.txt", "", hm450)
hm450 <- gsub("jhu-usc.edu_LAML.HumanMethylation450.2.lvl-3.", "", hm450)

all.hm450.methyl <- NULL
for(i in 1:length(hm450)){
  filename <- dir("/Users/peng/Documents/AML/450K/JHU_USC__HumanMethylation450/Level_3")[i]
  methyl <- read.table(paste0("/Users/peng/Documents/AML/450K/JHU_USC__HumanMethylation450/Level_3/", filename), header = T, sep = "\t", skip = 1)
  all.hm450.methyl <- cbind(all.hm450.methyl, methyl[, 2])
}

all.methyl.data <- data.frame(methyl[, c(1, 3, 4, 5)], all.hm450.methyl)
colnames(all.methyl.data)[-c(1:4)] <- hm450

na.flag <- apply(all.methyl.data[, -c(1:4)], 1, function(x){
  any(is.na(x))
})

sd.methyl <- apply(all.methyl.data[, -c(1:4)], 1, function(x){
  sd(x)
})

color.vector <- rep("white", length(hm450)) 
color.vector[hm450 %in% TET2.samples] <- "black"
color.vector[hm450 %in% IDH.samples] <- "green"
color.vector[hm450 %in% DNMT3A.samples] <- "yellow"
color.vector[hm450 %in% CEBPA.samples] <- "purple"
color.vector[hm450 %in% IDHDNMT3A.samples] <- "blue"
color.vector[hm450 %in% DNMT3ATET2.samples] <- "pink"

library(gplots)
palette.gr.marray <- colorRampPalette(c("blue", "white", "red"))(56)


selected.methyl <- all.methyl.data[!na.flag & sd.methyl > 0.27, ] 

heatmap.2(as.matrix(selected.methyl[, -c(1:4)]), trace = "none", col = palette.gr.marray, symbreaks = F, labRow = NA, labCol = NA, dendrogram = "column", ColSideColors = color.vector)  

legend("bottomleft", c("IDH", "TET2", "CEBPA","DNMT3A","IDH-DNMT3A","DNMT3A-TET2"), fill = c("green", "black", "purple","yellow","blue","pink"), bty = "n", cex = 0.6)