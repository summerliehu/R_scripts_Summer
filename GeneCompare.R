#需要在R中复制粘贴运行
#编辑以下进行配置
gene1="MAP2K7";
gene2="TET2";
#----------------------
library(cgdsr);
mycgds = CGDS("http://www.cbioportal.org/");
#test(mycgds);
getCancerStudies(mycgds);
#----------------------
cancerstudyID=36;
mycancerstudy = getCancerStudies(mycgds)[cancerstudyID,1];
mycaselist = getCaseLists(mycgds, mycancerstudy)[,1];
mycaselist;
#---------------------
#填写数据集
mycaselistID=1;
mycaselist = getCaseLists(mycgds, mycancerstudy)[mycaselistID,1];
mygeneticprofile = getGeneticProfiles(mycgds, mycancerstudy)[,1];
mygeneticprofile;
#---------------------

#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH下面的代码为 突变-表达 关系，做出箱图的代码HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
#选择数据类型
type1ID=1;
type2ID=4;
mutationprofile = getGeneticProfiles(mycgds, mycancerstudy)[type1ID,1];
mutation <- getProfileData(mycgds, c(gene1,gene2), mutationprofile,mycaselist);
rnaprofile=getGeneticProfiles(mycgds, mycancerstudy)[type2ID,1];
rna=getProfileData(mycgds, c(gene1,gene2), rnaprofile,mycaselist);
mutation_rna <- cbind(mutation, rna);
mutation_rna;
#输出为csv
write.csv(mutation_rna, file = "/Users/xiayukun/Desktop/mutation_rna.csv");
#------------------------------------------------------------------------------------------
#清洗csv数据
mutation_rna <- read.csv("/Users/xiayukun/Desktop/mutation_rna.csv");
#mutation_rna;
#----------------------
#画图
library(ggplot2)
qplot(x=MAP2K7, y=TET2.1, data=mutation_rna, fill=MAP2K7, geom=c("boxplot"),
	 xlab="KRAS mutation type", ylab="TET2 mRNA level",main="Effect of MAP2K7 mutation on mRNA level of TET2 in coadread_tcga_pub");
#----------------------
#计算p-value
t.test(TET2.1 ~ MAP2K7, mutation_rna)
#将p-value添加到图中
qplot(x=MAP2K7, y=TET2.1, data=mutation_rna, fill=MAP2K7, geom=c("boxplot", "jitter"), xlab="MAP2K7 mutation type", 
	ylab="TET2 mRNA level",main="MAP2K7 mutation regulates TET2 RNA in coadread_tcga_pub")+
	annotate("text",x =1.5 , y = 10, label = "p-value: 0.50", size=5)

#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH下面的代码为 表达-表达 关系，做出点图的代码HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
#选择数据类型
type1ID=4;
type2ID=4;
rnaprofile1 = getGeneticProfiles(mycgds, mycancerstudy)[type1ID,1];
rna1 <- getProfileData(mycgds, c(gene1,gene2), mutationprofile,mycaselist);
rnaprofile2=getGeneticProfiles(mycgds, mycancerstudy)[type2ID,1];
rna2=getProfileData(mycgds, c(gene1,gene2), rnaprofile,mycaselist);
rna1_rna2 <- cbind(rna1, rna2);
rna1_rna2;
#输出为csv
write.csv(rna1_rna2, file = "/Users/xiayukun/Desktop/mutation_rna.csv");
#------------------------------------------------------------------------------------------
#清洗csv数据
rna1_rna2 <- read.csv("/Users/xiayukun/Desktop/mutation_rna.csv");
#mutation_rna;
#----------------------
#画图
library(ggplot2)
qplot(x=MAP2K7, y=TET2, data=rna1_rna2, 
	main="correlation betewwn MAP2K7 and TET2 expression", 
	xlab="MAP2K7 expression", ylab="TET2 expression")+ geom_smooth(method = "lm", se = FALSE)

#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH