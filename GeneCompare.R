#需要在R中复制粘贴运行
#编辑以下进行配置
gene1="BRAF";
gene2="DNMT3A";
#----------------------
library(cgdsr);
mycgds = CGDS("http://www.cbioportal.org/");
test(mycgds);
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
#----------------------
#清洗csv数据
mutation_rna <- read.csv("/Users/xiayukun/Desktop/mutation_rna.csv");
mutation_rna;
#----------------------
#画图
library(ggplot2)
qplot(x=BRAF, y=DNMT3A.1, data=mutation_rna, fill=BRAF, geom=c("boxplot"),
	 xlab="KRAS mutation type", ylab="DNMT1 mRNA level",main="Effect of BRAF mutation on mRNA level of DNMT3A in coadread_tcga_pub");
#----------------------
#计算p-value
t.test(DNMT3A.1 ~ BRAF, mutation_rna)；
#将p-value添加到图中
qplot(x=BRAF, y=DNMT3A.1, data=mutation_rna, fill=BRAF, geom=c("boxplot", "jitter"), xlab="KRAS mutation type", 
	ylab="DNMT1 mRNA level",main="Effect of BRAF mutation on mRNA level of DNMT3A in coadread_tcga_pub")+
	annotate("text",x =1.5 , y = 10, label = "p-value: 0.0032", size=5)