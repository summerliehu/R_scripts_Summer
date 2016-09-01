#需要在R中复制粘贴运行
#编辑以下进行配置
gene1="DNMT1";
gene2="MAP2K7";
#----------------------
library(cgdsr);
mycgds = CGDS("http://www.cbioportal.org/");
#test(mycgds);
getCancerStudies(mycgds);
#----------------------
cancerstudyID=99;
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

ggplot(mutation_rna, aes(DNMT1, MAP2K7.1)) + geom_boxplot()+geom_point(color="darkblue", alpha = 0.2)+
labs(x="DNMT1 genotype",y="MAP2K7 expression",title = "Correlation of DNMT1 genotype and MAP2K7 mRNA");

#----------------------
#计算p-value
t.test(MAP2K7.1 ~ DNMT1, mutation_rna)
#将p-value添加到图中
qplot(x=DNMT1, y=MAP2K7.1, data=mutation_rna, fill=DNMT1, geom=c("boxplot", "jitter"), xlab="DNMT1 mutation type", 
	ylab="MAP2K7 mRNA level",main="DNMT1 mutation regulates MAP2K7 RNA in coadread_tcga_pub")+
	annotate("text",x =1.5 , y = 10, label = "p-value: 0.50", size=5)

#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH下面的代码为 表达-表达 关系，做出点图的代码HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
#选择数据类型
type1ID=4;
rnaprofile1 = getGeneticProfiles(mycgds, mycancerstudy)[type1ID,1];
rna1 <- getProfileData(mycgds, c(gene1,gene2), rnaprofile1,mycaselist);
#输出为csv
write.csv(rna1, file = "/Users/xiayukun/Desktop/mutation_rna.csv");
#------------------------------------------------------------------------------------------
#清洗csv数据
rna1 <- read.csv("/Users/xiayukun/Desktop/mutation_rna.csv");
#mutation_rna;
#----------------------
#画图
library(ggplot2)
qplot(x=DNMT1, y=MAP2K7, data=rna1, 
	main="Correlation betewwn DNMT1 and MAP2K7 expression in NEPC", 
	xlab="DNMT1 expression", ylab="MAP2K7 expression")+ geom_smooth(method = "lm", se = FALSE)

#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
#下面代码使用ggplot作出带有回归方程的两个基因表达关系的图：
type1ID=3;
rnaprofile1 = getGeneticProfiles(mycgds, mycancerstudy)[type1ID,1];
rna1 <- getProfileData(mycgds, c(gene1,gene2), rnaprofile1,mycaselist);
#输出为csv
#write.csv(rna1, file = "/Users/xiayukun/Desktop/mutation_rna.csv");
#清洗csv数据
#rna1 <- read.csv("/Users/xiayukun/Desktop/mutation_rna.csv");

library(gridExtra);

lm_eqn = function(m) {

  l <- list(a = format(coef(m)[1], digits = 2),
      b = format(abs(coef(m)[2]), digits = 2),
      r2 = format(summary(m)$r.squared, digits = 3));

  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }

  as.character(as.expression(eq));                 
}

p <- ggplot(rna1, aes(DNMT1, MAP2K7)) + geom_point() + 
geom_smooth(method="lm", colour="darkblue", size=1) +
geom_text(aes(x = 100, y = 3, label = lm_eqn(lm(MAP2K7 ~ DNMT1, rna1))), parse = TRUE)+
labs(x="DNMT1 expression",y="MAP2K7 expression",title = "Correlation of DNMT1 and MAP2K7 mRNA in GBM(TCGA)");