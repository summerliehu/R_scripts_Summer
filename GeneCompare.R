#需要在R中复制粘贴运行
#编辑以下进行配置
gene1="KRAS";
gene2="DNMT1";
#----------------------
library(cgdsr);
mycgds = CGDS("http://www.cbioportal.org/");
test(mycgds);
getCancerStudies(mycgds);
#----------------------
cancerstudyID=104;
mycancerstudy = getCancerStudies(mycgds)[cancerstudyID,1];
mycaselist = getCaseLists(mycgds, mycancerstudy)[,];
mycaselist;
#---------------------
#填写数据集
mycaselist = getCaseLists(mycgds, mycancerstudy)[,1];
#---------------------
#查看数据都包括哪些类型
mygeneticprofile = getGeneticProfiles(mycgds, mycancerstudy)[,1];
#---------------------
#选择数据类型
type1ID=8;
type2ID=4;
#----------------------
mygeneticprofile = getGeneticProfiles(mycgds, mycancerstudy)[type1ID,1];
mutation <- getProfileData(mycgds, c(gene1,gene2), mygeneticprofile,mycaselist);
rnaprofile=getGeneticProfiles(mycgds, mycancerstudy)[type2ID,1];
rna=getProfileData(mycgds, c(gene1,gene2), rnaprofile,mycaselist);
#输出为csv
mutation_rna <- cbind(mutation, rna);
write.csv(mutation_rna, file = "/Users/xiayukun/Desktop/mutation_rna.csv");
#----------------------
#清洗csv数据
mutation_rna <- read.csv("/Users/xiayukun/Desktop/mutation_rna.csv");
mutation_rna;
qplot(x=KRAS, y=DNMT1.1, data=mutation_rna, fill=KRAS, geom=c("boxplot"),
	 xlab="KRAS mutation type", ylab="DNMT1 mRNA level",main="Effect of KRAS mutation on mRNA level of DNMT1");
