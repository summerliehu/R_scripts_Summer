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
cancerstudyID=104
mycancerstudy = getCancerStudies(mycgds)[cancerstudyID,1];
mycaselist = getCaseLists(mycgds, mycancerstudy)[1,1];
mycaselist;
#---------------------
#查看数据都包括哪些类型
mygeneticprofile = getGeneticProfiles(mycgds, mycancerstudy)[,1];
#---------------------
#选择数据类型