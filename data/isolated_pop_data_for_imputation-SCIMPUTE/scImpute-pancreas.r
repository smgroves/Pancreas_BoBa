#Install packages
#install.packages("devtools")
library(devtools)
#install_github("Vivianstats/scImpute")
setwd('~/Documents/GitHub/Pancreas_BoBa/data/isolated_pop_data_for_imputation-SCIMPUTE/')
library(scImpute)

for (fil in list.files(path='.',pattern='*.csv')){
  print(fil)
  scimpute(count_path = fil,ncores=1,Kcluster=1,drop_thre=0.5,out_dir=paste0('./scimpute-ouput_',strsplit(fil,".csv")[1],'/'))
}
