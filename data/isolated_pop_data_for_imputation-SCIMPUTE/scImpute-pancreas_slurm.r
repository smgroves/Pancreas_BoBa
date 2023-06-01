#Load packages
.libPaths("~/R/rlib-3.4.3")
library(devtools)
library(scImpute)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
#Get list of csv files in current directory
fil = list.files(path='.',pattern='*_scimpute.csv')
#Loop through all files and run scimpute on the correct file
for (f in seq_along(fil)){
  print(fil[f])
  if (f==args[1]){
    scimpute(count_path = fil[f],ncores=1,Kcluster=1,drop_thre=0.5,out_dir=paste0('./scimpute-ouput_',strsplit(fil[f],".csv")[1],'/'))
  }
}
print('Success')
