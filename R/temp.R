
##--------------------------------------
#read paper table
ptable<-function(table_name="table2"){
  mr<-read.xlsx("E:\\weiyunSync\\博士-赵青松\\博士论文\\大豆异黄酮检测\\Zhao et al. -2022-01-22-tables-isoflavone.xlsx",table_name,startRow=3)
  mr<-modify_if(mr, ~is.numeric(.), ~round(.,2))
  return(mr)
}
##--------------------------------------
sname<-function(){
  mr<-getSheetNames("E:\\weiyunSync\\博士-赵青松\\博士论文\\大豆异黄酮检测\\Zhao et al. -2022-01-22-tables-isoflavone.xlsx")
  return(mr)
}






