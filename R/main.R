# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'



library("stringr")
library("ggplot2")
###f1 get marker with Mb
get_marker<-function(){
  mr<-SoyBase_marker
  return(mr)
}

###f2 get_qtl
get_qtl<-function(){
  mr<-SoyBase_qtl
  return(mr)
}



##f3get
get_marker_with_qtl<-function(){
  smarker<-get_marker()
  smarker$phy.loc<-paste("Con",str_remove(smarker$Chromosome,"Gm"),"_",smarker$Start,sep="")
  sqtl<-get_qtl()
  msqtl<-merge(sqtl,smarker,by.x="Object.Name",by.y="Marker.Name",all.x=TRUE)
  #msqtl$Object.Type[!duplicated(msqtl$Object.Type)]
  nse<-c( "Isozyme","PCR",  "RAPD",   "Marker", "RFLP",  "AFLP")
  mr<-subset(msqtl,!msqtl$Object.Type%in%nse)
  ch.lg<-ch_lg[,1:2]
  mr<-merge(mr,ch.lg,by.x="LG",by.y="Linkage_Group",all.x=TRUE)
  mr<-mr[order(mr$LG,mr$Start.cM,mr$Stop.cM),]
  return(mr)
}


#get consensus markers only
get_consensus_marker_only<-function(){
  mq<-get_marker_with_qtl()
  markers<-subset(mq,mq$Start>=0)
  return(markers)
}
#


#f4get some traits with consensus markers
get_traits_with_consensus_markers<-function(traits=c("oil","protein","isoflavone")){
  traits<-str_to_lower(traits)
  mq<-get_marker_with_qtl()

    allqtls<-NULL
  for(iname in traits){
    qtls<-mq[str_detect(str_to_lower(mq$Object.Name),iname),]
    allqtls<-rbind(allqtls,qtls)
  }
  markers<-subset(mq,mq$Start>=0)
  mrmq<-rbind(allqtls,markers)
  mrmq<-mrmq[order(mrmq$Chromosome_Number,mrmq$Start.cM,mrmq$Stop.cM),]
  return(mrmq)
}


#get selected traits only
get_traits_only<-function(traits=c("oil","protein","isoflavone")){
  traits<-str_to_lower(traits)
  mq<-get_marker_with_qtl()
  allqtls<-NULL
  for(iname in traits){
    qtls<-mq[str_detect(str_to_lower(mq$Object.Name),iname),]
    allqtls<-rbind(allqtls,qtls)
  }
  return(allqtls)
}
### get all genes in composite map
get_gene_only<-function(){
  mq<-get_marker_with_qtl()
  mr<-subset(mq,mq$Object.Type=="Gene")
  return(mr)
  }

#####################
get_jy_marker<-function(mymarker="06_16192915"){
  mk<-jymap3943
  smk<-mk[stringr::str_detect(mk$marker,mymarker),]
  return(smk)
}
##############
compare_v2v1<-function(chr="Gm06"){
  v1<-SoyBase_marker_v1
  v2<-SoyBase_marker
  #head(v2)
  mv2<-merge(v2,v1,by.x="Marker.Name",by.y="Marker.Name")
  mv2<-mv2[order(mv2$Chromosome.x,mv2$Start.x),]
  mv2$diff_mb<-(mv2$Start.x-mv2$Start.y)/1000000
  mr<-mv2[mv2$Chromosome.x==chr,]
  plot(mr$Start.x, mr$diff_mb)
  return(mr)
}

##########################################
#repair the jy map,fill the null value
#if up 5 and down 5 are the same ,then fill
#######################################
if_fill_A<-function(vac=c(rep("A",5),"-",rep("A",5))){
  if_has_A<-any(vac=="A")
  mr<-any(vac=="B")
  mr<-!mr&&if_has_A
  #if don't detect B but has A,then return TRUE,that mean you can fill A
  return(mr)
}
#
if_fill_B<-function(vac=c(rep("B",5),"-",rep("B",5))){
  if_has_B<-any(vac=="B")
  mr<-any(vac=="A")
  mr<-!mr&&if_has_B
  return(mr)
}

#fill one line
fill_one_line<-function(oneline=jymarker3943$JY135,mylen=3943){
  myloc<-oneline
  loci<-which(myloc=="-")
  myjA<-NULL
  myjB<-NULL
  for(i in loci){
    sti<-i-5
    eni<-i+5
    if(sti<0) sti<-0
    if(eni>mylen) eni<-mylen
    myjA<-c(myjA,if_fill_A(myloc[sti:eni]))
    myjB<-c(myjB,if_fill_B(myloc[sti:eni]))
  }
  myloc[loci][myjA]<-"A"
  myloc[loci][myjB]<-"B"
  return(myloc)
}
#
jyline_name<-function(){
  n1<-paste("JY00",1:9,sep="")
  n10<-paste("JY0",10:99,sep="")
  n100<-paste("JY",100:185,sep="")
  mr<-c(n1,n10,n100)
  return(mr)
}

############
jymarker3943_filled<-function(mydata=jymarker3943){
  for(iname in jyline_name()){
    mydata[,iname]<-fill_one_line(mydata[,iname])
  }
  return(mydata)
}

##----------------------------------------------------------------
##extract data from new near-inferred spectrometer.
extract_data_new_NI<-function(myd){
  md1<-dcast(myd,样品名~组分,mean,na.rm = TRUE,value.var = "预测")
  md1$'样品名'<-str_split_fixed(md1$'样品名',";",3)[,1]
  names(md1)<-str_replace(names(md1)," ","-")
  return(md1)
}
######----------------------------------------------------------------



