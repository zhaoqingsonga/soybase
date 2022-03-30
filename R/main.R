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



#f4get some traits with consensus markers
get_tratis_with_consensus_markers<-function(traits=c("oil","protein","isoflavone")){
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

