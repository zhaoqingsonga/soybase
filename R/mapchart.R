##get the necessary marker, marker with QTLs
get_locus_marker<-function(mychoosed_data=mychoosed){
  mym<-unlist(str_split(mychoosed_data$Locus,"-"))
  mym<-mym[mym!=""]
  return(mym)
}
## get the reduced marker To draw finely
get_reduced_marker<-function(mymap=jymap3943,distance=20){
  #reduce marker
  choosedmarker2<-mymap$marker[floor((mymap$cm/distance))-c(floor(mymap$cm/distance)[-1],NA)!=0]
  #0 location
  choosedmarker3<-mymap$marker[mymap$cm==0]
  myr01<-as.character(choosedmarker2)
  myr02<-as.character(choosedmarker3)
  myr<-c(myr01,myr02)
  myr<-myr[!is.na(myr)]
  return(myr)
}
### merge select marker and bolt the locus marker
###
merge_selected_marker<-function(mychoosed_data=mychoosed,mymap=jymap3943,distance=30){
  marker1<-get_locus_marker(mychoosed_data)
  marker2<-get_reduced_marker(mymap,distance)
  mym<-c(marker1,marker2)
  mym<-mym[!duplicated(mym)]
  mr<-subset(mymap,mymap$marker%in%mym)
  mr$i<-NA
  mr$b<-NA
  mr$c<-NA
  mr$b[mr$marker%in%marker1]<-"b"
  return(mr)
}

deal_qtl<-function(mychoosed_data=mychoosed){
  mychoosed_data$i<-"i"
  mychoosed_data$b<-"b"
  mychoosed_data$c<-"c2"
  mychoosed_data$c[mychoosed_data$LOD<3]<-"c1"
  mychoosed_data$c[mychoosed_data$LOD>5]<-"c3"
  return(mychoosed_data)
}
############
mapchart_format<-function(mychoosed_data=mychoosed,mymap=jymap3943,distance=100,filename="mymapchart.mct"){
  markers<-merge_selected_marker(mychoosed_data,mymap,distance)
  qtls<-deal_qtl(mychoosed_data)
  sm<-split(markers,markers$group)
  write("\t",filename)
  for(i in 1:20){
    #marker
    write(paste("group",i),filename,append = TRUE)
    write.table(sm[[i]][c("marker","cm","i","b","c","h")],filename,
                quote=FALSE,col.names=FALSE,na="",
                row.names=FALSE,append=TRUE)

    iqtl<-qtls[qtls$Group==i,]
    iqtl<-iqtl[order(iqtl$Position),]
    if(nrow(iqtl)>0){
      #qtls
      write("qtls",filename,append=TRUE)
      write.table(iqtl[c("traits","Position","Position","Position","Position","i","b","c")],filename,
                  quote=FALSE,col.names=FALSE,na="",
                  row.names=FALSE,append=TRUE)
      #segments
      write("segments",filename,append=TRUE)
      write.table(iqtl[c("Position","Position","c")],filename,
                  quote=FALSE,col.names=FALSE,na="",
                  row.names=FALSE,append=TRUE)

    }
  }
}
###############
#
####################

#SoyBase_GWAS_Positions
#library(stringr)
gwas_mapchart<-function(traits=c("isoflavone","protein","oil"),filename="gwaschart.mct"){
  gwas<-SoyBase_GWAS_Positions
  gwas$traits_loci<-paste(gwas$traits,"_",gwas$loci,sep="")
  gqtls<-gwas
  #排序
  gqtls<-gqtls[order(gqtls$chr,gqtls$Mb),]
  gqtls$i<-"i"
  gqtls$b<-"b"
  gqtls$c<-NA
  gqtls$c[str_detect(str_to_lower(gqtls$traits),traits[1])]<-"c2"
  gqtls$c[str_detect(str_to_lower(gqtls$traits),traits[2])]<-"c3"
  gqtls$c[str_detect(str_to_lower(gqtls$traits),traits[3])]<-"c4"

  splitG<-split(gqtls,gqtls$chr)
  write("\t",filename)
  for(iname in names(splitG)){
    #
    write(paste("group",iname),filename,append=TRUE)
    #
    write.table(splitG[[iname[1]]][c("traits_loci","Mb","b","c")],
                filename,col.names=FALSE,row.names=FALSE,
                na="",
                quote=FALSE,append=TRUE)

  }
}


################
#
######################
gwas_mapchart2<-function(traits=c("protein","oil"),filename="gwaschart.mct"){
  gwas<-SoyBase_GWAS_Positions
  gwas$chr_loci<-paste(gwas$chr,"_",gwas$loci,sep="")

  if(length(traits)>0){
    traits<-str_to_lower(traits)
    gqtls<-NULL
    for(iname in traits){
      igqtls<-gwas[str_detect(str_to_lower(gwas$traits),iname),]
      gqtls<-rbind(gqtls,igqtls)
    }
  }else{
    gqtls<-gwas
  }
  #排序
  gqtls<-gqtls[order(gqtls$chr,gqtls$Mb),]

  gqtls$i<-"i"
  gqtls$b<-"b"
  gqtls$c<-"c1"
  gqtls$c[str_detect(str_to_lower(gqtls$traits),traits[2])]<-"c2"
  gqtls$c[str_detect(str_to_lower(gqtls$traits),traits[3])]<-"c3"
  splitG<-split(gqtls,gqtls$chr)
  write("\t",filename)
  iname<-NULL
  for(iname in names(splitG)){
      mynrow<-nrow(splitG[[iname]])
     if(mynrow>0){
       #write marker
       write(paste("group",iname),filename,append=TRUE)
       #
       write("s 0",filename,append=TRUE)
       write.table(splitG[[iname]][c("chr_loci","Mb","b","c")],filename,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
       endMb<-splitG[[iname]]$Mb[mynrow]+1
       write(paste("e",endMb),filename,append=TRUE)
       #qtls

       write("qtls",filename,append=TRUE)
       write.table(splitG[[iname]][c("traits","Mb","Mb","Mb","Mb","i","b","c")],filename,
                   quote=FALSE,col.names=FALSE,na="",
                   row.names=FALSE,append=TRUE)
       #segments
       write("segments",filename,append=TRUE)
       write.table(splitG[[iname]][c("Mb","Mb","c")],filename,
                   quote=FALSE,col.names=FALSE,na="",
                   row.names=FALSE,append=TRUE)
     }
  }
}




# #add h jymap3943
# conm<-get_consensus_marker_only()
# jym<-jymap3943
#
# jym$h<-NA
# for(i in 1:nrow(jym)){
#   ibp<-jym$bp[i]
#   igroup<-jym$group[i]
#   #notice use Gm column,first replace "Gm",then as.numeric.
#   #the program become slow.
#   iconm<-conm[as.numeric(str_replace(conm$Chromosome,"Gm",""))==igroup,]
#   im<-iconm$Start
#   #the nearest marker
#   nearm<-im[which.min(abs(ibp-im))]
#   #the corresponding h
#   ih<-as.character(iconm$h[iconm$Start==nearm])
#   jym$h[i]<-ih[1]
# }
#
# write.csv(jym,"jym.csv")








