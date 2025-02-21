mapqtl_format_mb<-function(marker=jymarker3943,map=jymap3943){
  head(map)
  map<-merge(map,marker)
  map$mb<-map$bp/1000000
  map<-map[order(map$group,map$mb),]
  smap<-split(map,map$group)

  mapfile<-"jy185_chr_mb.map"

  locfile<-"jy185_chr_mb.loc"

  #########write mapfile
  group<-paste("group",1:20)
  write(";20220401",mapfile)
  write(";ngrp=20",mapfile,append=TRUE)
  write(paste(";nloc=",nrow(map)),mapfile,append=TRUE)
  for(i in 1:20){
    write("\n",mapfile,append=TRUE)
    write(group[i],mapfile,append=TRUE)
    write("\n",mapfile,append=TRUE)
    write.table(smap[[i]][c("marker","mb")],mapfile,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

  }
  #########write loc file
  jyname<-c("marker",paste("JY","00",1:9,sep=""),paste("JY","0",10:99,sep=""),paste("JY",100:185,sep=""))
  ######################################
  write(";",locfile)
  write("name=jy_chr_185_mb",locfile,append=TRUE)
  write("popt=RI8",locfile,append=TRUE)
  write(paste("nloc=",nrow(map),sep=""),locfile,append=TRUE)
  write("nind=185",locfile,append=TRUE)
  write("\n",locfile,append=TRUE)
  write.table(map[jyname],locfile,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  #################
}
