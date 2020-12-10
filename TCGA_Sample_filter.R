# setwd('./data1/2020Year/drug_repositioning/')
argv1<-c('./Cancer_Exp/', 'BRCA', '01,06', './BRCA_all.txt')
# argv[1]= help, All, ACC,BRCA ..  argv1[3]= All, 01, 02
library(stringr)
argv1<-commandArgs(trailingOnly = T)
dir_path=argv1[1] # TCGA dir 

#dir_path='./data1/2020Year/drug_repositioning/Cancer_Exp/'



if(argv1[2]=='help'|argv1[2]=='h'){
  print('Rscript TCGA_Sample_filter.R TCGA_Sample_folder_path Cancer_Type Sample_Type Output')
  print('For Example:')
  print('TCGA_Sample_folder_path (Cancer_Exp)-- ACC')
  print('                                     | BLCA')
  print('                                     | BRCA')
  print('                                     | CESC')
  print('                                     | CHOL')
  print('                                     L COAD')
  print('I want BRCA Cancer type and Sample type 01,06 output file name is BRCA_all.txt.gz')
  print('Rscript TCGA_Sample_filter.R ./Cancer_Exp/ BRCA 01,06 ./BRCA_all.txt # BRCA Sample, Sample type : 01,06 , output : BRCA_all.txt.gz')
  print('Another')
  print('I want All TCGA Cancer type and Sample type All(Tumor,Normal) output file name is all.txt.gz')
  print('Rscript TCGA_Sample_filter.R ./Cancer_Exp/ All All ./all.txt # All TCGA, Sample type : all , output : all.txt.gz')
  
}else{

######### data load
if(argv1[2]=='All'|argv1[2]=='ALL'|argv1[2]=='all'){
  all_list<-list()
  list_d<-list.files(dir_path)
  for(i in 1:length(list_d)){
    data_mat<-gzfile(paste0(dir_path,list_d[i],'/','HiSeqV2.gz'))
    data_mat<-read.table(data_mat,sep='\t',header = T,stringsAsFactors = F,check.names = F)
    rown<-data_mat[,1]
    all_list[[i]]=data_mat[,2:ncol(data_mat)]
  }
  all_data<-cbind.data.frame(sample=rown,all_list)
}

if(argv1[2]!='help'&argv1[2]!='All'&argv1[2]!='ALL'){
  if(argv1[2]%in%list.files(dir_path)){
    data_mat<-gzfile(paste0(dir_path,argv1[2],'/','HiSeqV2.gz'))
    data_mat<-read.table(data_mat,sep='\t',header = T,stringsAsFactors = F,check.names = F)
  }else{
    print('error')
  }
  all_data=data_mat
}
######### data load
## output all_data

if(argv1[3]=='All'|argv1[3]=='ALL'|argv1[3]=='all'){
  col_n<-colnames(all_data)
  col_n<-str_split_fixed(col_n,'-',4)[,4]
  col_n<-paste0('_',col_n)
  col_n[1]='Sample_type'
  all_data<-rbind.data.frame(col_n,all_data)
}else{
  re_list<-list()
  ty=strsplit(argv1[3],',')[[1]]
  for(i in ty){
    TN<-TN[[i]]
    col_n<-colnames(all_data)
    col_n<-str_split_fixed(col_n,'-',4)[,4]
    col_n<-paste0('_',col_n)
    col_n[1]='Sample_type'
    all_data<-rbind.data.frame(col_n,all_data)
    row_n<-all_data[,1]
    all_data1<-all_data[,grep(i,x = all_data[1,])]
  re_list[[i]]=all_data1[,2:ncol(all_data1)]
  }
  all_data<-cbind.data.frame(sample=row_n,re_list)
}

file1<-trimws(argv1[4])
gz_f<-gzfile(paste0(file1,'.gz'),'w')
write.table(all_data,gz_f,sep='\t',quote = F,row.names = F)
}













