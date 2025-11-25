###########################################demographic information for original data##########################################
summary_information=data.frame(matrix(ncol = 3, nrow =30))
colnames(summary_information)=c('Whole population','No MASLD','MASLD')
row.names(summary_information)=c('Total N','Age group','Age<65 years','>=65 yeas',
                                 'Sex','Females','Males','Race','Ethnic minorities','White',
                                 'Deprivation','Quarter 1','Quarter 2','Quarter 3','Quarter 4',
                                 'Education level','Above Colleage','Less than college',
                                 'Activity','Regularly active','Inactive','Smoking','Never','Ever','Alcohol','Never_','Ever_',
                                 'Metabolic syndrome','No','Yes')
for (masld_status in c(10,0,1)) {
  if (masld_status==10){
    index_masld=c(1:nrow(datas))#10=whole samples
    column_index=1
  }else{
    index_masld=which(datas$masld==masld_status)
    column_index=masld_status+2
  }
  temp_datas=datas[index_masld,]
  rm(index_masld)
  #N
  summary_information[1,column_index]=nrow(temp_datas)
  #age
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$age==i-1))
    summary_information[i+2,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #gender
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$gender==i-1))
    summary_information[i+5,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #race
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$race==i-1))
    summary_information[i+8,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #deprivation
  for (i in c(1:4)) {
    n_temp=length(which(temp_datas$deprivation==i-1))
    summary_information[i+11,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #education
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$university_education==i-1))
    summary_information[i+16,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #activity
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$activity==i-1))
    summary_information[i+19,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #smoking
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$smoking==i-1))
    summary_information[i+22,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #alcohol
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$alcohol==i-1))
    summary_information[i+25,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #metabolic_syndrome
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$metabolic_syndrome==i-1))
    summary_information[i+28,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp,i)
  rm(temp_datas)
}
rm(masld_status,column_index)
#write.csv(summary_information,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/demongraphic_masld_depression.csv')



