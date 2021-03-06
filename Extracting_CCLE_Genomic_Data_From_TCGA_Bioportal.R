#This script describes the process I used to obtain gene copy number data from
#the tcga cbioportal using the cgdsr R package

#I generated a list of HUGO human gene symbols by going here: http://www.genenames.org/cgi-bin/download
#and downloading the information using the default settings.  I then placed the file in my home directory
#and performed the following commands to get all of the gene symbols, remove withdrawn symbols and get 
#rid of duplicates

hugo<-read.csv(file.choose(), header=T)
not_Withdrawn<-hugo$Approved.Name!="entry withdrawn"
hugo1<-hugo[not_Withdrawn,]
withdrawn<-grep("^symbol withdrawn", hugo$Approved.Name)
hugo_final<-hugo[-(withdrawn),]
all_gene_symbols<-hugo_final[[2]]
all_gene_symbols<-as.character(all_gene_symbols)
all_gene_symbols<-unique(all_gene_symbols)

#I next open up the cgdsr library and set up the connection with the server that will be used in subsequent 
#cgdsr commands

library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

#In the next section of code I first initialize a dataframe and then use a for loop to populate
#that dataframe with gene copy number data for each gene symbol in all_gene_symbols.  I am using 
#the getProfileData function from the cgdsr library to pull the data.

Amp_Data<-data.frame(matrix(nrow=1019, ncol=length(all_gene_symbols)))

for(i in 1:length(all_gene_symbols)){
  Amp_Data[i] <- getProfileData(mycgds,all_gene_symbols[i],"cellline_ccle_broad_CNA", "cellline_ccle_broad_all")
}

#I add the appropriate gene symbols as column names

names(Amp_Data)<-all_gene_symbols[1:length(all_gene_symbols)]

#I add the appropriate cell line names as row names
x<-getProfileData(mycgds,all_gene_symbols[1],"cellline_ccle_broad_CNA", "cellline_ccle_broad_all")
rownames(Amp_Data)<-rownames(x)

#We'll next get rid of the tissue type from the end of the cell line name and also change the two cell lines
#that are both named "TT" into "TT1" and "TT2".  

new_names<-sub("_.*", "", row.names(Amp_Data))
new_names[906]<-"TT1"
new_names[907]<-"TT2"
row.names(Amp_Data)<-new_names

#Getting rid of columns that contain only NAs

Amp_Data_NA_Cols_Removed<-Filter(function(x)!all(is.na(x)), Amp_Data)

#Getting rid of rows that contain only NAs

Amp_Data_NA_Cols_Rows_Removed<-Amp_Data_NA_Cols_Removed[rowSums(is.na(Amp_Data_NA_Cols_Removed)) != ncol(Amp_Data_NA_Cols_Removed),]
  
#Let's write the file to the Git repo 

write.csv(Amp_Data_NA_Cols_Rows_Removed, file="CCLE_All_Copy_Number_Data_TCGA.csv")

#In the next section of code I first initialize a dataframe and then use a for loop to populate that dataframe
#with gene expression data for each gene symbol in all_gene_symbols.  I am using the getProfileData function
#from the cgdsr library to pull the data.

Exp_Data<-data.frame(matrix(nrow=1019, ncol=length(all_gene_symbols)))

for(i in 1:length(all_gene_symbols)) {
  Exp_Data[i] <- getProfileData(mycgds,all_gene_symbols[i],"cellline_ccle_broad_mrna_median_Zscores", "cellline_ccle_broad_all")
  }

#I add the appropriate gene symbols as column names

names(Exp_Data)<-all_gene_symbols[1:length(all_gene_symbols)]

#I add the appropriate cell line names as row names
x<-getProfileData(mycgds,all_gene_symbols[1],"cellline_ccle_broad_mrna_median_Zscores", "cellline_ccle_broad_all")
rownames(Exp_Data)<-rownames(x)

#We'll next get rid of the tissue type from the end of the cell line name and also change the two cell lines
#that are both named "TT" into "TT1" and "TT2".  

new_names<-sub("_.*", "", row.names(Exp_Data))
new_names[906]<-"TT1"
new_names[907]<-"TT2"
row.names(Exp_Data)<-new_names

#Getting rid of columns that contain only NAs

Exp_Data_NA_Cols_Removed<-Filter(function(x)!all(is.na(x)), Exp_Data)

#Getting rid of rows that contain only NAs

Exp_Data_NA_Cols_Rows_Removed<-Exp_Data_NA_Cols_Removed[rowSums(is.na(Exp_Data_NA_Cols_Removed)) != ncol(Exp_Data_NA_Cols_Removed),]

#Let's write the file to the Git repo 

write.csv(Exp_Data_NA_Cols_Rows_Removed, file="CCLE_All__Exp_Data_TCGA.csv")

#Now on to the mutation data.  We'll again use the cgdsr package and a for loop to extract mutation data for each gene symbol
#and store it in a dataframe

Mut_Data<-data.frame(matrix(nrow=1019, ncol=length(all_gene_symbols)))

for(i in 1:length(all_gene_symbols)) {
  Mut_Data[i] <- getProfileData(mycgds,all_gene_symbols[i],"cellline_ccle_broad_mutations", "cellline_ccle_broad_all")
}

#I add the appropriate gene symbols as column names

names(Mut_Data)<-all_gene_symbols[1:length(all_gene_symbols)]

#I add the appropriate cell line names as row names
x<-getProfileData(mycgds,all_gene_symbols[1],"cellline_ccle_broad_mutations", "cellline_ccle_broad_all")
rownames(Mut_Data)<-rownames(x)

#We'll next get rid of the tissue type from the end of the cell line name and also change the two cell lines
#that are both named "TT" into "TT1" and "TT2".  

new_names<-sub("_.*", "", row.names(Mut_Data))
new_names[906]<-"TT1"
new_names[907]<-"TT2"
row.names(Mut_Data)<-new_names

#Getting rid of columns that contain only NAs

Mut_Data_NA_Cols_Removed<-Filter(function(x)!all(is.na(x)), Mut_Data)

#Getting rid of rows that contain only NAs

Mut_Data_NA_Cols_Rows_Removed<-Mut_Data_NA_Cols_Removed[rowSums(is.na(Mut_Data_NA_Cols_Removed)) != ncol(Mut_Data_NA_Cols_Removed),]

#Let's write the file to the Git repo 

write.csv(Mut_Data_NA_Cols_Rows_Removed, file="CCLE_All__Mut_Data_TCGA.csv")
