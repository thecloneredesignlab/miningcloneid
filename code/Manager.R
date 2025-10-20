# To specify a different output file name:
setwd("/Users/4470246/Projects/Teaching/IMOworkshops/IMO-workshop-2025/data/DrugResponseData/")
source("../../code/ModelingDrug-inducedSelection.R")

################
###get ploidy###
if(file.exists("../ploidyPerCL_HeteroAndHomotypic.txt")){
  ploidy=read.table("../ploidyPerCL_HeteroAndHomotypic.txt",sep="\t")
}else{
  # Retrieve all rows where cellLine contains a slash
  mydb = cloneid::connect2DB()
  g <- dbGetQuery(mydb, "SELECT * FROM Passaging WHERE cellLine LIKE '%/%'")
  
  # Split each cellLine value by '/'
  split_lines <- strsplit(g$cellLine, "/")
  
  # Optionally, convert to a data frame with separate columns
  split_df <- data.frame(
    original = g$cellLine,
    left = sapply(split_lines, `[`, 1),
    right = sapply(split_lines, `[`, 2),
    stringsAsFactors = FALSE
  )
  
  # View result
  head(split_df)
  ucls=union(split_df$left,split_df$right)
  
  dat=read.xlsx("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/OtherCellLines/FlowPloidyQuantification/PloidySummary_merged.xlsx")
  # Apply replacements with gsub to standardize
  dat$X1 <- gsub("MDA[- ]?MB[- ]?231", "MDAMB231", dat$X1, ignore.case = TRUE)
  dat$X1 <- gsub("SUM[- ]?159", "SUM-159", dat$X1, ignore.case = TRUE)
  dat$X1 <- gsub("HS?578[ -]?T", "HS578T", dat$X1, ignore.case = TRUE)
  dat$X1 <- gsub("MCF10[ -]?DCIS", "MCFdcis", dat$X1, ignore.case = TRUE)
  dat$X1 <- gsub("SUM[- ]?149", "SUM-149", dat$X1, ignore.case = TRUE)
  rownames(dat)=dat$X1
  dat$ploidy_flow=as.numeric(dat$ploidy_flow)
  adj=dat["MDAMB231",'ploidy_flow']/dat["MDAMB231",'ploidy_WES']
  dat$ploidy_flow=dat$ploidy_flow/adj
  
  dat$ploidy_karyo=NA
  dat["SUM-159",'ploidy_karyo']=2.009975 
  dat["MDAMB231",'ploidy_karyo']=4.15675859608933 
  
  dat[c('MCFdcis','SUM-149'),'ploidy_flow']=2
  dat$ploidy = dat$ploidy_karyo
  dat$ploidy[is.na(dat$ploidy)]=dat$ploidy_flow[is.na(dat$ploidy)]
  split_df$ploidy=dat[split_df$right,'ploidy']+dat[split_df$left,'ploidy']
  split_df[nrow(split_df)+1,c('original','ploidy')]=c('SUM-159 (4N)',3.49593528659252)
  split_df[nrow(split_df)+1,c('original','ploidy')]=c('SUM-159 (2N)',dat["SUM-159",'ploidy_karyo'])
  write.table(split_df,"../ploidyPerCL_HeteroAndHomotypic.txt",sep='\t',quote=F)
  ploidy=split_df
}


############################################
##########fit drug response models##########
rownames(ploidy)=ploidy$original

# Rmax=1841919
# analyze_drug_response("Ceralasertib_ATRi.txt", output_xlsx_path = "~/Downloads/My_Drug_Analysis.xlsx",Rmax_optional = Rmax)

f=list.files(pattern = ".txt")
for(x in f){
  tab=read.table(x,header = T,check.names = F,sep="\t")
  tab <- standardize_colnames(tab)
  ii=intersect(colnames(tab), rownames(ploidy))
  df_subset <- ploidy[ii, 'ploidy', drop = FALSE]
  ploidy_map <- setNames(as.numeric(df_subset$ploidy), rownames(df_subset))
  
  try(analyze_drug_response(x, ploidy_map =ploidy_map , output_xlsx_path = "~/Downloads/My_Drug_Analysis.xlsx"));#,Rmax_optional = Rmax))
}

