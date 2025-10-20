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
  
  dat['SUM-149','ploidy_flow']=2.01
  dat['MCFdcis','ploidy_flow']=2
  dat$ploidy = dat$ploidy_karyo
  dat$ploidy[is.na(dat$ploidy)]=dat$ploidy_flow[is.na(dat$ploidy)]
  split_df$ploidy=dat[split_df$right,'ploidy']+dat[split_df$left,'ploidy']
  split_df[nrow(split_df)+1,c('original','ploidy')]=c('SUM-159 (4N)',3.49593528659252)
  split_df[nrow(split_df)+1,c('original','ploidy')]=c('SUM-159 (2N)',dat["SUM-159",'ploidy_karyo'])
  split_df[nrow(split_df)+1,c('original','ploidy')]=c('MDAMB231 (4N)',6.5)
  split_df[nrow(split_df)+1,c('original','ploidy')]=c('MDAMB231 (2N)',dat["MDAMB231",'ploidy_karyo'])
  split_df[nrow(split_df)+1,c('original','ploidy')]=c('MCF10A (2N)',2.00123)
  split_df[nrow(split_df)+1,c('original','ploidy')]=c('MCF10A (4N)',4.00123)
  write.table(split_df,"../ploidyPerCL_HeteroAndHomotypic.txt",sep='\t',quote=F)
  ploidy=split_df
}

############################################
##########fit drug response models##########
rownames(ploidy)=ploidy$original

# Rmax=1841919
# analyze_drug_response("Ceralasertib_ATRi.txt", output_xlsx_path = "~/Downloads/My_Drug_Analysis.xlsx",Rmax_optional = Rmax)

f=list.files(pattern = ".txt")
f=c(grep("emcitab",f, value=T), grep("emcitab",f, value=T, invert = T))
fit_results=list()
for(x in f){
  tab=read.table(x,header = T,check.names = F,sep="\t")
  tab <- standardize_colnames(tab)
  ii=intersect(colnames(tab), rownames(ploidy))
  df_subset <- ploidy[ii, 'ploidy', drop = FALSE]
  ploidy_map <- setNames(as.numeric(df_subset$ploidy), rownames(df_subset))
  
  source("../../code/ModelingDrug-inducedSelection.R")
  # try(analyze_drug_response(x, ploidy_map =ploidy_map , output_xlsx_path = "~/Downloads/My_Drug_Analysis.xlsx"));#,Rmax_optional = Rmax))
  fit_results[[fileparts(x)$name]]=try(analyze_drug_response_global(x, ploidy_map =ploidy_map , output_xlsx_path = "~/Downloads/My_Drug_Analysis.xlsx"));#,Rmax_optional = Rmax))
}

# Now, call the comparison function with these results
lpModels=c("Capecitabine", "Ceralasertib", "Adavociclib", "Hybrids_Gemcitabine", "Cytarabine")
hpModels=c("ABT199", "Ispinesib", "Volasertib", "Tegafur", "TAS2")

hpModels = unlist(sapply(hpModels, function(x) grep(x, names(fit_results), value=T)))
hpModels = sapply(hpModels, function(x) coef(fit_results[[x]]), simplify = F)
lpModels = unlist(sapply(lpModels, function(x) grep(x, names(fit_results), value=T)))
lpModels = sapply(lpModels, function(x) coef(fit_results[[x]]), simplify = F)
compare_model_fits(hpModels,  output_filename = "~/Downloads/HPregimes_Comparison.png",legloc="bottomleft", ploidy_range=seq(2.5,8,by=0.31))
compare_model_fits(lpModels,  output_filename = "~/Downloads/LPregimes_Comparison.png", ploidy_range=seq(1.5,2.5,by=0.1))

## compare EC50 to ploidy across cell lines
path="~/Downloads/My_Drug_Analysis.xlsx"
sheets <- readxl::excel_sheets(path)
dfs <- lapply(sheets, function(sh) readxl::read_excel(path, sheet = sh))
names(dfs)=sheets
tab=rbind(dfs$SUMPlusMDA_Gemcitabine3, dfs$HybridsPlusMDA_Gemcitabine)
# Define nice colors for the plot
# plot_data=dfs$SUMPlusMDA_Gemcitabine3
plot_data=dfs$MCFPlusSUMPlusMDA_Gemcitabine3
plot_data$cell_line=sapply(strsplit(plot_data$lineage," "),"[[",1)
lineage_colors <- c("SUM-159" = "#1E88E5", "MDAMB231" = "#D81B60","MCF10A"="purple")
# Create the ggplot object
ploidy_ec50_plot <- ggplot(plot_data, aes(x = Ploidy, y = EC50, color = cell_line)) +
  # Add points with some transparency
  geom_point(size = 5, alpha = 0.8) +
  # Add a linear trend line (on the log-transformed data) for each cell line
  geom_smooth(method = "lm", se = FALSE, aes(group = cell_line), linetype = "dashed", size = 1) +
  # Use a logarithmic scale for the y-axis because the EC50 values are far apart
  scale_y_log10(
    breaks = c(10, 20, 50, 100, 200, 500), # Custom breaks for clarity
    labels = scales::label_number()
  ) +
  # Manually set the colors for consistency
  scale_color_manual(values = lineage_colors) +
  # Add informative labels and a title
  labs(    title = "Drug Resistance (EC50) Increases with Ploidy",    subtitle = "Comparison of SUM-159 and MDAMB231 cell lines",
    x = "Ploidy (N)", y = "EC50 (Log Scale)", color = "Cell Line"  ) +
  # Apply a clean theme
  theme_minimal(base_size = 15) +
  theme(    legend.position = "bottom",    plot.title = element_text(face = "bold") )

# Save the plot to a file
ggsave("~/Downloads/Ploidy_vs_EC50_Comparison.png", plot = ploidy_ec50_plot, width = 7, height = 5, dpi = 300)
