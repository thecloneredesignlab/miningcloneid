library(dplyr)
library(gplots)
mydb = cloneid::connect2DB()
q <- "SELECT id, event,media, passaged_from_id1,cellLine, correctedCount,passage, date,lastModified,owner from Passaging"
pass <- dbGetQuery(mydb,q)
rownames(pass) <- pass$id
pass=pass[order(pass$date,decreasing=T),]
pass=pass[!is.na(pass$media),]
pass=pass[pass$media==133 | pass$media==134,]
pass_=pass[grep("Gem",pass$id, invert = T),]; ## exclude treatment start events


## gather harvests only
pass_[pass_$event=="harvest","pass_aged_from_id1"]
alive = setdiff(pass_$id[pass_$event=="seeding"], pass_[pass_$event=="harvest","passaged_from_id1"])
saced=intersect(pass_$id[pass_$event=="seeding"], pass_[pass_$event=="harvest","passaged_from_id1"])

## compute Kaplan Meier curves
library(survminer)
library(survival)
library(xlsx)
ii=match(pass_$passaged_from_id1, pass_$id)
dt=cbind(pass_[,c('id','passaged_from_id1','date')],pass_[ii,c('id','date')])
colnames(dt)[3]="date2"
dt=dt[!is.na(dt$date),]
dt$date2=as.POSIXct(dt$date2, format = "%Y-%m-%d %H:%M:%S")
dt$date=as.POSIXct(dt$date, format = "%Y-%m-%d %H:%M:%S")
dt$deltaDays=dt$date2-dt$date
write.xlsx(dt, "data/dt_Gem.xlsx")

# Fit a Kaplan-Meier survival model
my_data <- data.frame(
  time = dt$deltaDays,  # Time differences
  event = rep(1,nrow(dt)),      # Event occurrence (1 = occurred, 0 = censored)
  group= as.numeric(grepl('4N', dt$id))
)
# Create a survival object
surv_object <- Surv(my_data$time, my_data$event)
# Fit a Kaplan-Meier survival model, stratified by group
km_fit <- survfit(surv_object ~ group, data = my_data)
# Plot the Kaplan-Meier curves
ggsurvplot(
  km_fit,
  conf.int = TRUE,           # Show confidence intervals
  risk.table = TRUE,         # Include a risk table
  conf.int.style = "step",  
  xlab = "Time", 
  ylab = "Survival Probability",
  title = "Kaplan-Meier Survival Curve by Group",
  legend.title = "Group",
  legend.labs = c("2N", "4N") # Rename legend labels
)


## Read scRNAseq data
##f=list.files("../data/S3Buckets/scRNAseq_InVivo/A04_CLONEID_input",pattern = ".cbs",recursive = T, full.names = T)
##sc=read.table(f[1])

# install.packages("aws.s3")  # once
library(aws.s3)

# list objects under the prefix
objs <- get_bucket_df("cloneid4mysql8", prefix="scRNAseq_InVivo/A04_CLONEID_input/", max=Inf)
cbs_keys <- subset(objs, grepl("\\.cbs$", Key, ignore.case=TRUE))$Key

sc <- aws.s3::s3read_using(
  FUN = read.table,
  object = cbs_keys[1], bucket = "cloneid4mysql8",
  sep = "", header = TRUE, fill = TRUE,
  quote = "", comment.char = "", stringsAsFactors = FALSE
)



heatmap.2(as.matrix(sc),trace="n",Colv = F,margins = c(12,12))

# ---  Calculate Treatment Durations ---
# Filter for 'on treatment' events (media == 134).
# Group by cell line, order by date, and calculate the time difference
# between consecutive events for that same line.
pass_=pass[pass$event=='harvest',]; 
treatment_durations <- pass_ %>%
  filter(media == 134) %>%
  mutate(date = as.POSIXct(date, format = "%Y-%m-%d %H:%M:%S")) %>%
  arrange(passaged_from_id1, date) %>%
  group_by(passaged_from_id1) %>%
  mutate(duration_days = as.numeric(difftime(date, lag(date), units = "days"))) %>%
  ungroup() %>%
  filter(!is.na(duration_days)) # Remove first event for each group (NA duration)

# Print a summary of the calculated durations
print("--- Summary of Treatment Durations (days) ---")
print(summary(treatment_durations$duration_days))
print("---------------------------------------------")

# --- Visualize Durations as a Histogram ---
# Define the output path
# output_path <- "/path/to/your/treatment_duration_histogram.png"

# Create the histogram plot
duration_plot <- ggplot(treatment_durations, aes(x = duration_days)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Distribution of 'On Treatment' Durations",
       x = "Duration of Treatment Period (days)",
       y = "Frequency (Number of Periods)") +
  theme_minimal()

print(duration_plot)
# # Save the plot
# ggsave(output_path, plot = duration_plot, width = 8, height = 6, dpi = 300)
# 
# print(paste("Histogram saved to", output_path))

