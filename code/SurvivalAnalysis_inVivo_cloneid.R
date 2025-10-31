# SurvivalAnalysis_inVivo_cloneid.R
# This script focuses on performing survival and time-course analyses on data
# from in vivo experiments. Its primary goal is to compare the survival outcomes
# between different experimental groups (2N vs. 4N cells) by calculating the
# time between "seeding" and "harvest" events from a database. It then uses
# this information to fit and plot Kaplan-Meier survival curves, complete with
# risk tables and confidence intervals. The script also includes secondary analyses,
# such as calculating and visualizing the distribution of treatment durations.

library(dplyr)
library(gplots)

# --- 1. Database Query and Setup ---
# Connect to the database.
mydb = cloneid::connect2DB()
# Query the Passaging table for relevant metadata.
q <- "SELECT id, event,media, passaged_from_id1,cellLine, correctedCount,passage, date,lastModified,owner from Passaging"
pass <- dbGetQuery(mydb,q)
rownames(pass) <- pass$id
# Order data by date, descending.
pass=pass[order(pass$date,decreasing=T),]
# Filter for specific media types used in this experiment.
pass=pass[!is.na(pass$media),]
pass=pass[pass$media==133 | pass$media==134,]
# Exclude events that mark the start of Gemcitabine treatment.
pass_=pass[grep("Gem",pass$id, invert = T),];


## --- 2. Prepare Data for Survival Analysis ---
## The goal is to determine the time between a "seeding" event and a corresponding "harvest" event.

# (This section appears to be an initial attempt, the main analysis follows below)
# `alive` would be seedings that have not yet been harvested.
alive = setdiff(pass_$id[pass_$event=="seeding"], pass_[pass_$event=="harvest","passaged_from_id1"])
# `saced` (sacrificed) would be seedings that HAVE been harvested.
saced=intersect(pass_$id[pass_$event=="seeding"], pass_[pass_$event=="harvest","passaged_from_id1"])


## --- 3. Compute Kaplan Meier curves ---
library(survminer)
library(survival)
library(xlsx)

# To calculate survival time, we need to link each event to its parent event.
# `match` finds the row index of the parent (`passaged_from_id1`) in the `id` column.
ii=match(pass_$passaged_from_id1, pass_$id)
# Create a new dataframe linking child event date (`date`) with parent event date (`date2`).
dt=cbind(pass_[,c('id','passaged_from_id1','date')],pass_[ii,c('id','date')])
colnames(dt)[3]="date2"
dt=dt[!is.na(dt$date),]
# Convert date strings to datetime objects.
dt$date2=as.POSIXct(dt$date2, format = "%Y-%m-%d %H:%M:%S")
dt$date=as.POSIXct(dt$date, format = "%Y-%m-%d %H:%M:%S")
# The survival time is the difference between these two dates.
dt$deltaDays=dt$date2-dt$date
# Write the intermediate data to an Excel file.
write.xlsx(dt, "data/dt_Gem.xlsx")

# Structure the data for the `survival` package.
my_data <- data.frame(
  time = dt$deltaDays,  # The duration of survival.
  # The event status. Here, all events are considered "occurred" (1), meaning
  # no data is censored (e.g., from an animal being lost to follow-up).
  event = rep(1,nrow(dt)),
  # Define the groups to compare: 2N vs 4N, determined by searching the ID string.
  group= as.numeric(grepl('4N', dt$id))
)

# Create a survival object, which is the standard input for survival models in R.
surv_object <- Surv(my_data$time, my_data$event)
# Fit the Kaplan-Meier model, stratified by the 'group' variable.
km_fit <- survfit(surv_object ~ group, data = my_data)
# Plot the Kaplan-Meier curves using the `ggsurvplot` function.
ggsurvplot(
  km_fit,
  conf.int = TRUE,           # Show confidence intervals.
  risk.table = TRUE,         # Include a table of numbers at risk below the plot.
  conf.int.style = "step",
  xlab = "Time",
  ylab = "Survival Probability",
  title = "Kaplan-Meier Survival Curve by Group",
  legend.title = "Group",
  legend.labs = c("2N", "4N") # Provide clearer names for the groups in the legend.
)


## --- 4. scRNAseq Data ---
# This section reads single-cell RNAseq data from an S3 bucket for one mouse as example and creates a simple heatmap.
library(aws.s3)
# list objects under the prefix in the S3 bucket.
objs <- get_bucket_df("cloneid4mysql8", prefix="scRNAseq_InVivo/A04_CLONEID_input/", max=Inf)
# Filter for files ending in .cbs.
cbs_keys <- subset(objs, grepl("\\.cbs$", Key, ignore.case=TRUE))$Key

# Read the first .cbs file directly from S3 into a data frame.
sc <- aws.s3::s3read_using(
  FUN = read.table,
  object = cbs_keys[1], bucket = "cloneid4mysql8",
  sep = "", header = TRUE, fill = TRUE,
  quote = "", comment.char = "", stringsAsFactors = FALSE
)

# Generate a simple heatmap of the scRNAseq data.
heatmap.2(as.matrix(sc),trace="n",Colv = F,margins = c(12,12))

# --- 5. Calculate Treatment Durations ---
# This analysis calculates the time between consecutive "on treatment" harvest events
# for the same animal/lineage, providing a distribution of treatment period lengths.
pass_=pass[pass$event=='harvest',];
treatment_durations <- pass_ %>%
  # Filter for events where the treatment media (134) was used.
  filter(media == 134) %>%
  mutate(date = as.POSIXct(date, format = "%Y-%m-%d %H:%M:%S")) %>%
  # IMPORTANT: Arrange by the parent ID and then date to ensure correct ordering.
  arrange(passaged_from_id1, date) %>%
  # Group by the parent ID, so calculations are done within each animal/lineage.
  group_by(passaged_from_id1) %>%
  # `lag(date)` gets the date from the *previous* row within the same group.
  # `difftime` then calculates the duration between the current and previous event.
  mutate(duration_days = as.numeric(difftime(date, lag(date), units = "days"))) %>%
  ungroup() %>%
  # The first event for each group will have an NA duration, so remove it.
  filter(!is.na(duration_days))

# Print a summary of the calculated durations.
print("--- Summary of Treatment Durations (days) ---")
print(summary(treatment_durations$duration_days))
print("---------------------------------------------")

# --- 6. Visualize Durations as a Histogram ---
# Create a histogram to show the distribution of treatment period lengths.
duration_plot <- ggplot(treatment_durations, aes(x = duration_days)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Distribution of 'On Treatment' Durations",
       x = "Duration of Treatment Period (days)",
       y = "Frequency (Number of Periods)") +
  theme_minimal()

# Display the plot.
print(duration_plot)