# Export selected archived arrays only; zero model evaluations.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1]);entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"))
contract <- readRDS(file.path(study,"manifests/contract.rds"))
points <- list();labels <- list()
if(file.exists(file.path(study,"tables/baseline_policy_gates.tsv"))) {
 gates <- read_tsv(file.path(study,"tables/baseline_policy_gates.tsv"))
 for(i in seq_len(nrow(gates))) {
  policy <- gates$policy[i];name <- paste0("baseline_",policy)
  points[[name]] <- readRDS(file.path(study,"tables",paste0(name,".rds")))
  labels[[name]] <- data.frame(point=name,kind="unrestricted_candidate",policy=policy,density_index=0L,
   checks_passed=gates$pass[i],release_status=if(gates$pass[i]) "accepted_baseline" else "unresolved_candidate")
 }
}
if(file.exists(file.path(study,"tables/profile_point_review.tsv"))) {
 review <- read_tsv(file.path(study,"tables/profile_point_review.tsv"))
 for(i in seq_len(nrow(review))) {
  z <- review[i,];name <- paste0("profile_",z$policy,"_",z$density_index)
  points[[name]] <- readRDS(file.path(study,"tables",paste0(name,"_verified.rds")))
  labels[[name]] <- data.frame(point=name,kind="density_anchor",policy=z$policy,density_index=z$density_index,
   checks_passed=z$pass,release_status=if(z$pass) "anchor_checks_passed" else "unresolved_candidate")
 }
}
if(!length(points)) {cat("No closed selected candidate table; no export.\n");quit(status=0)}
streams <- list();chromosomes <- list();states <- list();components <- list();parameters <- list();observations <- list()
for(name in names(points)) {
 a <- points[[name]];r <- a$record;label <- labels[[name]]
 stopifnot(identical(r$contract_sha256,hash_file(file.path(study,"manifests/contract.rds"))),
  identical(r$source_manifest_sha256,hash_file(file.path(study,"manifests/input_manifest.tsv"))))
 components[[name]] <- cbind(label,source_task=a$source_task,as.data.frame(r$components))
 p <- r$parameters;nat <- p;use <- startsWith(names(p),"log10_");nat[use] <- 10^p[use]
 position <- (p-contract$lower)/(contract$upper-contract$lower)
 parameters[[name]] <- cbind(label[rep(1L,length(p)),],data.frame(parameter=names(p),transformed=unname(p),natural=unname(nat),
  normalized_position=unname(position),near_lower_bound=unname(position)<=1e-6,near_upper_bound=unname(position)>=1-1e-6,
  bound_flag_normalized_tolerance=1e-6))
 for(stream in names(r$predictions$tables)) {
  tab <- r$predictions$tables[[stream]]
  observations[[stream]][[name]] <- cbind(label[rep(1L,nrow(tab)),],tab)
 }
 for(i in seq_along(r$predictions$native)) {
  n <- r$predictions$native[[i]]$result;arg <- r$predictions$native[[i]]$arguments
  scenario <- contract$scenarios[[i]];days <- arg$obs_steps*arg$DT
  stopifnot(length(days)==length(n$Vmm3_total_obs),length(n$frac_N_live)==arg$N1max-arg$N1min+1L)
  key <- paste(name,i,sep="_")
  streams[[key]] <- cbind(label[rep(1L,length(days)),],data.frame(harvest=scenario$harvest,day=days,
   live=n$Ntot_live_obs,retained_hypoxia_dead=n$Ntot_dead_hypoxia_obs,retained_buffer_dead=n$Ntot_dead_buffer_obs,
   total_stock=n$Ntot_total_obs,total_volume_mm3=n$Vmm3_total_obs,live_volume_mm3=n$Vmm3_live_obs,
   hypoxia_dead_volume_mm3=n$Vmm3_dead_hypoxia_obs,buffer_dead_volume_mm3=n$Vmm3_dead_buffer_obs,
   oxygen_effective=n$O2_eff_obs,oxygen_target=n$O2_target_obs))
  chr <- seq.int(arg$N1min,arg$N1max)
  chromosomes[[key]] <- cbind(label[rep(1L,length(chr)),],data.frame(harvest=scenario$harvest,
   terminal_day=arg$sim_end_step*arg$DT,chromosome=chr,live_fraction=n$frac_N_live,
   live_count=n$frac_N_live*n$Ntot_live_terminal))
  state <- expand.grid(day=days,chromosome=chr)
  stopifnot(nrow(state)==length(n$live_state_obs),nrow(n$live_state_obs)==length(days),ncol(n$live_state_obs)==length(chr))
  states[[key]] <- cbind(label[rep(1L,nrow(state)),],harvest=scenario$harvest,state,
   live=as.vector(n$live_state_obs),retained_hypoxia_dead=as.vector(n$dead_hypoxia_state_obs),
   retained_buffer_dead=as.vector(n$dead_buffer_state_obs))
 }
}
output <- list(selected_candidate_components=components,selected_candidate_parameters=parameters,
 selected_candidate_prediction_streams=streams,selected_candidate_terminal_chromosomes=chromosomes,
 selected_candidate_chromosome_stocks=states)
for(name in names(output)) write_tsv(do.call(rbind,output[[name]]),file.path(study,"tables",paste0(name,".tsv")))
for(stream in names(observations)) write_tsv(do.call(rbind,observations[[stream]]),file.path(study,"tables",paste0("selected_candidate_",stream,"_predictions.tsv")))
assert_namespaces();cat("Exported",length(points),"archived candidates with explicit acceptance labels; zero model evaluations.\n")
