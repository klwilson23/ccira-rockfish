library(dplyr)
max_only <- function(x, na.rm = TRUE) {
  return((x- 0) /(max(x, na.rm = TRUE)-0))
}

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")
coral_spp <- c("calcigorgia","primnoa","stylaster","chrysopathes","paragorgia","swiftia")
coral_spp <- sort(coral_spp)


coral_scores <- read.csv("Data/coral species weights v2.csv")
coral_scores$score <- coral_scores$height/max(coral_scores$height)
coral_ncc <- read.csv("Data/coral normalized cpue by 4km puid and species.csv",header=TRUE)
coral_ncc2 <- read.csv("Data/coral normalized cpue by 1km puid and species.csv",header=TRUE)

coral_ncc$score <- coral_scores$score[match(coral_ncc$species,coral_scores$name)]
coral_ncc$raw_score <- coral_ncc$normalized_lambda*coral_ncc$score
coral_ncc$obs_score <- coral_ncc$normalized_cpue*coral_ncc$score

coral_ncc2$score <- coral_scores$score[match(coral_ncc2$species,coral_scores$name)]
coral_ncc2$raw_score <- coral_ncc2$normalized_lambda*coral_ncc2$score
coral_ncc2$obs_score <- coral_ncc2$normalized_cpue*coral_ncc2$score

hotspots_coast <- aggregate(raw_score~PU_4Km_ID+PU_1Km_ID,data=coral_ncc2,function(x){sum(x)})
hotspots_coast <- aggregate(raw_score~PU_4Km_ID,data=hotspots_coast,function(x){mean(max(x),ifelse(sum(x[-which.max(x)])>length(which.max(x)),mean(x[-which.max(x)]),mean(x)))})
hotspots_coast$obs_score <- aggregate(obs_score~PU_4Km_ID,data=coral_ncc,function(x){sum(x)})$obs_score
hotspots_coast$UpperOceanSR <- coral_ncc$UpperOceanSR[match(hotspots_coast$PU_4Km_ID,coral_ncc$PU_4Km_ID)]
hotspots_coast$depth <- aggregate(depth~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$depth
hotspots_coast$max_depth <- aggregate(max_depth~PU_4Km_ID,data=coral_ncc,max)$max_depth
hotspots_coast$dive_samps <- aggregate(dive_samps~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$dive_samps
hotspots_coast$hook_samps <- aggregate(hook_samps~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$hook_samps
hotspots_coast$mid_samps <- aggregate(mid_samps~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$mid_samps
hotspots_coast$deep_samps <- aggregate(deep_samps~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$deep_samps
hotspots_coast$normalized_obs_score <- max_only(hotspots_coast$obs_score)

hotspots_coast$normalized_score <- max_only(hotspots_coast$raw_score)
hotspots_coast <- hotspots_coast %>%
  mutate(coral_rank = ntile(normalized_score, 10))
hotspots_coast$coral_hotspot <- ifelse(hotspots_coast$coral_rank>=10,1,0)
hotspots_coast$coral_hotspot2 <- ifelse(hotspots_coast$normalized_score>=0.01,1,0)
write.csv(hotspots_coast,"Data/coral hotspots 4km PUID.csv")
