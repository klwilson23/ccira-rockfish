library(dplyr)
max_only <- function(x, na.rm = TRUE) {
  return((x- 0) /(max(x, na.rm = TRUE)-0))
}

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

sponge_ncc <- read.csv("Data/sponge normalized cpue by 4km puid and species.csv",header=TRUE)
sponge_ncc2 <- read.csv("Data/sponge normalized cpue by 1km puid and species.csv",header=TRUE)

sponge_ncc$score <- 1
sponge_ncc$raw_score <- sponge_ncc$normalized_lambda*sponge_ncc$score
sponge_ncc$obs_score <- sponge_ncc$normalized_cpue*sponge_ncc$score

sponge_ncc2$score <- 1
sponge_ncc2$raw_score <- sponge_ncc2$normalized_lambda*sponge_ncc2$score
sponge_ncc2$obs_score <- sponge_ncc2$normalized_cpue*sponge_ncc2$score

hotspots_coast <- aggregate(raw_score~PU_4Km_ID+PU_1Km_ID,data=sponge_ncc2,function(x){sum(x)})
hotspots_coast <- aggregate(raw_score~PU_4Km_ID,data=hotspots_coast,function(x){mean(max(x),ifelse(sum(x[-which.max(x)])>length(which.max(x)),mean(x[-which.max(x)]),mean(x)))})
hotspots_coast$obs_score <- aggregate(obs_score~PU_4Km_ID,data=sponge_ncc,function(x){sum(x)})$obs_score
hotspots_coast$UpperOceanSR <- sponge_ncc$UpperOceanSR[match(hotspots_coast$PU_4Km_ID,sponge_ncc$PU_4Km_ID)]
hotspots_coast$depth <- aggregate(depth~PU_4Km_ID,data=sponge_ncc,function(x){mean(x)})$depth
hotspots_coast$max_depth <- aggregate(max_depth~PU_4Km_ID,data=sponge_ncc,max)$max_depth
hotspots_coast$dive_samps <- aggregate(dive_samps~PU_4Km_ID,data=sponge_ncc,function(x){mean(x)})$dive_samps
hotspots_coast$hook_samps <- aggregate(hook_samps~PU_4Km_ID,data=sponge_ncc,function(x){mean(x)})$hook_samps
hotspots_coast$mid_samps <- aggregate(mid_samps~PU_4Km_ID,data=sponge_ncc,function(x){mean(x)})$mid_samps
hotspots_coast$deep_samps <- aggregate(deep_samps~PU_4Km_ID,data=sponge_ncc,function(x){mean(x)})$deep_samps
hotspots_coast$normalized_obs_score <- max_only(hotspots_coast$obs_score)

hotspots_coast$normalized_score <- max_only(hotspots_coast$raw_score)
hotspots_coast <- hotspots_coast %>%
  mutate(sponge_rank = ntile(normalized_score, 10))
hotspots_coast$sponge_hotspot <- ifelse(hotspots_coast$sponge_rank>=9,1,0)
hotspots_coast$sponge_hotspot2 <- ifelse(hotspots_coast$normalized_score>=0.01,1,0)
write.csv(hotspots_coast,"Data/sponge hotspots 4km PUID.csv")

hotspots_coast <- aggregate(cbind(normalized_lambda,normalized_cpue)~PU_4Km_ID,data=sponge_ncc,FUN=function(x){sum(x)})
hotspots_coast$sample_sizes <- aggregate(sample_sizes~PU_4Km_ID,data=sponge_ncc,FUN=sum)$sample_sizes
hotspots_coast$spp_richness <- aggregate(normalized_cpue~PU_4Km_ID,data=sponge_ncc,FUN=function(x){sum(x>0)})$normalized_cpue

hotspots_coast <- hotspots_coast %>%
  mutate(quantilelambda = ntile(normalized_lambda, 10),quantileCPUE = ntile(normalized_cpue, 10))

plot(quantilelambda~quantileCPUE,data=hotspots_coast)
summary(lm(quantilelambda~quantileCPUE,data=hotspots_coast))