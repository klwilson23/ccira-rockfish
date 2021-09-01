library(dplyr)
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}
rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")

new_df <- read.csv("Data/Rockfish counts by sample id.csv")
hotspots_df <- aggregate(depth~PU_4Km_ID+gear+species,data=new_df,FUN=mean)
hotspots_df$max_depth <- aggregate(max_depth~PU_4Km_ID+gear+species,data=new_df,FUN=max)$max_depth
hotspots_df$UpperOceanSR <- new_df$UpperOceanSR[match(hotspots_df$PU_4Km_ID,new_df$PU_4Km_ID)]

hotspots_df$effort <- aggregate(effort~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$effort
hotspots_df$counts <- aggregate(counts~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$counts
hotspots_df$sample_sizes <- aggregate(sample_size~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$sample_size
#hotspots_df$p_counts <- predict(m3fTMB,newdata = hotspots_df,type="r")
hotspots_df$p_counts <- aggregate(predicted_counts_TMB2~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$predicted_counts_TMB2
hotspots_df$lambda <- hotspots_df$p_counts/hotspots_df$effort
hotspots_df$cpue <- hotspots_df$counts/hotspots_df$effort

hotspots_df <- hotspots_df %>%
  group_by(species,gear) %>% 
  mutate(normalized_lambda = normalize(lambda),normalized_cpue = normalize(cpue)) %>%
  ungroup(species,gear)
hotspots_df <- hotspots_df[complete.cases(hotspots_df),]
hotspots_agg <- aggregate(cbind(normalized_lambda,normalized_cpue)~PU_4Km_ID+species,data=hotspots_df,FUN=mean)
hotspots_agg$depth <- aggregate(depth~PU_4Km_ID+species,data=hotspots_df,FUN=mean)$depth
hotspots_agg$max_depth <- aggregate(max_depth~PU_4Km_ID+species,data=hotspots_df,FUN=max)$max_depth
hotspots_agg$UpperOceanSR <- hotspots_df$UpperOceanSR[match(hotspots_agg$PU_4Km_ID,hotspots_df$PU_4Km_ID)]
hotspots_agg$sample_sizes <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df,FUN=sum)$sample_sizes
dive_samps <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df[hotspots_df$gear=="dive",],FUN=sum)
hotspots_agg$dive_samps <- pmax(0,dive_samps$sample_sizes[match(hotspots_agg$PU_4Km_ID,dive_samps$PU_4Km_ID)],na.rm=TRUE)
hook_samps <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df[hotspots_df$gear=="hook_line",],FUN=sum)
hotspots_agg$hook_samps <- pmax(0,hook_samps$sample_sizes[match(hotspots_agg$PU_4Km_ID,hook_samps$PU_4Km_ID)],na.rm=TRUE)
mid_samps <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df[hotspots_df$gear=="mid_video",],FUN=sum)
hotspots_agg$mid_samps <- pmax(0,mid_samps$sample_sizes[match(hotspots_agg$PU_4Km_ID,mid_samps$PU_4Km_ID)],na.rm=TRUE)
deep_samps <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df[hotspots_df$gear=="deep_video",],FUN=sum)
hotspots_agg$deep_samps <- pmax(0,deep_samps$sample_sizes[match(hotspots_agg$PU_4Km_ID,deep_samps$PU_4Km_ID)],na.rm=TRUE)
plot(normalized_lambda~depth,hotspots_df[hotspots_df$species=="quillback",])
plot(normalized_cpue~depth,hotspots_df[hotspots_df$species=="quillback",])

write.csv(hotspots_agg,"Data/rockfish normalized cpue by site and species.csv")

plot(normalized_lambda~normalized_cpue,data=hotspots_df)
summary(lm(normalized_lambda~normalized_cpue,data=hotspots_df))
plot(normalized_lambda~normalized_cpue,data=hotspots_agg)
summary(lm(normalized_lambda~normalized_cpue,data=hotspots_agg))

hotspots_coast <- aggregate(cbind(normalized_lambda,normalized_cpue)~PU_4Km_ID,data=hotspots_agg,FUN=function(x){sum(x)/length(rockfish_spp)})
hotspots_coast$sample_sizes <- aggregate(sample_sizes~PU_4Km_ID,data=hotspots_agg,FUN=sum)$sample_sizes
hotspots_coast$spp_richness <- aggregate(normalized_cpue~PU_4Km_ID,data=hotspots_agg,FUN=function(x){sum(x>0)})$normalized_cpue

hotspots_coast <- hotspots_coast %>%
  mutate(quantilelambda = ntile(normalized_lambda, 10),quantileCPUE = ntile(normalized_cpue, 10))

plot(quantilelambda~quantileCPUE,data=hotspots_coast)
summary(lm(quantilelambda~quantileCPUE,data=hotspots_coast))
hotspots_coast$hotspot <- ifelse(hotspots_coast$quantilelambda>=9,1,0)
hotspots_coast$hotspot_old <- ifelse(hotspots_coast$quantileCPUE>=9,1,0)
sum((hotspots_coast$hotspot-hotspots_coast$hotspot_old)==1)
sum((hotspots_coast$hotspot-hotspots_coast$hotspot_old)==-1)
sum((hotspots_coast$hotspot-hotspots_coast$hotspot_old)==0)
sum((hotspots_coast$hotspot+hotspots_coast$hotspot_old)==2)
