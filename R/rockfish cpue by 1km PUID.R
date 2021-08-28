library(dplyr)
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}
rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")

new_df <- read.csv("Data/Rockfish counts by sample id.csv")
hotspots_df <- aggregate(depth~PU_1Km_ID+gear+species,data=new_df,FUN=mean)
hotspots_df$effort <- aggregate(effort~PU_1Km_ID+gear+species,data=new_df,FUN=sum)$effort
hotspots_df$counts <- aggregate(counts~PU_1Km_ID+gear+species,data=new_df,FUN=sum)$counts
hotspots_df$sample_sizes <- aggregate(sample_size~PU_1Km_ID+gear+species,data=new_df,FUN=sum)$sample_size
hotspots_df$PU_4Km_ID <- new_df$PU_4Km_ID[match(hotspots_df$PU_1Km_ID,new_df$PU_1Km_ID)]
#hotspots_df$p_counts <- predict(m3fTMB,newdata = hotspots_df,type="r")
hotspots_df$p_counts <- aggregate(predicted_counts_TMB2~PU_1Km_ID+gear+species,data=new_df,FUN=sum)$predicted_counts_TMB2
hotspots_df$lambda <- hotspots_df$p_counts/hotspots_df$effort
hotspots_df$cpue <- hotspots_df$counts/hotspots_df$effort

hotspots_df <- hotspots_df %>%
  group_by(species,gear) %>% 
  mutate(normalized_lambda = normalize(lambda),normalized_cpue = normalize(cpue)) %>%
  ungroup(species,gear)
hotspots_df <- hotspots_df[complete.cases(hotspots_df),]
hotspots_agg <- aggregate(cbind(normalized_lambda,normalized_cpue)~PU_1Km_ID+species,data=hotspots_df,FUN=mean)
hotspots_agg$depth <- aggregate(depth~PU_1Km_ID+species,data=hotspots_df,FUN=mean)$depth
hotspots_agg$sample_sizes <- aggregate(sample_sizes~PU_1Km_ID+species,data=hotspots_df,FUN=sum)$sample_sizes
hotspots_agg$sample_sizes <- aggregate(sample_sizes~PU_1Km_ID+species,data=hotspots_df,FUN=sum)$sample_sizes
hotspots_agg$sample_sizes <- aggregate(sample_sizes~PU_1Km_ID+species,data=hotspots_df,FUN=sum)$sample_sizes
dive_samps <- aggregate(sample_sizes~PU_1Km_ID+species,data=hotspots_df[hotspots_df$gear=="dive",],FUN=sum)
hotspots_agg$dive_samps <- pmax(0,dive_samps$sample_sizes[match(hotspots_agg$PU_1Km_ID,dive_samps$PU_1Km_ID)],na.rm=TRUE)
hook_samps <- aggregate(sample_sizes~PU_1Km_ID+species,data=hotspots_df[hotspots_df$gear=="hook_line",],FUN=sum)
hotspots_agg$hook_samps <- pmax(0,hook_samps$sample_sizes[match(hotspots_agg$PU_1Km_ID,hook_samps$PU_1Km_ID)],na.rm=TRUE)
mid_samps <- aggregate(sample_sizes~PU_1Km_ID+species,data=hotspots_df[hotspots_df$gear=="mid_video",],FUN=sum)
hotspots_agg$mid_samps <- pmax(0,mid_samps$sample_sizes[match(hotspots_agg$PU_1Km_ID,mid_samps$PU_1Km_ID)],na.rm=TRUE)
deep_samps <- aggregate(sample_sizes~PU_1Km_ID+species,data=hotspots_df[hotspots_df$gear=="deep_video",],FUN=sum)
hotspots_agg$deep_samps <- pmax(0,deep_samps$sample_sizes[match(hotspots_agg$PU_1Km_ID,deep_samps$PU_1Km_ID)],na.rm=TRUE)
plot(normalized_lambda~depth,hotspots_df[hotspots_df$species=="quillback",])
plot(normalized_cpue~depth,hotspots_df[hotspots_df$species=="quillback",])

write.csv(hotspots_agg,"Data/rockfish normalized cpue by 1km puid and species.csv")
