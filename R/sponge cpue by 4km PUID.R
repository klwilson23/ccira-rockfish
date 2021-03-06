library(dplyr)
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}
rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")
coral_spp <- c("calcigorgia","primnoa","stylaster","chrysopathes","paragorgia","swiftia")
coral_spp <- sort(coral_spp)

new_df <- read.csv("Data/sponge counts by sample id.csv")
hotspots_df <- aggregate(depth~PU_4Km_ID+gear+species,data=new_df,FUN=mean)
hotspots_df$max_depth <- aggregate(max_depth~PU_4Km_ID+gear+species,data=new_df,FUN=max)$max_depth
hotspots_df$UpperOceanSR <- new_df$UpperOceanSR[match(hotspots_df$PU_4Km_ID,new_df$PU_4Km_ID)]

hotspots_df$effort <- aggregate(effort~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$effort
hotspots_df$counts <- aggregate(counts~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$counts
hotspots_df$sample_sizes <- aggregate(sample_size~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$sample_size
hotspots_df$p_counts <- ifelse(hotspots_df$gear=="deep_video",
                               aggregate(predicted_counts~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$predicted_counts,
                               aggregate(predicted_counts~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$predicted_counts)
hotspots_df$lambda <- ifelse(hotspots_df$gear=="deep_video",hotspots_df$p_counts/hotspots_df$effort,hotspots_df$p_counts/hotspots_df$effort)
hotspots_df$cpue <- ifelse(hotspots_df$gear=="deep_video",hotspots_df$counts/hotspots_df$effort,hotspots_df$counts/hotspots_df$effort)
hotspots_df <- hotspots_df %>%
  group_by(species,gear) %>% 
  mutate(normalized_lambda = normalize(lambda),normalized_cpue = normalize(cpue)) %>%
  ungroup(species,gear)
hotspots_df <- hotspots_df[complete.cases(hotspots_df),]
hotspots_agg <- aggregate(cbind(normalized_lambda,normalized_cpue)~PU_4Km_ID+species,data=hotspots_df,FUN=mean)
hotspots_agg$PU_4Km_ID <- hotspots_df$PU_4Km_ID[match(hotspots_agg$PU_4Km_ID,hotspots_df$PU_4Km_ID)]
hotspots_agg$depth <- aggregate(depth~PU_4Km_ID+species,data=hotspots_df,FUN=mean)$depth
hotspots_agg$max_depth <- aggregate(max_depth~PU_4Km_ID+species,data=hotspots_df,FUN=max)$max_depth
hotspots_agg$UpperOceanSR <- hotspots_df$UpperOceanSR[match(hotspots_agg$PU_4Km_ID,hotspots_df$PU_4Km_ID)]
hotspots_agg$sample_sizes <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df,FUN=sum)$sample_sizes
hotspots_agg$sample_sizes <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df,FUN=sum)$sample_sizes
hotspots_agg$sample_sizes <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df,FUN=sum)$sample_sizes
dive_samps <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df[hotspots_df$gear=="dive",],FUN=sum)
hotspots_agg$dive_samps <- pmax(0,dive_samps$sample_sizes[match(hotspots_agg$PU_4Km_ID,dive_samps$PU_4Km_ID)],na.rm=TRUE)
hotspots_agg$hook_samps <- 0
mid_samps <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df[hotspots_df$gear=="mid_video",],FUN=sum)
hotspots_agg$mid_samps <- pmax(0,mid_samps$sample_sizes[match(hotspots_agg$PU_4Km_ID,mid_samps$PU_4Km_ID)],na.rm=TRUE)
deep_samps <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df[hotspots_df$gear=="deep_video",],FUN=sum)
hotspots_agg$deep_samps <- pmax(0,deep_samps$sample_sizes[match(hotspots_agg$PU_4Km_ID,deep_samps$PU_4Km_ID)],na.rm=TRUE)
plot(normalized_lambda~depth,hotspots_df)
plot(normalized_cpue~depth,hotspots_df)

write.csv(hotspots_agg,"Data/sponge normalized cpue by 4km puid and species.csv")
plot(normalized_lambda~normalized_cpue,data=hotspots_df)
summary(lm(normalized_lambda~normalized_cpue,data=hotspots_df))
plot(normalized_lambda~normalized_cpue,data=hotspots_agg)
summary(lm(normalized_lambda~normalized_cpue,data=hotspots_agg))