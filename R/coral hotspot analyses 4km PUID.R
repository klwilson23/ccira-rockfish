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
coral_ncc <- read.csv("Data/coral normalized cpue by site and species.csv",header=TRUE)

coral_ncc$score <- coral_scores$score[match(coral_ncc$species,coral_scores$name)]
coral_ncc$raw_score <- coral_ncc$normalized_lambda*coral_ncc$score

hotspots_coast <- aggregate(raw_score~PU_4Km_ID,data=coral_ncc,function(x){sum(x)})
hotspots_coast$UpperOceanSR <- coral_ncc$UpperOceanSR[match(hotspots_coast$PU_4Km_ID,coral_ncc$PU_4Km_ID)]
hotspots_coast$depth <- aggregate(depth~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$depth
hotspots_coast$max_depth <- aggregate(max_depth~PU_4Km_ID,data=coral_ncc,max)$max_depth
hotspots_coast$dive_samps <- aggregate(dive_samps~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$dive_samps
hotspots_coast$hook_samps <- aggregate(hook_samps~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$hook_samps
hotspots_coast$mid_samps <- aggregate(mid_samps~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$mid_samps
hotspots_coast$deep_samps <- aggregate(deep_samps~PU_4Km_ID,data=coral_ncc,function(x){mean(x)})$deep_samps

hotspots_coast$normalized_score <- max_only(hotspots_coast$raw_score)
hotspots_coast <- hotspots_coast %>%
  mutate(coral_rank = ntile(normalized_score, 10))
hotspots_coast$coral_hotspot <- ifelse(hotspots_coast$coral_rank>=9,1,0)

table(hotspots_coast$coral_rank)
hotspots_coast[hotspots_coast$PU_4Km_ID=="6267",]

oldhotties <- c(6267, 6425, 6582, 6736, 7054, 7524, 7529, 7530, 7847, 7997, 8003, 8157, 8159, 8625, 8628, 8629, 8783, 8786, 8792, 8950, 9102, 9108, 9259, 9265, 9418, 9730, 9731, 9733, 9734, 9889, 9906, 10367, 10524, 10674, 10675, 11331, 11461, 11629, 11774, 11789, 11944, 12093, 12242, 12400, 12410, 12566, 12741, 12759, 12895, 12897, 13034, 13055, 13194, 13212, 13213, 14312)

sum(hotspots_coast[hotspots_coast$PU_4Km_ID%in%oldhotties,"coral_rank"]>=9)/sum(hotspots_coast$coral_rank>=9)
sum(hotspots_coast[hotspots_coast$PU_4Km_ID%in%oldhotties,"coral_rank"]>=7)/length(oldhotties)
plot(coral_rank~depth,data=hotspots_coast)
plot(coral_hotspot~depth,data=hotspots_coast)
mDepth <- glm(coral_hotspot~depth,data=hotspots_coast,family=binomial)
mDepth2 <- glm(coral_hotspot~poly(depth,2),data=hotspots_coast,family=binomial)
mDepth_poly <- glm(coral_hotspot~poly(depth,3),data=hotspots_coast,family=binomial)

AIC(mDepth,mDepth2,mDepth_poly)
depth_seq <- data.frame("depth"=seq(0,400,by=25))
pNew <- predict(mDepth,newdata = depth_seq,type="response")
lines(depth_seq$depth,pNew,lwd=2,col="black")
write.csv(hotspots_coast,"Data/coral hotspots 4km PUID.csv")
