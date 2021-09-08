library(dplyr)
max_only <- function(x, na.rm = TRUE) {
  return((x- 0) /(max(x, na.rm = TRUE)-0))
}

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

sponge_ncc <- read.csv("Data/sponge normalized cpue by 1km puid and species.csv",header=TRUE)

sponge_ncc$score <- 1
sponge_ncc$raw_score <- sponge_ncc$normalized_lambda*sponge_ncc$score
sponge_ncc$obs_score <- sponge_ncc$normalized_cpue*sponge_ncc$score

hotspots_coast <- aggregate(raw_score~PU_1Km_ID,data=sponge_ncc,function(x){sum(x)})
hotspots_coast$obs_score <- aggregate(obs_score~PU_1Km_ID,data=sponge_ncc,function(x){sum(x)})$obs_score
hotspots_coast$UpperOceanSR <- sponge_ncc$UpperOceanSR[match(hotspots_coast$PU_1Km_ID,sponge_ncc$PU_1Km_ID)]
hotspots_coast$PU_4Km_ID <- sponge_ncc$PU_4Km_ID[match(hotspots_coast$PU_1Km_ID,sponge_ncc$PU_1Km_ID)]
hotspots_coast$depth <- aggregate(depth~PU_1Km_ID,data=sponge_ncc,function(x){mean(x)})$depth
hotspots_coast$max_depth <- aggregate(max_depth~PU_1Km_ID,data=sponge_ncc,max)$max_depth
hotspots_coast$dive_samps <- aggregate(dive_samps~PU_1Km_ID,data=sponge_ncc,function(x){mean(x)})$dive_samps
hotspots_coast$hook_samps <- aggregate(hook_samps~PU_1Km_ID,data=sponge_ncc,function(x){mean(x)})$hook_samps
hotspots_coast$mid_samps <- aggregate(mid_samps~PU_1Km_ID,data=sponge_ncc,function(x){mean(x)})$mid_samps
hotspots_coast$deep_samps <- aggregate(deep_samps~PU_1Km_ID,data=sponge_ncc,function(x){mean(x)})$deep_samps
hotspots_coast$normalized_obs_score <- max_only(hotspots_coast$obs_score)
hotspots_coast$normalized_score <- max_only(hotspots_coast$raw_score)
hotspots_coast <- hotspots_coast %>%
  mutate(sponge_rank = ntile(normalized_score, 10))
hotspots_coast$sponge_hotspot <- ifelse(hotspots_coast$sponge_rank>=10,1,0)
hotspots_coast$sponge_hotspot2 <- ifelse(hotspots_coast$normalized_score>=0.01,1,0)
write.csv(hotspots_coast,"Data/sponge hotspots 1km PUID.csv")

head(sponge_ncc,15)
puid_4km <- aggregate(normalized_score~PU_4Km_ID,data=hotspots_coast,function(x){(max(x))})
puid_4km <- puid_4km %>%
  mutate(sponge_hotspot = ntile(normalized_score, 10))
puid_4km[puid_4km$sponge_hotspot==10,]

puid_4km <- aggregate(normalized_score~PU_4Km_ID,data=hotspots_coast,function(x){mean(max(x),ifelse(length(x)>1,mean(x[-which(x==max(x))]),mean(x)))})
puid_4km <- puid_4km %>%
  mutate(sponge_hotspot = ntile(normalized_score, 10))
puid_4km$samps <- as.vector(table(hotspots_coast$PU_4Km_ID))
plot(puid_4km$samps,puid_4km$sponge_hotspot,xlab="Number of 1km planning units sampled",ylab="Rank at 4km PUID scale")
sum(puid_4km$sponge_hotspot==10 & puid_4km$samps==1)
puid_4km[puid_4km$sponge_hotspot==10,]
hist(table(hotspots_coast$PU_4Km_ID))
sum(hotspots_coast$sponge_hotspot2)
table(hotspots_coast$sponge_rank)
hotspots_coast[hotspots_coast$PU_1Km_ID=="6267",]

oldhotties <- c(6267, 6425, 6582, 6736, 7054, 7524, 7529, 7530, 7847, 7997, 8003, 8157, 8159, 8625, 8628, 8629, 8783, 8786, 8792, 8950, 9102, 9108, 9259, 9265, 9418, 9730, 9731, 9733, 9734, 9889, 9906, 10367, 10524, 10674, 10675, 11331, 11461, 11629, 11774, 11789, 11944, 12093, 12242, 12400, 12410, 12566, 12741, 12759, 12895, 12897, 13034, 13055, 13194, 13212, 13213, 14312)

sum(hotspots_coast[hotspots_coast$PU_1Km_ID%in%oldhotties,"sponge_rank"]>=9)/sum(hotspots_coast$sponge_rank>=9)
sum(hotspots_coast[hotspots_coast$PU_1Km_ID%in%oldhotties,"sponge_rank"]>=7)/length(oldhotties)
plot(sponge_rank~depth,data=hotspots_coast)
mDepth <- glm(sponge_hotspot~depth,data=hotspots_coast,family=binomial)
mDepth2 <- glm(sponge_hotspot~poly(depth,2),data=hotspots_coast,family=binomial)
mDepth_poly <- glm(sponge_hotspot~poly(depth,3),data=hotspots_coast,family=binomial)

AIC(mDepth,mDepth2,mDepth_poly)
depth_seq <- data.frame("depth"=seq(0,400,by=25))
pNew <- predict(mDepth_poly,newdata = depth_seq,type="response")
plot(sponge_hotspot~depth,data=hotspots_coast)
lines(depth_seq$depth,pNew,lwd=2,col="black")
plot(raw_score~obs_score,hotspots_coast)

hotspots_coast <- aggregate(cbind(normalized_lambda,normalized_cpue)~PU_1Km_ID,data=sponge_ncc,FUN=function(x){sum(x)})
hotspots_coast$sample_sizes <- aggregate(sample_sizes~PU_1Km_ID,data=sponge_ncc,FUN=sum)$sample_sizes
hotspots_coast$spp_richness <- aggregate(normalized_cpue~PU_1Km_ID,data=sponge_ncc,FUN=function(x){sum(x>0)})$normalized_cpue

hotspots_coast <- hotspots_coast %>%
  mutate(quantilelambda = ntile(normalized_lambda, 10),quantileCPUE = ntile(normalized_cpue, 10))

plot(quantilelambda~quantileCPUE,data=hotspots_coast)
summary(lm(quantilelambda~quantileCPUE,data=hotspots_coast))