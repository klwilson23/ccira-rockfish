normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

coral_spp <- c("calcigorgia","primnoa","stylaster","chrysopathes","paragorgia","swiftia")
coral_spp <- sort(coral_spp)
sponge_spp <- c("aphrocallistidae","mycale","boot","farrea")

rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")
library(dplyr)
library(glmmTMB)
library(lme4)
library(brms)

new_df <- read.csv("Data/sponge counts PU4km.csv")
deep_df <- new_df[new_df$gear=="deep_video",]
pac_df <- new_df[new_df$gear!="deep_video",]
dive_df <- new_df[new_df$gear=="dive",]
mid_df <- new_df[new_df$gear=="mid_video",]

mSpTMBrand <- glmmTMB(counts~poly(depth,2)+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=mid_df,family=beta_family())
mSpTMBrand2 <- glmmTMB(counts~poly(depth,2)+UpperOceanSR+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=mid_df,family=beta_family())
mSpTMBrand3 <- glmmTMB(counts~depth+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=mid_df,family=beta_family())
mSpTMBrand4 <- glmmTMB(counts~depth+UpperOceanSR+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=mid_df,family=beta_family())
mSpTMBrand5 <- glmmTMB(counts~depth*UpperOceanSR+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=mid_df,family=beta_family())
mSpTMBrand6 <- glmmTMB(counts~depth+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=mid_df,family=beta_family())

AIC(mSpTMBrand,mSpTMBrand2,mSpTMBrand3,mSpTMBrand4,mSpTMBrand5,mSpTMBrand6)
mid_df$predicted_counts<- predict(mSpTMBrand6,type='r')
plot(mid_df$counts,mid_df$predicted_counts)

mDvTMBrand <- glmmTMB(counts~poly(depth,2)+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=dive_df,family=beta_family())
mDvTMBrand2 <- glmmTMB(counts~poly(depth,2)+UpperOceanSR+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=dive_df,family=beta_family())
mDvTMBrand3 <- glmmTMB(counts~poly(depth,2)*UpperOceanSR+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=dive_df,family=beta_family())
mDvTMBrand4 <- glmmTMB(counts~poly(depth,1)+(1|PU_1Km_ID),ziformula=~poly(depth,2)+(1|PU_1Km_ID),data=dive_df,family=beta_family())

AIC(mDvTMBrand,mDvTMBrand2,mDvTMBrand3,mDvTMBrand4)
dive_df$predicted_counts <- predict(mDvTMBrand,type='r')
plot(dive_df$counts,dive_df$predicted_counts)

mDpTMBrand <- glmmTMB(counts~depth+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID),dispformula = ~1,ziformula = ~1,data=deep_df,family=nbinom2)
mDpTMBrand2 <- glmmTMB(counts~depth+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID),dispformula = ~1,ziformula = ~poly(depth,2)+(1|PU_1Km_ID),data=deep_df,family=nbinom2)
mDpTMBrand3 <- glmmTMB(counts~depth+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID),dispformula = ~1,data=deep_df,family=nbinom2)
AIC(mDpTMBrand,mDpTMBrand2,mDpTMBrand3)
deep_df$predicted_counts <- predict(mDpTMBrand,type='r')
plot(deep_df$counts,deep_df$predicted_counts)
new_df2 <- rbind(dive_df,mid_df,deep_df)
new_df2$lambda <- ifelse(new_df2$gear=="deep_video",new_df2$predicted_counts/new_df2$effort,new_df2$predicted_counts)
new_df2$cpue <- ifelse(new_df2$gear=="deep_video",new_df2$counts/new_df2$effort,new_df2$counts)
write.csv(new_df2,"Data/sponge counts by sample id.csv")

plot(predicted_counts~counts,data=new_df2)
abline(b=1,a=0)
plot(lambda~cpue,data=new_df2)
abline(b=1,a=0)

hotspots_df <- aggregate(depth~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df2,FUN=mean)
hotspots_df$effort <- aggregate(effort~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df2,FUN=sum)$effort
hotspots_df$counts <- aggregate(counts~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df2,FUN=sum)$counts
hotspots_df$sample_sizes <- aggregate(sample_size~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df2,FUN=sum)$sample_size
hotspots_df$p_counts <- aggregate(predicted_counts~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df2,FUN=sum)$predicted_counts
hotspots_df$lambda <- hotspots_df$p_counts/hotspots_df$effort
hotspots_df$cpue <- hotspots_df$counts/hotspots_df$effort

hotspots_df <- hotspots_df %>%
  group_by(species,gear) %>% 
  mutate(normalized_lambda = normalize(lambda),normalized_cpue = normalize(cpue)) %>%
  ungroup(species,gear)
hotspots_df <- hotspots_df[complete.cases(hotspots_df),]
hotspots_agg <- aggregate(cbind(normalized_lambda,normalized_cpue)~PU_4Km_ID+species,data=hotspots_df,FUN=mean)
hotspots_agg$depth <- aggregate(depth~PU_4Km_ID+species,data=hotspots_df,FUN=mean)$depth
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

#write.csv(hotspots_agg,"Data/sponge normalized cpue by site.csv")

plot(normalized_lambda~normalized_cpue,data=hotspots_df)
summary(lm(normalized_lambda~normalized_cpue,data=hotspots_df))
plot(normalized_lambda~normalized_cpue,data=hotspots_agg)
summary(lm(normalized_lambda~normalized_cpue,data=hotspots_agg))

hotspots_coast <- aggregate(cbind(normalized_lambda,normalized_cpue)~PU_4Km_ID,data=hotspots_agg,FUN=function(x){sum(x)/length(coral_spp)})
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

simulationOutput <- DHARMa::simulateResiduals(fittedModel = m3fTMBrand7,plot = F,n=1000)
plot(simulationOutput)
DHARMa::testZeroInflation(simulationOutput)
DHARMa::testDispersion(simulationOutput)
DHARMa::testQuantiles(simulationOutput)
DHARMa::testOutliers(simulationOutput)
