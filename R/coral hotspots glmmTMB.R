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

new_df <- read.csv("Data/coral counts PU4km.csv")
m3fTMB <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))

m3fTMB2 <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_1Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMB3 <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula =~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMB4 <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula =~1,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMBrand <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/gear/species),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMBrand2 <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/gear/species),dispformula = ~gear,ziformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMBrand3 <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/species/gear),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMBrand4 <- glmmTMB(counts~depth+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/gear/species),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMBrand5 <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_4Km_ID/gear/species),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMBrand6 <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/species),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMBrand7 <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))
m3fTMBrand8 <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e6,eval.max=5e6)))

AIC(m3fTMB,m3fTMB2,m3fTMB3,m3fTMB4,m3fTMBrand,m3fTMBrand2,m3fTMBrand3,m3fTMBrand4,m3fTMBrand5,m3fTMBrand6,m3fTMBrand7,m3fTMBrand7zi,m3fTMBrand8)

new_df$predicted_counts <- predict(m3fTMBrand7,type='r')
new_df$lambda <- new_df$predicted_counts/new_df$effort
new_df$cpue <- new_df$counts/new_df$effort
write.csv(new_df,"Data/coral counts by sample id.csv")

plot(predicted_counts~counts,data=new_df)
abline(b=1,a=0)
plot(lambda~cpue,data=new_df)
abline(b=1,a=0)

hotspots_df <- aggregate(depth~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df,FUN=mean)
hotspots_df$effort <- aggregate(effort~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df,FUN=sum)$effort
hotspots_df$counts <- aggregate(counts~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df,FUN=sum)$counts
hotspots_df$sample_sizes <- aggregate(sample_size~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df,FUN=sum)$sample_size
hotspots_df$p_counts_2 <- aggregate(predicted_counts~PU_4Km_ID+PU_1Km_ID+gear+species,data=new_df,FUN=sum)$predicted_counts

hotspots_df$p_counts <- predict(m3fTMBrand6,newdata = hotspots_df,type="r")
hotspots_df$lambda <- hotspots_df$p_counts_2/hotspots_df$effort
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
plot(normalized_lambda~depth,hotspots_df[hotspots_df$species=="paragorgia",])
plot(normalized_cpue~depth,hotspots_df[hotspots_df$species=="paragorgia",])

write.csv(hotspots_agg,"Data/coral normalized cpue by site and species.csv")

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
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m3fTMBrand7,plot = F,n=1000)
plot(simulationOutput)
DHARMa::testZeroInflation(simulationOutput)
DHARMa::testDispersion(simulationOutput)
DHARMa::testQuantiles(simulationOutput)
DHARMa::testOutliers(simulationOutput)
