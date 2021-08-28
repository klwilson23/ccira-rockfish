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
#new_df <- new_df[!complete.cases(new_df),]

#new_df <- new_df[!complete.cases(new_df),]
m1TMB <- glmmTMB(counts~depth+species+offset(log(effort)),data=new_df,family=poisson)
m2TMB <- glmmTMB(counts~depth+species+offset(log(effort)),data=new_df,family=nbinom2)
m2aTMB <- glmmTMB(counts~depth+species+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m2bTMB <- glmmTMB(counts~depth+species+offset(log(effort)),dispformula = ~species,data=new_df,family=nbinom2)
m2cTMB <- glmmTMB(counts~depth+species+offset(log(effort)),dispformula = ~species*gear,data=new_df,family=nbinom2)
m3aTMB <- glmmTMB(counts~depth+I(depth^2)+species+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3bTMB <- glmmTMB(counts~depth+I(depth^2)+species+gear+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3cTMB <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3dTMB <- glmmTMB(counts~depth+depth:species+species+gear+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3eTMB <- glmmTMB(counts~depth+depth:species+I(depth^2)+I(depth^2):species+species+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3fTMB <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3gTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3hTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula = ~1,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3iTMB <- glmmTMB(counts~species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3jTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula = ~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3kTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(survey_id)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))

#m3eTMB_rand <- glmmTMB(counts~depth+depth:species+I(depth^2)+I(depth^2):species+species+offset(log(effort)) + 1|PU_4Km_ID,dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3cTMBzi <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+offset(log(effort)),dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3cTMBzi2 <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+offset(log(effort)),dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3cTMBzi3 <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3cTMBzi4 <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3cTMBzi5 <- glmmTMB(counts~depth+depth:species+species+gear+offset(log(effort)),dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3cTMBzi6 <- glmmTMB(counts~depth+depth:species+I(depth^2)+I(depth^2):species+species+gear+offset(log(effort)),dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3cTMBzi7 <- glmmTMB(counts~depth+depth:species+species+gear+offset(log(effort)),dispformula = ~gear,ziformula=~as.factor(PU_4Km_ID),data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3cTMBzi8 <- glmmTMB(counts~as.factor(PU_4Km_ID)+offset(log(effort))+species+gear,dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3cTMBzi9 <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+offset(log(effort))+as.factor(PU_4Km_ID),dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e5,eval.max=5e5)))
m3cTMBzi10 <- glmmTMB(counts~poly(depth,2)*species+gear+offset(log(effort)),dispformula=~gear*species,ziformula=~as.factor(PU_4Km_ID),data=new_df,family=truncated_nbinom2(),control=glmmTMBControl(optCtrl=list(iter.max=5e5,eval.max=5e5)))
m3cTMBzi11 <- glmmTMB(counts~depth+depth:species+species+gear+offset(log(effort)),dispformula = ~gear,ziformula=~gear*effort+species*poly(depth,2),data=new_df,family=truncated_nbinom2(),control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))

#m3f_rand_brms <- brm(bf(counts~depth+depth:species+I(depth^2)+species+gear+offset(log(effort))+(1|PU_1Km_ID),shape~gear,zi~species),data=new_df,family=zero_inflated_negbinomial(),iter=2000,chains=4,cores=4)

m1 <- glm(counts~depth+depth:species+I(depth^2)+I(depth^2):species+species*gear+as.factor(PU_4Km_ID)+offset(log(effort)),data=new_df,family=poisson(link="log"),weights=sqrt(new_df$sample_size))

#m3fbrms <- brm(bf(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),shape~gear),data=new_df,family=negbinomial())

#m3fTMB_rand <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+offset(log(effort))+1|PU_4Km_ID/gear,dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))

# C=cpue*effort
# cpue = exp(coefficients*x)
# coefficinets{depth,species,gear,location}
m2 <- glm(counts~depth+I(depth^2)+species*gear+PU_4Km_ID,data=new_df,family=poisson(link="log"),weights=sqrt(new_df$sample_size))
m1a <- glmer(counts~depth+species+gear+offset(log(effort))+(1|PU_4Km_ID),data=new_df,family=poisson(link="log"),lme4::glmerControl(optCtrl=list(method='nlminb')),weights=sqrt(new_df$sample_size))
m1b <- glmer(counts~depth+I(depth^2)+gear+PU_4Km_ID+offset(log(effort))+(1|species),data=new_df,family=poisson(link="log"),lme4::glmerControl(optCtrl=list(method='nlminb')),weights=sqrt(new_df$sample_size))

AIC(m1,m1a,m1b,m2,m1TMB,m2TMB,m2aTMB,m2bTMB,m2cTMB,m3aTMB,m3bTMB,m3cTMB,m3dTMB,m3eTMB,m3fTMB,m3gTMB,m3hTMB,m3iTMB,m3jTMB,m3kTMB,m3cTMBzi,m3cTMBzi2,m3cTMBzi3,m3cTMBzi4,m3cTMBzi5,m3cTMBzi6,m3cTMBzi7,m3cTMBzi8)

new_df$predicted_counts <- predict(m1,type='r')
new_df$predicted_counts_TMB_space <- predict(m3cTMBzi7,type='r')
new_df$predicted_counts_TMB2 <- predict(m3cTMBzi10,type='r')
new_df$predicted_counts_TMB <- predict(m3cTMBzi ,type='r')

write.csv(new_df,"Data/coral counts by sample id.csv")

plot(predicted_counts_TMB_space~predicted_counts,data=new_df)
abline(b=1,a=0)
plot(predicted_counts_TMB~predicted_counts_TMB_space,data=new_df)
abline(b=1,a=0)
plot(predicted_counts_TMB~predicted_counts_TMB2,data=new_df)
abline(b=1,a=0)
new_df$lambda <- new_df$predicted_counts_TMB2/new_df$effort
new_df$cpue <- new_df$counts/new_df$effort
plot(predicted_counts_TMB2~counts,data=new_df)
abline(b=1,a=0)
plot(lambda~cpue,data=new_df)
abline(b=1,a=0)

hotspots_df <- aggregate(depth~PU_4Km_ID+gear+species,data=new_df,FUN=mean)
hotspots_df$effort <- aggregate(effort~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$effort
hotspots_df$counts <- aggregate(counts~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$counts
hotspots_df$sample_sizes <- aggregate(sample_size~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$sample_size
hotspots_df$p_counts_2 <- aggregate(predicted_counts_TMB2~PU_4Km_ID+gear+species,data=new_df,FUN=sum)$predicted_counts_TMB2

hotspots_df$p_counts <- predict(m3fTMB,newdata = hotspots_df,type="r")
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
hotspots_agg$sample_sizes <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df,FUN=sum)$sample_sizes
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
hotspots_coast$hotspot <- ifelse(hotspots_coast$quantilelambda>=9,1,0)
hotspots_coast$hotspot_old <- ifelse(hotspots_coast$quantileCPUE>=9,1,0)
sum((hotspots_coast$hotspot-hotspots_coast$hotspot_old)==1)
sum((hotspots_coast$hotspot-hotspots_coast$hotspot_old)==-1)
sum((hotspots_coast$hotspot-hotspots_coast$hotspot_old)==0)
sum((hotspots_coast$hotspot+hotspots_coast$hotspot_old)==2)

simulationOutput <- DHARMa::simulateResiduals(fittedModel = m3cTMBzi7,plot = F)
plot(simulationOutput)
DHARMa::testZeroInflation(simulationOutput)
DHARMa::testDispersion(simulationOutput)
DHARMa::testQuantiles(simulationOutput)
DHARMa::testOutliers(simulationOutput)
