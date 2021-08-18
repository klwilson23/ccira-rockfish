normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")


library(dplyr)
library(glmmTMB)
library(lme4)
library(brms)

new_df <- read.csv("Data/Rockfish counts PU4km v2.csv")
#new_df <- new_df[!complete.cases(new_df),]
m1TMB <- glmmTMB(counts~depth+species+offset(log(effort)),data=new_df,family=poisson)
m2TMB <- glmmTMB(counts~depth+species+offset(log(effort)),data=new_df,family=nbinom2)
m2aTMB <- glmmTMB(counts~depth+species+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m2bTMB <- glmmTMB(counts~depth+species+offset(log(effort)),dispformula = ~species,data=new_df,family=nbinom2)
m2cTMB <- glmmTMB(counts~depth+species+offset(log(effort)),dispformula = ~species*gear,data=new_df,family=nbinom2)
m3aTMB <- glmmTMB(counts~depth+I(depth^2)+species+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3bTMB <- glmmTMB(counts~depth+I(depth^2)+species+gear+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3cTMB <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3dTMB <- glmmTMB(counts~depth+depth:species+species+gear+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3eTMB <- glmmTMB(counts~depth+depth:species+I(depth^2)+I(depth^2):species+species+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3fTMB <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3gTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3hTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula = ~1,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3iTMB <- glmmTMB(counts~species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3jTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula = ~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3kTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(survey_id)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))

m1 <- glm(counts~depth+depth:species+I(depth^2)+I(depth^2):species+species*gear+as.factor(PU_4Km_ID)+offset(log(effort)),data=new_df,family=poisson(link="log"),weights=sqrt(new_df$sample_size))

#m3fbrms <- brm(bf(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),shape~gear),data=new_df,family=negbinomial())

#m3f_rand_brms <- brm(bf(counts~depth+depth:species+I(depth^2)+species+gear+offset(log(effort))+(1|as.factor(PU_4Km_ID)),shape~gear),data=new_df,family=negbinomial())

curve(exp(-3.33e+00+-5.46e-01-2.51e-02*x+7.49e-02*x-2.18e-04*x^2+log(120)),lwd=2,col="blue",from=0,to=500,ylab="Expected counts per transect") # quillback
curve(exp(-3.33e+00+-3.33e+00-2.51e-02*x+9.16e-02*x-2.18e-04*x^2+log(120)),lwd=2,col="tomato",from=0,to=500,add=TRUE) # yelloweye
curve(exp(-3.33e+00+-1.99e+01-2.51e-02*x+1.56e-01*x-2.18e-04*x^2+log(120)),lwd=2,col="orange",from=0,to=500,ylab="Expected counts per transect",add=TRUE) # redbanded

# C=cpue*effort
# cpue = exp(coefficients*x)
# coefficinets{depth,species,gear,location}
m2 <- glm(counts~depth+I(depth^2)+species*gear+PU_4Km_ID,data=new_df,family=poisson(link="log"),weights=sqrt(new_df$sample_size))
m1a <- glmer(counts~depth+species+gear+offset(log(effort))+(1|PU_4Km_ID),data=new_df,family=poisson(link="log"),lme4::glmerControl(optCtrl=list(method='nlminb')),weights=sqrt(new_df$sample_size))
m1b <- glmer(counts~depth+I(depth^2)+gear+PU_4Km_ID+offset(log(effort))+(1|species),data=new_df,family=poisson(link="log"),lme4::glmerControl(optCtrl=list(method='nlminb')),weights=sqrt(new_df$sample_size))

AIC(m1,m1a,m1b,m2,m1TMB,m2TMB,m2aTMB,m2bTMB,m2cTMB,m3aTMB,m3bTMB,m3cTMB,m3dTMB,m3eTMB,m3fTMB,m3gTMB,m3hTMB,m3iTMB,m3jTMB,m3kTMB)

new_df$predicted_counts <- predict(m1,type='r')
new_df$predicted_counts_mm <- predict(m1a,type='r')
new_df$predicted_counts_TMB <- predict(m3gTMB,type='r')
new_df$predicted_counts_TMB2 <- predict(m3fTMB ,type='r')

plot(predicted_counts_mm~predicted_counts,data=new_df)
abline(b=1,a=0)
plot(predicted_counts_TMB~predicted_counts,data=new_df)
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

hotspots_df$p_counts <- predict(m3fTMB,newdata = hotspots_df,type="r")
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
hotspots_agg$sample_sizes <- aggregate(sample_sizes~PU_4Km_ID+species,data=hotspots_df,FUN=sum)$sample_sizes
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

simulationOutput <- DHARMa::simulateResiduals(fittedModel = m3fTMB,plot = F)
plot(simulationOutput)
DHARMa::testZeroInflation(simulationOutput)
DHARMa::testDispersion(simulationOutput)
DHARMa::testQuantiles(simulationOutput)
DHARMa::testOutliers(simulationOutput)
