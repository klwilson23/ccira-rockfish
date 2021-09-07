normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")


library(dplyr)
library(glmmTMB)
library(lme4)
library(brms)

new_df <- read.csv("Data/Rockfish counts PU4km v2.csv")
km1_puid <- read.csv("Data/1kmPUID.csv")
km1_puid$X1km.puid[!km1_puid$X1km.puid%in%new_df$PU_1Km_ID]

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
m3fTMB <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e5,eval.max=5e5)))
m3fTMBrand <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/species),dispformula = ~gear,data=new_df,family=nbinom2)
m3fTMBrand2 <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,data=new_df,family=nbinom2)
m3fTMBrand3 <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
m3fTMBrand4 <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/species/gear),dispformula = ~gear,data=new_df,family=nbinom2)
m3fTMBrand5 <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2)
m3fTMBrandzi <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~(1|PU_1Km_ID/gear:species),data=new_df,family=nbinom2)
m3fTMBrandzi2 <- glmmTMB(counts~poly(depth,2)+species+gear+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~(1|PU_1Km_ID/species),data=new_df,family=nbinom2)

m3gTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e5,eval.max=5e5)))
m3hTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula = ~1,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3iTMB <- glmmTMB(counts~species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2)
m3jTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,ziformula = ~species,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3kTMB <- glmmTMB(counts~depth+depth:species+species+gear+as.factor(survey_id)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))
m3lTMB <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+UpperOceanSR+gear+as.factor(PU_4Km_ID)+offset(log(effort)),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))


m1 <- glm(counts~depth+depth:species+I(depth^2)+I(depth^2):species+species*gear+as.factor(PU_4Km_ID)+offset(log(effort)),data=new_df,family=poisson(link="log"))

#m3fbrms <- brm(bf(counts~depth+depth:species+I(depth^2)+species+gear+as.factor(PU_4Km_ID)+offset(log(effort)),shape~gear),data=new_df,family=negbinomial())

#m3f_rand_brms <- brm(bf(counts~depth+depth:species+I(depth^2)+species+gear+offset(log(effort))+(1|as.factor(PU_4Km_ID)),shape~gear),data=new_df,family=negbinomial())

#m3fTMB_rand <- glmmTMB(counts~depth+depth:species+I(depth^2)+species+offset(log(effort))+1|PU_4Km_ID/gear,dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=5e4,eval.max=5e4)))

curve(exp(-3.33e+00+-5.46e-01-2.51e-02*x+7.49e-02*x-2.18e-04*x^2+log(120)),lwd=2,col="blue",from=0,to=500,ylab="Expected counts per transect") # quillback
curve(exp(-3.33e+00+-3.33e+00-2.51e-02*x+9.16e-02*x-2.18e-04*x^2+log(120)),lwd=2,col="tomato",from=0,to=500,add=TRUE) # yelloweye
curve(exp(-3.33e+00+-1.99e+01-2.51e-02*x+1.56e-01*x-2.18e-04*x^2+log(120)),lwd=2,col="orange",from=0,to=500,ylab="Expected counts per transect",add=TRUE) # redbanded

# C=cpue*effort
# cpue = exp(coefficients*x)
# coefficinets{depth,species,gear,location}
m2 <- glm(counts~depth+I(depth^2)+species*gear+PU_4Km_ID,data=new_df,family=poisson(link="log"),weights=sqrt(new_df$sample_size))
m1a <- glmer(counts~depth+species+gear+offset(log(effort))+(1|PU_4Km_ID),data=new_df,family=poisson(link="log"),lme4::glmerControl(optCtrl=list(method='nlminb')),weights=sqrt(new_df$sample_size))
m1b <- glmer(counts~depth+I(depth^2)+gear+PU_4Km_ID+offset(log(effort))+(1|species),data=new_df,family=poisson(link="log"),lme4::glmerControl(optCtrl=list(method='nlminb')),weights=sqrt(new_df$sample_size))

AIC(m1a,m1b,m2,m1TMB,m2TMB,m2aTMB,m2bTMB,m2cTMB,m3aTMB,m3bTMB,m3cTMB,m3dTMB,m3eTMB,m3fTMB,m3fTMBrand,m3fTMBrand2,m3fTMBrand3,m3fTMBrand4,m3gTMB,m3hTMB,m3iTMB,m3jTMB,m3kTMB,m3lTMB)

new_df$predicted_counts <- predict(m3fTMBrand5,type='r')
new_df$lambda <- new_df$predicted_counts/new_df$effort
new_df$cpue <- new_df$counts/new_df$effort
write.csv(new_df,"Data/Rockfish counts by sample id.csv")

plot(predicted_counts~counts,data=new_df)
abline(b=1,a=0)
plot(lambda~cpue,data=new_df)
abline(b=1,a=0)

simulationOutput <- DHARMa::simulateResiduals(fittedModel = m3fTMBrand5,plot = F,n=1000)
plot(simulationOutput)
DHARMa::testZeroInflation(simulationOutput)
DHARMa::testDispersion(simulationOutput)
DHARMa::testQuantiles(simulationOutput)
DHARMa::testOutliers(simulationOutput)
