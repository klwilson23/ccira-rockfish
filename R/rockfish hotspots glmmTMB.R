normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail","lingcod")


library(dplyr)
library(glmmTMB)
library(lme4)
library(brms)

new_df <- read.csv("Data/Rockfish counts PU4km v2.csv")
km1_puid <- read.csv("Data/1kmPUID.csv")
km1_puid$X1km.puid[!km1_puid$X1km.puid%in%new_df$PU_1Km_ID]

#new_df <- new_df[!complete.cases(new_df),]
m3fTMBrand <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/species),dispformula = ~gear,data=new_df,family=nbinom2)
m3fTMBrand2 <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,data=new_df,family=nbinom2)
m3fTMBrand3 <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
m3fTMBrand4 <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/species/gear),dispformula = ~gear,data=new_df,family=nbinom2)
m3fTMBrand5 <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~species,data=new_df,family=nbinom2)
m3fTMBrandzi <- glmmTMB(counts~depth+depth:species+I(depth^2)+gear+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~(1|PU_1Km_ID/gear:species),data=new_df,family=nbinom2)
m3fTMBrandzi2 <- glmmTMB(counts~poly(depth,2)+species+gear+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~(1|PU_1Km_ID/species),data=new_df,family=nbinom2)
m3fTMBrandzi3 <- glmmTMB(counts~poly(depth,2)+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~(1|PU_1Km_ID/species),data=new_df,family=nbinom2)
m3fTMBrandzi4 <- glmmTMB(counts~depth+depth:species+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~(1|PU_1Km_ID/species),data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
m3fTMBrandzi5 <- glmmTMB(counts~depth+depth:species+I(depth^2)+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~(1|PU_1Km_ID/species),data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
m3fTMBrandzi6 <- glmmTMB(counts~poly(depth,2)+depth:species+offset(log(effort))+(1|PU_1Km_ID/gear:species),dispformula = ~gear,ziformula=~(1|PU_1Km_ID/species),data=new_df,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
AIC(m3fTMBrand,m3fTMBrand2,m3fTMBrand3,m3fTMBrand4,m3fTMBrand5,m3fTMBrandzi,m3fTMBrandzi2,m3fTMBrandzi3,m3fTMBrandzi4,m3fTMBrandzi5)

new_df$predicted_counts <- predict(m3fTMBrandzi4,type='r')
new_df$lambda <- new_df$predicted_counts/new_df$effort
new_df$cpue <- new_df$counts/new_df$effort
write.csv(new_df,"Data/Rockfish counts by sample id.csv")

plot(predicted_counts~counts,data=new_df)
abline(b=1,a=0)
plot(lambda~cpue,data=new_df)
abline(b=1,a=0)

simulationOutput <- DHARMa::simulateResiduals(fittedModel = m3fTMBrandzi4,plot = F,n=1000)
plot(simulationOutput)
DHARMa::testZeroInflation(simulationOutput)
DHARMa::testDispersion(simulationOutput)
DHARMa::testQuantiles(simulationOutput)
DHARMa::testOutliers(simulationOutput)
