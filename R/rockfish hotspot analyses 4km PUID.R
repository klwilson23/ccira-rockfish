#devtools::install_github("James-Thorson-NOAA/FishLife")
library(FishLife)
library(TMB)
library(dplyr)

max_only <- function(x, na.rm = TRUE) {
  return((x- 0) /(max(x, na.rm = TRUE)-0))
}

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail","lingcod")

rockfish_vuln <- read.csv("Data/Fish life histories for scoring v2.csv")
rockfish_vuln <- rockfish_vuln[,-ncol(rockfish_vuln)]
rockfish_vuln$Common.name[rockfish_vuln$Common.name=="silvergrey"] <- "silvergray"
rockfish_ncc <- read.csv("Data/rockfish normalized cpue by 4km puid and species.csv",header=TRUE)
rockfish_ncc <- rockfish_ncc[!rockfish_ncc%in%c("splitnose","stripetail"),]
rockfish_ncc2 <- read.csv("Data/rockfish normalized cpue by 1km puid and species.csv",header=TRUE)
rockfish_ncc2 <- rockfish_ncc2[!rockfish_ncc2%in%c("splitnose","stripetail"),]

m1 <- glm(Loo_fem~log(maxTL),family=gaussian(link="log"),data=rockfish_vuln)
with(summary(m1), 1 - deviance/null.deviance)
coef(m1)
rockfish_vuln$Loo <- rockfish_vuln$Loo_fem
rockfish_vuln$Loo[is.na(rockfish_vuln$Loo_fem)] <- predict(m1,newdata=rockfish_vuln[is.na(rockfish_vuln$Loo_fem),],type="response")
plot(Loo_fem~maxTL,data=rockfish_vuln)
thorson_traits <- as.data.frame(FishLife::FishBase_and_RAM$beta_gv)
thorson_covar <- as.data.frame(FishLife::FishBase_and_RAM$Cov_gvv)
thorson_covar <- thorson_covar[,intersect(colnames(thorson_covar),paste(colnames(thorson_traits),".",colnames(thorson_traits),sep=""))]
colnames(thorson_covar) <- colnames(thorson_traits)
leading_traits <- c("Loo","K","Winfinity","tmax","tm","M","Lm","ln_var","ln_MASPS","ln_margsd","ln_Fmsy_over_M","ln_Fmsy")
thorson_traits[,leading_traits] <- exp(thorson_traits[,leading_traits]+0.5*thorson_covar[,leading_traits])
taxa <- "Actinopterygii_Scorpaeniformes_Sebastidae"
rockfish_rows <- row.names(thorson_traits)[grep(taxa,row.names(thorson_traits))]
rockfish <- thorson_traits[rockfish_rows,]
rockfish[,c("Loo","K","tmax","tm","M","h","r","ln_Fmsy_over_M","ln_Fmsy","G")]

plot(r~tmax,data=rockfish,xlab="Maximum age (all Sebastidae)",ylab="Intrinsic rate of growth",pch=21,bg="grey90",ylim=c(0,1))
m1 <- glm(r~log(tmax)+log(Loo),family=gaussian(link="log"),data=rockfish)
summary(m1)
with(summary(m1), 1 - deviance/null.deviance)
curve(exp(coef(m1)[1] + coef(m1)[2]*log(x)+coef(m1)[3]*log(mean(rockfish$Loo))),add=TRUE)
summary(lm(log(r)~log(tmax),data=rockfish))
summary(lm(log(r)~log(tmax)+log(K),data=rockfish))

rockfish_vuln$tmax <- rockfish_vuln$Max.Age
rockfish_vuln$r <- predict(m1,newdata = rockfish_vuln,type="r")
rockfish_vuln$vuln <- 1*max_only(1/rockfish_vuln$r)
plot(vuln~Vulnerability.Cheung,data=rockfish_vuln)
plot(vuln~Vulnerability.Magnuson,data=rockfish_vuln)
summary(lm(vuln~Vulnerability.Magnuson,data=rockfish_vuln))

plot(1/r~tmax,data=rockfish_vuln,xlab="Maximum age (all Sebastidae)",ylab="Inverse of Rmax",pch=21,bg="grey90")
plot(1/r~Loo,data=rockfish_vuln,xlab="Maximum age (all Sebastidae)",ylab="Inverse of Rmax",pch=21,bg="grey90")

rockfish_ncc$vuln <- rockfish_vuln$vuln[match(rockfish_ncc$species,rockfish_vuln$Common.name)]
rockfish_ncc$trophic <- 1*normalize(rockfish_vuln$Trophic.level[match(rockfish_ncc$species,rockfish_vuln$Common.name)])
rockfish_ncc$depletion<- 1*max_only(1-rockfish_vuln$By.B0[match(rockfish_ncc$species,rockfish_vuln$Common.name)])
rockfish_ncc$evo_dist<- 1*max_only(rockfish_vuln$Evolutionary.Distinctiveness[match(rockfish_ncc$species,rockfish_vuln$Common.name)])
weights <- c(4,2,4,2) # lambda, vulnerability, trophic level, depletion, evolutionary distinctiveness

spp_scores <- t(weights*t(cbind(rockfish_vuln$vuln,normalize(rockfish_vuln$Trophic.level),max_only(1-rockfish_vuln$By.B0),max_only(rockfish_vuln$Evolutionary.Distinctiveness))))
colnames(spp_scores) <- c("weighted_vuln","weighted_trophic","weighted_depleted","weighted_evol")

rockfish_vuln <- cbind(rockfish_vuln,spp_scores)
rockfish_vuln$scores <- rowSums(rockfish_vuln[,c("weighted_vuln","weighted_trophic","weighted_depleted","weighted_evol")],na.rm=TRUE)
rockfish_vuln$total <- rowSums(t(weights*t(!is.na(rockfish_vuln[,c("weighted_vuln","weighted_trophic","weighted_depleted","weighted_evol")]))),na.rm=TRUE)
rockfish_vuln$weighted_score <- rockfish_vuln$scores/rockfish_vuln$total
rockfish_vuln <- rockfish_vuln[order(rockfish_vuln$weighted_score,decreasing = TRUE),]

rockfish_ncc$scores <- rowSums(t(weights*t(cbind(rockfish_ncc$vuln,rockfish_ncc$trophic,rockfish_ncc$depletion,rockfish_ncc$evo_dist))),na.rm=TRUE)
rockfish_ncc$total <- rowSums(t(weights*t(!is.na(cbind(rockfish_ncc$vuln,rockfish_ncc$trophic,rockfish_ncc$depletion,rockfish_ncc$evo_dist)))),na.rm=TRUE)
rockfish_ncc$weighted_score <- rockfish_ncc$scores/rockfish_ncc$total
rockfish_ncc$raw_score <- rockfish_ncc$normalized_lambda*rockfish_ncc$weighted_score
rockfish_ncc$obs_score <- rockfish_ncc$normalized_cpue*rockfish_ncc$weighted_score

rockfish_ncc2$vuln <- rockfish_vuln$vuln[match(rockfish_ncc2$species,rockfish_vuln$Common.name)]
rockfish_ncc2$trophic <- 1*normalize(rockfish_vuln$Trophic.level[match(rockfish_ncc2$species,rockfish_vuln$Common.name)])
rockfish_ncc2$depletion<- 1*max_only(1-rockfish_vuln$By.B0[match(rockfish_ncc2$species,rockfish_vuln$Common.name)])
rockfish_ncc2$evo_dist<- 1*max_only(rockfish_vuln$Evolutionary.Distinctiveness[match(rockfish_ncc2$species,rockfish_vuln$Common.name)])
rockfish_ncc2$scores <- rowSums(t(weights*t(cbind(rockfish_ncc2$vuln,rockfish_ncc2$trophic,rockfish_ncc2$depletion,rockfish_ncc2$evo_dist))),na.rm=TRUE)
rockfish_ncc2$total <- rowSums(t(weights*t(!is.na(cbind(rockfish_ncc2$vuln,rockfish_ncc2$trophic,rockfish_ncc2$depletion,rockfish_ncc2$evo_dist)))),na.rm=TRUE)
rockfish_ncc2$weighted_score <- rockfish_ncc2$scores/rockfish_ncc2$total
rockfish_ncc2$raw_score <- rockfish_ncc2$normalized_lambda*rockfish_ncc2$weighted_score
rockfish_ncc2$obs_score <- rockfish_ncc2$normalized_cpue*rockfish_ncc2$weighted_score

hotspots_coast <- aggregate(raw_score~PU_4Km_ID+PU_1Km_ID,data=rockfish_ncc2,function(x){sum(x)})
hotspots_coast <- aggregate(raw_score~PU_4Km_ID,data=hotspots_coast,function(x){mean(max(x),ifelse(sum(x[-which.max(x)])>length(which.max(x)),mean(x[-which.max(x)]),mean(x)))})
hotspots_coast$obs_score <- aggregate(obs_score~PU_4Km_ID,data=rockfish_ncc,function(x){sum(x)})$obs_score
hotspots_coast$UpperOceanSR <- rockfish_ncc$UpperOceanSR[match(hotspots_coast$PU_4Km_ID,rockfish_ncc$PU_4Km_ID)]
hotspots_coast$depth <- aggregate(depth~PU_4Km_ID,data=rockfish_ncc,function(x){mean(x)})$depth
hotspots_coast$max_depth <- aggregate(max_depth~PU_4Km_ID,data=rockfish_ncc,max)$max_depth
hotspots_coast$dive_samps <- aggregate(dive_samps~PU_4Km_ID,data=rockfish_ncc,function(x){mean(x)})$dive_samps
hotspots_coast$hook_samps <- aggregate(hook_samps~PU_4Km_ID,data=rockfish_ncc,function(x){mean(x)})$hook_samps
hotspots_coast$mid_samps <- aggregate(mid_samps~PU_4Km_ID,data=rockfish_ncc,function(x){mean(x)})$mid_samps
hotspots_coast$deep_samps <- aggregate(deep_samps~PU_4Km_ID,data=rockfish_ncc,function(x){mean(x)})$deep_samps
hotspots_coast$normalized_obs_score <- max_only(hotspots_coast$obs_score)

hotspots_coast$normalized_score <- max_only(hotspots_coast$raw_score)
hotspots_coast <- hotspots_coast %>%
  mutate(rockfish_rank = ntile(normalized_score, 10))
hotspots_coast$rockfish_hotspot <- ifelse(hotspots_coast$rockfish_rank>=9,1,0)
hotspots_coast$rockfish_hotspot2 <- ifelse(hotspots_coast$normalized_score>=0.01,1,0)

oldhotties <- c(6267, 6425, 6582, 6736, 7054, 7524, 7529, 7530, 7847, 7997, 8003, 8157, 8159, 8625, 8628, 8629, 8783, 8786, 8792, 8950, 9102, 9108, 9259, 9265, 9418, 9730, 9731, 9733, 9734, 9889, 9906, 10367, 10524, 10674, 10675, 11331, 11461, 11629, 11774, 11789, 11944, 12093, 12242, 12400, 12410, 12566, 12741, 12759, 12895, 12897, 13034, 13055, 13194, 13212, 13213, 14312)

sum(hotspots_coast[hotspots_coast$PU_4Km_ID%in%oldhotties,"rockfish_rank"]>=9)/sum(hotspots_coast$rockfish_rank>=9)
write.csv(hotspots_coast,"Data/Rockfish hotspots 4km PUID.csv")
write.csv(rockfish_vuln,"Data/Rockfish species scores.csv")
