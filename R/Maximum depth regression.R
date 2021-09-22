library(MASS)
library(ggplot2)
library(tidyr)
library(ggExtra)
library(cowplot)

##### data ####
data = read.csv("Data/NCC hotspots 1km puid.csv")
data$PU_1Km_ID <- data$PUID
data$SebHotspot <- data$rf_hotspot
data$CorHotspot <- data$cor_hotspot
data$SpHotspot <- data$sponge_hotspot
data$OvHotspot <- data$all_hotspots

data$OceanSR <- data$UpperOceanSR
data$maxDepth <- data$max_depth
data$OceanSR = factor(data$OceanSR,levels = c("(11) Mainland Fjords","(10 Aristazabal Banks Upwelling","(13) Eastern Queen Charlotte Sound"))
data2 = subset(data, data$OceanSR!="(10 Aristazabal Banks Upwelling")
data2$OceanSR = factor(data2$OceanSR,levels = c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound"))

#I re-ordered the factors because I wanted to compare MF to to other two.
#for corals and spongens I had to drop aristazabal entirely - see data2

#binomial regression, but more complicated ####
m1 = glm(SebHotspot~maxDepth+OceanSR,data=data,family=binomial)
m1poly = glm(SebHotspot~poly(maxDepth,2)+OceanSR,data=data,family=binomial)
m1polint = glm(SebHotspot~poly(maxDepth,2)*OceanSR,data=data,family=binomial)

AIC(m1,m1poly,m1polint)
m2 = glm(CorHotspot~maxDepth+OceanSR,data=data2,family=binomial,na.action = na.omit)
m2poly = glm(CorHotspot~poly(maxDepth,2)+OceanSR,data=data2,family=binomial,na.action = na.omit)
AIC(m2,m2poly)
m3 = glm(SpHotspot~maxDepth+OceanSR,data=data2,family=binomial,na.action = na.omit)
m3poly = glm(SpHotspot~poly(maxDepth,2)+OceanSR,data=data2,family=binomial,na.action = na.omit)
AIC(m3,m3poly)

m4 = glm(OvHotspot~maxDepth+OceanSR,data=data2,family=binomial,na.action = na.omit)
m4poly = glm(OvHotspot~poly(maxDepth,2)+OceanSR,data=data2,family=binomial,na.action = na.omit)
AIC(m4,m4poly)

summary(m1poly)
summary(m2poly)
summary(m3poly)
summary(m4poly)

# plots
newdat = data.frame(
  puid = data$PU_1Km_ID,
  SebHotspot = data$SebHotspot,
  CorHotspot = data$CorHotspot,
  SpHotspot = data$SpHotspot,
  OvHotspot = data$OvHotspot,
  maxDepth = data$maxDepth,
  OceanSR = data$OceanSR
)

test= predict(m1poly, newdat,type = "response",se.fit = TRUE)
newdat$prediction = test$fit
newdat$upr = pmin(1,test$fit+test$se.fit)
newdat$lwr = pmax(0,test$fit-test$se.fit)

test= predict(m2poly, subset(newdat,newdat$OceanSR!="(10 Aristazabal Banks Upwelling"),type = "response",se.fit = TRUE)
newdat$Corprediction[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = test$fit
newdat$Corupr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = pmin(1,test$fit+test$se.fit)
newdat$Corlwr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = pmax(0,test$fit-test$se.fit)

test= predict(m3poly, subset(newdat,newdat$OceanSR!="(10 Aristazabal Banks Upwelling"),type = "response",se.fit = TRUE)
newdat$Spprediction[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = test$fit
newdat$Spupr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = pmin(1,test$fit+test$se.fit)
newdat$Splwr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = pmax(0,test$fit-test$se.fit)

test= predict(m4poly, subset(newdat,newdat$OceanSR!="(10 Aristazabal Banks Upwelling"),type = "response",se.fit = TRUE)
newdat$Ovprediction[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = test$fit
newdat$Ovupr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = pmin(1,test$fit+test$se.fit)
newdat$Ovlwr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = pmax(0,test$fit-test$se.fit)


myPlot = function(newdat,OSR,ylabels=c("","","",""),xlabels=""){
  p1 = ggplot(data = subset(newdat,newdat$OceanSR==OSR),aes(x=maxDepth,y=SebHotspot))+geom_point()+
    geom_line(aes(y=prediction))+geom_ribbon(aes(ymin = lwr,ymax=upr),alpha=0.25)+
    theme_classic()+ theme(legend.position='none')+ylim(0,1)+ylab(ylabels[1])+xlab("")
  
  g1 = ggExtra::ggMarginal(p1, type = "histogram")
  
  p2 = ggplot(data = subset(newdat,newdat$OceanSR==OSR),aes(x=maxDepth,y=CorHotspot))+geom_point()+
    geom_line(aes(y=Corprediction))+geom_ribbon(aes(ymin=Corlwr, ymax=Corupr),alpha=0.25)+
    theme_classic()+ theme(legend.position='none')+ylim(0,1)+ylab(ylabels[2])+xlab("")
  
  g2 = ggExtra::ggMarginal(p2, type = "histogram")
  
  p3 = ggplot(data = subset(newdat,newdat$OceanSR==OSR),aes(x=maxDepth,y=SpHotspot))+geom_point()+
    geom_line(aes(y=Spprediction))+geom_ribbon(aes(ymin=Splwr, ymax=Spupr),alpha=0.25)+
    theme_classic()+ theme(legend.position='bottom')+ylim(0,1)+ylab(ylabels[3])+xlab("")
  
  g3 = ggExtra::ggMarginal(p3, type = "histogram")
  
  p4 = ggplot(data = subset(newdat,newdat$OceanSR==OSR),aes(x=maxDepth,y=OvHotspot))+geom_point()+
    geom_line(aes(y=Ovprediction))+geom_ribbon(aes(ymin=Ovlwr, ymax=Ovupr),alpha=0.25)+
    theme_classic()+ theme(legend.position='bottom')+ylim(0,1)+ylab(ylabels[4])+xlab(xlabels)
  
  g4 = ggExtra::ggMarginal(p4, type = "histogram")
  
  FullPlot = plot_grid(
    g1, g2, g3, g4, ncol = 1
  )
  
  return(FullPlot)
}

a = myPlot(newdat, "(10 Aristazabal Banks Upwelling",ylabels=c("Sebastidae","Corals","Sponges","Overall"),xlabels="Maximum depth (m)")
b = myPlot(newdat, "(11) Mainland Fjords",xlabels="Maximum depth (m)")
c = myPlot(newdat, "(13) Eastern Queen Charlotte Sound",xlabels="Maximum depth (m)")

x = plot_grid(a,b,c,labels=c("AU","MF","EQCS"),nrow=1)

jpeg("Figures/2021_Sep_biohotspot.jpeg", width = 11, height = 9, units = 'in', res = 300)
x
dev.off()

pdf("Figures/2021_Sep_biohotspot.pdf", width = 11, height = 11)
x
dev.off()
