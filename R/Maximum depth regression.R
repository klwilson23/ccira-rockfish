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
data$OceanSR <- data$UpperOceanSR
data$maxDepth <- data$max_depth
data$OceanSR = factor(data$OceanSR,levels = c("(11) Mainland Fjords","(10 Aristazabal Banks Upwelling","(13) Eastern Queen Charlotte Sound"))
data2 = subset(data, data$OceanSR!="(10 Aristazabal Banks Upwelling")
data2$OceanSR = factor(data2$OceanSR,levels = c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound"))

#I re-ordered the factors because I wanted to compare MF to to other two.
#for corals and spongens I had to drop aristazabal entirely - see data2

#binomial regression, but more complicated ####
m1 = glm(SebHotspot~maxDepth+OceanSR,data=data,family=binomial)
m2 = glm(CorHotspot~maxDepth+OceanSR,data=data2,family=binomial,na.action = na.omit)
m3 = glm(SpHotspot~maxDepth+OceanSR,data=data2,family=binomial,na.action = na.omit)

summary(m1)
summary(m2)
summary(m3)

# plots
newdat = data.frame(
  puid = data$PU_1Km_ID,
  SebHotspot = data$SebHotspot,
  CorHotspot = data$CorHotspot,
  SpHotspot = data$SpHotspot,
  maxDepth = data$maxDepth,
  OceanSR = data$OceanSR
)

test= predict(m1, newdat,type = "response",se.fit = TRUE)
newdat$prediction = test$fit
newdat$upr = newdat$prediction+test$se.fit
newdat$lwr = newdat$prediction-test$se.fit

test= predict(m2, subset(newdat,newdat$OceanSR!="(10 Aristazabal Banks Upwelling"),type = "response",se.fit = TRUE)
newdat$Corprediction[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = test$fit
newdat$Corupr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = test$fit+test$se.fit
newdat$Corlwr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = test$fit-test$se.fit

test= predict(m3, subset(newdat,newdat$OceanSR!="(10 Aristazabal Banks Upwelling"),type = "response",se.fit = TRUE)
newdat$Spprediction[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = test$fit
newdat$Spupr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = test$fit+test$se.fit
newdat$Splwr[which(newdat$OceanSR!="(10 Aristazabal Banks Upwelling")] = test$fit-test$se.fit


myPlot = function(newdat,OSR,ylabels=c("","",""),xlabels=""){
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
    theme_classic()+ theme(legend.position='bottom')+ylim(0,1)+ylab(ylabels[3])+xlab(xlabels)
  
  g3 = ggExtra::ggMarginal(p3, type = "histogram")
  
  FullPlot = plot_grid(
    g1, g2, g3, ncol = 1
  )
  
  return(FullPlot)
}

a = myPlot(newdat, "(10 Aristazabal Banks Upwelling",ylabels=c("Sebastidae","Corals","Sponges"),xlabels="Maximum depth (m)")
b = myPlot(newdat, "(11) Mainland Fjords",xlabels="Maximum depth (m)")
c = myPlot(newdat, "(13) Eastern Queen Charlotte Sound",xlabels="Maximum depth (m)")

x = plot_grid(a,b,c,labels=c("AU","MF","EQCS"),nrow=1)

jpeg("Figures/2021_Sep_biohotspot.jpeg", width = 11, height = 9, units = 'in', res = 300)
x
dev.off()

pdf("Figures/2021_Sep_biohotspot.pdf", width = 11, height = 9)
x
dev.off()
