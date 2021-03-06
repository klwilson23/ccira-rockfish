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
m1poly = glm(SebHotspot~poly(maxDepth,2,raw=TRUE)+OceanSR,data=data,family=binomial)
m1polint = glm(SebHotspot~poly(maxDepth,2,raw=TRUE)*OceanSR,data=data,family=binomial)

AIC(m1,m1poly,m1polint)
m2 = glm(CorHotspot~maxDepth+OceanSR,data=data2,family=binomial,na.action = na.omit)
m2poly = glm(CorHotspot~poly(maxDepth,2,raw=TRUE)+OceanSR,data=data2,family=binomial,na.action = na.omit)
m2polyint = glm(CorHotspot~poly(maxDepth,2,raw=TRUE)*OceanSR,data=data2,family=binomial,na.action = na.omit)
AIC(m2,m2poly,m2polyint)
m3 = glm(SpHotspot~maxDepth+OceanSR,data=data2,family=binomial,na.action = na.omit)
m3poly = glm(SpHotspot~poly(maxDepth,2,raw=TRUE)+OceanSR,data=data2,family=binomial,na.action = na.omit)
m3polyint = glm(CorHotspot~poly(maxDepth,2,raw=TRUE)*OceanSR,data=data2,family=binomial,na.action = na.omit)
AIC(m3,m3poly,m3polyint)

m4 = glm(OvHotspot~maxDepth+OceanSR,data=data2,family=binomial,na.action = na.omit)
m4poly = glm(OvHotspot~poly(maxDepth,2,raw=TRUE)+OceanSR,data=data2,family=binomial,na.action = na.omit)
m4polyint = glm(CorHotspot~poly(maxDepth,2,raw=TRUE)*OceanSR,data=data2,family=binomial,na.action = na.omit)
AIC(m4,m4poly,m4polyint)

rockfish_peak <- -coef(m1poly)[2]/(2*coef(m1poly)[3])
coral_peak <- -coef(m2poly)[2]/(2*coef(m2poly)[3])
sponge_peak <- -coef(m3poly)[2]/(2*coef(m3poly)[3])
overall_peak <- -coef(m4poly)[2]/(2*coef(m4poly)[3])

summary(m1poly)
summary(m2poly)
summary(m3poly)
summary(m4poly)
depth_SR <- aggregate(max_depth~OceanSR,data=data,quantile,probs=(0.9))
depth_SR2 <- aggregate(max_depth~OceanSR,data=data2,quantile,probs=(0.9))

rf <- predict(m1poly,type="r",newdata=data.frame("maxDepth"=depth_SR$max_depth[2],"OceanSR"=c("(11) Mainland Fjords","(10 Aristazabal Banks Upwelling","(13) Eastern Queen Charlotte Sound")))
cor <- predict(m2poly,type="r",newdata=data.frame("maxDepth"=depth_SR$max_depth[2],"OceanSR"=c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound")))
sp <- predict(m3poly,type="r",newdata=data.frame("maxDepth"=depth_SR$max_depth[2],"OceanSR"=c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound")))
ov <- predict(m4poly,type="r",newdata=data.frame("maxDepth"=depth_SR$max_depth[2],"OceanSR"=c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound")))
((1-rf)/rf)/((1-rf)/rf)[1]
((1-cor)/cor)/((1-cor)/cor)[1]
((1-sp)/sp)/((1-sp)/sp)[1]
((1-ov)/ov)/((1-ov)/ov)[1]


rf <- predict(m1poly,type="r",newdata=data.frame("maxDepth"=rockfish_peak,"OceanSR"=c("(11) Mainland Fjords","(10 Aristazabal Banks Upwelling","(13) Eastern Queen Charlotte Sound")))
cor <- predict(m2poly,type="r",newdata=data.frame("maxDepth"=coral_peak,"OceanSR"=c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound")))
sp <- predict(m3poly,type="r",newdata=data.frame("maxDepth"=sponge_peak,"OceanSR"=c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound")))
ov <- predict(m4poly,type="r",newdata=data.frame("maxDepth"=overall_peak,"OceanSR"=c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound")))
((rf-1)/rf)/((rf-1)/rf)[1]
((cor-1)/cor)/((cor-1)/cor)[1]
((sp-1)/sp)/((sp-1)/sp)[1]
((ov-1)/ov)/((ov-1)/ov)[1]

cor <- predict(m2poly,type="r",newdata=data.frame("maxDepth"=depth_SR$max_depth[3],"OceanSR"=c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound")))
sp <- predict(m3poly,type="r",newdata=data.frame("maxDepth"=depth_SR$max_depth[3],"OceanSR"=c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound")))
ov <- predict(m4poly,type="r",newdata=data.frame("maxDepth"=depth_SR$max_depth[3],"OceanSR"=c("(11) Mainland Fjords","(13) Eastern Queen Charlotte Sound")))
((cor-1)/cor)/((cor-1)/cor)[1]
((sp-1)/sp)/((sp-1)/sp)[1]
((ov-1)/ov)/((ov-1)/ov)[1]


depth_seq <- seq(0,400,by=25)
rf <- predict(m1poly,type="r",newdata=data.frame("maxDepth"=c(0.5*rockfish_peak,rockfish_peak,1.5*rockfish_peak),"OceanSR"=c("(11) Mainland Fjords")))
cor <- predict(m2poly,type="r",newdata=data.frame("maxDepth"=c(0.5*coral_peak,coral_peak,1.5*coral_peak),"OceanSR"=c("(11) Mainland Fjords")))
sp <- predict(m3poly,type="r",newdata=data.frame("maxDepth"=c(0.5*sponge_peak,sponge_peak,1.5*sponge_peak),"OceanSR"=c("(11) Mainland Fjords")))
ov <- predict(m4poly,type="r",newdata=data.frame("maxDepth"=c(0.5*overall_peak,overall_peak,1.5*overall_peak),"OceanSR"=c("(11) Mainland Fjords")))
((rf-1)/rf)/((rf-1)/rf)[1]
((cor-1)/cor)/((cor-1)/cor)[1]
((sp-1)/sp)/((sp-1)/sp)[1]
((ov-1)/ov)/((ov-1)/ov)[1]

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
