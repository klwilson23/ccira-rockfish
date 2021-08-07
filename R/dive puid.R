#Organize Dive Data

#load libraries ####
library(plyr)

#load data ####
data = read.csv("Data/DiveData.csv")
data$transectID = paste(data$Descriptive.ID,data$Transect.Transect.Number,sep="_")
sebastidae = read.csv("Data/Sebastidae list.csv")
#I need to make sure the species names match between these two files. 
setdiff(sebastidae[,1],levels(data$Species))
setdiff(levels(data$Species),sebastidae[,1])

# fix the species names to match sebastidae
data$Species = as.character(data$Species)
data$Species[which(data$Species%in%c("dusky","dark"))]="Dusky-Dark"
data$Species=as.factor(data$Species)

#re-organize glass sponge and coral data #
data$avgGlassSponge = apply(data[,c("GlassSponge1","GlassSponge2","GlassSponge3")],1,function(x)mean(x, na.rm=T))

GSna = which(is.na(data$avgGlassSponge))
data$BootCloud = as.numeric(data$Boot.Sponge==1|data$Cloud.Sponge==1)
data$avgGlassSponge[GSna] = data$BootCloud[GSna]


# summarize at the transect-level ####
transectData = ddply(data,"transectID",summarize,
                     Descriptive.ID = Descriptive.ID[1],
                     Date = Date[1],
                     PU_4Km_ID = PU_4Km_ID[1],
                     #RCA_DistKM = RCA_DistKM[1],
                     #RCA_Name = RCA_Name[1],
                     # mean.Lat = mean(Lat),
                     # mean.Long = mean(Long),
                     mean.Depth = mean(Depth),
                     min.Depth = min(Depth),
                     max.Depth = max(Depth),
                     
                     #sponge data
                     Sponge = mean(avgGlassSponge,na.rm=T),
                     
                     
                     #number of species
                     #Only sebastes, not lingcod
                     #all sizes
                     num.sp.Sebastes = length(unique(Species[which(Species%in%sebastidae[,1])])),
                     
                     #individual species 10P
                     Bocaccio10P = sum(Multiplier[which(Species=="bocaccio"&TL>=10)]),
                     Canary10P = sum(Multiplier[which(Species=="canary"&TL>=10)]),
                     China10P = sum(Multiplier[which(Species=="china"&TL>=10)]),
                     Quillback10P = sum(Multiplier[which(Species=="quillback"&TL>=10)]),
                     Silvergrey10P = sum(Multiplier[which(Species=="silvergray"&TL>=10)]),
                     Tiger10P = sum(Multiplier[which(Species=="tiger"&TL>=10)]),
                     Yelloweye10P = sum(Multiplier[which(Species=="yelloweye"&TL>=10)]),
                     Black10P = sum(Multiplier[which(Species=="black"&TL>=10)]),
                     Brown10P = sum(Multiplier[which(Species=="brown"&TL>=10)]),
                     Copper10P = sum(Multiplier[which(Species=="copper"&TL>=10)]),
                     Deacon10P = sum(Multiplier[which(Species=="deacon"&TL>=10)]),
                     DuskyDark10P = sum(Multiplier[which(Species=="Dusky-Dark"&TL>=10)]),
                     Greenstriped10P = sum(Multiplier[which(Species=="greenstriped"&TL>=10)]),
                     PugetSound10P = sum(Multiplier[which(Species=="puget sound"&TL>=10)]),
                     Pygmy10P = sum(Multiplier[which(Species=="pygmy"&TL>=10)]),
                     Redbanded10P = sum(Multiplier[which(Species=="redbanded"&TL>=10)]),
                     Redstripe10P = sum(Multiplier[which(Species=="redstripe"&TL>=10)]),
                     Rosethorn10P = sum(Multiplier[which(Species=="rosethorn"&TL>=10)]),
                     REBS10P = sum(Multiplier[which(Species=="REBS"&TL>=10)]),
                     Sharpchin10P = sum(Multiplier[which(Species=="sharpchin"&TL>=10)]),
                     Shortbelly10P = sum(Multiplier[which(Species=="shortbelly"&TL>=10)]),           
                     Shortraker10P = sum(Multiplier[which(Species=="shortraker"&TL>=10)]),
                     Splitnose10P = sum(Multiplier[which(Species=="splitnose"&TL>=10)]),
                     Vermillion10P = sum(Multiplier[which(Species=="vermillion"&TL>=10)]),
                     Widow10P = sum(Multiplier[which(Species=="widow"&TL>=10)]),
                     ShortspineThornyhead10P = sum(Multiplier[which(Species=="shortspine thornyhead"&TL>=10)]),
                     Yellowtail10P = sum(Multiplier[which(Species=="yellowtail"&TL>=10)])
)


# summarize habitat data by PUID ####

PUID.data = ddply(transectData,"PU_4Km_ID",summarize,
                  num.transects = length(unique(transectID)),
                  area.surveyed = length(unique(transectID))*120,
                  #mean.Lat = mean(mean.Lat),
                  #mean.Long = mean(mean.Long),
                  mean.Depth = mean(mean.Depth),
                  min.Depth = min(min.Depth),
                  max.Depth = max(max.Depth),
                  #mean.Complex = mean(mean.Complexity),
                  
                  #sponges
                  sponge = sum(Sponge)/length(unique(transectID)),
                  
                  #richness
                  sp.density = mean(num.sp.Sebastes))

index = which(names(transectData)=="Bocaccio10P")

a = transectData[,index:length(transectData)]
x = aggregate(a,list(transectData$PU_4Km_ID),sum) 
PUID.data = merge(PUID.data,x,by.x = "PU_4Km_ID",by.y = "Group.1")

PUID.data[,which(names(PUID.data)=="Bocaccio10P"):length(PUID.data)] = 
  PUID.data[,which(names(PUID.data)=="Bocaccio10P"):length(PUID.data)]/(PUID.data$num.transects*120)


#normalize data ####
myFun = function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
b = PUID.data[,c(which(names(PUID.data)=="sp.density"),which(names(PUID.data)=="sponge"),which(names(PUID.data)=="Bocaccio10P"):length(PUID.data))]
row.names(b)=PUID.data$PU_4Km_ID
y = apply(b,2,myFun)
colnames(y)=paste(colnames(y),".z",sep="")
y=as.data.frame(y)
y=tibble::rownames_to_column(y,"PU_4Km_ID")

PUID.data.norm = merge(PUID.data,y,by.x="PU_4Km_ID",by.y="PU_4Km_ID")


#we also wanted to know the 90% min and max depths across the whole data set
PUID.data.norm$min90 = quantile(data$Depth,0.1,na.rm=T)
PUID.data.norm$max90 = quantile(data$Depth,0.9,na.rm=T)
