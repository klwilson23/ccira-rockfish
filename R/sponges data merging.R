normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

CategoryToPercent = function(x){
  ifelse(is.na(x),NA,((x*25)-ifelse(x<1,x*12.5,12.5))*ceiling(x/10))
}
spongeCountToArea <- 0.04

library(glmmTMB)
coral_spp <- c("calcigorgia","primnoa","stylaster","chrysopathes","paragorgia","swiftia")
coral_spp <- sort(coral_spp)
sponge_names <- c("aphrocallistidae","mycale","boot","farrea")

rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")

deep_vid <- read.csv("Data/DeepVideo v3.csv")
dive <- read.csv("Data/dive data sponge.csv")
dive_old <- read.csv("Data/DiveData v5.csv")
mid_vid <- read.csv("Data/MidDepthVideo v3.csv")
selectivity <- read.csv("Data/sponge gear selectivity.csv")

col_names <- colnames(deep_vid)[c(5,6,8,9)]
deep_vid <- deep_vid[deep_vid$binArea>=75 & deep_vid$binArea<=130,]
deep_vid <- deep_vid[,-c(7,10:23)]
deep_vid_wide <- tidyr::pivot_longer(deep_vid,cols=col_names)
deep_vid_wide$effort <- deep_vid_wide$binArea
deep_vid_wide <- deep_vid_wide[deep_vid_wide$effort>=75 & deep_vid_wide$effort<=130,] # subset based on bins between 75-130
deep_vid_wide$sample_size <- 1
deep_vid_wide$name <- "sponge"
deep_vid_agg <- aggregate(cbind(value,effort,AvgDepth)~Dive_Bin+name,data=deep_vid_wide,function(x){c(sum(x),mean(x,na.rm=TRUE))})
deep_vid$effort <- deep_vid$binArea
deep_vid$sample_size <- 1
effort <- aggregate(effort~Dive_Bin,data=deep_vid[deep_vid$effort>=75 & deep_vid$effort<=130,],FUN=sum,na.rm=TRUE)
sample_size <- aggregate(sample_size~Dive_Bin,data=deep_vid[deep_vid$effort>=75 & deep_vid$effort<=130,],FUN=sum)
depth <- aggregate(AvgDepth~Dive_Bin,data=deep_vid[deep_vid$effort>=75 & deep_vid$effort<=130,],FUN=max,na.rm=TRUE)

deep_vid_agg$X4km2grid <- deep_vid_wide$X4km2grid[match(deep_vid_agg$Dive_Bin,deep_vid_wide$Dive_Bin)]
deep_vid_agg$X1km2grid <- deep_vid_wide$X1km2grid[match(deep_vid_agg$Dive_Bin,deep_vid_wide$Dive_Bin)]
deep_vid_agg$UpperOceanSR <- deep_vid_wide$UpperOceanSR[match(deep_vid_agg$Dive_Bin,deep_vid_wide$Dive_Bin)]
deep_vid_agg$effort <- effort$effort[match(deep_vid_agg$Dive_Bin,effort$Dive_Bin)]
deep_vid_agg$sample_size <- sample_size$sample_size[match(deep_vid_agg$Dive_Bin,sample_size$Dive_Bin)]
deep_vid_agg$counts <- deep_vid_agg$value[,1]
deep_vid_agg$depth <- deep_vid_agg$AvgDepth[,2]
deep_vid_agg$max_depth <- depth$AvgDepth[match(deep_vid_agg$Dive_Bin,depth$Dive_Bin)]
deep_vid_agg$species <- unlist(lapply(strsplit(deep_vid_agg$name,'\\.'),function(x){x[1]}))
deep_vid_agg$species <- tolower(deep_vid_agg$species)
deep_vid_agg <- deep_vid_agg[,c("X4km2grid","X1km2grid","UpperOceanSR","Dive_Bin","name","counts","effort","sample_size","species","depth","max_depth")]
colnames(deep_vid_agg) <- c("PU_4Km_ID","PU_1Km_ID","UpperOceanSR","survey_id","name","counts","effort","sample_size","species","depth","max_depth")
#deep_vid_agg$counts <- (deep_vid_agg$counts*spongeCountToArea)/deep_vid_agg$effort

mid_vid$sponge <- CategoryToPercent(mid_vid$Glass.sponge.cover)/100
col_names <- "sponge"
mid_vid <- mid_vid[mid_vid$area.m2>=75 & mid_vid$area.m2<=130,]
mid_vid <- mid_vid[,-(3:28)]
mid_vid_wide <- tidyr::pivot_longer(mid_vid,cols=col_names)
mid_vid_wide$effort <- mid_vid_wide$area.m2
mid_vid_wide$sample_size <- 1
mid_vid_wide <- mid_vid_wide[mid_vid_wide$effort>=75 & mid_vid_wide$effort<=130,] # subset based on bins between 75-130
mid_vid_agg <- aggregate(cbind(value,depth)~binGrp+name,data=mid_vid_wide,function(x){c(sum(x),mean(x,na.rm=TRUE))})
mid_vid$effort <- mid_vid$area.m2
mid_vid$sample_size <- 1
effort <- aggregate(effort~binGrp,data=mid_vid[mid_vid$effort>=75 & mid_vid$effort<=130,],function(x){sum(x,na.rm=TRUE)})
sample_size <- aggregate(sample_size~binGrp,data=mid_vid_wide,FUN=function(x){length(unique(x))})
depth <- aggregate(depth~binGrp,data=mid_vid[mid_vid$effort>=75 & mid_vid$effort<=130,],FUN=max,na.rm=TRUE)

mid_vid_agg$PU_4Km_ID <- mid_vid_wide$PU_4Km_ID[match(mid_vid_agg$binGrp,mid_vid_wide$binGrp)]
mid_vid_agg$PU_1Km_ID <- mid_vid_wide$PU_1Km_ID[match(mid_vid_agg$binGrp,mid_vid_wide$binGrp)]
mid_vid_agg$UpperOceanSR <- mid_vid_wide$UpperOceanSR[match(mid_vid_agg$binGrp,mid_vid_wide$binGrp)]
mid_vid_agg$effort <- effort$effort[match(mid_vid_agg$binGrp,effort$binGrp)]
mid_vid_agg$sample_size <- sample_size$sample_size[match(mid_vid_agg$binGrp,sample_size$binGrp)]
mid_vid_agg$counts <- mid_vid_agg$value[,1]
mid_vid_agg$depth <- mid_vid_agg$depth[,2]
mid_vid_agg$max_depth <- depth$depth[match(mid_vid_agg$binGrp,depth$binGrp)]

mid_vid_agg$species <- unlist(lapply(strsplit(mid_vid_agg$name,'\\_'),function(x){x[1]}))
mid_vid_agg$species <- tolower(mid_vid_agg$species)
mid_vid_agg <- mid_vid_agg[,c("PU_4Km_ID","PU_1Km_ID","UpperOceanSR","binGrp","name","counts","effort","sample_size","species","depth","max_depth")]
colnames(mid_vid_agg) <- c("PU_4Km_ID","PU_1Km_ID","UpperOceanSR","survey_id","name","counts","effort","sample_size","species","depth","max_depth")

dive$sponge <- apply(dive[,c("GlassSponge1","GlassSponge2","GlassSponge3")],1,function(x)mean(CategoryToPercent(x),na.rm=TRUE))/100
dive$species <- "sponge"
dive$UpperOceanSR <- dive_old$UpperOceanSR[match(dive$PU_1Km_ID,dive_old$PU_1Km_ID)]
dive <- dive[!is.na(dive$PU_4Km_ID),]
dive <- dive[!is.na(dive$Transect.Number),]
dive$transects <- paste(dive$DescriptiveID,dive$Transect.Number,sep="_")
effort <- rep(0,length(dive$transects))
effort[match(unique(dive$transects),dive$transects)] <- 1
dive$effort <- effort
dive$sample_size <- 1
dive_agg <- aggregate(cbind(sponge,Depth)~transects+species,data=dive,function(x){c(sum(x),mean(x,na.rm=TRUE))})
effort <- aggregate(sample_size~transects,data=dive,function(x){length(unique(x))})
depth <- aggregate(Depth~transects,data=dive,FUN=max,na.rm=TRUE)
dive_agg$effort <- effort$sample_size[match(dive_agg$transects,effort$transects)]
dive_agg$PU_4Km_ID <- dive$PU_4Km_ID[match(dive_agg$transects,dive$transects)]
dive_agg$PU_1Km_ID <- dive$PU_1Km_ID[match(dive_agg$transects,dive$transects)]
dive_agg$UpperOceanSR <- dive$UpperOceanSR[match(dive_agg$transects,dive$transects)]
dive_agg$max_depth <- depth$Depth[match(dive_agg$transects,depth$transects)]

dive_counts <- dive_agg[,c("PU_4Km_ID","PU_1Km_ID","UpperOceanSR","transects","species")]
colnames(dive_counts) <- c("PU_4Km_ID","PU_1Km_ID","UpperOceanSR","survey_id","species")
dive_counts$counts <- dive_agg$sponge[,1]
dive_counts$effort <- dive_agg$effort
dive_counts$sample_size <- dive_agg$effort
dive_counts$depth <- dive_agg$Depth[,2]
dive_counts$max_depth <- dive_agg$max_depth

dive_counts$effort <- dive_counts$effort*120
dive_counts$gear <- "dive"
deep_vid_agg$gear <- "deep_video"
mid_vid_agg$gear <- "mid_video"
merged_data <- merge(merge(deep_vid_agg[,-5],mid_vid_agg[,-5],all=TRUE),dive_counts,all=TRUE)
site_ids <- sort(unique(merged_data$PU_4Km_ID))
surv_ids <- sort(unique(merged_data$survey_id))
merged_data$PU_4Km_ID <- factor(merged_data$PU_4Km_ID,levels=site_ids)
merged_data$survey_id <- factor(merged_data$survey_id,levels=surv_ids)

merged_data$counts[merged_data$counts==0 & merged_data$depth<15] <- NA

merged_data <- merged_data[!is.na(merged_data$counts),]
write.csv(merged_data,"Data/sponge counts PU4km.csv")