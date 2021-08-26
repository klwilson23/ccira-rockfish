normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

library(glmmTMB)
coral_spp <- c("calcigorgia","primnoa","stylaster","chrysopathes","paragorgia","swiftia")
coral_spp <- sort(coral_spp)
sponge_names <- c("aphrocallistidae","mycale","boot","farrea")

rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")

deep_vid <- read.csv("Data/DeepVideo v2.csv")
dive <- read.csv("Data/DiveData v4.csv")
mid_vid <- read.csv("Data/MidDepthVideo v2.csv")
dive_coral <- read.csv("Data/StructuralInverts.csv")
selectivity <- read.csv("Data/coral gear selectivity.csv")

col_names <- colnames(deep_vid)[-c(1:4,24)]
deep_vid <- deep_vid[deep_vid$binArea>=75 & deep_vid$binArea<=130,]
deep_vid_wide <- tidyr::pivot_longer(deep_vid,cols=col_names)
deep_vid_wide$effort <- deep_vid_wide$binArea
deep_vid_wide <- deep_vid_wide[deep_vid_wide$effort>=75 & deep_vid_wide$effort<=130,] # subset based on bins between 75-130
deep_vid_wide$sample_size <- 1
deep_vid_agg <- aggregate(cbind(value,effort,AvgDepth)~Dive_Bin+name,data=deep_vid_wide,function(x){c(sum(x,na.rm=TRUE),mean(x,na.rm=TRUE))})
deep_vid$effort <- deep_vid$binArea
deep_vid$sample_size <- 1
effort <- aggregate(effort~Dive_Bin,data=deep_vid[deep_vid$effort>=75 & deep_vid$effort<=130,],FUN=sum,na.rm=TRUE)
sample_size <- aggregate(sample_size~Dive_Bin,data=deep_vid[deep_vid$effort>=75 & deep_vid$effort<=130,],FUN=sum)

deep_vid_agg$X4km2grid <- deep_vid_wide$X4km2grid[match(deep_vid_agg$Dive_Bin,deep_vid_wide$Dive_Bin)]
deep_vid_agg$X1km2grid <- deep_vid_wide$X1km2grid[match(deep_vid_agg$Dive_Bin,deep_vid_wide$Dive_Bin)]
deep_vid_agg$effort <- effort$effort[match(deep_vid_agg$Dive_Bin,effort$Dive_Bin)]
deep_vid_agg$sample_size <- sample_size$sample_size[match(deep_vid_agg$Dive_Bin,sample_size$Dive_Bin)]
deep_vid_agg$counts <- deep_vid_agg$value[,1]
deep_vid_agg$depth <- deep_vid_agg$AvgDepth[,2]
deep_vid_agg$species <- unlist(lapply(strsplit(deep_vid_agg$name,'\\.'),function(x){x[1]}))
deep_vid_agg$species <- tolower(deep_vid_agg$species)
deep_vid_agg <- deep_vid_agg[,c("X4km2grid","X1km2grid","Dive_Bin","name","counts","effort","sample_size","species","depth")]
colnames(deep_vid_agg) <- c("PU_4Km_ID","PU_1Km_ID","survey_id","name","counts","effort","sample_size","species","depth")

col_names <- colnames(mid_vid)[-c(1,2,29,30,33,34)]
mid_vid <- mid_vid[mid_vid$area.m2>=75 & mid_vid$area.m2<=130,]
mid_vid_wide <- tidyr::pivot_longer(mid_vid,cols=col_names)
mid_vid_wide$effort <- mid_vid_wide$area.m2
mid_vid_wide$sample_size <- 1
mid_vid_wide <- mid_vid_wide[mid_vid_wide$effort>=75 & mid_vid_wide$effort<=130,] # subset based on bins between 75-130
mid_vid_agg <- aggregate(cbind(value,depth)~binGrp+name,data=mid_vid_wide,function(x){c(sum(x,na.rm=TRUE),mean(x,na.rm=TRUE))})
mid_vid$effort <- mid_vid$area.m2
mid_vid$sample_size <- 1
effort <- aggregate(effort~binGrp,data=mid_vid[mid_vid$effort>=75 & mid_vid$effort<=130,],function(x){sum(x,na.rm=TRUE)})
sample_size <- aggregate(sample_size~binGrp,data=mid_vid_wide,FUN=function(x){length(unique(x))})

mid_vid_agg$PU_4Km_ID <- mid_vid_wide$PU_4Km_ID[match(mid_vid_agg$binGrp,mid_vid_wide$binGrp)]
mid_vid_agg$PU_1Km_ID <- mid_vid_wide$PU_1Km_ID[match(mid_vid_agg$binGrp,mid_vid_wide$binGrp)]
mid_vid_agg$effort <- effort$effort[match(mid_vid_agg$binGrp,effort$binGrp)]
mid_vid_agg$sample_size <- sample_size$sample_size[match(mid_vid_agg$binGrp,sample_size$binGrp)]
mid_vid_agg$counts <- mid_vid_agg$value[,1]
mid_vid_agg$depth <- mid_vid_agg$depth[,2]
mid_vid_agg$species <- unlist(lapply(strsplit(mid_vid_agg$name,'\\_'),function(x){x[1]}))
mid_vid_agg$species <- tolower(mid_vid_agg$species)
mid_vid_agg <- mid_vid_agg[,c("PU_4Km_ID","PU_1Km_ID","binGrp","name","counts","effort","sample_size","species","depth")]
colnames(mid_vid_agg) <- c("PU_4Km_ID","PU_1Km_ID","survey_id","name","counts","effort","sample_size","species","depth")

dive_old <- dive
dive_coral$survey_id <- paste(dive_coral$DescriptiveSurveyID,dive_coral$TransectNumber,sep="_")
dive_coral$counts <- dive_coral$Multiplier
dive_coral$species <- tolower(dive_coral$Species)
dive_coral$species[which(dive_coral$species%in%c("primnoa pacifica"))]="primnoa"
dive <- dive[!is.na(dive$PU_4Km_ID),]
dive <- dive[!is.na(dive$Transect.Transect.Number),]
dive$counts <- dive$Multiplier
dive$transects <- paste(dive$Descriptive.ID,dive$Transect.Transect.Number,sep="_")
effort <- rep(0,length(dive$transects))
effort[match(unique(dive$transects),dive$transects)] <- 1
dive$effort <- effort
dive$Species[which(dive$Species%in%c("dusky","dark"))]="dusky-dark"
dive$Species[which(dive$Species%in%c("rougheye"))]="blackspotted"

dive <- dive[!dive$Species=="black & yellowtail",]
dive$species <- dive$Species
dive$counts[dive$TL>=10] <- 0 # set the small fish counts to 0 for the purposes of summing at the transect level
dive$sample_size <- 1
dive_agg <- aggregate(cbind(counts,Depth)~transects+species,data=dive,function(x){c(sum(x,na.rm=TRUE),mean(x,na.rm=TRUE))})
effort <- aggregate(sample_size~transects,data=dive,function(x){length(unique(x))})
dive_agg$effort <- effort$sample_size[match(dive_agg$transects,effort$transects)]

dive_agg$PU_4Km_ID <- dive$PU_4Km_ID[match(dive_agg$transects,dive$transects)]
dive_agg$PU_1Km_ID <- dive$PU_1Km_ID[match(dive_agg$transects,dive$transects)]

dive_counts <- dive_agg[,c("PU_4Km_ID","PU_1Km_ID","transects","species")]
colnames(dive_counts) <- c("PU_4Km_ID","PU_1Km_ID","survey_id","species")
dive_counts$counts <- dive_agg$count[,1]
dive_counts$effort <- dive_agg$effort
dive_counts$sample_size <- dive_agg$effort
dive_counts$depth <- dive_agg$Depth[,2]
dive_counts$effort <- dive_counts$effort*120
dive_counts$gear <- "dive"
deep_vid_agg$gear <- "deep_video"
mid_vid_agg$gear <- "mid_video"
merged_data <- merge(merge(deep_vid_agg[,-4],mid_vid_agg[,-4],all=TRUE),dive_counts,all=TRUE)
merged_data$species[merged_data$species=="puget.sound"]="puget sound"
merged_data$species[merged_data$species=="puget"]="puget sound"
merged_data$species[merged_data$species=="red stripe"]="redstriped"
merged_data$species[merged_data$species=="paragorgia.pac"]="paragorgia"
merged_data$species[merged_data$species=="greenstriped"]="greenstripe"
merged_data$species[merged_data$species=="redstriped"]="redstripe"
merged_data$species[merged_data$species=="dusky"]="dusky-dark"

site_ids <- sort(unique(merged_data$PU_4Km_ID))
surv_ids <- sort(unique(merged_data$survey_id))
merged_data$PU_4Km_ID <- factor(merged_data$PU_4Km_ID,levels=site_ids)
merged_data$survey_id <- factor(merged_data$survey_id,levels=surv_ids)

new_data <- expand.grid("survey_id"=surv_ids,"species"=unique(c(merged_data$species,coral_spp,sponge_names))) # create a combination of all 4x4 sites, all species, and all gear types
new_data$gear <- merged_data$gear[match(new_data$survey_id,merged_data$survey_id)]
new_data$effort <- merged_data$effort[match(paste(new_data$survey_id,new_data$gear),paste(merged_data$survey_id,merged_data$gear),nomatch=NA)]
new_data$sample_size <- merged_data$sample_size[match(paste(new_data$survey_id,new_data$gear),paste(merged_data$survey_id,merged_data$gear),nomatch=NA)]
new_data$depth <- merged_data$depth[match(paste(new_data$survey_id,new_data$gear),paste(merged_data$survey_id,merged_data$gear),nomatch=NA)]
new_data$counts <- 0
new_data$counts[match(paste(merged_data$survey_id,merged_data$gear,merged_data$species),paste(new_data$survey_id,new_data$gear,new_data$species))] <- merged_data$counts
new_df <- new_data
new_df$PU_4Km_ID <- merged_data$PU_4Km_ID[match(new_df$survey_id,merged_data$survey_id)]
new_df$PU_1Km_ID <- merged_data$PU_1Km_ID[match(new_df$survey_id,merged_data$survey_id)]

new_dat <- NULL
for(i in 1:length(unique(new_df$survey_id)))
{
  #i <- which(unique(new_df$survey_id)=="wds-07-03-2021_2")
  sub_dat <- new_df[new_df$survey_id==unique(new_df$survey_id)[i],]
  sub_dive <- dive_coral[dive_coral$survey_id==unique(new_df$survey_id)[i],]
  for(k in 1:length(coral_spp))
  {
    sub_dat2 <- sub_dat[sub_dat$species==coral_spp[k],]
    sub_dive_coral <- sub_dive[sub_dive$species==coral_spp[k],]
    tru_dat <- sub_dat2
    if(any(sub_dat2$species==sub_dive_coral$species))
    {
      tru_dat$counts <- sub_dive_coral$counts
    }
    if(any(selectivity$Species==coral_spp[k]))
    {
      sub_sel <- selectivity[selectivity$Species==coral_spp[k],]
      zeros <- sub_sel[,c(1:3,which(colnames(sub_sel)==tru_dat$gear))]
      omit <- ifelse(zeros[,4]==0 & tru_dat$counts==0,1,0)
      if(omit==1)
      {
        tru_dat$counts <- NA
      }else{
        tru_dat$counts <- tru_dat$counts
      }
    }
    new_dat <- rbind(new_dat,tru_dat)
  }
}
full_dat <- new_dat[!is.na(new_dat$counts),]
nrow(full_dat[full_dat$species=="calcigorgia" & full_dat$counts>0,])

write.csv(full_dat,"Data/coral counts PU4km.csv")
