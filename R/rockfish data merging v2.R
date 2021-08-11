normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

library(glmmTMB)
rockfish_spp <- c("black","blackspotted","bocaccio","brown","canary","china","copper","deacon","dusky-dark","greenstripe","puget sound","pygmy","quillback","redbanded","redstripe","rosethorn","sebastolobus","sharpchin","shortbelly","shortraker","silvergray","splitnose","stripetail","tiger","vermillion","widow","yelloweye","yellowtail")

hook_line <- read.csv("Data/All_HLData v2.csv")
deep_vid <- read.csv("Data/DeepVideo v2.csv")
dive <- read.csv("Data/DiveData v4.csv")
mid_vid <- read.csv("Data/MidDepthVideo v2.csv")
selectivity <- read.csv("Data/rockfish gear selectivity v2.csv")

hook_line <- hook_line[hook_line$Survey.Type=="hl-syst",]
hook_line <- hook_line[complete.cases(hook_line[,-which(colnames(hook_line)=="TL")]),]
hook_line$count <- 1
hook_line$effort <- hook_line$Bottom.Time..min.
hook_line$sample_size <- 1
#hook_line$species[which(hook_line$species%in%c("REBS"))]="blackspotted"

hook_agg <- aggregate(cbind(count,depth)~survey_id+species,data=hook_line,function(x){c(sum(x,na.rm=TRUE),mean(x,na.rm=TRUE))})
effort <- aggregate(cbind(effort,sample_size)~survey_id,data=hook_line,function(x){sum(x,na.rm=TRUE)})

hook_agg$PU_4Km_ID <- hook_line$PU_4Km_ID[match(hook_agg$survey_id,hook_line$survey_id)]
hook_agg$effort <- effort$effort[match(hook_agg$survey_id,effort$survey_id)]
hook_agg$sample_size <- effort$sample_size[match(hook_agg$survey_id,effort$survey_id)]
hook_agg <- hook_agg[complete.cases(hook_agg),]
hook_counts <- hook_agg[,c("PU_4Km_ID","survey_id","species")]
hook_counts$counts <- hook_agg$count[,1]
hook_counts$depth <- hook_agg$depth[,2]
hook_counts$effort <- hook_agg$effort
hook_counts$sample_size <- hook_agg$sample_size
zero_counts <- hook_counts[!hook_counts$species%in%rockfish_spp,]
zero_count_df <- expand.grid("survey_id"=unique(zero_counts$survey_id),"species"=rockfish_spp)
zero_count_df$effort <- zero_counts$effort[match(zero_count_df$survey_id,zero_counts$survey_id)]
zero_count_df$sample_size <- zero_counts$sample_size[match(zero_count_df$survey_id,zero_counts$survey_id)]
zero_count_df$depth <- zero_counts$depth[match(zero_count_df$survey_id,zero_counts$survey_id)]
zero_count_df$counts <- 0
zero_count_df$PU_4Km_ID <- hook_agg$PU_4Km_ID[match(zero_count_df$survey_id,hook_agg$survey_id)]
hook_counts <- rbind(hook_counts[hook_counts$species%in%rockfish_spp,],zero_count_df)
#hook_counts <- hook_counts[hook_counts$species!="none",]

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
effort$X4km2grid <- deep_vid_wide$X4km2grid[match(effort$Dive_Bin,deep_vid_wide$Dive_Bin)]
sample_size$X4km2grid <- deep_vid_wide$X4km2grid[match(sample_size$Dive_Bin,deep_vid_wide$Dive_Bin)]
deep_vid_agg$effort <- effort$effort[match(deep_vid_agg$Dive_Bin,effort$Dive_Bin)]
deep_vid_agg$sample_size <- sample_size$sample_size[match(deep_vid_agg$Dive_Bin,sample_size$Dive_Bin)]
deep_vid_agg$counts <- deep_vid_agg$value[,1]
deep_vid_agg$depth <- deep_vid_agg$AvgDepth[,2]
deep_vid_agg$species <- unlist(lapply(strsplit(deep_vid_agg$name,'\\.'),function(x){x[1]}))
deep_vid_agg$species <- tolower(deep_vid_agg$species)
deep_vid_agg <- deep_vid_agg[,c("X4km2grid","Dive_Bin","name","counts","effort","sample_size","species","depth")]
colnames(deep_vid_agg) <- c("PU_4Km_ID","survey_id","name","counts","effort","sample_size","species","depth")

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
effort$PU_4Km_ID <- mid_vid_wide$PU_4Km_ID[match(effort$binGrp,mid_vid_wide$binGrp)]
sample_size$PU_4Km_ID <- mid_vid_wide$PU_4Km_ID[match(sample_size$binGrp,mid_vid_wide$binGrp)]
mid_vid_agg$effort <- effort$effort[match(mid_vid_agg$binGrp,effort$binGrp)]
mid_vid_agg$sample_size <- sample_size$sample_size[match(mid_vid_agg$binGrp,sample_size$binGrp)]
mid_vid_agg$counts <- mid_vid_agg$value[,1]
mid_vid_agg$depth <- mid_vid_agg$depth[,2]
mid_vid_agg$species <- unlist(lapply(strsplit(mid_vid_agg$name,'\\_'),function(x){x[1]}))
mid_vid_agg$species <- tolower(mid_vid_agg$species)
mid_vid_agg <- mid_vid_agg[,c("PU_4Km_ID","binGrp","name","counts","effort","sample_size","species","depth")]
colnames(mid_vid_agg) <- c("PU_4Km_ID","survey_id","name","counts","effort","sample_size","species","depth")

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
dive <- dive[dive$TL>=10,]
dive$sample_size <- 1
dive_agg <- aggregate(cbind(counts,Depth)~transects+species,data=dive,function(x){c(sum(x,na.rm=TRUE),mean(x,na.rm=TRUE))})
effort <- aggregate(sample_size~transects,data=dive,function(x){length(unique(x))})
dive_agg$effort <- effort$sample_size[match(dive_agg$transects,effort$transects)]

dive_agg$PU_4Km_ID <- dive$PU_4Km_ID[match(dive_agg$transects,dive$transects)]
dive_counts <- dive_agg[,c("PU_4Km_ID","transects","species")]
colnames(dive_counts) <- c("PU_4Km_ID","survey_id","species")
dive_counts$counts <- dive_agg$count[,1]
dive_counts$effort <- dive_agg$effort
dive_counts$sample_size <- dive_agg$effort
dive_counts$depth <- dive_agg$Depth[,2]
dive_counts$effort <- dive_counts$effort*120
dive_counts$gear <- "dive"
deep_vid_agg$gear <- "deep_video"
hook_counts$gear <- "hook_line"
mid_vid_agg$gear <- "mid_video"
merged_data <- merge(merge(merge(deep_vid_agg[,-3],hook_counts,all=TRUE),mid_vid_agg[,-3],all=TRUE),dive_counts,all=TRUE)
merged_data$species[merged_data$species=="puget.sound"]="puget sound"
merged_data$species[merged_data$species=="puget"]="puget sound"
merged_data$species[merged_data$species=="red stripe"]="redstriped"
merged_data$species[merged_data$species=="paragorgia.pac"]="paragorgia"
merged_data$species[merged_data$species=="greenstriped"]="greenstripe"
merged_data$species[merged_data$species=="redstriped"]="redstripe"
merged_data$species[merged_data$species=="dusky"]="dusky-dark"
sponge_names <- c("aphrocallistidae","mycale","chrysopathes","paragorgia","primnoa","stylaster","swiftia","boot","farrea")
site_ids <- sort(unique(merged_data$PU_4Km_ID))
surv_ids <- sort(unique(merged_data$survey_id))
#merged_data <- merged_data[!merged_data$species%in%sponge_names,]
merged_data$PU_4Km_ID <- factor(merged_data$PU_4Km_ID,levels=site_ids)
merged_data$survey_id <- factor(merged_data$survey_id,levels=surv_ids)
#merged_data <- merged_data[merged_data$species!="none",]
#merged_data <- merged_data[merged_data$species%in%rockfish_spp,]
new_data <- expand.grid("survey_id"=surv_ids,"species"=unique(merged_data$species)) # create a combination of all 4x4 sites, all species, and all gear types
new_data$gear <- merged_data$gear[match(new_data$survey_id,merged_data$survey_id)]
new_data$effort <- merged_data$effort[match(paste(new_data$survey_id,new_data$gear),paste(merged_data$survey_id,merged_data$gear),nomatch=NA)]
new_data$sample_size <- merged_data$sample_size[match(paste(new_data$survey_id,new_data$gear),paste(merged_data$survey_id,merged_data$gear),nomatch=NA)]
new_data$depth <- merged_data$depth[match(paste(new_data$survey_id,new_data$gear),paste(merged_data$survey_id,merged_data$gear),nomatch=NA)]
new_data$counts <- 0
new_data$counts[match(paste(merged_data$survey_id,merged_data$gear,merged_data$species),paste(new_data$survey_id,new_data$gear,new_data$species))] <- merged_data$counts
#new_df <- new_data[which(match(paste(new_data$survey_id,new_data$gear),paste(merged_data$survey_id,merged_data$gear),nomatch=0)>0),]
new_df <- new_data
#new_df$depth[match(paste(merged_data$survey_id,merged_data$gear,merged_data$species),paste(new_df$survey_id,new_df$gear,new_df$species))] <- merged_data$depth

new_df$PU_4Km_ID <- merged_data$PU_4Km_ID[match(new_df$survey_id,merged_data$survey_id)]

new_dat <- NULL
for(i in 1:length(unique(new_df$survey_id)))
{
  sub_dat <- new_df[new_df$survey_id==unique(new_df$survey_id)[i],]
  for(k in 1:length(rockfish_spp))
  {
    sub_dat2 <- sub_dat[sub_dat$species==rockfish_spp[k],]
    
    if(any(selectivity$Common.name==rockfish_spp[k]))
    {
      sub_sel <- selectivity[selectivity$Common.name==rockfish_spp[k],]
      zeros <- sub_sel[,c(1:3,which(colnames(sub_sel)==sub_dat2$gear))]
      omit <- ifelse(zeros[,4]==0 & sub_dat2$counts==0,1,0)
      if(omit==1)
      {
        tru_dat <- sub_dat2
        tru_dat$counts <- NA
      }else{
        tru_dat <- sub_dat2
      }
    }else{
      sub_sel <- selectivity[selectivity$Common.name=="REBS",]
      zeros <- sub_sel[,c(1:3,which(colnames(sub_sel)==sub_dat2$gear))]
      omit <- ifelse(zeros[,4]==0 & sub_dat2$counts==0,1,0)
      if(omit==1)
      {
        tru_dat <- sub_dat2
        tru_dat$counts <- NA
      }else{
        tru_dat <- sub_dat2
      }
    }
    new_dat <- rbind(new_dat,tru_dat)
  }
}

full_dat <- new_dat[!is.na(new_dat$counts),]
write.csv(full_dat,"Data/Rockfish counts PU4km v2.csv")

merged_data[merged_data$PU_4Km_ID=="6265",]
full_dat[full_dat$PU_4Km_ID=="6265",]
unique(merged_data$PU_4Km_ID)
merged_data[merged_data$PU_4Km_ID=="8944",]
full_dat[full_dat$PU_4Km_ID=="8944",]
mid_vid[mid_vid$PU_4Km_ID=="8944",]
sum(mid_vid[mid_vid$PU_4Km_ID=="8944","area.m2"])
full_dat[full_dat$PU_4Km_ID=="8473",]
full_dat[full_dat$PU_4Km_ID=="11310",]
