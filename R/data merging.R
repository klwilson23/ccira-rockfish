library(openxlsx)
library(lubridate)

setwd("~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/CCIRA Tow Video/Full Transcription/2020 Tow Video Transcription")
ndir <- length(grep('Drop',dir()))
file_name <- dir()[grep('Drop',dir())[1]]
file_path <- paste("/",file_name,sep="")
merge_data <- read.xlsx(paste(dir()[grep('Drop',dir())[1]],"/",dir(paste(getwd(),file_path,sep=""))[match(paste(file_name,'.xlsx',sep=""),dir(paste(getwd(),file_path,sep="")))],sep=""))
if(is.na(match('Time',colnames(merge_data))))
{
  colnames(merge_data)[grep("\\X",colnames(merge_data))[1]] <- "Time" 
}
merge_data$Date <- strftime(convertToDateTime(merge_data$Date),format="%Y-%m-%d")
merge_data$Time <- convertToDateTime(merge_data$Time,origin=merge_data$Date)

merge_data$duration <- max(merge_data$Time)-merge_data$Time[1]

for(i in 2:ndir)
{
  file_name <- dir()[grep('Drop',dir())[i]]
  file_path <- paste("/",file_name,sep="")
  temp_data <- read.xlsx(paste(dir()[grep('Drop',dir())[i]],"/",dir(paste(getwd(),file_path,sep=""))[match(paste(file_name,'.xlsx',sep=""),dir(paste(getwd(),file_path,sep="")))],sep=""))  
  if(is.na(match('Time',colnames(temp_data))))
  {
    colnames(temp_data)[grep("\\X",colnames(temp_data))[1]] <- "Time" 
  }
  colnames(temp_data)[grep("Time",colnames(temp_data))] <- "Time"
  colnames(temp_data)[grep("Temperature",colnames(temp_data))] <- "Temperature"
  colnames(temp_data)[grep("Depth",colnames(temp_data))] <- "Depth"
  colnames(temp_data)[colnames(temp_data)=='Lat'] <- "Latitude"
  colnames(temp_data)[colnames(temp_data)=='Long'] <- "Longitude"
  colnames(temp_data)[colnames(temp_data)=='drop.ID'] <- "Drop.ID"
  temp_data <- temp_data[,match(colnames(merge_data),colnames(temp_data),nomatch=0)]
  temp_data$Date <- strftime(convertToDateTime(temp_data$Date),format="%Y-%m-%d")
  temp_data$Time <- convertToDateTime(temp_data$Time,origin=temp_data$Date)
  temp_data$duration <- max(temp_data$Time)-temp_data$Time[1]
  merge_data <- merge(merge_data,temp_data,by = intersect(names(merge_data), names(temp_data)),all=TRUE)
}
write.csv(merge_data,"~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/CCIRA Tow Video/Full Transcription/2020_tow_video_merger.csv")

rockfish_cpue <- aggregate(count~Drop.ID+duration+fish.species,data=merge_data[merge_data$sp.id.certainty>=3,],sum,na.rm=TRUE)
rockfish_cpue$cpue <- rockfish_cpue$count/as.numeric(rockfish_cpue$duration)
rockfish_cpue <- rockfish_cpue[order(factor(rockfish_cpue$Drop.ID)),]
rockfish_cpue$Latitude <- merge_data$Latitude[match(rockfish_cpue$Drop.ID,merge_data$Drop.ID)]
rockfish_cpue$Longitude <- merge_data$Longitude[match(rockfish_cpue$Drop.ID,merge_data$Drop.ID)]

write.csv(rockfish_cpue,"~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/CCIRA Tow Video/Full Transcription/2020_rockfish_site_cpue.csv",row.names = FALSE)
