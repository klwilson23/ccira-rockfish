library(dplyr)
library(tidyr)
rockfish <- read.csv("Data/Rockfish hotspots 4km PUID.csv")
corals <- read.csv("Data/coral hotspots 4km PUID.csv")
sponges <- read.csv("Data/sponge hotspots 4km PUID.csv")
sponges$UpperOceanSR[sponges$PU_4Km_ID=="6896"] <- "(11) Mainland Fjords"

rockfish_old <- rockfish
corals_old <- corals
sponges_old <- sponges

rockfish <- rockfish[,-((ncol(rockfish)-2):ncol(rockfish))]
corals <- corals[,-((ncol(corals)-2):ncol(corals))]
sponges <- sponges[,-((ncol(sponges)-2):ncol(sponges))]

nrow(corals);nrow(sponges)
sponges$PU_4Km_ID[!sponges$PU_4Km_ID%in%corals$PU_4Km_ID]
sponges[which(!sponges$PU_4Km_ID%in%corals$PU_4Km_ID),]
rockfish$taxa <- "rockfish"
corals$taxa <- "corals"
sponges$taxa <- "sponges"
ncc <- rbind(sponges,rbind(rockfish,corals))
ncc <- ncc %>%
  group_by(taxa) %>% 
  mutate(rank = ntile(normalized_score,10)) %>%
  ungroup(taxa)
ncc_wide <- pivot_wider(ncc,names_from=taxa,id_cols=c(PU_4Km_ID),values_from=c(normalized_score))
ncc_agg <- data.frame("PUID"=ncc_wide$PU_4Km_ID,ncc_wide[,-1])
ncc_agg$group_score <- rowSums(ncc_agg[,2:4])
ncc_agg <- data.frame(ncc_agg,ncc_agg[,2:4]/ncc_agg$group_score)
ncc_agg$even_score <- apply(ncc_agg[,5:8],1,function(x){x[1]+-sum(x[2:4]*log(x[2:4]))}) # evenness

ncc_agg$UpperOceanSR <- rockfish$UpperOceanSR[match(ncc_agg$PUID,rockfish$PU_4Km_ID)]
ncc_agg$max_depth <- rockfish$max_depth[match(ncc_agg$PUID,rockfish$PU_4Km_ID)]
ncc_agg$mean_depth <- rockfish$mean_depth[match(ncc_agg$PUID,rockfish$PU_4Km_ID)]
ncc_agg$rf_hotspot <- rockfish_old$rockfish_hotspot[match(ncc_agg$PUID,rockfish_old$PU_4Km_ID,nomatch = NA)]
ncc_agg$cor_hotspot <- corals_old$coral_hotspot[match(ncc_agg$PUID,corals_old$PU_4Km_ID,nomatch = NA)]
ncc_agg$sponge_hotspot <- sponges_old$sponge_hotspot[match(ncc_agg$PUID,sponges_old$PU_4Km_ID,nomatch = NA)]
ncc_agg$rf_cor <- ifelse(ncc_agg$rf_hotspot==1 & ncc_agg$cor_hotspot==1, 1, 0)
ncc_agg$rf_sp <- ifelse(ncc_agg$rf_hotspot==1 & ncc_agg$sponge_hotspot==1, 1, 0)
ncc_agg$cor_sp <- ifelse(ncc_agg$cor_hotspot==1 & ncc_agg$sponge_hotspot==1, 1, 0)
ncc_agg$rf_cor_sp <- ifelse(ncc_agg$rf_hotspot==1 & ncc_agg$cor_hotspot==1 & ncc_agg$sponge_hotspot==1, 1, 0)
ncc_agg <- ncc_agg %>%
  mutate(puid_rank = ntile(even_score, 10))
ncc_agg$all_hotspots <- ifelse(ncc_agg$puid_rank==10,1,0)
write.csv(ncc_agg,"Data/NCC hotspots 4km puid.csv")
plot(ncc_agg$max_depth,ncc_agg$puid_rank)
ncc_agg[is.na(ncc_agg$group_score),]
table(ncc_agg$rf_cor)
table(ncc_agg$rf_sp)
table(ncc_agg$cor_sp)
table(ncc_agg$rf_cor_sp)



rockfish <- read.csv("Data/Rockfish hotspots 1km PUID.csv")
corals <- read.csv("Data/coral hotspots 1km PUID.csv")
sponges <- read.csv("Data/sponge hotspots 1km PUID.csv")

rockfish_old <- rockfish
corals_old <- corals
sponges_old <- sponges

rockfish <- rockfish[,-((ncol(rockfish)-2):ncol(rockfish))]
corals <- corals[,-((ncol(corals)-2):ncol(corals))]
sponges <- sponges[,-((ncol(sponges)-2):ncol(sponges))]

nrow(corals);nrow(sponges)
sponges$PU_1Km_ID[!sponges$PU_1Km_ID%in%corals$PU_1Km_ID]
sponges[which(!sponges$PU_1Km_ID%in%corals$PU_1Km_ID),]
rockfish$taxa <- "rockfish"
corals$taxa <- "corals"
sponges$taxa <- "sponges"
ncc <- rbind(sponges,rbind(rockfish,corals))
ncc <- ncc %>%
  group_by(taxa) %>% 
  mutate(rank = ntile(normalized_score,10)) %>%
  ungroup(taxa)
ncc_wide <- pivot_wider(ncc,names_from=taxa,id_cols=c(PU_1Km_ID),values_from=c(normalized_score))
ncc_agg <- data.frame("PUID"=ncc_wide$PU_1Km_ID,ncc_wide[,-1])
ncc_agg$group_score <- rowSums(ncc_agg[,2:4])
ncc_agg <- data.frame(ncc_agg,ncc_agg[,2:4]/ncc_agg$group_score)
ncc_agg$even_score <- apply(ncc_agg[,5:8],1,function(x){x[1]+-sum(x[2:4]*log(x[2:4]))}) # evenness
ncc_agg$PUID_4km <- rockfish$PU_4Km_ID[match(ncc_agg$PUID,rockfish$PU_1Km_ID)]
ncc_agg$UpperOceanSR <- rockfish$UpperOceanSR[match(ncc_agg$PUID,rockfish$PU_1Km_ID)]
ncc_agg$max_depth <- rockfish$max_depth[match(ncc_agg$PUID,rockfish$PU_1Km_ID)]
ncc_agg$mean_depth <- rockfish$mean_depth[match(ncc_agg$PUID,rockfish$PU_1Km_ID)]
ncc_agg$rf_hotspot <- rockfish_old$rockfish_hotspot[match(ncc_agg$PUID,rockfish_old$PU_1Km_ID,nomatch = NA)]
ncc_agg$cor_hotspot <- corals_old$coral_hotspot[match(ncc_agg$PUID,corals_old$PU_1Km_ID,nomatch = NA)]
ncc_agg$sponge_hotspot <- sponges_old$sponge_hotspot[match(ncc_agg$PUID,sponges_old$PU_1Km_ID,nomatch = NA)]
ncc_agg$rf_cor <- ifelse(ncc_agg$rf_hotspot==1 & ncc_agg$cor_hotspot==1, 1, 0)
ncc_agg$rf_sp <- ifelse(ncc_agg$rf_hotspot==1 & ncc_agg$sponge_hotspot==1, 1, 0)
ncc_agg$cor_sp <- ifelse(ncc_agg$cor_hotspot==1 & ncc_agg$sponge_hotspot==1, 1, 0)
ncc_agg$rf_cor_sp <- ifelse(ncc_agg$rf_hotspot==1 & ncc_agg$cor_hotspot==1 & ncc_agg$sponge_hotspot==1, 1, 0)
ncc_agg <- ncc_agg %>%
  mutate(puid_rank = ntile(even_score, 10))
ncc_agg$all_hotspots <- ifelse(ncc_agg$puid_rank==10,1,0)
write.csv(ncc_agg,"Data/NCC hotspots 1km puid.csv")

plot(ncc_agg$max_depth,ncc_agg$puid_rank)
ncc_agg[is.na(ncc_agg$group_score),]
table(ncc_agg$rf_cor)
table(ncc_agg$rf_sp)
table(ncc_agg$cor_sp)
table(ncc_agg$rf_cor_sp)
ncc_agg[ncc_agg$all_hotspots==1,]
