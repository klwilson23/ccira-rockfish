#devtools::install_github("James-Thorson-NOAA/FishLife")
library(FishLife)
library(TMB)
thorson_traits <- as.data.frame(FishLife::FishBase_and_RAM$beta_gv)
thorson_covar <- as.data.frame(FishLife::FishBase_and_RAM$Cov_gvv)
thorson_covar <- thorson_covar[,intersect(colnames(thorson_covar),paste(colnames(thorson_traits),".",colnames(thorson_traits),sep=""))]
colnames(thorson_covar) <- colnames(thorson_traits)
leading_traits <- c("Loo","K","Winfinity","tmax","tm","M","Lm","ln_var","ln_MASPS","ln_margsd","ln_Fmsy_over_M","ln_Fmsy")
thorson_traits[,leading_traits] <- exp(thorson_traits[,leading_traits]+0.5*thorson_covar[,leading_traits])
row.names(thorson_traits)[grep("Sebastes",row.names(thorson_traits))]
rockfish <- rbind(thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_predictive",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_ruberrimus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_borealis",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_aleutianus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_melanostictus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_paucispinis",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_pinniger",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_maliger",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_babcocki",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_alascanus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_melanops",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_nebulosus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_caurinus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_helvomaculatus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_proriger",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_brevispinis",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_miniatus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_auriculatus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_variabilis ",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_miniatus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_ciliatus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_elongatus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_zacentrus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_entomelas",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_flavidus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_diaconus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_wilsoni",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_emphaeus",],
                  thorson_traits["Actinopterygii_Scorpaeniformes_Sebastidae_Sebastes_jordani",])
rockfish[,c("Loo","K","tmax","tm","M","h","r","ln_Fmsy_over_M","ln_Fmsy","G")]
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}
rockfish_norm <- apply(rockfish,2,normalize)
plot(ln_Fmsy_over_M~tmax,data=rockfish_norm)
plot(r~tmax,data=rockfish_norm)
plot(h~tmax,data=rockfish_norm)

vignette("tutorial","FishLife")

# Re-run results with a different model configuration
Ynew_ij = matrix( c("Loo"=NA,"K"=NA,"Winfinity"=NA,"tmax"=119,"tm"=NA,"M"=NA,"Lm"=NA,"Temperature"=NA), nrow=1)
Update = Update_prediction( Taxon=Search_species(Genus="Sebastes",Species="ruberrimus",add_ancestors=FALSE)$match_taxonomy, Ynew_ij=Ynew_ij)
Update$updateMean_j