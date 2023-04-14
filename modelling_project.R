### 
# modelling project
# Diana, Ismael, Benjamin

setwd("C:/Users/diana/OneDrive/Área de Trabalho/ENM course Évora 2022/modelling_project/GBIF/")

library(readxl)
data<-read_xlsx("alldata.xlsx")
head(data)

library(dplyr)
unique(data$occurrenceStatus)
data <- data %>% filter(occurrenceStatus=="PRESENT")
unique(data$occurrenceStatus)

setwd("C:/Users/diana/OneDrive/Área de Trabalho/ENM course Évora 2022/modelling_project/")
getwd()

#crop to the atlantic
library(raster)
n1 <- list.files(pattern="1976-2005.tif$") # variables of the present time
preds_pres <- stack(n1)
plot(preds_pres[[2]]) # plot first layer
n2 <- list.files(pattern="rcp26_2046-2055.tif$") # variables of the mid century
preds_mid26 <- stack(n2)
plot(preds_mid26[[1]])
n3 <- list.files(pattern="rcp45_2046-2055.tif$") # variables of the mid century
preds_mid45 <- stack(n3)
plot(preds_mid45[[1]])
n4 <- list.files(pattern="rcp85_2046-2055.tif$") # variables of the mid century
preds_mid85 <- stack(n4)
plot(preds_mid85[[1]])

e<-extent(c(-70,0,-60,80)) #look at the map and check the coordinates
bioc <- crop(preds_pres,e) # xmin, xmax, ymin, ymax
plot(bioc[[1]])

# load species points
table(data$basisOfRecord)
data <- data %>% filter(basisOfRecord=="HUMAN_OBSERVATION")
nrow(data)
sp<- data %>% dplyr::select(decimalLongitude,decimalLatitude, species)
head(sp)

sp$decimalLatitude<- as.numeric(sp$decimalLatitude)
sp$decimalLongitude<- as.numeric(sp$decimalLongitude)
sp<- sp %>% tidyr::drop_na()
table(sp$species)

sp2<-sp %>% filter(!species=="Balaenoptera_brydei") # ! is for selecting the opposite
sp2<-sp2 %>% filter(!species=="Tursiops_truncatus________") 
sp2<-sp2 %>% filter(!species=="Balaenoptera_acutorostrata") 
sp2<-sp2 %>% filter(!species=="Balaenoptera_acutorostrata") 
sp2<-sp2 %>% filter(!species=="Balaenoptera_physalus_") 
sp2<-sp2 %>% filter(!species=="Balaenoptera_spp.") 
sp2<-sp2 %>% filter(!species=="Delphinus_delphis_")
sp2<-sp2 %>% filter(!species=="Globicephala_sp.")
sp2<-sp2 %>% filter(!species=="Mesoplodon_mirus")
sp2<-sp2 %>% filter(!species=="Mesoplodon_sp")
sp2<-sp2 %>% filter(!species=="Mesoplodon_sp.")
sp2<-sp2 %>% filter(!species=="Stenella_Frontalis")
sp2<-sp2 %>% filter(!species=="Tursiops_truncatus________")
sp2 <- sp2 %>% filter(!species=="Tursiops_truncatus")
sp2 <- sp2 %>% filter(!species=="Stenella_frontalis")
sp2 <- sp2 %>% filter(!species=="Stenella_coeruleoalba")
sp2 <- sp2 %>% filter(!species=="Physeter_macrocephalus")
sp2 <- sp2 %>% filter(!species=="Mesoplodon_bidens")
sp2 <- sp2 %>% filter(!species=="Globicephala_macrorhynchus")
sp2 <- sp2 %>% filter(!species=="Balaenoptera_borealis")
sp2 <- sp2 %>% filter(!species=="Balaenoptera_edeni")
table(sp2$species)

dups2 <-duplicated(sp2[,c("decimalLatitude","decimalLongitude")])
sp3<- sp2[!dups2,]
points(sp3,col="red")
table(sp3$species)

# to select only the points inside the Atlantic area #bioc
library(terra)
sp4<- vect(sp3,geom=c("decimalLatitude", "decimalLongitude"),crs=projection(bioc))
sp5<-crop(sp4,e)
plot(bioc[[3]])
points(sp5,col="red")
nrow(sp5)
table(sp5$species) 
points(sp6$species, add=T)
na.omit(sp5)

sp6<-as(sp5,"Spatial")
points(sp6,col="red")
sp7 <- as.data.frame(sp6)
unique(sp6$species)
sp7.1<- sp7 %>% filter(species=="Balaenoptera_musculus")
sp7.2<- sp7 %>% filter(species=="Pseudorca_crassidens")
points(sp7.2$x,sp7.2$y, col="red", pch= 8, cex=0.5)

unique(sp7$species)

##

#######
## correlation
library(usdm)
v <- vifstep(bioc) # just for the present variables
vif(bioc) # tos is correlated with dissolved oxygen so is dropped
v
bioc2 <- exclude(bioc,v) #exclude the correlated variables
bioc2


### modelling
library(sdm)

spn <-unique(sp6$species)
n <- spn[1]
biocterra <- rast(bioc2)
biocdf <- as.data.frame(biocterra,cell=T)
head(biocdf)

for (n in spn) {
  sp1<- sp6[sp6$species == n,]
  sp1$species <- 1
  d<-sdmData(species~., sp1, predictors = bioc2, bg = list(n=500))
  
  m1<-sdm(species~., d,methods=c('glm','rf','maxent','gam'),
          replication='sub',test.p=30,n=2,parallelSettings=list(method="parallel",ncore=6))
  write.sdm(m1,paste0("MODEL_",n,".sdm"))
  # en <- ensemble(m1, bioc2,filename = paste0("Ens_",n,".tif"),
  #                setting = list(method="weighted",stat="auc",power=2))
  rr <- rast(biocterra[[1]])
  en <- ensemble(m1, biocdf,
                 setting = list(method="weighted",stat="auc",power=2))
  rr[biocdf$cell] <- en
  writeRaster(rr,filename = paste0("Ens_",n,".tif"))
  
  cat("\n The Model for species",n,"is done!\n ********************")
  
  
}



lst <- list.files(pattern = "Ens_")
r <- rast(lst)
plot(r)

rich <- app(r, sum)
plot(rich)

for (n in spn) {
  en <- raster(paste0("Ens_",n,".tif"))
  m <- read.sdm(paste0("MODEL_",n,".sdm"))
  th <- evaluates(m@data,en)@threshold_based$threshold[2]
  pa <- en
  pa[] <- ifelse(en[] >= th, 1,0)
  writeRaster(rast(pa),paste0("PA_",n,".tif"))
  cat(".")
  
}

#----------------

# calculate species richness
list.files(pattern = "PA_")
lst <- list.files(pattern = "PA_")
r <- rast(lst)
plot(r)


rich <- app(r, sum) / 7
plot(rich)
writeRaster(rich, "Richness_7species.tif")


d<-sdmData(species~., sp6, predictors = bioc2, bg = list(n=500))
class(sp6)
class(bioc2)
m1<-sdm(species~., d,methods=c('glm','rf','maxent','gam'),
       replication='sub',test.p=30,n=2)
# 
# .ens_terra <- function(m,pr)  {
#   ensemble(m,pr,setting = list(method="weighted",stat="auc",power=2))
# }
# plot(en)
# n
# .r <- predict(biocterra,m,fun=.ens_terra)
# plot(.r)


#################################

## cropping environmental variables of the future
bioc_mid26 <- crop(preds_mid26,e) # xmin, xmax, ymin, ymax
plot(bioc_mid26[[1]])
bioc_mid45 <- crop(preds_mid45,e) # xmin, xmax, ymin, ymax
plot(bioc_mid45[[1]])
bioc_mid85 <- crop(preds_mid85,e) # xmin, xmax, ymin, ymax
plot(bioc_mid85[[1]])

#transform the names of the future variables as the same as present
names(bioc_mid26)=names(bioc2)
names(bioc_mid45)=names(bioc2)
names(bioc_mid85)=names(bioc2)
# this must be changed according with the future 26 45 85
plot(bioc_mid26[[1]])
biocterra <- rast(bioc_mid26)
biocdf <- as.data.frame(biocterra,cell=T)
head(biocdf)

for (n in spn) {
  sp1<- sp6[sp6$species == n,]
  sp1$species <- 1
  d<-sdmData(species~., sp1, predictors = bioc_mid26, bg = list(n=500))
  
  m1<-sdm(species~., d,methods=c('glm','rf','maxent','gam'),
          replication='sub',test.p=30,n=2,parallelSettings=list(method="parallel",ncore=6))
  write.sdm(m1,paste0("MODEL_future_mid26",n,".sdm"))
  #en <- ensemble(m1, bioc_mid26,filename = paste0("Ens_future_mid26",n,".tif"),
  #               setting = list(method="weighted",stat="auc",power=2))
  rr <- rast(biocterra[[1]])
  en <- ensemble(m1, biocdf,
                 setting = list(method="weighted",stat="auc",power=2))
  rr[biocdf$cell] <- en
  writeRaster(rr,filename = paste0("Ens_future_mid26",n,".tif"))
  
  cat("\n The Model for species26",n,"is done!\n ********************")
  
} #26

plot(bioc_mid45[[1]])
biocterra <- rast(bioc_mid45)
biocdf <- as.data.frame(biocterra,cell=T)
head(biocdf)
for (n in spn) {
  sp1<- sp6[sp6$species == n,]
  sp1$species <- 1
  d<-sdmData(species~., sp1, predictors = bioc_mid45, bg = list(n=500))
  
  m1<-sdm(species~., d,methods=c('glm','rf','maxent','gam'),
          replication='sub',test.p=30,n=2,parallelSettings=list(method="parallel",ncore=6))
  write.sdm(m1,paste0("MODEL_future_mid45",n,".sdm"))
  #en <- ensemble(m1, bioc_mid26,filename = paste0("Ens_future_mid26",n,".tif"),
  #               setting = list(method="weighted",stat="auc",power=2))
  rr <- rast(biocterra[[1]])
  en <- ensemble(m1, biocdf,
                 setting = list(method="weighted",stat="auc",power=2))
  rr[biocdf$cell] <- en
  writeRaster(rr,filename = paste0("Ens_future_mid45",n,".tif"))
  
  cat("\n The Model for species26",n,"is done!\n ********************")
  
} #45


plot(bioc_mid85[[1]])
biocterra <- rast(bioc_mid85)
biocdf <- as.data.frame(biocterra,cell=T)
head(biocdf)
for (n in spn) {
  sp1<- sp6[sp6$species == n,]
  sp1$species <- 1
  d<-sdmData(species~., sp1, predictors = bioc_mid85, bg = list(n=500))
  
  m1<-sdm(species~., d,methods=c('glm','rf','maxent','gam'),
          replication='sub',test.p=30,n=2,parallelSettings=list(method="parallel",ncore=6))
  write.sdm(m1,paste0("MODEL_future_mid85",n,".sdm"))
  #en <- ensemble(m1, bioc_mid26,filename = paste0("Ens_future_mid26",n,".tif"),
  #               setting = list(method="weighted",stat="auc",power=2))
  rr <- rast(biocterra[[1]])
  en <- ensemble(m1, biocdf,
                 setting = list(method="weighted",stat="auc",power=2))
  rr[biocdf$cell] <- en
  writeRaster(rr,filename = paste0("Ens_future_mid85",n,".tif"))
  
  cat("\n The Model for species26",n,"is done!\n ********************")
  
} #85

m1 <- raster("Ens_mid26Ziphius_cavirostris.tif")
plot(m1)
read.sdm("MODEL_future_mid26Pseudorca_crassidens.sdm")



for (n in spn) {
  en <- raster(paste0("Ens_future_mid26",n,".tif"))
  m <- read.sdm(paste0("MODEL_future_mid26",n,".sdm"))
  th <- evaluates(m@data,en)@threshold_based$threshold[2]
  pa <- en
  pa[] <- ifelse(en[] >= th, 1,0)
  writeRaster(rast(pa),paste0("PA_26",n,".tif"))
  cat(".")
  
}

for (n in spn) {
  en <- raster(paste0("Ens_future_mid45",n,".tif"))
  m <- read.sdm(paste0("MODEL_future_mid45",n,".sdm"))
  th <- evaluates(m@data,en)@threshold_based$threshold[2]
  pa <- en
  pa[] <- ifelse(en[] >= th, 1,0)
  writeRaster(rast(pa),paste0("PA_45",n,".tif"))
  cat(".")
  
}


for (n in spn) {
  en <- raster(paste0("Ens_future_mid85",n,".tif"))
  m <- read.sdm(paste0("MODEL_future_mid85",n,".sdm"))
  th <- evaluates(m@data,en)@threshold_based$threshold[2]
  pa <- en
  pa[] <- ifelse(en[] >= th, 1,0)
  writeRaster(rast(pa),paste0("PA_85",n,".tif"))
  cat(".")
  
}

for (n in spn) {
  en <- raster(paste0("Ens_future_mid26",n,".tif"))
  m <- read.sdm(paste0("MODEL_future_mid26",n,".sdm"))
  th <- evaluates(m@data,en)@threshold_based$threshold[2]
  pa <- en
  pa[] <- ifelse(en[] >= th, 3,0)
  writeRaster(rast(pa),paste0("PA_26_colonization",n,".tif"))
  cat(".")
  
}

for (n in spn) {
  en <- raster(paste0("Ens_future_mid26",n,".tif"))
  m <- read.sdm(paste0("MODEL_future_mid26",n,".sdm"))
  th <- evaluates(m@data,en)@threshold_based$threshold[2]
  pa <- en
  pa[] <- ifelse(en[] >= th, 0,1)
  writeRaster(rast(pa),paste0("PA_26_extintion",n,".tif"))
  cat(".")
  
}


for (n in spn) {
  en <- raster(paste0("Ens_future_mid85",n,".tif"))
  m <- read.sdm(paste0("MODEL_future_mid85",n,".sdm"))
  th <- evaluates(m@data,en)@threshold_based$threshold[2]
  pa <- en
  pa[] <- ifelse(en[] >= th, 0,1)#extinti
  writeRaster(rast(pa),paste0("PA_85_extintion",n,".tif"))
  cat(".")
  
}

for (n in spn) {
  en <- raster(paste0("Ens_future_mid85",n,".tif"))
  m <- read.sdm(paste0("MODEL_future_mid85",n,".sdm"))
  th <- evaluates(m@data,en)@threshold_based$threshold[2]
  pa <- en
  pa[] <- ifelse(en[] >= th, 3,0)# colonization
  writeRaster(rast(pa),paste0("PA_85_colonization",n,".tif"))
  cat(".")}


t11 <- raster("Ens_Balaenoptera_musculus.tif")
plot(t11)
points(sp7.1$x,sp7.1$y, col="red", pch= 16, cex=0.5)
head(sp7.1)
coordinates(sp7.1) = c("x","y")
mapview(t11)+ sp7.1

t12 <- raster("Ens_future_mid26Balaenoptera_musculus.tif")
plot(t12)

t13 <- raster("Ens_future_mid45Balaenoptera_musculus.tif")
plot(t13)

t14 <- raster("Ens_future_mid85Balaenoptera_musculus.tif")
plot(t14)

plot(stack(t11,t12,t13,t14))

t15<- raster("Ens_Pseudorca_crassidens.tif")
plot(t15)

points(sp7.2$x,sp7.2$y, col="red", pch= 16, cex=0.5)
head(sp7.2)
coordinates(sp7.2) = c("x","y")
mapview(t15)+ sp7.2

#variable importance 
library(sdm)
a<-read.sdm("MODEL_Pseudorca_crassidens.sdm")
getVarImp(a)
