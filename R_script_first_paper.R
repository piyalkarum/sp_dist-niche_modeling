
# loading the collection data
FCM <- read.table("R_fcm28112016.txt", header = T)
head(FCM, 6) # take a look at the data


### Spatial analysis ###
library (spatstat)

# create a point pattern dataset from the collection points
FCM_ppp <- ppp(FCM$Longitude, FCM$Latitude, c(-65, -55), c(-33, -23)) # coordinates are the window values
marks (FCM_ppp) <- factor(FCM$Ploidy1, labels = c("Diploid", "Tetraploids") # adding mark labels to the ppp object

#quadrat test
QT <- quadrat.test(FCM_ppp, nx=10, ny=10) 
#Performs a test of Complete Spatial Randomness for a given point pattern, based on #quadrat counts. Alternatively performs a goodness-of-fit test of a fitted #inhomogeneous Poisson model. By default performs chi-squared tests; can also perform #Monte Carlo based tests
QT
plot(FCM_ppp, pch=19, cols = c(3,2), cex = 0.75, main = "Chi-squared test using quadrat counts for cytotype distribution")
plot(QT, add=T, cex=0.5) # adding quadrat test values to the plot

# K means test
Kest(FCM_ppp)
#Estimates Ripley's reduced second moment function K(r) from a point pattern in a #window of arbitrary shape
plot(Kest(FCM_ppp)) # plot K means test


#### Environmental data ######

# raster brick of evironmental variables (the stack is prepared separately: see the script data.prep.R)
ENV.26 <- brick("two.six.variables.grd")

# extracting environmental data from raster files
Dip <- subset(FCM, FCM$Ploidy1==2) #get the diploids
Tet <- subset(FCM, FCM$Ploidy1==4) #get the tetraploids
Bclimdata2 <- extract(ENV.26, Dip[,c("Longitude", "Latitude")]) #extract climate data for diploids
Bclimdata4 <- extract(ENV.26, Tet[,c("Longitude", "Latitude")]) #extract climate data for tetraploids
Bclimdata2 <- data.frame(lon=Dip$Longitude, lat=Dip$Latitude, cyt=Dip$Ploidy1, Bclimdata2) # make a data frame with coordinates and ploidy levels for diploids
Bclimdata4 <- data.frame(lon=Tet$Longitude, lat=Tet$Latitude, cyt=Tet$Ploidy1, Bclimdata4) # make a data frame with coordinates and ploidy levels for tetraploids 
allclim <- merge(Bclimdata2, Bclimdata4, all=T) # merge the data frames to make one
allclim$cyt <- factor(allclim$cyt, labels = c("Dip", "Tet")) # factorize ploidy level

# saving new data set allclim
write.table(allclim,"allclim.26.var.txt")

## GLM on the climate data 
ml <- allclim
ml$pl <- relevel(ml$cyt, ref = "Dip")
Multinom <- multinom(pl ~ bio1+bio2+bio3+bio4+bio5+bio6+bio7+bio15+bio16+bio18+bio19+uvb+elevation+soiltype+PAR+cloud+frost+vapor, data = ml)
Z <- summary(Multinom)$coefficients/summary(Multinom)$standard.errors
P <- (1-pnorm(abs(Z), 0,1))*2

predict(Multinome, ml)
table (predict(Multinom, ml), ml$cyt)

#logistic regression output
p.value <- c()
for (i in 4:length(allclim[1,])){
  
  GLMtest <- summary(glm(allclim$cyt~allclim[,i], family = binomial))
  p.value[i-3] <- GLMtest$coefficients[8]
  
}
Log.reg <- data.frame(cbind(names(allclim[4:29])), p.value) # this table gives the significance of each climate variable for the geographical distribution of cytotypes in the studied area. The variables with higher significance are retained for further analysis i.e. niche assessment, distribution modeling.


########## PCA with ade4 and factoextra  ######

library(ade4)
library(factoextra)

ADEpca <- dudi.pca(allclim[,-c(1:3)], scale = T) # PCA with scaling

#simple visualization of the PCA
s.corcircle(ADEpca$co)
fviz_pca_var(ADEpca)

ADEpca$co # get PC values

INt <- inertia.dudi(ADEpca, row.inertia = T, col.inertia = T) # get intertia values of the PCA
varcos2 <- abs(INt$row.rel/10000) # calculate cosine of the variable inertia
VAR <- get_pca_var(ADEpca) # get PCA variables

# plot the contribution of each variable to the PCA
fviz_pca_var(ADEpca, col.var="contrib")+scale_color_gradient2(low="white", mid="blue", 
high="red", midpoint=3.2)
fviz_contrib(ADEpca, choice = "var", axes = 1) # variable contribution for PCs
fviz_contrib(ADEpca, choice = "ind", axes = 1) # individual contribution for PCs
# variables with higher contribution (>50%) are retained for further analysis

s.label(ADEpca$li, xax = 1, yax = 2, label = allclim$cyt, boxes = F) # plot points
s.arrow(7*ADEpca$c1, add.plot = T) # add variables
s.class(ADEpca$li, fac = allclim$cyt, xax = 1, yax = 2, col=c("red", "blue") # group individuals and make ellipses
# modification to get only the ellipse for both cytotypes
s.class(ADEpca$li, fac = factor(rep(1,68)), xax = 1, yax = 2, cpoint = 0, cstar = 0, col = "black", add.plot = T, cellipse = 1.5)
#get only the ellepses for separate cytotypes
s.class(ADEpca$li, fac = allclim$cyt, addaxes = T, xax = 1, yax = 2, cpoint = 0, cstar = 0, col = c("blue", "red"))

# biplot with a different perspective
fviz_pca_biplot(ADEpca, habillage = allclim$cyt, addEllipses =F,col.var = "red", alpha.var ="cos2",label = "var") +scale_color_brewer(palette="Dark2")

###boxplots of env variables (selected)###
par(mfrow=c(1,7))
plot(allclim$cyt,allclim$bio1, main="Bio1",frame=F, ylab="Temperature x 10 °C" )
plot(allclim$cyt,allclim$uvb, main="UV-B",frame=F, ylab="Radiation J/m^2/day")
plot(allclim$cyt,allclim$PAR, main="PAR",frame=F, ylab="Radiation J/m^2/day")
plot(allclim$cyt,allclim$bio3, main="Bio3",frame=F,ylab="Temperature x 10 °C %" )
plot(allclim$cyt,allclim$bio4, main="Bio4",frame=F, ylab="Temperature x 10 °C x 100")
plot(allclim$cyt,allclim$bio7, main="Bio7",frame=F, ylab="Temperature annual range x10°C")
plot(allclim$cyt,allclim$bio6, main="Bio6",frame=F, ylab="Temperature x 10 °C")
par(mfrow=c(1,7))
plot(allclim$cyt,allclim$bio8, main="Bio8",frame=F, ylab="Temperature x 10 °C")
plot(allclim$cyt,allclim$bio9, main="Bio9",frame=F, ylab="Temperature x 10 °C")
plot(allclim$cyt,allclim$bio10, main="Bio10",frame=F, ylab="Temperature x 10 °C")
plot(allclim$cyt,allclim$bio11, main="Bio11",frame=F, ylab="Temperature x 10 °C")
plot(allclim$cyt,allclim$elevation, main="Elevation",frame=F, ylab="elevation/m")
plot(allclim$cyt,allclim$bio15, main="Precipitation seasanality",frame=F, ylab="Coefficient of variation")
plot(allclim$cyt, allclim$vapor, main="Vapor pressure", ylab="Pressure hPa", frame=F)
#plot(allclim$cyt, allclim$cloud/10, main="Cloud cover", ylab="Cover % ", frame=F)
plot(allclim$cyt, allclim$frost/10, main="Frost day frequency", ylab="No. of days", frame=F)
#plot(allclim$cyt,allclim$bio1, main="Bio1",frame=F, ylab="Temperature x 10 °C" )


#### plotting niche density with pie charts and boxplots ####

nichall <- allenvdata
nichedip <- subset(nichall, nichall$cyt=="Dip" )
nichedip <- nichedip[,-c(1:3,5,7,8,10,15:22,25)]
nichetet <- subset(nichall, nichall$cyt=="Tet" )
nichetet <- nichetet[,-c(1:3,5,7,8,10,15:22,25)]

nichedipdf <- data.frame(nichedip)
nichetetdf <- data.frame(nichetet)


d2 <- sapply(nichedip, log)
d2 <- rowSums(d2)
#d2 <- density(d2, bw="SJ", adjust=1.5) # an alternative kernel density
d2 <- density(d2, bw="nrd", adjust=1, cut=c(47.21239,48.52767))

d4 <- sapply(nichetet, log)
d4 <- rowSums(d4)
#d4 <- density(d4, bw="nrd", adjust=1) # an alternative kernel density
d4 <- density(d4, bw="nrd", adjust=1, cut=c(45.67400,48.61277))

# making a separate dataframe for reference
D2 <- sapply(nichedip, log)
D2 <- rowSums(D2)
D2 <- data.frame(D2)
D4 <- sapply(nichetet, log)
D4 <- rowSums(D4)
D4 <- data.frame(D4)
cyt2 <- rep("Dip", 23)
D2 <- cbind(D2, cyt2)
cyt4 <- rep("Tet", 41)
D4 <- cbind(D4, cyt4)
colnames(D2) <- c("sum", "cyt")
colnames(D4) <- c("sum", "cyt")
Dall <- merge(D2,D4, all=T)

par(fig=c(0,1,0,0.8))
#par(mar=c(1,1,1,1))
plot(d2, type="n",ylim=c(0,1), xlim=c(44.5, 49.6), xlab=NA, xaxt="n", ylab="Probability of occurence", main=NA, frame=F, srt=45)
# make transparent colors
mycol1 <- rgb(255,0,0, max=255, alpha = 125, names = "redt")
mycol2 <- rgb(0,0,255, max=255, alpha = 125, names = "bluet")
polygon(d2, col=mycol1)
polygon(d4, col=mycol2)

# adding pie charts of mixed ploidy populations (proportions) to the density plot
# data
mixedpop <- read.table("mixpop.txt", header = T)
colnames(mixedpop) <- c("2Q", "2Y", "2T", "2W")
twoQ <- mixedpop$2Q
twoY <- mixedpop$2Y
twoT <- mixedpop$2T
twoW <- mixedpop$2W

library(plotrix)
floating.pie(47.86037,0.15, x=twoQ, radius = 0.09, col = c("blue", "red"))
floating.pie(47.44111,0.25, x=twoY, radius = 0.09, col = c("blue", "red"))
floating.pie(47.83318,0.25, x=twoT, radius = 0.09, col = c("blue", "red"))
floating.pie(47.49500,0.15, x=twoW, radius = 0.09, col = c("blue", "red"))

#text(47.86037,0, "2Q", cex=0.7)
#text(47.44111,0.2, "2Y", cex=0.7)
#text(47.83318,0.2, "2T", cex=0.7)
#text(47.49500,0, "2W", cex=0.7)

# add the boxplots of the density of cytotypes 
par(fig=c(0,1,0.6,1), new=TRUE)
BxPlot <- plot(Dall$cyt, Dall$sum, horizontal=T, ylim=c(44.5, 49.6), xaxt="n", frame=F, range=0, col=c(mycol1,mycol2))


## Niche breadth and overlap analysis ###


library(ecospat)
library(dismo)
library(maptools)

#data for ecospat niche analysis
niche_dip_data <- subset(alenv, alenv$cyt=="Dip")
niche_tet_data <- subset(alenv, alenv$cyt=="Tet")

#extracting absence points for diploid
niche_dip_circles <- circles(niche_dip_data[,c("lon", "lat")], d=65000, lonlat=T)
niche_dip_ab_samples <- spsample(niche_dip_circles@polygons, 500, type = "random", iter=5000)
#extracting pseudo-absence points for tetraploid
niche_tet_circles <- circles(niche_tet_data[,c("lon", "lat")], d=100000, lonlat=T)
niche_tet_ab_samples <- spsample(niche_tet_circles@polygons, 800, type = "random", iter=5000)

#500 and 800 absence points produced good results

#extracting data for diploid
envdmdip_env <- extract(ENV, envdmdip[,c("lon", "lat")]) #already there in alenv
niche_dip_bgpoints <- extract(envar, niche_dip_ab_samples) #background data
#extracting data for tetraploid
envdmtet_env <- extract(ENV, envdmtet[,c("lon", "lat")]) #already there in alenv
niche_tet_bgpoints <- extract(envar, niche_tet_ab_samples) #background data

#adding coordinates to background points
dip_bgpoints <- niche_dip_ab_samples@coords
colnames(dip_bgpoints) <- c("lon", "lat")
niche_dip_bgpoints <- data.frame(cbind(dip_bgpoints, niche_dip_bgpoints))
length(which(is.na(niche_dip_bgpoints$bio2)))

tet_bgpoints <- niche_tet_ab_samples@coords
colnames(tet_bgpoints) <- c("lon", "lat")
niche_tet_bgpoints <- data.frame(cbind(tet_bgpoints, niche_tet_bgpoints))
length(which(is.na(niche_tet_bgpoints$bio6)))

#dip and tet data combined
niche_dip_and_tet <- rbind(niche_dip_data, niche_tet_data)

#ploidy levels
niche_pl.levels <- levels(niche_dip_and_tet[,3])

#number of interation for the tests of equivalency and similarity
iterations <- 100
#resolution of the gridding of the climate space
R <- 100

###pca environment

#combining all data to one data frame
niche_all_data <- rbind(niche_dip_data[,c(4:26)], niche_tet_data[,c(4:26)], niche_dip_bgpoints[,c(3:25)], niche_tet_bgpoints[,c(3:25)])
ones_and_zeros <- c(rep(0, nrow(rbind(niche_dip_data, niche_tet_data))), rep(1, nrow(rbind(niche_dip_bgpoints[-8,], niche_tet_bgpoints))))

#there is a NA in the 8th row of soiltype of niche_dip_bgpoints data frame
niche_pca <- dudi.pca(niche_all_data, row.w = ones_and_zeros, center = T, scale = T, scannf = F, nf=2)

# selection of ploidies
ploidy_combn <- combn(1:2,2)

#extracting scores from the PCA
for(i in 1:ncol(ploidy_combn)){
  row.dip <- which(niche_dip_and_tet[,3] == niche_pl.levels[ploidy_combn[1,i]])
  row.tet <- which(niche_dip_and_tet[,3] == niche_pl.levels[ploidy_combn[2,i]])
  name.dip <- niche_pl.levels[ploidy_combn[1,i]]
  name.tet <- niche_pl.levels[ploidy_combn[2,i]]
  #predict the scores on the axes
  niche_scores.env <- niche_pca$li[(nrow(niche_dip_and_tet[,-c(1:2)])+1):nrow(niche_all_data),]
  scores.dip <- niche_pca$li[row.dip,]
  scores.tet <- niche_pca$li[row.tet,]
}

#calculation of occurence density and test of niche equivalencey and similarity
niche.dip.grid <- ecospat.grid.clim.dyn(niche_scores.env, niche_scores.env, scores.dip, R=100)
niche.tet.grid <- ecospat.grid.clim.dyn(niche_scores.env, niche_scores.env, scores.tet, R=100)

#niche equivalencey test
niche.eq <- ecospat.niche.equivalency.test(niche.dip.grid,niche.tet.grid, alternative = "lower",rep=100)
#plotting niche overlap under niche equivalency
ecospat.plot.overlap.test(niche.eq, "D", "niche overlap") #instead use the modified function below

#niche similarity test
niche.sm <- ecospat.niche.similarity.test (niche.dip.grid,niche.tet.grid, rep=100)
#plotting niche overlap under niche similarity
ecospat.plot.overlap.test(niche.sm, "D", "Niche overlap under similarity")#use the modified funtion

#better plotting
par(mfrow=c(1,2))
Nicheoverlap.plot(niche.eq, "D", "Niche overlap under equivalency")
Nicheoverlap.plot(niche.sm, "D", "Niche overlap under similarity")


### species distribution modeling #####

#extracting data for multinomial analysis and dis.mod.
Dip <- subset(FCM, FCM$Ploidy1==2) #get the diploids
Tet <- subset(FCM, FCM$Ploidy1==4) #get the tetraploids
envdata2 <- extract(ENV.26, Dip[,c("Longitude", "Latitude")]) #extract climate data for diploids
envdata4 <- extract(ENV.26, Tet[,c("Longitude", "Latitude")]) #extract climate data for tetraploids
envdata2 <- data.frame(lon=Dip$Longitude, lat=Dip$Latitude, cyt=Dip$Ploidy1, envdata2) # make a data frame with coordinates and ploidy levels for diploids
envdata4 <- data.frame(lon=Tet$Longitude, lat=Tet$Latitude, cyt=Tet$Ploidy1, envdata4) # make a data frame with coordinates and ploidy levels for tetraploids 
envall <- merge(envdata2, envdata4, all=T) # merge the data frames to make one
envall$cyt <- factor(envall$cyt, labels = c("Dip", "Tet")) # factorize ploidy level

envdmdip <- data.frame(cbind(envdata2$lon, envdata2$lat, factor(envdata2$cyt)))
colnames(envdmdip) <- c("lon", "lat", "cyt")
envdmtet <- data.frame(cbind(envdata4$lon, envdata4$lat, factor(envdata4$cyt)))
colnames(envdmtet) <- c("lon", "lat", "cyt")

#extracting pseudo-absence points for diploid
envX2 <- circles(envdmdip[,c("lon", "lat")], d=65000, lonlat=T)
envbg2 <- spsample(envX2@polygons, 500, type = "random", iter=5000)
#extracting pseudo-absence points for tetraploid
envX4 <- circles(envdmtet[,c("lon", "lat")], d=100000, lonlat=T)
envbg4 <- spsample(envX4@polygons, 800, type = "random", iter=5000)


#extracting data for diploid
envdmdip_env <- extract(ENV.26, envdmdip[,c("lon", "lat")]) #diploid point data
envbg2_env <- extract(ENV.26, envbg2) #background data
#extracting data for tetraploid
envdmtet_env <- extract(ENV.26, envdmtet[,c("lon", "lat")]) #diploid point data
envbg4_env <- extract(ENV.26, envbg4) #background data


envdmdip_env <- data.frame(lon=envdmdip$lon, lat=envdmdip$lat, envdmdip_env)
envbgpoints2 <- envbg2@coords
colnames(envbgpoints2) <- c("lon", "lat")
envbg2_env <- data.frame(cbind(envbgpoints2, envbg2_env))
length(which(is.na(envbg2_env$bio2)))

envdmtet_env <- data.frame(lon=envdmtet$lon, lat=envdmtet$lat, envdmtet_env)
envbgpoints4 <- envbg4@coords
colnames(envbgpoints4) <- c("lon", "lat")
envbg4_env <- data.frame(cbind(envbgpoints4, envbg4_env))
length(which(is.na(envbg4_env$bio2)))



#modeling for diploid
envme_dip <- maxent(ENV, p=cbind (envdmdip_env$lon, envdmdip_env$lat), a=cbind(envbg2_env$lon, envbg2_env$lat))
#prediction
env_pred_me_dip <- predict(envme_dip, ENV)
#for the larger range
env_pred_me1_dip <- predict(envme_dip, ENV1)

#modeling for tetraploid
envme_tet <- maxent(ENV, p=cbind (envdmtet_env$lon, envdmtet_env$lat), a=cbind(envbg4_env$lon, envbg4_env$lat))
#prediction
env_pred_me_tet <- predict(envme_tet, ENV)
#for the larger range
env_pred_me1_tet <- predict(envme_tet, ENV1)

###testing (evaluation of the model)
e2 <- evaluate(p=cbind (envdmdip_env$lon, envdmdip_env$lat), a=cbind(envbg2_env$lon, envbg2_env$lat), envme_dip, ENV)
e4 <- evaluate(p=cbind (envdmtet_env$lon, envdmtet_env$lat), a=cbind(envbg4_env$lon, envbg4_env$lat), envme_tet, ENV)


#Plotting
plot(env_pred_me_dip, 1, cex=0.5, legend=T, mar=par("mar"), xaxt="n", yaxt="n", main="Predicted presence of Diploids (variables)")
plot(cnt1, fill=F, add=T, border="cornsilk")
#for the larger range data
plot(env_pred_me1_dip, 1, cex=0.5, legend=T, mar=par("mar"), xaxt="n", yaxt="n", main="Predicted presence of Diploids (variables)")
plot(cnt2, fill=F, add=T, border="cornsilk", box=F)

#tetraploid
plot(env_pred_me_tet, 1, cex=0.5, legend=T, mar=par("mar"), xaxt="n", yaxt="n", main="Predicted presence of Tetraploids (variables)")
plot(cnt1, fill=F, add=T, border="cornsilk")
#for the larger range data
plot(env_pred_me1_tet, 1, cex=0.5, legend=T, mar=par("mar"), xaxt="n", yaxt="n", main="Predicted presence of Tetraploids (variables)")
plot(cnt2, fill=F, add=T, border="cornsilk", box=F)

#response curves for each variable
response(envme_dip) #for diploids
response (envme_tet) # for tetraploids
