#---------------------#
# Fit age depth model #
#---------------------#

library(Bchron)

# Import data
dates <- read.csv("LakeHindDates.csv", header=TRUE, 
                  stringsAsFactors=FALSE)
#Depths and names of unit boundaries
UnitBounds <- c(0, 22, 30, 33, 43, 68)
BoundaryNames <- c("B2 to B3", "B1 to B2", "Sand", "YDB proxies",
                   "A2 to B1", "A1 to A2")

# Fit age-depth model
AgeDepthModel <- Bchronology(ages=dates$Mu,
                             ageSds=dates$SE,
                             positions=dates$Depth,
                             positionThicknesses=dates$Thickness,
                             calCurves=rep("intcal13", nrow(dates)),
                             ids=dates$LabNo,
                             outlierProbs=dates$OutlierProb,
                             predictPositions=seq(0, 70, by=0.5),
                             jitterPositions=TRUE)


# Extract age samples at unit boundaries
BoundaryPositions <- which(AgeDepthModel$predictPositions %in% UnitBounds)
BoundarySamples <- AgeDepthModel$thetaPredict[,BoundaryPositions]

# Find maximum age densities for both modelled boundaries and radiocarbon
# dates. This will be used to scale densities for plotting
maxbound <- max(sapply(1:ncol(BoundarySamples), function(i){
  max(density(BoundarySamples[,i])$y)
}))
maxdate <- max(sapply(1:nrow(dates), function(i){
  max(density(AgeDepthModel$theta[,i])$y)
}))
densscale <- max(c(maxbound, maxdate))/3
rm(maxbound, maxdate)

# Convert unit bound samples to long format densities for plotting
BoundarySamples <- lapply(1:ncol(BoundarySamples), function(i){
  d <- density(BoundarySamples[,i])
  ddf <- data.frame(dens=UnitBounds[i]-(d$y/densscale),
                    ymin=rep(UnitBounds[i], length(d$x)),
                    date=d$x, id=rep(BoundaryNames[i], length(d$x)),
                    stringsAsFactors=FALSE)
  return(ddf)
})
BoundarySamples <- do.call("rbind", BoundarySamples)

# Convert calibrated date samples to long format densities for plotting
DateSamples <- lapply(1:nrow(dates), function(i){
  
  d <- AgeDepthModel$calAges[[i]]
  wa <- which(d$densities>1e-5)
  
  y <- d$densities[wa]
  x <- d$ageGrid[wa]
  
  ddf <- data.frame(dens=dates$Depth[i]-(y/densscale),
                    ymin=rep(dates$Depth[i], length(x)),
                    date=x, id=rep(dates$LabNo[i], length(x)),
                    stringsAsFactors=FALSE)
  return(ddf)
})
DateSamples <- do.call("rbind", DateSamples)

# Convert modelled date samples to long format densities for plotting
DateSamples2 <- lapply(which(dates$OutlierProb<0.5), function(i){
  d <- density(AgeDepthModel$theta[,i])
  ddf <- data.frame(dens=dates$Depth[i]-(d$y/densscale),
                    ymin=rep(dates$Depth[i], length(d$x)),
                    date=d$x, id=rep(dates$LabNo[i], length(d$x)),
                    stringsAsFactors=FALSE)
  return(ddf)
})
DateSamples2 <- do.call("rbind", DateSamples2)


# Extract quantiles of modelled ages through the strata
q <- sapply(seq(0.98, 0.52, -0.02), function(x) c(x+0.01, 1-x-0.01))
PredInts <- lapply(1:ncol(q), function(x){
 xu <- sapply(1:ncol(AgeDepthModel$thetaPredict), function(z){
   quantile(AgeDepthModel$thetaPredict[,z], probs=q[1,x])
 })
 xl <- rev(sapply(1:ncol(AgeDepthModel$thetaPredict), function(z){
   quantile(AgeDepthModel$thetaPredict[,z], probs=q[2,x])
 }))
  
 df <- data.frame(Qu=rep(q[1,x]-q[2,x], length(xu)*2),
                  x=c(xu, xl),
                  y=c(AgeDepthModel$predictPositions,
                      rev(AgeDepthModel$predictPositions)))
 return(df)
})

PredInts <- do.call("rbind", PredInts)
PredInts <- PredInts[which(PredInts$y<max(UnitBounds)+1),]


# Obtain 95% Intervals for each boundary
BoundaryInts <- sapply(1:ncol(AgeDepthModel$thetaPredict), function(x){
  q <- quantile(AgeDepthModel$thetaPredict[,x], probs=c(0.975, 0.025))
  q <- round(q, 0)
  return(paste0(q[1],"-",q[2]))
})[BoundaryPositions]

# Obtain 95% Intervals for modelled dates
ModelledInts <- sapply(1:ncol(AgeDepthModel$theta), function(x){
  q <- quantile(AgeDepthModel$theta[,x], probs=c(0.975, 0.025))
  q <- round(q, 0)
  return(paste0(q[1],"-",q[2]))
})

  
#---------------------#
#---- Plot results ---#
#---------------------#
library(ggplot2)
library(patchwork)


p1 <- ggplot()+
  annotate("rect", xmin=rep(-Inf, 3), xmax=rep(Inf, 3),
           ymin=UnitBounds[c(2, 5, 6)], ymax=UnitBounds[c(1, 2, 5)], 
           alpha=0.2,
           fill=c("grey", "deepskyblue2", "palegoldenrod"))+
  annotate("text", label="a", x=max(DateSamples$date)-100, y=-1,
           hjust=0, vjust=1, size=5)+
  annotate("segment", x=rep(-Inf, 2), xend=rep(Inf, 2),
           y=UnitBounds[c(3,4)], yend=UnitBounds[c(3,4)],
           size=1, color="deepskyblue2", linetype="dashed")+
  annotate("text", label=paste(BoundaryInts, "BP (95% intervla)"), 
           y=UnitBounds+c(0,0,-1,1,0,0), x=rep(19000,6), 
           hjust=0, size=2.5, vjust=c(0.5, 0.5, 0, 1, 0.5, 0.5),
           color=c("blue4", "blue4", "darkred", "darkred",
                   "blue4", "blue4"))+
  geom_polygon(data=PredInts, aes(x=x, y=y, group=Qu), alpha=0.04)+
  geom_ribbon(data=DateSamples, aes(group=id, ymin=ymin, 
                                    ymax=dens, x=date),
              alpha=0.2, fill="purple")+
  geom_ribbon(data=DateSamples2, aes(group=id, ymin=ymin, 
                                    ymax=dens, x=date),
              alpha=0.8, fill="purple")+
  annotate("text", label=c("Subunit B2", "Upper Subunit B1", 
                           "YDB", "Lower Subunit B1", 
                           "Unit A"), 
           x=rep(20600, 5), y=c(20, 24, 31.5, 41, 45),
           color=c("blue4", "blue4", "darkred",
                   "blue4", "blue4"),
           hjust=0)+
  scale_x_reverse(expand=c(0,0))+
  scale_y_reverse(expand=c(0.01, 0), breaks=UnitBounds)+
  labs(x="Cal yr BP", y="Depth (cm)")+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(), 
        axis.ticks=element_line(color="grey"),
        axis.line=element_line(color="grey"))


# Subsample boundary date samples for 3-cm interval
Boundary_Dens <- BoundarySamples[which(BoundarySamples$id %in%
                                         BoundaryNames[3:4]),]
Boundary_Dens$dens <- abs(Boundary_Dens$dens - Boundary_Dens$ymin)
Boundary_Dens$ymin <- 0

# Subsample radiocarbon date samples for 3-cm interval
Date_Dens <- DateSamples[which(DateSamples$id=="UCIAMS-29317"),]
Date_Dens$dens <- abs(Date_Dens$dens - Date_Dens$ymin)
Date_Dens$ymin <- 0

Date_Dens2 <- DateSamples2[which(DateSamples2$ymin > 30 &
                                   DateSamples2$ymin < 33),]
Date_Dens2$dens <- abs(Date_Dens2$dens - Date_Dens2$ymin)
Date_Dens2$ymin <- 0

ymax <- max(Date_Dens$dens)

p2 <- ggplot(data=NULL, aes(ymin=0, ymax=dens, x=date))+
  annotate("rect", ymin=0, ymax=ymax*1.2, xmin=12735, xmax=12835,
           fill="orange", alpha=0.4)+
  annotate("text", label="b", x=12855, y=ymax*1.15, 
           hjust=0, vjust=1, size=5)+
  annotate("text", label="Hypothesized", color="orange", 
           x=12825, y=ymax*0.6, angle=90)+
  annotate("text", label="Younger Dryas", color="orange", 
           x=12805, y=ymax*0.6, angle=90)+
  annotate("text", label="Boundary", color="orange", 
           x=12785, y=ymax*0.6, angle=90)+
  annotate("text", label="12835-12735", color="orange", 
           x=12765, y=ymax*0.6, angle=90)+
  annotate("text", label="cal yr BP", color="orange", 
           x=12745, y=ymax*0.6, angle=90)+
  annotate("text", label=paste("YDB bounds:", BoundaryInts[4], "to",
                                BoundaryInts[3], "cal yr BP (95% intervals)"),
           color="blue", x=12500, y=ymax*1.05, vjust=0, size=3)+
  annotate("text", label="Modelled ages (95% intervals)", x=12320, 
           y=ymax*0.9, size=3, color="purple")+
  annotate("text", label=paste(dates$LabNo[13], ModelledInts[13], 
                               "cal yr BP"), x=12320, 
           y=ymax*0.74, size=3, color="purple", vjust=0)+
  annotate("text", label=paste(dates$LabNo[14], ModelledInts[14], 
                               "cal yr BP"), x=12320, 
           y=ymax*0.72, size=3, color="purple", vjust=1)+
  geom_ribbon(data=Boundary_Dens, aes(group=id), fill="blue", alpha=0.5)+
  geom_ribbon(data=Date_Dens2, aes(group=id), fill="purple", alpha=0.5)+
  geom_ribbon(data=Date_Dens, aes(group=id), fill="purple", alpha=0.2)+
  scale_x_reverse(expand=c(0,0))+
  scale_y_continuous(expand=c(0.01, 0))+
  labs(x="Cal yr BP")+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(), 
        axis.ticks.x=element_line(color="grey"),
        axis.line.x=element_line(color="grey"),
        axis.line.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank())

p1 + p2 + plot_layout(ncol=1, heights=c(1, 0.3))
