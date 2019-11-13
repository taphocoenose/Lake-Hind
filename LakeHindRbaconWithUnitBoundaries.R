
# Lasted Edited by Ryan Breslawski on Nov 11, 2019
# email: rbreslawski@smu.edu
# 
# This script was last edited in R 3.6.1 on a Windows 10 machine

library(rbacon)
library(ggplot2)
library(patchwork)

# Fix randomization for reproducibility
set.seed(42334532)

#---------------------#
# Fit age depth model #
#---------------------#

#Depths and names of unit boundaries
UnitBounds <- c(0, 22, 30, 33, 43, 68)
BoundaryNames <- c("B2 to B3", "B1 to B2", "Sand", "YDB proxies",
                   "A2 to B1", "A1 to A2")

# Prep data for Bacon
d <- read.csv("LakeHindDates.csv", header=TRUE, stringsAsFactors=FALSE)
dir.create("Bacon_runs")
dir.create("Bacon_runs/LakeHindDates")
write.csv(d[,1:(ncol(d)-2)], file="Bacon_runs/LakeHindDates/LakeHindDates.csv",
          row.names=FALSE)

# Fit age-depth model. Increase iterations above 1000 default
# to improve stability of modelled ages across simulation runs.
Bacon(core="LakeHindDates", d.min=0, d.max=70, d.by=0.5, cc=1, 
      MaxAge=40000, ssize=10000, rotate.axes=TRUE, 
      boundary=UnitBounds[2:5], thick=1)

# Extract age samples at unit boundaries
BoundarySamples <- sapply(UnitBounds, function(x) Bacon.Age.d(x))

# Find maximum age densities for both modelled dates and unmodelled
# radiocarbon dates. This will be used to scale densities for plotting
densscale <- max(sapply(1:length(info$calib$probs), function(i){
  s <- sample(info$calib$probs[[i]][,1], size=1e5, 
              prob=info$calib$probs[[i]][,2], replace=TRUE)
  max(density(s, adjust=1.5)$y)
}))/10

# Convert unit bound samples to long format densities for plotting
BoundarySamples <- lapply(1:ncol(BoundarySamples), function(i){
  d <- density(BoundarySamples[,i], adjust=1.5)
  ddf <- data.frame(dens=UnitBounds[i]-(d$y/densscale),
                    ymin=rep(UnitBounds[i], length(d$x)),
                    date=d$x, id=rep(BoundaryNames[i], length(d$x)),
                    stringsAsFactors=FALSE)
  return(ddf)
})
BoundarySamples <- do.call("rbind", BoundarySamples)

# Convert calibrated date samples to long format densities for plotting.
DateSamples <- lapply(1:nrow(info$dets), function(i){
  
  d <- sample(info$calib$probs[[i]][,1], size=1e5, 
              prob=info$calib$probs[[i]][,2], replace=TRUE)
  wa <- density(d, adjust=1.5)
  y <- wa$y
  x <- wa$x
  
  ddf <- data.frame(dens=info$dets$depth[i]-(y/densscale),
                    ymin=rep(info$dets$depth[i], length(x)),
                    date=x, id=rep(info$dets$labID[i], length(x)),
                    stringsAsFactors=FALSE)
  return(ddf)
})
DateSamples <- do.call("rbind", DateSamples)


# Extract quantiles of calibrated ages through the strata
q <- sapply(seq(0.98, 0.52, -0.02), function(x) c(x+0.01, 1-x-0.01))
s <- sapply(info$depths, function(z) Bacon.Age.d(z))
PredInts <- lapply(1:ncol(q), function(x){
  
 xu <- sapply(1:ncol(s), function(z) quantile(s[,z], probs=q[1,x]))
 xl <- rev(sapply(1:ncol(s), function(z) quantile(s[,z], probs=q[2,x])))
  
 df <- data.frame(Qu=rep(q[1,x]-q[2,x], length(xu)*2),
                  x=c(xu, xl),
                  y=c(info$depths, rev(info$depths)))
 return(df)
})
rm(s)

PredInts <- do.call("rbind", PredInts)
PredInts <- PredInts[which(PredInts$y<max(UnitBounds)+1),]


# Obtain 95% Intervals for each boundary
BoundaryInts <- sapply(UnitBounds, function(x){
  q <- quantile(Bacon.Age.d(x), probs=c(0.975, 0.025))
  q <- round(q, 0)
  return(paste0(q[1],"-",q[2]))
})

# Obtain 95% Intervals for modelled dates
ModelledInts <- sapply(1:nrow(info$dets), function(i){
  s <- sample(info$calib$probs[[i]][,1], size=1e5, 
              prob=info$calib$probs[[i]][,2], replace=TRUE)
  q <- quantile(s, probs=c(0.975, 0.025))
  q <- round(q, 0)
  return(paste0(q[1],"-",q[2]))
})

  
#---------------------#
#---- Plot results ---#
#---------------------#

# Age-depth model plot
p1 <- ggplot()+
  annotate("rect", xmin=rep(-Inf, 3), xmax=rep(Inf, 3),
           ymin=UnitBounds[c(2, 5, 6)], ymax=UnitBounds[c(1, 2, 5)], 
           alpha=0.2,
           fill=c("grey", "deepskyblue2", "palegoldenrod"))+
  annotate("text", label="a", x=max(DateSamples$date)-100, y=-1,
           hjust=0, vjust=1, size=5)+
  annotate("segment", x=rep(-Inf, 2), xend=rep(Inf, 2),
           y=UnitBounds[c(3,4)], yend=UnitBounds[c(3,4)],
           size=0.7, color="deepskyblue2", linetype="dashed")+
  annotate("text", label="Unit boundary ages (95% intervals)", 
           x=18100, y=-2, vjust=0, size=3, color="blue4")+
  annotate("text", label=paste(BoundaryInts, "cal yr BP"), 
           y=UnitBounds+c(0,0,-.5,.5,0,0), x=rep(18100,6), 
           size=2.7, vjust=c(0.5, 0.5, 0, 1, 0.5, 0.5),
           color=rep("blue4", 6))+
  geom_polygon(data=PredInts, aes(x=x, y=y, group=Qu), alpha=0.04)+
  geom_ribbon(data=DateSamples, aes(group=id, ymin=ymin, 
                                    ymax=dens, x=date),
              alpha=0.75, fill="purple")+
  annotate("text", label=c("Subunit B2", "Upper Subunit B1", 
                           "Middle Subunit B1", "Lower Subunit B1", 
                           "Unit A"), 
           x=rep(21000, 5), y=c(21, 23, 31.5, 42, 44),
           color=rep("blue4", 5),
           hjust=0, size=3.2)+
  scale_x_reverse(expand=c(0,0))+
  scale_y_reverse(expand=c(0.01, 0), breaks=UnitBounds, limits=c(69, -3))+
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
Date_Dens <- DateSamples[which(DateSamples$id %in% c("UCIAMS-88701",
                                                     "PSUAMS-1572",
                                                     "UCIAMS-29317")),]
Date_Dens$dens <- abs(Date_Dens$dens - Date_Dens$ymin)
Date_Dens$ymin <- 0

# Truncate Date_Dens density for plotting
Date_Dens <- Date_Dens[which(Date_Dens$date>12220),]

# Y-axis max for plot 2
ymax <- max(c(Date_Dens$dens, Boundary_Dens$dens))

# Middle Subunit B1 plot
p2 <- ggplot(data=NULL, aes(ymin=0, ymax=dens, x=date))+
  annotate("rect", ymin=0, ymax=Inf, xmin=12735, xmax=12835,
           fill="orange", alpha=0.4)+
  annotate("text", label="b", x=12855, y=ymax*1.15, 
           hjust=0, vjust=1, size=5)+
  annotate("text", label="Hypothesized", color="orange3", 
           x=12833, y=ymax, size=3.5, hjust=0)+
  annotate("text", label="Younger Dryas", color="orange3", 
           x=12833, y=ymax*0.8, size=3.5, hjust=0)+
  annotate("text", label="Boundary", color="orange3", 
           x=12833, y=ymax*0.6, size=3.5, hjust=0)+
  annotate("text", label="12835-12735", color="orange3", 
           x=12833, y=ymax*0.4, size=3.5, hjust=0)+
  annotate("text", label="cal yr BP", color="orange3", 
           x=12833, y=ymax*0.2, size=3.5, hjust=0)+
  annotate("text", label=paste("Middle Subunit B1 bounds:", 
                               BoundaryInts[4], "to",
                                BoundaryInts[3], 
                               "cal yr BP (95% intervals)"),
           color="blue", x=12480, y=ymax*1.05, vjust=0, size=3)+
  annotate("text", label=paste(info$dets$labID[15], ModelledInts[15], 
                               "cal yr BP (95% interval)"), x=12520, 
           y=ymax*(-0.05), size=2.2, color="purple", vjust=1, hjust=0)+
  annotate("text", label=paste(info$dets$labID[13], ModelledInts[13], 
                               "cal yr BP (95% interval)"), x=12470, 
           y=ymax*(-0.18), size=2.2, color="purple", vjust=1, hjust=0)+
  annotate("text", label=paste(info$dets$labID[14], ModelledInts[14], 
                               "cal yr BP (95% interval)"), x=12800, 
           y=ymax*(-0.05), size=2.2, color="purple", vjust=1, hjust=0)+
  geom_ribbon(data=Boundary_Dens, aes(group=id), fill="blue", alpha=0.5)+
  geom_ribbon(data=Date_Dens, aes(group=id), fill="purple", alpha=0.5)+
  scale_x_reverse(expand=c(0,0))+
  scale_y_continuous(expand=c(0.01, 0), limits=c(ymax*(-0.29), ymax*1.2))+
  labs(x="Cal yr BP")+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(), 
        axis.ticks.x=element_line(color="grey"),
        axis.line.x=element_line(color="grey"),
        axis.line.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank())

# Plot combined panels
p <- p1 + p2 + plot_layout(ncol=1, heights=c(1, 0.3))

ggsave(filename="FigS2.pdf", plot=p, device="pdf", units="in",
       width=7, height=6)
