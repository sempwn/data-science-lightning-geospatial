############################################################################################
###########    Fitting of presence/absence data with                      ##################
###########                            R-INLA package                     ##################
###########    Ref: Wnv Kaggle competition                                ##################
###########                                                               ##################
###########    Written by: Mike Irvine                                    ##################
###########                                                               ##################
###########    Last updated: 25/11/17                                     ##################
############################################################################################
#if running for first time, uncomment following to install relevent packages
#install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/stable"); install.packages("gridExtra"); install.packages("sp"); install.packages("geoR");
#install.packages("ggmap");
#load libraries
#install.packages(c("ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library(gridExtra); library(ggplot2); library(lattice); library(INLA); library(splancs); library(fields); library(gdata);
library(MASS);
library(tidyverse);
library(lubridate);
theme_set(theme_bw())
library(sf)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir) #change working directory to model file path.


data <- read_csv("../data/seals.csv")
data <- data %>% 
    mutate(Count = as.integer(Count)) %>%
    filter(Longitude < 0.) %>%
    drop_na() %>%
    mutate(year = year(Date))


p <- ggplot() + geom_point(aes(x=data$Longitude,y=data$Latitude,colour=data$Count), size=2,
                           alpha=0.5) + scale_colour_gradientn(colours=tim.colors(100),trans="log10") +
  labs(list(color = "Seals", x = expression("Longitude ("*degree*")"), y = expression("Latitude ("*degree*")")))
print(p)
#ggsave(file="villages-mf_ml.pdf")

# load nautral Earth
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

plot_seal_data <- function(data){
  p <- ggplot(data = world) +
    geom_sf() +
    geom_point(data = data, aes(x=Longitude,y=Latitude,colour=Count), size=2,
                 alpha=0.5) + scale_colour_gradientn(colours=tim.colors(100),trans="log10") +
    labs(list(color = "Seals", x = expression("Longitude ("*degree*")"), y = expression("Latitude ("*degree*")"))) +
    coord_sf(xlim = c(-126, -122), ylim = c(48., 51.), expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.4) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering) 
  return(p)

}
p <- plot_seal_data(data)
print(p)
# coords <- as.matrix(dat[,8:9]) 

# just use data from 2014
data <- data %>% filter(year == 2014)

coords <- data %>% select(c(Longitude,Latitude)) %>%
  as.matrix

# 
Y <- data$Count
year <- data$year

# #create mesh with less triangles
prmesh <- inla.mesh.create(coords)
plot(prmesh, asp=1, main="")
points(coords[,1], coords[,2], pch=19, cex=.5, col="red")

# #create observation matrix based on mesh
A <- inla.spde.make.A(prmesh, loc=coords) 
spde <- inla.spde2.matern(prmesh, alpha=2)
# 
# #INLA stack
mesh.index <- inla.spde.make.index(name="field", n.spde=spde$n.spde)
effects <- list(
  c(mesh.index,list(Intercept=1)),
  covar=list(
    long=inla.group(coords[,1]),
    lat=inla.group(coords[,2])
  )
)
stk.dat <- inla.stack(
  data=list(y=Y), A=list(A,1), tag="est",
  effects=effects
)

f.s <- y ~ -1 + Intercept + long + lat + f(field, model=spde)
r.s <- inla(f.s, family="Poisson", data=inla.stack.data(stk.dat), verbose=TRUE,
            control.predictor=list(A=inla.stack.A(stk.dat), compute=TRUE),
            control.fixed = list(expand.factor.strategy='inla'),
            control.compute = list(dic=T))
r.s$summary.fixed
r.s$summary.hyperpar


# mm <- max(r.s$summary.random[[1]]$mean)
# vv <- max(r.s$summary.random[[1]]$sd)
# pp <- 2*(1 - pnorm(abs(mm/vv)))
# print(cat("max breeding site coefficient mean : ", mm," sd: ", vv ," p-value : ",pp,"\n"))
# r.f <- inla.spde2.result(r.s, "field", spde, do.transf=TRUE)
# inla.emarginal(function(x) x, r.f$marginals.variance.nominal[[1]])
# plot.default(r.f$marginals.variance.nominal[[1]], type="l",
#              xlab=expression(sigma[x]^2), ylab="Density") 
# p <- ggplot() + geom_path(aes(x=r.f$marginals.range.nominal[[1]][,1],y=r.f$marginals.range.nominal[[1]][,2])) +
#   labs(list(x = "Practical range (km)", y = "Density"))
# print(p)
# #ggsave("practical_range_mf.pdf")
# 
# par(mfrow=c(2,3), mar=c(3,3.5,0,0), mgp=c(1.5, .5, 0), las=0) 
# plot(r.s$marginals.fix[[1]], type='l', xlab='Intercept', ylab='Density') 
# plot(r.s$summary.random[[1]][,1:2], type='l',
#      xlab='Distance to breed site (Km)', ylab='Coeficient'); abline(h=0, lty=3)
# for (i in c(4,6)) lines(r.s$summary.random[[1]][,c(1,i)], lty=2) 
# plot(r.s$marginals.hy[[1]], type='l', ylab='Density', xlab=expression(phi)) 
# plot.default(inla.tmarginal(function(x) 1/exp(x), r.s$marginals.hy[[2]]),type='l', xlab=expression(kappa), ylab='Density') 
# plot.default(r.f$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[x]^2), ylab='Density')
# plot.default(r.f$marginals.range.nominal[[1]], type='l', xlab='Practical range', ylab='Density')
# par(mfrow=c(1,1))
# 
# ##################################################################
# ###############       Kriging procedure    #######################
# ##################################################################
nxy <- c(100,100)
projgrid <- inla.mesh.projector(prmesh, dims=nxy)

thull <- chull(x=coords[,1],y=coords[,2])
thull <- c(thull, thull[1])
PRborder <- coords[thull,]

xy.in <- inout(projgrid$lattice$loc, cbind(PRborder[,1], PRborder[,2]))

coord.prd <- projgrid$lattice$loc[xy.in,]
plot(coord.prd, type="p", cex=.1); lines(PRborder)
points(coords[,1], coords[,2], pch=19, cex=.5, col="red")
A.prd <- projgrid$proj$A[xy.in, ]

ef.prd = list(c(mesh.index,list(Intercept=1)), list(long=inla.group(coord.prd[,1]),
                                                    lat=inla.group(coord.prd[,2])))
stk.prd <- inla.stack(data=list(y=NA), A=list(A.prd,1), tag="prd", effects=ef.prd)
stk.all <- inla.stack(stk.dat, stk.prd)
r2.s <- inla(f.s, family="Poisson", data=inla.stack.data(stk.all),
             control.predictor=list(A=inla.stack.A(stk.all), compute=TRUE,link = 1),
             quantiles=NULL, control.results=list(return.marginals.random=F, return.marginals.predictor=F),
             verbose=TRUE,
             control.inla = list(int.strategy = "eb"),
             control.fixed = list(expand.factor.strategy='inla'),
             control.compute = list(dic=T))
id.prd <- inla.stack.index(stk.all, "prd")$data
sd.prd <- m.prd <- matrix(NA, nxy[1], nxy[2])
m.prd[xy.in] <- r2.s$summary.fitted.values$mean[id.prd]
sd.prd[xy.in] <- r2.s$summary.fitted.values$sd[id.prd]
grid.arrange(levelplot(m.prd, col.regions=tim.colors(99), xlab="", ylab="", main="mean",
                       scales=list(draw=FALSE)), levelplot(sd.prd, col.regions=topo.colors(99),
                                                           xlab="", ylab="", scales=list(draw=FALSE), main="standard deviation"))

levelplot(m.prd[40:60,18:35], col.regions=tim.colors(99), xlab="", ylab="", main="mean",
          panel=function(...){
            panel.levelplot(...)
            grid.points(as.vector(20*(coords[,1]-142.679)/0.023), as.vector(20*(coords[,2]+3.581)/0.012), pch=2)
          }
)
levelplot(m.prd, col.regions=tim.colors(99), xlab="", ylab="", main="mean",
          panel=function(...){
            panel.levelplot(...)
            grid.points(as.vector(20*(coords[,1]-142.679)/0.023), as.vector(20*(coords[,2]+3.581)/0.012), pch=2)
          }
)

image.plot(projgrid$x,projgrid$y,m.prd,col=tim.colors(99),legend.lab='mean',xlab='Longitude',ylab='Latitude'); #points(coords[,1],coords[,2]);


x <- matrix(projgrid$x,nrow=100,ncol=length(projgrid$x),byrow=TRUE)
y <- matrix(projgrid$y,nrow=100,ncol=length(projgrid$y),byrow=FALSE)
df<- data.frame(x=as.vector(t(x)),y=as.vector(t(y)),m = as.vector(as.matrix(m.prd)), sd = as.vector(as.matrix(sd.prd)))
df[is.na(df$m),]$m <- 0.
df[is.na(df$sd),]$sd <- max(df$sd,na.rm=TRUE)+0.001

p<- ggplot(data = world) +
  geom_raster(data=df, aes(x, y, fill = m, alpha = m)) +
  geom_sf() +
  scale_fill_gradientn(name="Seals",colours = topo.colors(20),limits=c(0.0,400.0))+
  guides(fill = guide_colorbar()) +
  geom_tile(data=df, aes(x,y,alpha=sd)) +
  scale_alpha("uncertainty",range=c(0.0,1.0),limits=c(0.,120.)) + 
  #scale_alpha("uncertainty", limits=c(0.09,16.), guide = 'none',range=c(0.0,1.0),trans="reverse") +
  scale_x_continuous(name=expression(paste("Longitude (",degree,")")),
                     expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude (",degree,")")),
                     expand=c(0,0)) +
  coord_sf(xlim = c(-126, -122), ylim = c(48., 51.), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) 

print(p)
# ggsave('spatial_fit.pdf')