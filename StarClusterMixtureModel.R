# May 8, 2018
# 
# ================================================================================
# Based on code written for 2014ApJ...787..107K.
# Title: The Spatial Structure of Young Stellar Clusters. I. Subclusters
# Authors: Michael A. Kuhn, Eric D. Feigelson, Konstantin V. Getman, Adrian J. Baddeley, 
# Leisa K. Townsley, Patrick S. Broos, Alison Sills, Matthew R. Bate, Matthew S. Povich, 
# Kevin L. Luhman, Heather A. Busk, Tim Naylor, Robert R. King
# ================================================================================
# 
# Code names: StarClusterMixtureModel.R
# Language: R
# Code tested under the following compilers/operating systems: 
#      R version 3.3.1 (2018-04-29) -- "Bug in Your Hair" on Mac OS 10.9.5
# 
# Description of input and output: 
#    Input and output of the various R functions in this library are described below. Here, 
#    we give a general overview of the input and output of the fitting recipe, described
#    in the article, using these functions.
#    Input: Machine-readable-table versions of Tables 2 and 3 and machine-readable-table
#           "data behind the figure" for Figure 2
#    Output: Fitted subcluster parameters are saved in the varable "ellipse.out"
#            Surface density maps of clusters with subclusters ellipses overplotted as shown
#            in Figure 2. Colorscale units of observed stars per square parsec. Axis units of
#            parsecs.
#            Model residual maps with subcluster ellipses overplotted as shown in 
#            Figure 4 (right panel). Colorscale units of observed stars per square parsec.
#            Axis units of parsecs.
#
# System requirements: R version 2.15.2 or later. Spatstat version 1.13-1 or later.
# 
# Calls to external routines: 
#    This code makes use of spatstat functions owin, ppp, unique, as.im, integral.im,
#    eval.im, ppm, plot.im, plot.ppp, and diagnose.ppm.
#    This code makes use of plotrix function draw.ellipse.
#    This code makes use of numerous Base R functions (in particular the function optim).
# 
# List of Functions in this library:
#    ell.model
#    const.model
#    multi.model
#    model.lik
#    model.ppm
#    draw.model
#    param2ellipse
#    ellipse2param
#    mask.freeze
#    mask.model
#    mask.adjust_positions
#    mask.adjust_rotations
#    mask.adjust_size
#    mask.adjust_shape
#    mask.adjust_mix
#    ds9.colors
#    american.colors
#    ad2xy.tangent
#    xy2ad.tangent
#    star.ppp
#    make.fig2
#    make.fig4
# 
# Additional comments: The library R functions defined in sourcecode.R will be accessed by the
# finite mixture model fitting example in the article.
# 
# =====================================================================================================
# The AAS gives permission to anyone who wishes to use these subroutines to run their own calculations.
# Permission to republish or reuse these routines should be directed to permissions@aas.org.
#
# Note that the AAS does not take responsibility for the content of the source code. Potential users should
# be wary of applying the code to conditions that the code was not written to model and the accuracy of the
# code may be affected when compiled and executed on different systems.
# =====================================================================================================



# 
# ell.model (The isothermal ellipsoid surface density distribution)
# -----------------------------------------------------------------
# The inputs are x and y coordinates and a vector "param" containing the (x0,y0) coordinates of the subcluster
# center, the core radius of the isothermal ellipsoid, the orientation of the ellipse in radians, and the
# ratio of semi-major and semi-minor ellipse axes.
# The output is a value of the model (not normalized) at the location (x,y).
# Usage> ell.model(0,1,param=c(0.0, 0.0, 1.0, pi/4, 0.5))
#
ell.model <- function(x1,y1,param=param,param2=NULL) {
  x0=param[1]
  y0=param[2]
  core=param[3]
  theta=param[4]
  b=param[5]
  rr <- sqrt(((x1-x0)*cos(theta)+(y1-y0)*sin(theta))^2+(-(x1-x0)*sin(theta)+(y1-y0)*cos(theta))^2/b^2)
  1.0/(1.0+rr^2/core^2)
}

# 
# const.model (The constant surface density distribution)
# -------------------------------------------------------
# The inputs are x and y coordinates.
# The output is the value 1 regardless of the input coordinates, representing a constant surface density.
# Usage> const.model(0,1,param2=NULL)
#
const.model <- function(x1,y1,param=NULL,param2=param2) { 1.0 }

# 
# multi.model (The finite mixture model surface density distribution)
# -------------------------------------------------------------------
# The inputs are x and y coordinates, a vector of model parameters, and a vector describing the 
# component models to be used. The vector of model parameters "param" is the concatination of 
# the model parameters for each component model. After the second, third, fourth (and so on) 
# model parameters, there is a mixture coefficient indicating the weight of this model relative
# to the first model component. I.e. For a mixture component for the second component 
# log mix = -1, the peak surface density of component 2 would be one tenth the peak surface 
# density of component 1.
# The vector describing the component models "param2" starts with the total number of components, 
# followed by the function names and the number of parameters to pass to each function.
# The output is a value of the model (not normalized) at the location (x,y).
# Usage> multi.model(0,1,param2=c(0.0, 0.0, 1.0, pi/4, 0.5, 2.0, 2.0, 0.5, pi/6, 1.0, 0.0),
#        param2=c(2,ell.model,5,ell.model,5))
#
multi.model <- function(x1,y1,param=param,param2=param2) {
 n_models <- param2[[1]]
 p1 <- 1
 p2 <- 2
 n1 <- length(param)
 n2 <- length(param2)
 output <- 0.0

 i <- 1
 while (i < n_models+1) {

   par1 <- param[p1:n1]
   par2 <- param2[p2:n2]
   model <- par2[[1]]
   n_param <- par2[[2]]

   param_new=NULL
   if (n_param > 0) { param_new=par1[1:n_param] }
   scale <- 1.0
   if (i > 1) {
     scale <- 10.0^(par1[[n_param+1]])
     p1 <- p1 + 1
   }

   output <- output + scale*model(x1,y1,param=param_new)

   p2 <- p2 + 2
   p1 <- p1 + n_param
   i <- i+1

 }
 output
}



#
# model.lik (Returns the negative of the log likelihood)
# ------------------------------------------------------
# The inputs are the "param" and "param2" vectors (described in multi.model), the
# model to be used, and the spatial point pattern data as a ppp object "clust".
# The output is the negative of the log likelihood for the model. (A negative value 
# is given so that optim can find the maximum likelihood by minimization.)
# Usage> loglik <- -model.lik(param,model=model.multi,clust=clust,param2=param2)
#
model.lik <- function(param,model=model,clust=clust,param2=NULL) {
  A <- as.im(model,clust$win,param=param,param2=param2)
  norm <- integral.im(A)
  result <- sum(log(model(clust$x,clust$y,param=param,param2=param2)/norm*clust$n))-clust$n
  -result
}

#
# model.im (Creates a pixelated image of the model)
# ------------------------------------------------------
# The inputs are the "param" and "param2" vectors (described in multi.model), the
# model to be used, and the spatial point pattern data as a ppp object "clust".
# The ouput is a pixilated array of values of the model as a spatstat "image" 
# object.
# Usage> plot(model.im(param,model=model.multi,clust=clust,param2=param2))
#
model.im <- function(param,model=model,clust=clust,param2=NULL) {
  A <- as.im(model,clust$win,param=param,param2=param2)
  norm <- integral.im(A)
  nnn <- clust$n
  eval.im(A/norm*nnn)
}

#
# model.im (Runs the spatstat "ppm" function on a model to scale it to point process data)
# ------------------------------------------------------
# The inputs are the "param" and "param2" vectors (described in multi.model), the
# model to be used, and the spatial point pattern data as a ppp object "clust".
# The output is a point process model of spatstat class "ppm".
# Usage> fmm.fit <- model.ppm(param,model=model.multi,clust=clust,param2=param2)
#
model.ppm <- function(param,model=model,clust=clust,param2=NULL) {
  A <- as.im(model,clust$win,param=param,param2=param2)
  ppm(clust, ~offset(log(B)), covariates = list(B=A))
}


#
# draw.model (Overplots model ellipses on an image of )
# ------------------------------------------------------
# The inputs are the "param" and "param2" vectors (described in multi.model), 
# the size of the ellipse in core radii, the color of the ellipse, and the 
# width of the ellipse.
# There is no output. This function causes ellipses to be plotted in the graphics
# window.
# Usage> draw.ellipse(param=param,param2=param2,size=1,border=2,lwd=2)
#
draw.model <- function(param=param,param2=param2,size=1,border=2,lwd=2) {
 n_models <- param2[[1]]
 p1 <- 1
 p2 <- 2
 n1 <- length(param)
 n2 <- length(param2)
 i <- 1
 while (i < n_models+1) {
   par1 <- param[p1:n1]
   par2 <- param2[p2:n2]
   model <- par2[[1]]
   n_param <- par2[[2]]
   param_new=NULL
   if (n_param > 0) { param_new=par1[1:n_param] }
   scale <- 1.0
   if (i > 1) {
     scale <- 10.0^(par1[[n_param+1]])
     p1 <- p1 + 1
   }
   if (n_param == 5) {
      draw.ellipse(param_new[1],param_new[2],a=param_new[3]*size,b=param_new[3]*param_new[5]*size,
      angle=param_new[4]*180.0/pi %% 360.0,deg=T,border=border,lwd=lwd)   
   }
   if (n_param == 3) {
      draw.ellipse(param_new[1],param_new[2],a=param_new[3]*size,b=param_new[3]*size,angle=0.0,border=border,lwd=lwd)   
   }
   p2 <- p2 + 2
   p1 <- p1 + n_param
   i <- i+1
 }
 # No output
}

#
# param2ellipse (Converts an array of parameters into an R data frame)
# -------------------------------------------------------------------
# The input is an array of parameters in the variable "param" as described
# for the function multi.model. The assumed model form is
# for a model with "param2" of c(<n+1>,model.ell,5, ... , model.ell,5,model.const,0)
# The output is a data.frame object with columns, x, y, core, theta, b, and mix. 
# Usage> ellipse <- param2ellipse(param)
# 
param2ellipse <-function(param){
  n <- length(param)/6 # i.e. six parameters for ellipses, minus one for the first mix, plus one for constant component
  ellipsoid <- data.frame(x=rep(NA,n),y=rep(NA,n),core=rep(NA,n),theta=rep(NA,n),b=rep(NA,n),mix=rep(NA,n))
  adjust <- 0 # compensate for missing mix[1]
  for (i in 1:n){
    ellipsoid[i,1] <- param[6*(i-1)+1-adjust]
    ellipsoid[i,2] <- param[6*(i-1)+2-adjust]
    ellipsoid[i,3] <- param[6*(i-1)+3-adjust]
    ellipsoid[i,4] <- param[6*(i-1)+4-adjust]
    ellipsoid[i,5] <- param[6*(i-1)+5-adjust]
    if (i > 1) {ellipsoid[i,6] <- param[6*(i-1)+6-adjust]}	
    if (i == 1) {ellipsoid[i,6] <- 0.0}	
    adjust <- 1
  }
ellipsoid
}

#
# ellipse2param (Converts an R data frame into an array of parameters)
# -------------------------------------------------------------------
# The input is a data.frame object with columns, x, y, core, theta, b, 
# and mix, describing a collection of ellipsoids. The variable
# "const_mix" is the log mixing parameter for the constant component.
# The output is an array of parameters in the format of the "param" argument
# to the function multi.model. The assciated "param2" variable for this 
# "param" array has the form c(<n+1>,model.ell,5, ... , model.ell,5,model.const,0)
# Usage> param <- ellipse2param(ellipse)
# 
ellipse2param <-function(ellipse,const_mix=-1.5){
param <- c(t(as.matrix(ellipse)))
if (param[6] != 0.0) {print("Error!")}
param <- param[-6]
c(param,const_mix)
}


# 
# mask.freeze (optimization procedure will only allow certain parameters to be fit by optim, while the others stay fixed)
# -----------------------------------------------------------------------------------------------------------------------
# The input is the "param" array assuming the form for "param2" of c(<n+1>,model.ell,5, ... , model.ell,5,model.const,0),
# a vecor of 1s and 0s "mask" which indicates which parameters in "param" should be fit (1) or frozen (0), and the ppp object
# clust.
# The output is an optimized parameter array.
# Usage> param <- mask.freeze(param,mask.adjust_positions(n),clust=clust)
# 
mask.freeze <- function(param_input,mask,clust=clust){

  ind_masked <- which(mask == 0)
  ind_unmasked <- which(mask == 1)

  ocf <- optim(param_input[ind_unmasked],model.lik,model=mask.model,clust=clust,param2=c(param_input[ind_masked],ind_masked))
  param_out <- param_input
  param_out[ind_unmasked] <- ocf$par
  param_out
}

# 
# mask.model (model that is passed to optim by the mask.freeze function)
# ----------------------------------------------------------------------
# The input is the same as other models, x and y coordinates, and the model "param" parameters, and
# the model "param2" description.
# The output is the value of the model at the position (x,y).
# Usage> ocf <- optim(param_input[ind_unmasked],model.lik,model=mask.model,clust=clust,param2=c(param_input[ind_masked],ind_masked))
# 
mask.model <- function(x1,y1,param=NULL,param2=param2) {

  n_masked <- length(param2)/2
  n_unmasked <- length(param)
  n_par <- n_masked + n_unmasked
  mask <- rep(1,n_par)

  ind_masked=round(param2[(n_masked+1):(n_masked*2)])
  mask[ind_masked] <- 0 
  ind_unmasked=which(mask == 1)

  param_input <- rep(0.0,n_par)  
  param_input[ind_masked] <- param2[1:n_masked]
  param_input[ind_unmasked] <- param

  n <- length(param_input)/6
  param_models <- c(n+1)
  for (i in 1:n) { param_models=c(param_models,ell.model,5) }
  param_models <- c(param_models,const.model,0)
  multi.model(x1,y1,param=param_input,param2=param_models)
}

# 
# mask.adjust_positions (creats an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
# ----------------------------------------------------------------------
# The input is an integer "n" corresponding to the number of ellipsoids in a finite mixture model described by the 
# "param2" variable c(<n+1>,model.ell,5, ... , model.ell,5,model.const,0).
# The output is an array of 1s and 0s with length 6*n. Position and mixing parameters are thawed.
# Usage> param <- mask.freeze(param,mask.adjust_positions(n),clust=clust)
# 
mask.adjust_positions <- function(n) {
  mask <- rep(0,6*n)
  mask[1:2] <- 1
  if (n > 1) {
     indx <- (1:(n-1))*6
     indy <- 1+(1:(n-1))*6
     indmix <- 5+(1:(n-1))*6
     mask[indx] <- 1
     mask[indy] <- 1
     mask[indmix] <- 1
  }
  mask
}

# 
# mask.adjust_rotations (creats an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
# ----------------------------------------------------------------------
# The input is an integer "n" corresponding to the number of ellipsoids in a finite mixture model described by the 
# "param2" variable c(<n+1>,model.ell,5, ... , model.ell,5,model.const,0).
# The output is an array of 1s and 0s with length 6*n. Rotation parameters are thawed.
# Usage> param <- mask.freeze(param,mask.adjust_positions(n),clust=clust)
# 
mask.adjust_rotations <- function(n) {
  mask <- rep(0,6*n)
  mask[4] <- 1
  if (n > 1) {
     indrot <- (1:(n-1))*6+3
     mask[indrot] <- 1
  }
  mask
}

# 
# mask.adjust_size (creats an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
# ----------------------------------------------------------------------
# The input is an integer "n" corresponding to the number of ellipsoids in a finite mixture model described by the 
# "param2" variable c(<n+1>,model.ell,5, ... , model.ell,5,model.const,0).
# The output is an array of 1s and 0s with length 6*n. Core radius and mixing parameters are thawed.
# Usage> param <- mask.freeze(param,mask.adjust_positions(n),clust=clust)
# 
mask.adjust_size <- function(n) {
  mask <- rep(0,6*n)
  mask[3] <- 1
  if (n > 1) {
     indcore <- (1:(n-1))*6+2
     indmix <- 5+(1:(n-1))*6
     mask[indcore] <- 1
     mask[indmix] <- 1
  }
  mask
}

# 
# mask.adjust_shape (creats an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
# ----------------------------------------------------------------------
# The input is an integer "n" corresponding to the number of ellipsoids in a finite mixture model described by the 
# "param2" variable c(<n+1>,model.ell,5, ... , model.ell,5,model.const,0).
# The output is an array of 1s and 0s with length 6*n. Core radius and axis ratio parameters are thawed.
# Usage> param <- mask.freeze(param,mask.adjust_positions(n),clust=clust)
# 
mask.adjust_shape <- function(n) {
  mask <- rep(0,6*n)
  mask[3] <- 1
  mask[5] <- 1
  if (n > 1) {
     indcore <- (1:(n-1))*6+2
     indb <- (1:(n-1))*6+4
     indmix <- 5+(1:(n-1))*6
     mask[indcore] <- 1
     mask[indb] <- 1
     mask[indmix] <- 1
  }
  mask
}

# 
# mask.adjust_mix (creats an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
# ----------------------------------------------------------------------
# The input is an integer "n" corresponding to the number of ellipsoids in a finite mixture model described by the 
# "param2" variable c(<n+1>,model.ell,5, ... , model.ell,5,model.const,0).
# The output is an array of 1s and 0s with length 6*n. Mixture coefficient are thawed.
# Usage> param <- mask.freeze(param,mask.adjust_positions(n),clust=clust)
# 
mask.adjust_mix <- function(n) {
  mask <- rep(0,6*n)
  if (n > 1) {
     indmix <- 5+(1:(n-1))*6
     mask[indmix] <- 1
  }
  mask[length(mask)] <- 1
  mask
}

# 
# ds9.colors (creates an intensity colormap mimicking the ds9 hsv scheme)
# ----------------------------------------------------------------------
# The input "n" is the number of steps in the map.
# The output is an array of hexadecimal RGB colors.
# Usage> image(t(volcano)[ncol(volcano):1,],col=ds9.colors(512))
# 
ds9.colors <- function(n){
x <- (0:n)*1.0/n
h <- (270.0/360.0-x+1.0) %% 1.0
s <- 4.0*x*(1.0-x)
v <- x^0.33
color <- hsv(h=h,s=s,v=v)
color
}

# 
# american.colors (creates an intensity colormap with white as the middle color)
# ----------------------------------------------------------------------
# The input "n" is the number of steps in the map.
# The output is an array of hexadecimal RGB colors.
# Usage> image(t(volcano)[ncol(volcano):1,],col=american.colors(512))
# 
american.colors <- function(n){
x <- (0:n)*1.0/n
h <- round(1-x)*0.667
v <- pmin(1.0,(1.1-abs(1-2*x))*3)
s <- abs(2*x-1)
hsv(h=h,s=s,v=v)
}

#
# ad2xy.tangent (Gnomonic (Tangent Plane) Projection)
# ---------------------------------------------------
# Input: Vector of RAs (alpha), vector of DECs (delta), center 
# RA and Dec (alpha0,delta0), and conversion from degrees to pixel
# All units are in decimal degrees.
# The output is a vector of projected (x,y) coordinates
# Usage> ad2xy.tangent(alpha,delta,alpha0,delta0,scale=60.0)
#
ad2xy.tangent <- function(alpha,delta,alpha0,delta0,scale=1.0){

alph <- as.double(alpha)*pi/180.0
delt <- as.double(delta)*pi/180.0

alph0 <- as.double(alpha0)*pi/180.0
delt0 <- as.double(delta0)*pi/180.0

A <- cos(delt) * cos(alph - alph0)
F <- scale * (180.0/pi) / (sin(delt0) * sin(delt) + A*cos(delt0))

x <- -F * cos(delt) * sin(alph - alph0) # SAMPLE
y <- F * (cos(delt0) * sin(delt) - A*sin(delt0)) # LINE

data.frame(x=x,y=y)
}

#
# xy2ad.tangent (Inverse of ad2xy.tangent)
# ---------------------------------------------------
# Input: Vector of x's, vector of y's, center 
# RA and Dec (alpha0,delta0), and conversion from degrees to pixel
# All units are in decimal degrees.
# The output is a vector of projected (x,y) coordinates
# Usage> ad2xy.tangent(x,y,alpha0,delta0,scale=60.0)
#
# Corrected March 30, 2016
#
xy2ad.tangent <- function(x,y,alpha0,delta0,scale=1.0){
  x2 = x/(scale*180.0/pi)
  y2 = -y/(scale*180.0/pi)
  alph0 <- as.double(alpha0)*pi/180.0
  delt0 <- as.double(delta0)*pi/180.0
  D = atan(sqrt(x2^2+y2^2))
  B = atan2(-x2,y2)
  XX = sin(delt0)*sin(D)*cos(B) + cos(delt0)*cos(D)
  YY = sin(D)*sin(B)
  alpha = alph0 + atan2(YY,XX)
  delta = asin(sin(delt0)*cos(D) - cos(delt0)*sin(D)*cos(B))
  data.frame(alpha=alpha*180.0/pi,delta=delta*180.0/pi)
}


#
# star.ppp (Creates a spatstat point process object from machine readable tables)
# -------------------------------------------------------------------------------
# The input is the name of the target region and the distance to the target.
# The output is a "ppp" object.
# Usage> clust <- star.ppp(target='w40',distance=0.50)
#
star.ppp <- function(target=target,distance=distance){
# conversion between degrees and parsecs
deg2pc = distance*1000.0*3600.0/206264.8
# create window
dat.win <- read.table('fov_mrt.txt',skip=18,col.names=c('target','ra','dec','type'))
ra0 <- dat.win[which(dat.win$target == target & dat.win$type == "origin"),2]
de0 <- dat.win[which(dat.win$target == target & dat.win$type == "origin"),3]
ra.win <- dat.win[which(dat.win$target == target & dat.win$type == "vertex"),2]
de.win <- dat.win[which(dat.win$target == target & dat.win$type == "vertex"),3]
xy.win <- ad2xy.tangent(ra.win,de.win,ra0,de0,scale=deg2pc)
win <- owin(poly = xy.win)
# get points
dat.stars <- read.table('stars_mrt.txt',skip=18,col.names=c('target','desig','ra','dec'))
ra.stars <- dat.stars[which(dat.stars$target == target),3]
de.stars <- dat.stars[which(dat.stars$target == target),4]
xy.stars <- ad2xy.tangent(ra.stars,de.stars,ra0,de0,scale=deg2pc)
# output spatstat ppp object
ppp(xy.stars$x,xy.stars$y,window=win)
}

# 
# make.fig2 (Plotting commands to produce Figure 2)
# -------------------------------------------------
# Input: the param, param2, ppp object, pixelated surface density image, and min and max colorscale values
# Output: There is no output. Produces image in graphics device.
# Usage> make.fig2(param, param2, clust=clust, image=density, min.im=5.0,max.im=20000.0)
# 
make.fig2 <- function(param,param2,clust=clust,image=image, min.im=min.im,max.im=max.im){
im2 <- eval.im(5.0*(density < min.im) + density*(density > min.im)) # To make the plot look nice.
plot(clust$window,main="Adaptive Surface Density",lwd=10)
plot(eval.im(log10(im2)),add=T,col=ds9.colors(1000),transparent.color='transparent',zlim=log10(c(min.im,max.im)))
draw.model(param,param2,size=1)
plot(clust$window,add=T,lwd=10)
}

# 
# make.fig4 (Plotting commands to produce Figure 4 (left panel))
# --------------------------------------------------------------
# Input: param, model, ppp object, param2, smoothing bandwith
# Output: There is no output. Produces an image in graphics device.
# Usage> make.fig4(param,model=multi.model,clust=clust,param2=param2,bandwidth=0.38)
# 
make.fig4 <- function(param,model=model,clust=clust,param2=param2, bandwidth=bandwidth){
mlfit <- suppressWarnings(
model.ppm(param,model=multi.model,clust=clust,param2=param2))
dd <- diagnose.ppm(mlfit,which="smooth",sigma=bandwidth,plot.smooth='image',main='',plot.it=FALSE)
ddmax <- max(c(max(dd$smooth$Z),-min(dd$smooth$Z)))
plot(dd$smooth$Z,box=FALSE,main='',zlim=c(-ddmax,ddmax),col=american.colors(512),
ribsep=0.02,ribside='right',ribargs=)
draw.model(param,param2,size=1,border=1)
plot(clust$window,add=T,lwd=10)
}

#
# ds9.poly
# ---------------------------------------------------------------
# Input: ra, de, file
# Output: There is no output. Creates a DS9 region file.
# Usage> ds9.poly(ra.win,de.win,filename='ds9.reg')
#
ds9.poly <- function(ra,de,file='ds9.reg'){
  out <- paste('fk5;polygon(',toString(paste(ra.win,de.win,sep=','),sep=','),')',sep='')
  write.table(out,row.names=F,col.names=F,quote=F,file=file)
}





