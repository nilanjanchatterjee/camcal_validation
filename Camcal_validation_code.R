library(tidyverse, quietly = T)
library(MASS)
library(lattice)
library(pbapply)
library(lme4)

### R code to fit model for effective detection distance estimation
devtools::source_url("https://raw.githubusercontent.com/MarcusRowcliffe/distanceDF/master/distancedf.r")

##############################################################################
################## Data files for calibration validation 
predval <-read.csv("camera_deployment_validation_.csv")
sitecal <-read.csv("Site_digitisation_data.csv")
modcoef <- read.csv("pole_11_mod_param.csv")
posdat_mov <-read.csv("Speed_sequences_data.csv")

#########################################################################################
################### Simulation for Effective radius estimation
#Function:
# 1. samples from x,y co-ordinate distributions
# 2. calculates "true" radius from these using the calibration model
# 3. adds error to true radius, sampled from observed error distribution
# 4. fits detection functions to true and error radii using half normal and 
#    hazard rate keys
# 5. returns a vector of radius estimates and AICs for each of the 4 models, 
#    plus differences between true and error estimates for half normal and 
#    hazard rate models
#INPUTS
# points: number of points to sample
# b: calibration coeficients
# maxr: maximum radius to retain in samples
# mnx, mny, mnerr: means of x,y co-ordinate and error distributions
# sdx, sdy, sderr: sds of x,y co-ordinate and error distributions
sim_rep <- function(points, b, maxr, mnx, mny, mnerr, sdx, sdy, sderr){
  xsim <-(rnorm(points, mnx, sdx) / 2000) - 0.5
  ysim <-rnorm(points, mny, sdy) / 1500
  radius <- b[1] / (ysim - (b[2] + b[3] * xsim))
  radius <- radius[radius>0 & radius<maxr]
  err <- rnorm(length(radius), mnerr, sderr)
  pred <- radius + err
  dat <- data.frame(true=radius, err=pred)
#df models
  hn_true <- try(fitdf(true~1, key="hn", transect="point", order=0, data=dat))
  hr_true <- try(fitdf(true~1, key="hr", transect="point", order=0, data=dat))
  hn_err <- try(fitdf(err~1, key="hn", transect="point", order=0, data=dat))
  hr_err <- try(fitdf(err~1, key="hr", transect="point", order=0, data=dat))
#radius estimates  
  r_hn_true <- ifelse(class(hn_true)=="try-error", NA, hn_true$edd$estimate)
  r_hr_true <- ifelse(class(hr_true)=="try-error", NA, hr_true$edd$estimate)
  r_hn_err <- ifelse(class(hn_err)=="try-error", NA, hn_err$edd$estimate)
  r_hr_err <- ifelse(class(hr_err)=="try-error", NA, hr_err$edd$estimate)
#err-true differences
  diff_hn <- hn_err$edd$estimate - hn_true$edd$estimate
  diff_hr <- hr_err$edd$estimate - hr_true$edd$estimate
#model AICs
  aic_hn_true <- hn_true$ddf$criterion
  aic_hr_true <- hr_true$ddf$criterion
  aic_hn_err <- hn_err$ddf$criterion
  aic_hr_err <- hr_err$ddf$criterion
  
  c(r_hn_true=r_hn_true, r_hr_true=r_hr_true, 
    r_hn_err=r_hn_err, r_hr_err=r_hr_err,
    diff_hn=diff_hn, diff_hr=diff_hr,
    aic_hn_true=aic_hn_true, aic_hr_true=aic_hr_true,
    aic_hn_err=aic_hn_err, aic_hr_err=aic_hr_err)
}

mnx <- mean(sitecal$x)
mny <- mean(sitecal$y)
sdx <- sd(sitecal$x)
sdy <- sd(sitecal$y)
err <- predval$radius - predval$distance/100
mnerr <- mean(err, na.rm=T)
sderr <- sd(err, na.rm=T)
bflat <- subset(modcoef, location_type=="flat")[, 2:4] %>%
  apply(2, mean)

reps <- 250
flat50 <- data.frame(t(pbreplicate(reps, suppressMessages(
  sim_rep(points=50, b=bflat, maxr=25, mnx, mny, mnerr, sdx, sdy, sderr)
))))
flat100 <- data.frame(t(pbreplicate(reps, suppressMessages( 
  sim_rep(points=100, b=bflat, maxr=25, mnx, mny, mnerr, sdx, sdy, sderr)
))))
flat200 <- data.frame(t(pbreplicate(reps, suppressMessages(
  sim_rep(points=200, b=bflat, maxr=25, mnx, mny, mnerr, sdx, sdy, sderr)
))))

boxplot(flat50$r_hn_true, flat50$r_hn_err, flat50$r_hr_true, flat50$r_hr_err,
        names=rep(c("true", "error"), 2),
        col=rep(c("orange", "skyblue"), each=2),
        ylab="Effective radius (m)", main="Flat ground, 50 points")
legend("topright", c("Half normal", "Hazard rate"), fill=c("orange", "skyblue"))

boxplot(flat50$diff_hn, flat100$diff_hn, flat200$diff_hn, 
        flat50$diff_hr, flat100$diff_hr, flat200$diff_hr,
        names=rep(c(50,100,200), 2),
        col=rep(c("orange", "skyblue"), each=3),
        xlab="Points", ylab="Error (m)", main="Flat ground")
legend("topleft", c("Half normal", "Hazard rate"), fill=c("orange", "skyblue"))
lines(c(0,7), rep(mnerr,2), col="magenta")

#Hazard rate almost always preferred
mean(flat50$aic_hr_true < flat50$aic_hn_true)
mean(flat50$aic_hr_err < flat50$aic_hn_err)
mean(flat100$aic_hr_true < flat100$aic_hn_true)
mean(flat100$aic_hr_err < flat100$aic_hn_err)
mean(flat200$aic_hr_true < flat200$aic_hn_true)
mean(flat200$aic_hr_err < flat200$aic_hn_err)

#####################################################################################
############### Sloping surfaces
bslop <- subset(modcoef, location_type=="sloping")[, 2:4] %>%
  apply(2, mean)

reps <- 100
slop50 <- data.frame(t(pbreplicate(reps, suppressMessages(
  sim_rep(points=50, b=bslop, maxr=15, mnx, mny, mnerr, sdx, sdy, sderr)
))))
slop100 <- data.frame(t(pbreplicate(reps, suppressMessages( 
  sim_rep(points=100, b=bslop, maxr=15, mnx, mny, mnerr, sdx, sdy, sderr)
))))
slop200 <- data.frame(t(pbreplicate(reps, suppressMessages(
  sim_rep(points=200, b=bslop, maxr=15, mnx, mny, mnerr, sdx, sdy, sderr)
))))

boxplot(slop50$r_hn_true, slop50$r_hn_err, slop50$r_hr_true, slop50$r_hr_err,
        slop100$r_hn_true, slop100$r_hn_err, slop100$r_hr_true, slop100$r_hr_err,
        slop200$r_hn_true, slop200$r_hn_err, slop200$r_hr_true, slop200$r_hr_err,
        names=c("50_true", "50_error","50_true", "50_error",
                "100_true", "100_error","100_true", "100_error",
                "200_true", "200_error","200_true", "200_error"), 
        col=rep(c("orange", "skyblue"), each=2),
        ylab="Effective radius (m)", main="Sloping ground")
legend("topright", c("Half normal", "Hazard rate"), fill=c("orange", "skyblue"))

boxplot(slop50$diff_hn, slop100$diff_hn, slop200$diff_hn, 
        slop50$diff_hr, slop100$diff_hr, slop200$diff_hr,
        names=rep(c(50,100,200), 2),
        col=rep(c("orange", "skyblue"), each=3),
        xlab="Points", ylab="Error (m)", main="Sloping ground")
legend("topleft", c("Half normal", "Hazard rate"), fill=c("orange", "skyblue"))
lines(c(0,7), rep(mnerr,2), col="magenta")

#Hazard rate almost always preferred
mean(slop50$aic_hr_true < slop50$aic_hn_true)
mean(slop50$aic_hr_err < slop50$aic_hn_err)
mean(slop100$aic_hr_true < slop100$aic_hn_true)
mean(slop100$aic_hr_err < slop100$aic_hn_err)
mean(slop200$aic_hr_true < slop200$aic_hn_true)
mean(slop200$aic_hr_err < slop200$aic_hn_err)

####################################################################################
#################### Speed simulation function
devtools::source_url("https://raw.githubusercontent.com/MarcusRowcliffe/CTtracking/master/CTtracking.r")


### Mixed effects models for the error segregation
lmrsum <-summary(lmer(as.numeric(diff) ~ 1+ (1|deployment), data= predval))
lmrsum

depmn <-lmrsum$coefficients[1] 
depsd <-sqrt(lmrsum$varcor$deployment[1])
obssd <-lmrsum$sigma
deps <- unique(predval$deployment)
probdep <-as.numeric(xtabs(~deployment, predval)/nrow(predval))

### Load the speed sequences data
head(posdat_mov)
posdat_mov <- rename(posdat_mov, sequence_id_original=sequence_id)
posdat_mov$sequence_id <- paste(posdat_mov$siteid, posdat_mov$sequence_id_original, sep = "-")

### Adding a column for matching the folder name and prepare a new column
s1 <-as.numeric(xtabs(~sequence_id, posdat_mov))
posdat_mov$new_seq <-rep(sample(deps, 718,replace = T, prob = probdep), times= s1)
posdat_mov$new_seq1 <- paste(posdat_mov$sequence_id,posdat_mov$new_seq, sep = "-")
names(posdat_mov)[26] <-"new_seq_id"
names(posdat_mov)[28] <-"sequence_id"

seqdat <- seq.summary(posdat_mov, nframes=6) ## movement analysis with capture sequences
View(seqdat)
View(posdat_mov)
#timediff to drop some outliers and error in data
#seqdat1 <-subset(seqdat, seqdat$timediff<2000 & seqdat$species=="takin")
seqdat1 <-subset(seqdat, timediff<120 & dist>0.1 & dist<10 & species=="takin")
hist(log10(seqdat1$speed), breaks=40)
View(seqdat1)
1/mean(1/seqdat1$speed)
posdat1 <- subset(posdat_mov, sequence_id %in% seqdat1$sequence_id)

#avgpixdif <- 400 #average pixel difference between points used to generate error distributions... Actual around ~392
##

pixdiff <- seq.data(posdat1)$pixdiff
maxpixdif <- max(pixdiff, na.rm=T)
deps <- unique(posdat1$deployment)

speeds_err <- replicate(100, {
  dep_error <- rnorm(length(deps), depmn, depsd) * 0.01 #deployment-specific errors (m)
  obsd <- 0.01 * obssd * pixdiff/maxpixdif #observation specific standard deviation
  obsd[is.na(obsd)] <- 0
#  obsd[obsd>1] <- 1
  posdat_try <- posdat1
  posdat_try$radius <- posdat_try$radius + 
    dep_error[match(posdat_try$deployment, deps)] + #deployment-level error
    rnorm(nrow(posdat_try), sd=obsd) #observation-level error
  seqdat_try <- seq.summary(posdat_try)
  1/(mean(1/seqdat_try$speed))#harmonic mean speed
})

### Compare and plot the simulated speed distribution with error
truespeed <- 1/(mean(1/seqdat1$speed))
boxplot(speeds_err, ylim=c(0,max(speeds_err)), ylab= "Estimated speed (m/s)")
abline(h= truespeed, lwd=2, col="red")
median(speeds_err) / truespeed

reldist <- seqdat_try$dist / seqdat1$dist
hist(log(reldist))
sq <- seqdat1$sequence_id[order(reldist, decreasing=T)[1]]
x <- with(subset(posdat1, sequence_id==sq), radius*sin(angle))
y <- with(subset(posdat1, sequence_id==sq), radius*cos(angle))
xtry <- with(subset(posdat_try, sequence_id==sq), radius*sin(angle))
ytry <- with(subset(posdat_try, sequence_id==sq), radius*cos(angle))
plot(x, y, type="l", asp=1, xlim=range(x, xtry), ylim=range(y, ytry))
lines(xtry, ytry, col=2)

