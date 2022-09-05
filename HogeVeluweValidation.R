devtools::source_url("https://raw.githubusercontent.com/MarcusRowcliffe/CTtracking/master/CTtracking.r")

dat <- read.csv("HogeVeluweDigidat.csv")

# Function for one replicate test for deployment dep
oneRep <- function(dat, dep){
  sq <- 1:nrow(dat)
  i <- sq %in% sample(sq, round(0.5*nrow(dat)))
  dmod <- cal.dep(dat[i, ], flex=FALSE)
  names(dmod) <- dep
  dat[!i, ] %>% rename(x=xg, y=yg) %>% predict.pos(dmod, "folder")
}

# Function to train and test one deployment (dep), replicated reps times
oneDep <- function(dep, reps){
  dat <- filter(dat, folder==dep)
  bind_rows(replicate(reps, oneRep(dat, dep), simplify = FALSE))
}

#Replicated tests for all deployments
deps <- unique(dat$folder)
predictions <- lapply(deps, oneDep, 20)

#Plots
par(mfrow=c(5,4), mar=c(2,2,1,1))
lim <- c(0,20)
for(i in 1:length(predictions)){
  with(predictions[[i]], plot(distance, radius))
  legend("topleft", deps[i], bty="n")
  lines(lim, lim, col=2)
}
