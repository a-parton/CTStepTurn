obs <- read.table("observations.txt", header = T)
library(moveHMM)

obs <- prepData(obs, type = "UTM")
plot(obs)

fit <- fitHMM(obs, nbState = 2, 
              stepPar0 = c(100,500,100,100),
              anglePar0 = c(pi,0,5,5))
fit
plot(fit)
plotStates(fit)

states <- stateProbs(fit)
sim_behav <- read.table("behav_proc.txt", header = T)
plot(seq(0,1000,5),states[,2],type="l")
lines(sim_behav$times,sim_behav$states-1, type="s", col = "red")
