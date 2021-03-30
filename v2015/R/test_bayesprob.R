# Source the algorithms
source("binprobbayes.r")
source("rtnorm.r")

# Read in the data
df <- read.csv("probitdata.csv")

# Execute algorithm
tic <- proc.time()
out <- binprobbayes(df$y,as.matrix(df[,2:4]),c(0,0,0),diag(3)*1000,5000,1)
proc.time()-tic

# Check Bayes estimate
colMeans(out)

# Plot trace plots
matplot(out,typ="l",lty=1)
