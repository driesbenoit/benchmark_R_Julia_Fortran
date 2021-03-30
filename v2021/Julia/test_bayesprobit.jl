# Set wd
cd("./")

# Libraries
using CSV
using Tables 

# Source the algorithms
include("binprobbayes.jl")
include("rtnorm.jl")

# Read in csv file
println("Reading data from csv...")
dat = Tables.matrix(CSV.File("probitdata.csv"));
y =  dat[:,1].==1; # Make Boolean
X = dat[:,2:4];

# Run + time Bayesian binary probit model
@time out = binprobbayes(y,X,zeros(3),Symmetric(diagm(ones(3)*1000.0)),5000);

# Check Bayes estimates
mean(out,dims=1)

# Plot trace plots
using Plots
plot(out)
