# Source the algorithms
include("binprobbayes.jl")
include("rtnorm.jl")

# Read in the data
using DataFrames
df = readtable("probitdata.csv")

# Execute algorithm
tic();
out = binprobbayes(convert(Array,df[:y]),convert(Array,df[collect(2:4)]),[0.; 0.; 0.],eye(3)*1000.,5000,1);
toc();

# Check Bayes estimate
mean(out,1)

# Plot trace plots
using PyPlot
plot(out)
