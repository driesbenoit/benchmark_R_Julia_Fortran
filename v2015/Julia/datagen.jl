# Generate data

## Number of observations 
n = 5000;

## Parameters
beta = [-1.; 1.; .5];
sigsq = 1.;

## Number of regression coefs
k = size(beta,1);

## Generate independent vars
X = rand(n,k-1);
X = [ones(n) X]; # concat column of ones

## Generate dependent var
y_star = X*beta + rand(Normal(0.,sigsq),n);
y = map((x) -> ifelse(x<0.,0,1),y_star);


# Save as csv

## First create dataframe
using DataArrays, DataFrames
df = hcat(DataFrame(y=y),convert(DataFrame,X));

## Write to disk
writetable("probitdata.csv", df)
