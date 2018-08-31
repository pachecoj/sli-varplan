function p = mycdf(mean,sigma,x)
p = 0.5*erfc(-(x-mean)/sqrt(2)./sigma);