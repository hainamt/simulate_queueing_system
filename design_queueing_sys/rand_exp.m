function [x]=rand_exp(lambda)
y=rand;
if y==0
    y=rand;
end

x=-1/lambda*log(1-y);
