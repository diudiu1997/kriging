## this is the covariance function. 
## output is bigger when dist == zero 
## this adds uncorrelated noise
function out =  gauCov2(a,h,sil,dist)
  bigzero = (dist>0);
  iszero = (dist ==0);
  out = (a*exp(-dist.^2/h) + a*sil).*iszero+(a*exp(-dist.^2/h)).*bigzero;