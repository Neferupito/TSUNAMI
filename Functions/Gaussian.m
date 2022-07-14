function [y]=Gaussian(x,mu,sigma)

%Create a 1D gaussian signal given:
% x : vector
% mu : mean value
% sigma : standard deviation


y=exp(-(x-mu).^2/(2*sigma^2));

return