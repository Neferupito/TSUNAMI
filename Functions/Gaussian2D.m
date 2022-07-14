function [z]=Gaussian2D(x,y,mu_x,mu_y,sigma_x,sigma_y)

%Create a 2D gaussian signal given:
% x,y : vectors
% mu_x,mu_y : x,y mean values
% sigma_x,sigma_y : x,y standard deviations

[X,Y]=meshgrid(x,y);
z=exp(-((X-mu_x).^2/(2*sigma_x^2)+(Y-mu_y).^2/(2*sigma_y^2)));

return