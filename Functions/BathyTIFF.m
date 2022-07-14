function A=BathyTIFF(File,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output the smoothed Bathymetry matrix 
%from a GEOTIFF file


%File : File name/path of the Geotiff
%Sigma : Standard deviation of the gaussian filtering


A = readgeoraster(File); % Read georaster
A=imgaussfilt(A,sigma); % Smoothing

end
