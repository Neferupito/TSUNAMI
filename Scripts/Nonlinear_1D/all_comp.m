clear all;
close all;

%Read matrices
LW=readmatrix('Nonlinear_1D_LW.txt');
CN=readmatrix('Nonlinear_1D_CN.txt');
X=readmatrix('Linear_1D_x.txt');
T=readmatrix('Linear_1D_t.txt');
ANALYTICAL=readmatrix('Linear_1D_ANALYTICAL.txt');
EXPLICIT=readmatrix('Linear_1D_EXPLICIT.txt');
IMPLICIT=readmatrix('Linear_1D_IMPLICIT.txt');

ERR_EXPLICIT=vecnorm(EXPLICIT-ANALYTICAL);
ERR_IMPLICIT=vecnorm(IMPLICIT-ANALYTICAL);
ERR_LW=vecnorm(LW-ANALYTICAL);
ERR_CN=vecnorm(CN-ANALYTICAL);

%Plot
figure(1)
plot(T,ERR_EXPLICIT,T,ERR_IMPLICIT,T,ERR_LW,T,ERR_CN);
xlabel('Distance (m)');
ylabel('Amplitude');
title('Nonlinear propagation');