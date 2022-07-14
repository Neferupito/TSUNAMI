clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D LINEAR ERROR COMPARISON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROOTS

root1 ='/Users/francoisdallasta/Documents/Numerical_project/txt_files/'; % Txt files folder
root2='/Users/francoisdallasta/Documents/Numerical_project/Figures/Linear_1D/COMPARAISON/'; % Figures folder

%% Add txt_files path

addpath(root1);

%% Read matrices

ANALYTICAL=readmatrix('Linear_1D_ANALYTICAL.txt');
EXPLICIT=readmatrix('Linear_1D_EXPLICIT.txt');
IMPLICIT=readmatrix('Linear_1D_IMPLICIT.txt');
IMPLICIT_BIS=readmatrix('Linear_1D_IMPLICIT_bis.txt');
X=readmatrix('Linear_1D_x.txt');
T=readmatrix('Linear_1D_t.txt');

%% ERROR EVALUTATION

% Least square (L2 norm)

ERR_EXPLICIT=vecnorm(EXPLICIT-ANALYTICAL);
ERR_IMPLICIT=vecnorm(IMPLICIT-ANALYTICAL);
ERR_IMPLICIT_BIS=vecnorm(IMPLICIT_BIS-ANALYTICAL);

%% Plot

figure(1)
hold on;
plot(log10(T),log10(ERR_EXPLICIT),'DisplayName','Explicit scheme Error');
plot(log10(T),log10(ERR_IMPLICIT),'DisplayName','Implicit scheme Error');
plot(log10(T),log10(ERR_IMPLICIT_BIS),'DisplayName','Implicit Bis scheme Error');
xlabel('log10(Time (s))');
ylabel('log10(Error)');
title('Error comparison');
hold off;
legend;
saveas(gcf,strcat(root2,'Error_evolution'),'epsc')