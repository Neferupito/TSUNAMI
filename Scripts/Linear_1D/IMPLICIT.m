clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LINEAR 1D TSUNAMI MODEL IMPLICIT SCHEME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROOTS

root1 ='/Users/francoisdallasta/Documents/Numerical_project/Functions'; % Call functions
%root2 ='/Users/francoisdallasta/Documents/Numerical_project/Figures/Linear_1D/IMPLICIT/'; % Figures folder
%root3='/Users/francoisdallasta/Documents/Numerical_project/txt_files/'; % Txt files folder

%% Add function path

addpath(root1);

%% Boundaey condition type :

%BC = 0 : Null displacement 
%BC = 1 : Free displacement
%BC = 2 : Periodic

BC =1; % Boundary condition type


%% Parameters

% Space

dx=1250; % Spatial step (m)
L=1000000; % Model length (m)
x=0:dx:L; % Space vector
Nx=length(x); % Number of x samples

% Time

dt=10; % Time step (s)
T=8000; % Model duration (s)
t=0:dt:T; % Time vector

% Constants

g=9.81; % Gravity acceleration (N/kg)
H=1000; % (Constant!!!) Bathymetry (+ down) (m)
c=sqrt(g*H); % Propagation speed (m/s)
r=0.25*c^2*dt^2/dx^2; 

% Gaussian source
% 6 sigma = source width

width=100000; % Gaussian width (m)
sigma=width/6; % Gaussian standard deviation (m)
mean=round((x(end)-x(1))/2); % Gaussian mean (m)
h=1; % Gaussian height (m)
s=h*Gaussian(x,mean,sigma); % Gaussian source
S=fft(s); % Fourier transform of the source
freq=(-(Nx-1)/2:(Nx-1)/2)/L; % Frequency vector (m-1)
fmax=4*1.e-5; % Source maximum spatial frequency (m-1)
dmin=round(1/fmax); % Source minimum wavelength (m)


%% Conditions

% Shallow water approximation

if H/width <= 0.05 
    disp('Shallow water approximation OK');
else
    disp('Water depth / wavelength > 0.05 , increase wavelength or decrease water depth');
    return
end

% Numerical dispersion condition

if dmin/20 >= dx
    disp('Spatial step OK');
else
    disp('Spatial step must be < '+string(dmin/20)+' m');
    return
end

%% Model calculation

f = waitbar(0,'Please wait...');

A=-r*ones(length(x)-2,1);
B=(1+2*r)*ones(length(x)-2,1);
C=-r*ones(length(x)-2,1);
R=sparse(Tridiag_matrix_build(A,B,C,Nx-2)); % R : Tridiagonal matrix with constant values

% Boundary conditions

if BC==1
            
    % Free displacement 

    R(1,2)=R(1,2)*2;
    R(end,end-1)=R(end,end-1)*2;

elseif BC==2 
    
    % Periodic

    R(1,end)=-r;
    R(end,1)=-r;
    
end

n=zeros(length(x),length(t)); % Vertical displacement storage vector
n(:,1)=s; % Initial condition : Gaussian source in x at t=0

% First time step

V=2*r*n(1:end-2,1)+(2-4*r)*n(2:end-1,1)+2*r*n(3:end,1); % Values n at time step 1

n(2:end-1,2)=2*R\V; % Solve the linear system Rn=V to calculate n at time step 2

%Boundary condition

if BC==1
            
    % Free displacement 

    n(1,2)=n(3,2);
    n(end,2)=n(end-2,2);
    
elseif BC==2
    
    % Periodic     

    n(1,2)=n(end-1,2);
    n(end,2)=n(2,2);
            
end
        
% Next time steps

for j=2:length(t)-1

    waitbar(j/(length(t)-1),f,'Model calculation '+string(round(100*j/(length(t)-1)))+' %');

    V=2*r*n(1:end-2,j)+(2-4*r)*n(2:end-1,j)+2*r*n(3:end,j)+r*n(1:end-2,j-1)-(1+2*r)*n(2:end-1,j-1)+r*n(3:end,j-1); % Values n at time steps j-1 & j

    n(2:end-1,j+1)=R\V; % Solve the linear system Rn=V to calculate n at time step j+1
    
    % Boundary condition

    if BC==1
            
        % Free displacement 

        n(1,j+1)=n(3,j+1);
        n(end,j+1)=n(end-2,j+1);
            
    elseif BC==2
            
        % Periodic

        n(1,j+1)=n(end-1,j+1);
        n(end,j+1)=n(2,j+1);
            
    end
    
end

close(f);

%% Ploting

% Animation

figure(1)
for j=1:length(t)
    plot(n(:,j));
    title('Propagation');
    xlabel('Distance (km)');
    ylabel('Vertical displacement (m)');
    ylim([-1*h 1*h]);
    pause(0.01);
end

% Gaussian source

figure(2);
plot(x*1.e-3,s);
title('Initial Gaussian signal');
xlabel('Distance (km)');
ylabel('Vertical displacement (m)');
%saveas(gcf,strcat(root2,'Source'),'png') % Save figure

% Gaussian source frequency spectrum

figure(3);
plot(freq,fftshift(abs(S)));
title('Gaussian amplitude spectrum');
xlabel('Spatial frequency (m-1)');
ylabel('Amplitude');
%saveas(gcf,strcat(root2,'Spectrum'),'png') % Save figure

% Propagation - Vertical displacement

figure(4)
hold on;
for j=[1 16 61 141 261 381]
    plot(x*1.e-3,n(:,j),'DisplayName','t = '+string(t(j))+' s');
    title('Propagation');
    xlabel('Distance (km)');
    ylabel('Vertical displacement (m)');
    ylim([-1*h 1*h]);
end
hold off;
legend;
%saveas(gcf,strcat(root2,'Propagation_BC_',string(BC)),'png') % Save figure




%% DISPERSION ANALYSIS

figure(5)
P=2:0.01:20;
k=2*pi./(P*dx);
alpha=0.5:1:5;
dt=alpha*dx/c;
hold on;
for i=1:length(alpha)
    Ck=1./(dt(i)*k).*acos((2-alpha(i)^2+alpha(i)^2*cos(k*dx))./(2+alpha(i)^2-alpha(i)^2*cos(k*dx)));
    plot(P,Ck/c,'DisplayName','\beta = '+string(alpha(i)));
    title('Numerical dispersion');
    xlabel('Number of grid point per wavelength');
    ylabel('Ck/Co');
end
hold off;
legend('Location','Southeast');

%saveas(gcf,strcat(root2,'Dispersion'),'png') % Save figure

%% Save result as txt

%writematrix(n,strcat(root3,'Linear_1D_IMPLICIT.txt'));
