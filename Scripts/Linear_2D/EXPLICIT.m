clear all;
close all;

%% LINEAR 2D TSUNAMI MODEL W EXPLICIT SCHEME

%% ROOTS

root1 ='/Users/francoisdallasta/Documents/Numerical_project/Functions'; % Call functions
%root2 = '/Users/francoisdallasta/Documents/Numerical_project/Figures/Linear_2D/EXPLICIT/'; % Figures folder
%root3='/Users/francoisdallasta/Documents/Numerical_project/txt_files/'; % Txt files folder

%% Add functions

addpath(root1);

%% Parameters

% Space

dx=1250; % Spatial step x (m)
L=1000000; % Model length x (m)
x=0:dx:L; % Space vector x
Nx=length(x); % Number of x samples
dy=1250; % Spatial step y (m)
L=1000000; % Model length y (m)
y=0:dy:L; % Space vector y
Ny=length(y); % Number of y samples

% Time

dt=5; % Time step (s)
T=8000; % Model duration (s)
t=0:dt:T; % Time vector
Nt=length(t); % Number of t samples

% Constants

g=9.81; % Gravity acceleration (N/kg)
H=1000; % Bathymetry (+ down) (m)
c=sqrt(g*H); % Propagation speed (m/s)
rx=c*dt/dx;
ry=c*dt/dy;

% Gaussian source
% 6 sigma = source width
 
width=100000; % Gaussian width (m)
sigma=width/6; % Gaussian standard deviation (m)
mean=round((x(end)-x(1))/2); % Gaussian mean (m)
h=1; % Gaussian height (m)
s=h*Gaussian2D(x,y,mean,mean,sigma,sigma); % Gaussian source
S=fft2(s); % 2D Fourier transform of the source
freqx=(-(Nx-1)/2:(Nx-1)/2)/L; % Frequency vector x (m-1)
freqy=(-(Ny-1)/2:(Ny-1)/2)/L; % Frequency vector y (m-1)
fmax=4*1.e-5; % Source maximum spatial frequency (m-1)
dmin=round(1/fmax); % Source minimum wavelength (m)

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

%Numerical stability condition
if min(dx,dy)/(c*sqrt(2)) >= dt
    disp('Time step OK');
else
    disp('Time step must be < '+string(min(dx,dy)/(c*sqrt(2)))+' m');
    return
end

%% Model calculation

f = waitbar(0,'Please wait...');

snap=10;
n_final=zeros(length(y),length(x),round(Nt/snap)); % Vertical displacement storage vector
nnew=zeros(length(y),length(x)); % n+1 time step
n=zeros(length(y),length(x)); % n time step
nold=zeros(length(y),length(x)); % n-1 time step

% Initial condition

n_final(:,:,1)=s;
n=s; % Initial condition : Gaussian source in x at t=0

% First time step

k=2;

waitbar(k/length(t),f,'Model calculation '+string(round(100*k/length(t)))+' %');

nnew(2:end-1,2:end-1)=n(2:end-1,2:end-1)+ry^2*(n(3:end,2:end-1)-2*n(2:end-1,2:end-1)+n(1:end-2,2:end-1))/2+rx^2*(n(2:end-1,3:end)-2*n(2:end-1,2:end-1)+n(2:end-1,1:end-2))/2;
nold=n;
n=nnew;

% Save needed snapshots

if mod(k,snap)==0
    n_final(:,:,fix(k/snap))=n;
end

% Next time steps

for k=3:length(t)

    waitbar(k/length(t),f,'Model calculation '+string(round(100*k/length(t)))+' %');
    
    nnew(2:end-1,2:end-1)=2*n(2:end-1,2:end-1)-nold(2:end-1,2:end-1)+ry^2*(n(3:end,2:end-1)-2*n(2:end-1,2:end-1)+n(1:end-2,2:end-1))+rx^2*(n(2:end-1,3:end)-2*n(2:end-1,2:end-1)+n(2:end-1,1:end-2));
    
    nold=n;
    n=nnew;
    
    % Save needed snapshots

    if mod(k,snap)==0
        n_final(:,:,fix(k/snap))=n;
    end
end

close(f);


%% Ploting

% Animation

figure(1)
for j=1:round(Nt/snap)
    imagesc(x*1.e-3,y*1.e-3,n_final(:,:,j));
    title('Propagation');
    xlabel('X (km)');
    ylabel('Y (km)');
    c=colorbar();
    colormap('turbo')
    c.Label.String = 'Vertical displacement (m)';
    pause(0.01);
end

% Gaussian source

figure(2);
imagesc(x*1.e-3,y*1.e-3,s);
title('Gaussian source');
xlabel('X (km)');
ylabel('Y (km)');
colormap('turbo')
c=colorbar();
c.Label.String = 'Vertical displacement (m)';
%saveas(gcf,strcat(root2,'Source'),'png')

% Gaussian source frequency spectrum

figure(3);
imagesc(freqx,freqy,fftshift(abs(S)));
c=colorbar();
title('2D Gaussian amplitude spectrum');
xlabel('X Spatial frequency (m-1)');
ylabel('Y Spatial frequency (m-1)');
colormap('turbo')
c.Label.String = 'Amplitude';
%saveas(gcf,strcat(root2,'Spectrum'),'png')

% Propagation - Vertical displacement

figure(4)
i=1;
for j=[1 11 21 41 61 81]
    subplot(3,2,i);
    imagesc(x*1.e-3,y*1.e-3,n_final(:,:,j));
    set(gca,'DataAspectRatio',[1 1 1]);
    c=colorbar;
    c.Label.String = 'Vertical displacement (m)';
    title('t = '+string(t(j*snap))+' s');
    xlabel('X Distance (km)');
    ylabel('Y Distance (km)');
    colormap('turbo')
    i=i+1;
end
%saveas(gcf,strcat(root2,'Propagation'),'png')
