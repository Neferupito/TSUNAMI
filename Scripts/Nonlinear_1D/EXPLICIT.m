clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NONLINEAR 1D TSUNAMI MODEL EXPLICIT SCHEME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROOTS

root1 = '/Users/francoisdallasta/Documents/Numerical_project/Functions'; % Call functions
%root2 = '/Users/francoisdallasta/Documents/Numerical_project/Figures/Nonlinear_1D/EXPLICIT/'; % Figures folder
%root3='/Users/francoisdallasta/Documents/Numerical_project/txt_files/'; % Txt files folder
root4='/Users/francoisdallasta/Documents/Numerical_project/Videos/'; % mp4 file

%% Add function path

addpath(root1);

%% Parameters

%% Boundaey condition type : 

% BC = 0 : Null displacement Free velocity
% BC = 1 : Free displacement Null velocity

BC = 1; % Boundary condition type


%% Bathymétrie

B=double(-BathyTIFF('/Users/francoisdallasta/Documents/Numerical_project/Bathymetry/Bathy.tif',20)); % Read tiff
C=B(2200,1000:2500);
H=C;
xC=0:1000:1000*(length(H)-1);

%% Space

dx=100; % Spatial step (m)
L=1000*(length(C)-1); % Model length (m)
x=0:dx:L; % Space vector
H= interp1(xC,H,x,'spline');
H=H';
Nx=length(x); % Number of x samples


%% Time

dt=0.5; % Time step (s)
T=6000; % Model duration (s)
t=0:dt:T; % Time vector
Nt=length(t); % Number of time samples

%% Constants

g=9.81; % Gravity acceleration (N/kg)
m=0.1; % Bottom friction coefficient
c=sqrt(g*H(:)); % Propagation spead (m/s)
r=dt/dx;

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
    disp('Water depth/wavelength < 0.05 , increase wavelength or decrease water depth ');
    return
end

% Numerical dispersion condition

if dmin/20 >= dx
    disp('Spatial step OK');
else
    disp('Spatial step must be < '+string(dmin/20)+' m');
    return
end

% Numerical stability condition

if dx/max(c) >= dt
    disp('Time step OK');
else
    disp('Time step must be < '+string(dx/max(c))+' m');
    return
end

%% Model calculation

f = waitbar(0,'Please wait...');

snap=10;
n_final=zeros(length(x),round(Nt/snap)); % Vertical displacement storage vector
nnew=zeros(length(x),1); % n+1 time step
n=zeros(length(x),1); % n time step

M_final=zeros(length(x),round(Nt/snap)); % Vertical displacement storage vector
Mnew=zeros(length(x),1); % n+1 time step
M=zeros(length(x),1); % n time step

n_final(:,1)=s;
n(:)=s; % Initial condition
H1= 1/2*(H(3:end)+H(2:end-1));
H2= 1/2*(H(2:end-1)+H(1:end-2));

% Loop

for j=2:length(t)

    waitbar(j/length(t),f,'Model calculation '+string(round(100*j/length(t)))+' %');
 
    % STEP1

    n_1=1/2*(n(3:end)+n(2:end-1))-r/2*(M(3:end)-M(2:end-1));
    n_2=1/2*(n(2:end-1)+n(1:end-2))-r/2*(M(2:end-1)-M(1:end-2));
    M_1=1/2*(M(3:end)+M(2:end-1))-r/2*(M(3:end).^2./(H(3:end)+n(3:end))-M(2:end-1).^2./(H(2:end-1)+n(2:end-1)))-g*r/4*(H(3:end)+n(3:end)+H(2:end-1)+n(2:end-1)).*(n(3:end)-n(2:end-1))-g*m^2*dt/4*(M(3:end).^2./real((H(3:end)+n(3:end)).^(7/3))+M(2:end-1).^2./real((H(2:end-1)+n(2:end-1)).^(7/3)));
    M_2=1/2*(M(2:end-1)+M(1:end-2))-r/2*(M(2:end-1).^2./(H(2:end-1)+n(2:end-1))-M(1:end-2).^2./(H(1:end-2)+n(1:end-2)))-g*r/4*(H(2:end-1)+n(2:end-1)+H(1:end-2)+n(1:end-2)).*(n(2:end-1)-n(1:end-2))-g*m^2*dt/4*(M(2:end-1).^2./real((H(2:end-1)+n(2:end-1)).^(7/3))+M(1:end-2).^2./real((H(1:end-2)+n(1:end-2)).^(7/3)));
    
    
    % STEP2 

    nnew(2:end-1) = n(2:end-1)-r*(M_1-M_2);
    Mnew(2:end-1) = M(2:end-1)-r*(M_1.^2./(H1+n_1)-M_2.^2./(H2+n_2))-g*r/2*(H1+n_1+H2+n_2).*(n_1-n_2)-g*m^2*dt/2*(M_1.^2./real((H1+n_1).^(7/3))+M_2.^2./real((H2+n_2).^(7/3)));
    
    
    % Boundary condition

    if BC==0
            
        % Null displacement free velocity

        Mnew(end)=Mnew(end-2);
        Mnew(1)=Mnew(3);
           
    elseif BC==1
       
        % Free displacement null velocity

        nnew(end)=nnew(end-2);
        nnew(1)=nnew(3);
            
        
    end

    n=nnew;
    M=Mnew;

    % save needed snapshots

    if mod(j,snap)==0
        n_final(:,fix(j/snap))=n;
        M_final(:,fix(j/snap))=M;
    end
end

close(f);

%% Initialize video
myVideo = VideoWriter(strcat(root4,'Nonlinear_1D_EXPLICIT'), 'MPEG-4'); %open video file
myVideo.FrameRate = 12;  
open(myVideo)
figure()
for j=1:3:length(n_final(1,:))
    subplot(3,1,1)
    plot(x*1.e-3,n_final(:,j))
    title('1D Propagation');
    xlabel('X (km)');
    ylabel('Vertical displacement (m)','Interpreter','latex');
    legend(string(dt*j*snap)+' s')
    legend('Location','northwest')

    subplot(3,1,2)
    plot(x*1.e-3,M_final(:,j))
    title('1D Propagation');
    xlabel('X (km)');
    ylabel('Flux ($m^{2}/s$)','Interpreter','latex');
    legend(string(dt*j*snap)+' s')
    legend('Location','northwest')
    
    subplot(3,1,3)
    plot(x*1.e-3,-H)
    title('Bathymetry');
    xlabel('X (km)');
    ylabel('Elevation (m)','Interpreter','latex');
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

%% Ploting

% Gaussian source

figure(1);
plot(x,s);
title('Gaussian source');
xlabel('Distance (m)');
ylabel('Vertical displacement (m)');
%saveas(gcf,strcat(root2,'Source'),'png')

% Gaussian source frequency spectrum 

figure(2);
plot(freq,fftshift(abs(S)));
title('Gaussian amplitude spectrum');
xlabel('Spatial frequency (m-1)');
ylabel('Amplitude');
%saveas(gcf,strcat(root2,'Spectrum'),'png')

% Propagation - Vertical displacement

figure(3)
i=1;
for j=[1 150 300 450 600 750]
    subplot(3,2,i)
    plot(x*1.e-3,n_final(:,j))
    title('1D Propagation','Interpreter','latex');
    xlabel('X (km)','Interpreter','latex');
    ylabel('Vertical displacement (m)','Interpreter','latex');
    title('T = '+string(round(snap*dt*j))+' s');
    i=i+1;
end
%saveas(gcf,strcat(root2,'n'),'png')

% Propagation - Flux

figure(4)
i=1;
for j=[1 150 300 450 600 750]
    subplot(3,2,i)
    plot(x*1.e-3,M_final(:,j))
    title('T = '+string(round(snap*dt*j))+' s','Interpreter','latex');
    xlabel('X (km)','Interpreter','latex');
    ylabel('Flux ($m^{2}/s$)','Interpreter','latex');
    i=i+1;
end
%saveas(gcf,strcat(root2,'M'),'png')

% Bathymetry

figure(5)
plot(x*1.e-3,-H)
title('Bathymetry','Interpreter','latex');
xlabel('X (km)','Interpreter','latex');
ylabel('Elevation (m)','Interpreter','latex');

%saveas(gcf,strcat(root2,'Bathymetry'),'png')


%% Save matrix as txt

%writematrix(n,strcat(root3,'Nonlinear_1D_EXPLICIT.txt'));