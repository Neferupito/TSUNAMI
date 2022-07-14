clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NONLINEAR 1D TSUNAMI MODEL IMPLICIT SCHEME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROOTS

root1 = '/Users/francoisdallasta/Documents/Numerical_project/Functions'; % Call functions
%root2 = '/Users/francoisdallasta/Documents/Numerical_project/Figures/Nonlinear_1D/IMPLICIT/'; % Figures folder
%root3='/Users/francoisdallasta/Documents/Numerical_project/txt_files/'; % Txt files folder
root4='/Users/francoisdallasta/Documents/Numerical_project/Videos/'; % mp4 file
%% Add function path

addpath(root1);

%% Boundaey condition type : 

% BC = 0 : Null displacement Free velocity
% BC = 1 : Free displacement Null velocity

BC = 1; % Boundary condition type

%% Bathymetry

B=double(-BathyTIFF('/Users/francoisdallasta/Documents/Numerical_project/Bathymetry/Bathy.tif',20)); % Read tiff
C=B(2200,1000:2500);
H=C;
xC=0:1000:1000*(length(H)-1);

%% Parameters

% Space

dx=1000; % Spatial step (m)
L=1000*(length(C)-1); % Model length (m)
x=0:dx:L; % Space vector
H= interp1(xC,H,x,'spline');
H=H';
Nx=length(x); % Number of x samples

% Time vector

dt=10; % Time step (s)
T=6000; % Model duration (s)
t=0:dt:T; % Time vector (s)

% Constants

g=9.81; % Gravity acceleration (N/kg)
c=sqrt(g*H); % m/s
r=0.25*dt/dx;
m=0.1; % Bottom friction coefficient

% Gaussian source
% 6 sigma = source width

width=100000; % Gaussian width (m)
sigma=width/6; % Gaussian standard deviation (m)
mean=round((x(end)-x(1))/2); % Gaussian mean (m)
h=1; % Gaussian height (m)
s=h*Gaussian(x,mean,sigma); % Source
S=fft(s); % Fourier transform of the source
freq=(-(Nx-1)/2:(Nx-1)/2)/L;%m-1
fmax=4*1.e-5; %m-1
dmin=round(1/fmax); % Source minimum wavelength (m)

% Shallow water approximation

if H/width <= 0.05
    disp('Shallow water approximation OK');
else
    disp('Water depth/wavelength < 0.05 , increase wavelength or decrease water depth ');
    return
end

% Spatial sampling condition

if dmin/20 >= dx
    disp('Spatial sampling step OK');
else
    disp('Spatial sampling step must be < '+string(dmin/20)+' m');
    return
end

%% Model calculation

f = waitbar(0,'Please wait...');

z=zeros(length(x)-2,1);
o=ones(length(x)-2,1);
nf=zeros(length(x),1); % Vertical displacement
Mf=zeros(length(x),1); % Flux
nf(:)=s; % Initial condition
S=[]

for j=1:length(t)-1

    waitbar(j/(length(t)-1),f,'Model calculation '+string(round(100*j/(length(t)-1)))+' %');

    n=nf(:,end);
    M=Mf(:,end);
    nn=nf(:,end);
    MM=Mf(:,end);

    % Newton iteration

    convergence=1;
     
    while convergence>1.e-10
        
        %F1 F2
        F1=nn(2:end-1)-n(2:end-1)+r*(MM(3:end)-MM(1:end-2)+M(3:end)-M(1:end-2));
        F2=MM(2:end-1)-M(2:end-1)+r*(MM(3:end).^2./(H(3:end)+nn(3:end))-MM(1:end-2).^2./(H(1:end-2)+nn(1:end-2))+M(3:end).^2./(H(3:end)+n(3:end))-M(1:end-2).^2./(H(1:end-2)+n(1:end-2)))+r*g*((H(2:end-1)+nn(2:end-1)).*(nn(3:end)-nn(1:end-2))+(H(2:end-1)+n(2:end-1)).*(n(3:end)-n(1:end-2)))+dt*g*m^2/2*(MM(2:end-1).^2./real((H(2:end-1)+nn(2:end-1)).^(7/3))+M(2:end-1).^2./real((H(2:end-1)+n(2:end-1)).^(7/3)));
       
        
        % Jacobians

        [J1n,e1,e2]=Tridiag_matrix_build(z,o,z,Nx-2);
        [J1M,f1,f2]=Tridiag_matrix_build(-o*r,z,r*o,Nx-2);
        [J2n,g1,g2]=Tridiag_matrix_build(r*MM(1:end-2).^2.*(H(1:end-2)+nn(1:end-2)).^(-2)-r*g*(H(2:end-1)+nn(2:end-1)),r*g*(nn(3:end)-nn(1:end-2))-7*dt*g*m^2/6*MM(2:end-1).^2./real((H(2:end-1)+nn(2:end-1)).^(10/3)),-r*MM(3:end).^2.*(H(3:end)+nn(3:end)).^(-2)+r*g*(H(2:end-1)+nn(2:end-1)),Nx-2);
        [J2M,h1,h2]=Tridiag_matrix_build(-2*r*MM(1:end-2)./(H(1:end-2)+nn(1:end-2)),o+g*m^2*dt*MM(2:end-1)./real((H(2:end-1)+nn(2:end-1)).^(7/3)),2*r*MM(3:end)./(H(3:end)+nn(3:end)),Nx-2);
        
        % Boundary condition

        if BC==0 
            
            % Null displacement free velocity

            J1M(1,2)=J1M(1,2)*2;
            J1M(end,end-1)=J1M(end,end-1)*2;
            J2M(1,2)=J2M(1,2)*2;
            J2M(end,end-1)=J2M(end,end-1)*2;
         
        elseif BC==1
       
            % Free displacement null velocity

            J1n(1,2)=J1n(1,2)*2;
            J1n(end,end-1)=J1n(end,end-1)*2;
            J2n(1,2)=J2n(1,2)*2;
            J2n(end,end-1)=J2n(end,end-1)*2;
            

        end
        
        
        
        % Full Jacobian

        J1=[J1n,J1M];
        J2=[J2n,J2M];
        J=sparse([J1;J2]);

        F=[F1;F2];
        
        % Solve system of linear equations J x S = -F

        S=J\-F;
        convergence=norm(F);

        % Add increments of the variables

        nn(2:end-1)=nn(2:end-1)+S(1:Nx-2);
        MM(2:end-1)=MM(2:end-1)+S(Nx-1:end);
        
        % Boundary condition
        if BC==0
            
            % Null displacement free velocity

            MM(1)=MM(2);
            MM(end)=MM(end-1);
           
        elseif BC==1
       
            % Free displacement null velocity

            nn(1)=nn(2);
            nn(end)=nn(end-1);
         
        end
        
    end    

    % Solution n and M of the time step    

    n=nn;
    M=MM;
    nf=[nf,n];
    Mf=[Mf,M];
    
end

close(f);

%% Initialize video

myVideo = VideoWriter(strcat(root4,'Nonlinear_1D_EXPLICIT'), 'MPEG-4'); %open video file
myVideo.FrameRate = 12;  
open(myVideo)
figure()
for j=1:3:length(nf(1,:))
    subplot(3,1,1)
    plot(x*1.e-3,nf(:,j))
    title('1D Propagation');
    xlabel('X (km)');
    ylabel('Vertical displacement (m)','Interpreter','latex');
    legend(string(dt*j)+' s')
    legend('Location','northwest')

    subplot(3,1,2)
    plot(x*1.e-3,Mf(:,j))
    title('1D Propagation');
    xlabel('X (km)');
    ylabel('Flux ($m^{2}/s$)','Interpreter','latex');
    legend(string(dt*j)+' s')
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
xlabel('time (s)');
ylabel('Amplitude');
%saveas(gcf,strcat(root2,'Source'),'epsc')

% Gaussian source frequency spetcrum

figure(2);
plot(freq,abs(fft(s)));
title('Gaussian amplitude spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
%saveas(gcf,strcat(root2,'Spectrum'),'epsc')

% Propagation - Vertical displacement

figure(3)
i=1;
for j=[1 100 200 300 400 500]
    subplot(3,2,i)
    plot(x*1.e-3,nf(:,j))
    title('1D Propagation','Interpreter','latex');
    xlabel('X (km)','Interpreter','latex');
    ylabel('Displacement (m)','Interpreter','latex');
    title('T = '+string(round(dt*j))+' s');
    i=i+1;
end
%saveas(gcf,strcat(root2,'n'),'png')

% Propagation - Flux

figure(4)
i=1;
for j=[1 100 200 300 400 500]
    subplot(3,2,i)
    plot(x*1.e-3,Mf(:,j))
    title('T = '+string(round(dt*j))+' s','Interpreter','latex');
    xlabel('X (km)','Interpreter','latex');
    ylabel('X Flux ($m^{2}/s$)','Interpreter','latex');
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

%% Save txt files

%writematrix(nf,strcat(root3,'Nonlinear_1D_IMPLICIT.txt'));