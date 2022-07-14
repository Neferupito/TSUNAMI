clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LINEAR 1D TSUNAMI MODEL W IMPLICIT SCHEME (1ST ORSER COUPLED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROOTS

root1 = '../../../Functions/Functions'; % Call functions
%root2 = '/Users/francoisdallasta/Documents/Numerical_project/Figures/Linear_1D/IMPLICIT_BIS/'; % Save figures folder
%root3='/Users/francoisdallasta/Documents/Numerical_project/txt_files/'; % Save txt files folder

%% Add function path

addpath(root1);

%% Boundaey condition type : 

% BC = 0 : Null displacement Free velocity
% BC = 1 : Free displacement Null velocity
% BC = 2 : Periodic

BC = 1; % Boundary condition type

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


%% Model calculations

f = waitbar(0,'Please wait...');

z=zeros(length(x)-2,1);
o=ones(length(x)-2,1);
nf=zeros(length(x),length(t)); % Vertical displacement
Mf=zeros(length(x),length(t)); % Flux
nf(:,1)=s; % Initial condition

for j=1:length(t)-1

    waitbar(j/(length(t)-1),f,'Model calculation '+string(round(100*j/(length(t)-1)))+' %');
    
    n=nf(:,j);
    M=Mf(:,j);
    nn=nf(:,j);
    MM=Mf(:,j);

    % Newton iteration

    convergence=1;

    while convergence>1.e-10
        
        %F1 F2
        F1=nn(2:end-1)-n(2:end-1)+H*r/4*(MM(3:end)-MM(1:end-2)+M(3:end)-M(1:end-2));
        F2=MM(2:end-1)-M(2:end-1)+g*r/4*(nn(3:end)-nn(1:end-2)+n(3:end)-n(1:end-2));
        
        
        %Jacobians
        [J1n,e1,e2]=Tridiag_matrix_build(z,o,z,Nx-2);
        [J1M,f1,f2]=Tridiag_matrix_build(-o*H*r/4,z,r*H*o/4,Nx-2);
        [J2n,g1,g2]=Tridiag_matrix_build(-o*g*r/4,z,o*g*r/4,Nx-2);
        [J2M,h1,h2]=Tridiag_matrix_build(z,o,z,Nx-2);
        
        
        %Boundary condition
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
            
        else 
            
            %Periodic
            J1n(1,end)=e1;
            J1n(end,1)=e2;
            J1M(1,end)=f1;
            J1M(end,1)=f2;
            J2n(1,end)=g1;
            J2n(end,1)=g2;
            J2M(1,end)=h1;
            J2M(end,1)=h2;
            
        end
       
        %Full Jacobian
        J1=[J1n,J1M];
        J2=[J2n,J2M];
        J=sparse([J1;J2]);
        F=[F1;F2];
        
        %Solve system of linear equations J x S = -F
        S=J\-F;
        convergence=norm(F);
        
        %Add increments of the variables
        nn(2:end-1)=nn(2:end-1)+S(1:Nx-2);
        MM(2:end-1)=MM(2:end-1)+S(Nx-1:end);
        
        %Boundary condition
        if BC==0
            
            % Null displacement free velocity
            MM(1)=MM(2);
            MM(end)=MM(end-1);
           
        elseif BC==1
       
            % Free displacement null velocity
            nn(1)=nn(2);
            nn(end)=nn(end-1);
         
        else
            
            %Periodic
            nn(1)=nn(end-1);
            nn(end)=nn(2);
            MM(1)=MM(end-1);
            MM(end)=MM(2);
        end
        
        
    end

    % Solution n and M of the time step 

    n=nn;
    M=MM;
    nf(:,j+1)=n;
    Mf(:,j+1)=M;
    
end

%% Ploting

% Animation

figure(1)
for j=1:1:length(t)
    plot(x*1.e-3,nf(:,j))
    title('Propagation');
    xlabel('Distance (km)');
    ylabel('Vertical displacement (m)');
    ylim([-1*h 1*h])
    pause(0.01)
end

% Gaussian source

figure(2);
plot(x,s);
title('Gaussian source');
xlabel('time (s)');
ylabel('Vertical displacement (m)');
%saveas(gcf,strcat(root2,'Source'),'epsc')

% Gaussian source frequency spectrum

figure(3);
plot(freq,fftshift(abs(fft(s))));
title('Gaussian amplitude spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
%saveas(gcf,strcat(root2,'Spectrum'),'epsc')

% Propagation

figure(4)
hold on 
for j=[1 16 61 141 261 381]
    plot(x*1.e-3,nf(:,j),'DisplayName','t = '+string(round(t(j)))+' s')
    title('Wave propagation');
    xlabel('Distance (km)');
    ylabel('Vertical displacement (m)');
end
hold off
legend show
%saveas(gcf,strcat(root2,'Propagation'),'epsc')

%% Save result as txt

%writematrix(nf,strcat(root3,'Linear_1D_IMPLICIT_bis.txt'));
