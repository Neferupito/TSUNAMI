clear all,
close all,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPACE LOCAL TRUNCATION ERROR LINEAR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ROOTS
%%%%%%

root1 = '/Users/francoisdallasta/Documents/Numerical_project/Functions'; % Call functions
root2 = '/Users/francoisdallasta/Documents/Numerical_project/Figures/Linear_1D/COMPARISON/'; % Save figure


%Add function path
%%%%%%%%%%%%%%%%%%

addpath(root1);


%%%%%%%%%%%
%Parameters
%%%%%%%%%%%


%Spatial step 
dx=1250; % m


%Time step vector
dt=[0.0125,0.025,0.05,0.1,0.2]; % s

%Model length
L=1000000; % m

%X vector
x=0:dx:L;

%Constants
g=9.81; % Gravity acceleration (N/kg)
H=1000; % Bathymetry (+ down) (m)
c=sqrt(g*H); % Propagation speed (m/s)

%Gaussian source
%6 sigma = source width
width=100000; % Gaussian width (m)
sigma=width/6; % Gaussian standard deviation (m)
center=L/2; % Gaussian mean (m)


%%%%%%%%%%%%%%
% Calculations
%%%%%%%%%%%%%%

%EXPLICIT SCHEME

EULER_ERROR=[];

for j=1:length(dt)
    
    r=c^2*dt(j)^2/dx^2;
    D=[];
    
    for i=1:3
        
        analytical=1/2*(exp(-(x+c*(i-1)*dt(j)-center).^2/(2*sigma^2))+exp(-(x-c*(i-1)*dt(j)-center).^2/(2*sigma^2)));
        analytical=analytical';
        D=[D,analytical];
        
    end
    
    EULER=2*D(2:end-1,2)-D(2:end-1,1)+r*(D(3:end,2)-2*D(2:end-1,2)+D(1:end-2,2));
    EULER_ERROR=[EULER_ERROR;norm(D(2:end-1,3)-EULER)];

end


%IMPLICIT SCHEME (2nd order wave equation)

CN_ERROR=[];

for j=1:length(dt)

    r=0.25*c^2*dt(j)^2/dx^2; 
    Nx=length(x);
    D=[];
    
    for i=1:3
        
        analytical=1/2*(exp(-(x+c*(i-1)*dt(j)-center).^2/(2*sigma^2))+exp(-(x-c*(i-1)*dt(j)-center).^2/(2*sigma^2)));
        analytical=analytical';
        D=[D,analytical];
        
    end
    
    A=-r*ones(Nx-2,1);
    B=(1+2*r)*ones(Nx-2,1);
    C=-r*ones(Nx-2,1);
    R=sparse(Tridiag_matrix_build(A,B,C,Nx-2)); % R : Tridiagonal matrix with constant values

    
    % Vector of values n at time steps j-1 & j
    V=2*r*D(1:end-2,2)+(2-4*r)*D(2:end-1,2)+2*r*D(3:end,2)+r*D(1:end-2,1)-(1+2*r)*D(2:end-1,1)+r*D(3:end,1);
 
    
    %solve the linear system Rn=V to calculate n at time step j+1
    CN=R\V;
    CN_ERROR=[CN_ERROR;norm(D(2:end-1,3)-CN)];
end


%IMPLICIT SCHEME (1st order wave equation)


CN_BIS_ERROR=[];
for j=1:length(dt)

    r=dt(j)/dx;
    Nx=length(x);
    D=[];
    
    for i=1:2
        
        analytical=1/2*(exp(-(x+c*(i-1)*dt(j)-center).^2/(2*sigma^2))+exp(-(x-c*(i-1)*dt(j)-center).^2/(2*sigma^2)));
        analytical=analytical';
        D=[D,analytical];
        
    end
    
    z=zeros(length(x)-2,1);
    o=ones(length(x)-2,1);
    n=D(:,1);
    M=ones(size(D(:,1)))*0;
    nn=D(:,1);
    MM=ones(size(D(:,1)))*0;
    
    %Newton iteration
    convergence=1;
    while convergence>1.e-10
        %F1 F2
        F1=nn(2:end-1)-n(2:end-1)+H*r/4*(MM(3:end)-MM(1:end-2)+M(3:end)-M(1:end-2));
        F2=MM(2:end-1)-M(2:end-1)+g*r/4*(nn(3:end)-nn(1:end-2)+n(3:end)-n(1:end-2));
        
        %Jacobians
        J1n=Tridiag_matrix_build(z,o,z,Nx-2);
        J1M=Tridiag_matrix_build(-o*H*r/4,z,r*H*o/4,Nx-2);
        J2n=Tridiag_matrix_build(-o*g*r/4,z,o*g*r/4,Nx-2);
        J2M=Tridiag_matrix_build(z,o,z,Nx-2);
        
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
        
    end
      
    CN_BIS=nn;
    
    CN_BIS_ERROR=[CN_BIS_ERROR;norm(D(:,2)-CN_BIS)];
end


%%%%%
%PLOT
%%%%%

figure(1)
hold on;
plot(log10(dt),log10(EULER_ERROR),'DisplayName','EXPLICIT 2nd order eq')
plot(log10(dt),log10(CN_ERROR),'DisplayName','IMPLICIT 2nd order eq')
plot(log10(dt),log10(CN_BIS_ERROR),'DisplayName','IMPLICIT 1st order eq')
plot(log10(dt),2*log10(dt)-7,'DisplayName','Slope = 2')
xlabel('log10(dt (s))')
ylabel('log10(Error)')
title('Time local truncation error')
hold off;
legend('Location','Southeast');
saveas(gcf,strcat(root2,'Error_T'),'png')

%%%%%%%
%SLOPES
%%%%%%%

EULER_SLOPE=log10(EULER_ERROR(end)/EULER_ERROR(1))/log10(dt(end)/dt(1));
CN_SLOPE=log10(CN_ERROR(end)/CN_ERROR(1))/log10(dt(end)/dt(1));
CN_BIS_SLOPE=log10(CN_BIS_ERROR(end)/CN_BIS_ERROR(1))/log10(dt(end)/dt(1));





