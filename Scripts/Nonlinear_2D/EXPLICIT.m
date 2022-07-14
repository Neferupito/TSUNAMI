clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NONLINEAR 2D TSUNAMI MODEL EXPLICIT SCHEME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add function path

root1 ='/Users/francoisdallasta/Documents/Numerical_project/Functions'; % Call functions
root2='/Users/francoisdallasta/Documents/Numerical_project/Videos/'; % mp4 file folder
%root3='/Users/francoisdallasta/Documents/Numerical_project/txt_files/'; % Txt files folder
%root4='/Users/francoisdallasta/Documents/Numerical_project/Figures/Nonlinear_2D/EXPLICIT/'; % Png files folder
addpath(root1);

%% Boundaey condition type : 

% BC = 0 : Null displacement Free velocity
% BC = 1 : Free displacement Null velocity

BC=1; % Boundary condition type

%% Bathymetry

B=double(-BathyTIFF('/Users/francoisdallasta/Documents/Numerical_project/Bathymetry/Bathy.tif',20)); % Read tiff
H=B(1700:end-200,1700:end-200);

%% Parameters

% Space

dx=1000; % Spatial step x (m)
Lx=1000*(length(B(1,1700:end-200))-1); % Model length x (m)
x=0:dx:Lx; % Space vector x (m)
Nx=length(x); % Number of x samples
dy=1000; % Spatial step y (m)
Ly=1000*(length(B(1700:end-200,1))-1); % Model length y (m)
y=0:dy:Ly; % Space vector y (m)
Ny=length(y); % Number of y samples

% Time

dt=3; % Time step (s)
T=8000; % Model duration (s)
t=0:dt:T; % Time vector (s)
Nt=length(t); % Number of t samples

% Constants

g=9.81; % Gravity acceleration (N/kg)
c=sqrt(g*H); % m/s

% Gaussian source
% 6 sigma = source width

width=100000; % Gaussian width (m)
sigma=width/6; % Gaussian standard deviation (m)
mean=round((x(end)-x(1))/2); % Gaussian mean (m)
h=1; % Gaussian height (m)
s=h*Gaussian2D(x,y,mean,mean,sigma,sigma); %source
S=fft2(s);


%% Conditions

% Shallow water approximation

if max(H(:))/width <= 0.05
    disp('Shallow water approximation OK');
else
    disp('Water depth/wavelength < 0.05 , increase wavelength or decrease water depth ');
    return
end

% Numerical stability condition

if min(dx,dy)/(max(c(:))*sqrt(2)) >= dt
    disp('Time step OK');
else
    disp('Time step must be < '+string(dx/(max(c(:))*sqrt(2)))+' m');
    return
end

%% Model calculations

f = waitbar(0,'Please wait...');

snap=5;
n_final=zeros(length(y),length(x),round(Nt/snap));
nnew=zeros(length(y),length(x));
n=zeros(length(y),length(x));

M_final=zeros(length(y),length(x),round(Nt/snap));
Mnew=zeros(length(y),length(x));
M=zeros(length(y),length(x));

N_final=zeros(length(y),length(x),round(Nt/snap));
Nnew=zeros(length(y),length(x));
N=zeros(length(y),length(x));


% Initial condition

n_final(:,:,1)=s;
n=s;
alpha = dt/(2*dx);
beta = dt/(2*dy);

%loop

for k=2:length(t)

    waitbar(k/length(t),f,'Model calculation '+string(round(100*k/length(t)))+' %');

    n_1=1/4*(n(5:end,3:end-2)+n(3:end-2,3:end-2) +n(4:end-1,4:end-1) +n(4:end-1,2:end-3))-alpha/2*(M(4:end-1,4:end-1)-M(4:end-1,2:end-3))...
    -beta/2*(N(5:end,3:end-2)-N(3:end-2,3:end-2));

    n_2=1/4*(n(1:end-4,3:end-2)+n(3:end-2,3:end-2) +n(2:end-3,2:end-3) +n(2:end-3,4:end-1))-alpha/2*(M(2:end-3,4:end-1)-M(2:end-3,2:end-3))...
    -beta/2*(N(3:end-2,3:end-2)-N(1:end-4,3:end-2));

    n_3=1/4*(n(4:end-1,4:end-1)+n(2:end-3,4:end-1) +n(3:end-2,3:end-2) +n(3:end-2,5:end))-alpha/2*(M(3:end-2,5:end)-M(3:end-2,3:end-2))...
    -beta/2*(N(4:end-1,4:end-1)-N(2:end-3,4:end-1));

    n_4=1/4*(n(4:end-1,2:end-3)+n(2:end-3,2:end-3)+n(3:end-2,3:end-2)+n(3:end-2,1:end-4))-alpha/2*(M(3:end-2,3:end-2)-M(3:end-2,1:end-4))...
    -beta/2*(N(4:end-1,2:end-3)-N(2:end-3,2:end-3));
    
    M_1=1/4*(M(5:end,3:end-2)+M(3:end-2,3:end-2) +M(4:end-1,4:end-1) +M(4:end-1,2:end-3))-alpha/2*(M(4:end-1,4:end-1).^2./(H(4:end-1,4:end-1)+n(4:end-1,4:end-1))-M(4:end-1,2:end-3).^2./(H(4:end-1,2:end-3)+n(4:end-1,2:end-3)))...
    -beta/2*(M(5:end,3:end-2).*N(5:end,3:end-2)./(H(5:end,3:end-2)+n(5:end,3:end-2))-M(3:end-2,3:end-2).*N(3:end-2,3:end-2)./(H(3:end-2,3:end-2)+n(3:end-2,3:end-2)))...
    -g*alpha/2*(H(4:end-1,3:end-2)+n(4:end-1,3:end-2)).*(n(4:end-1,4:end-1)-n(4:end-1,2:end-3));

    M_2=1/4*(M(1:end-4,3:end-2)+M(3:end-2,3:end-2) +M(2:end-3,2:end-3) +M(2:end-3,4:end-1))-alpha/2*(M(2:end-3,4:end-1).^2./(H(2:end-3,4:end-1)+n(2:end-3,4:end-1))-M(2:end-3,2:end-3).^2./(H(2:end-3,2:end-3)+n(2:end-3,2:end-3)))...
    -beta/2*(M(3:end-2,3:end-2).*N(3:end-2,3:end-2)./(H(3:end-2,3:end-2)+n(3:end-2,3:end-2))-M(1:end-4,3:end-2).*N(1:end-4,3:end-2)./(H(1:end-4,3:end-2)+n(1:end-4,3:end-2)))...
    -g*alpha/2*(H(2:end-3,3:end-2)+n(2:end-3,3:end-2)).*(n(2:end-3,4:end-1)-n(2:end-3,2:end-3));

    M_3=1/4*(M(4:end-1,4:end-1)+M(2:end-3,4:end-1) +M(3:end-2,3:end-2) +M(3:end-2,5:end))-alpha/2*(M(3:end-2,5:end).^2./(H(3:end-2,5:end)+n(3:end-2,5:end))-M(3:end-2,3:end-2).^2./(H(3:end-2,3:end-2)+n(3:end-2,3:end-2)))...
    -beta/2*(M(4:end-1,4:end-1).*N(4:end-1,4:end-1)./(H(4:end-1,4:end-1)+n(4:end-1,4:end-1))-M(2:end-3,4:end-1).*N(2:end-3,4:end-1)./(H(2:end-3,4:end-1)+n(2:end-3,4:end-1)))...
    -g*alpha/2*(H(3:end-2,4:end-1)+n(3:end-2,4:end-1)).*(n(3:end-2,5:end)-n(3:end-2,3:end-2));

    M_4=1/4*(M(4:end-1,2:end-3)+M(2:end-3,2:end-3)+M(3:end-2,3:end-2)+M(3:end-2,1:end-4))-alpha/2*(M(3:end-2,3:end-2).^2./(H(3:end-2,3:end-2)+n(3:end-2,3:end-2))-M(3:end-2,1:end-4).^2./(H(3:end-2,1:end-4)+n(3:end-2,1:end-4)))...
    -beta/2*(M(4:end-1,2:end-3).*N(4:end-1,2:end-3)./(H(4:end-1,2:end-3)+n(4:end-1,2:end-3))-M(2:end-3,2:end-3).*N(2:end-3,2:end-3)./(H(2:end-3,2:end-3)+n(2:end-3,2:end-3)))...
    -g*alpha/2*(H(3:end-2,2:end-3)+n(3:end-2,2:end-3)).*(n(3:end-2,3:end-2)-n(3:end-2,1:end-4));
    
    N_1=1/4*(N(5:end,3:end-2)+N(3:end-2,3:end-2) +N(4:end-1,4:end-1) +N(4:end-1,2:end-3))-alpha/2*(M(4:end-1,4:end-1).*N(4:end-1,4:end-1)./(H(4:end-1,4:end-1)+n(4:end-1,4:end-1))-M(4:end-1,2:end-3).*N(4:end-1,2:end-3)./(H(4:end-1,2:end-3)+n(4:end-1,2:end-3)))...
    -beta/2*(N(5:end,3:end-2).^2./(H(5:end,3:end-2)+n(5:end,3:end-2))-N(3:end-2,3:end-2).^2./(H(3:end-2,3:end-2)+n(3:end-2,3:end-2)))...
    -g*beta/2*(H(4:end-1,3:end-2)+n(4:end-1,3:end-2)).*(n(5:end,3:end-2)-n(3:end-2,3:end-2));

    N_2=1/4*(N(1:end-4,3:end-2)+N(3:end-2,3:end-2) +N(2:end-3,2:end-3) +N(2:end-3,4:end-1))-alpha/2*(M(2:end-3,4:end-1).*N(2:end-3,4:end-1)./(H(2:end-3,4:end-1)+n(2:end-3,4:end-1))-M(2:end-3,2:end-3).*N(2:end-3,2:end-3)./(H(2:end-3,2:end-3)+n(2:end-3,2:end-3)))...
    -beta/2*(N(3:end-2,3:end-2).^2./(H(3:end-2,3:end-2)+n(3:end-2,3:end-2))-N(1:end-4,3:end-2).^2./(H(1:end-4,3:end-2)+n(1:end-4,3:end-2)))...
    -g*beta/2*(H(2:end-3,3:end-2)+n(2:end-3,3:end-2)).*(n(3:end-2,3:end-2)-n(1:end-4,3:end-2));

    N_3=1/4*(N(4:end-1,4:end-1)+N(2:end-3,4:end-1) +N(3:end-2,3:end-2) +N(3:end-2,5:end))-alpha/2*(M(3:end-2,5:end).*N(3:end-2,5:end)./(H(3:end-2,5:end)+n(3:end-2,5:end))-M(3:end-2,3:end-2).*N(3:end-2,3:end-2)./(H(3:end-2,3:end-2)+n(3:end-2,3:end-2)))...
    -beta/2*(N(4:end-1,4:end-1).^2./(H(4:end-1,4:end-1)+n(4:end-1,4:end-1))-N(2:end-3,4:end-1).^2./(H(2:end-3,4:end-1)+n(2:end-3,4:end-1)))...
    -g*beta/2*(H(3:end-2,4:end-1)+n(3:end-2,4:end-1)).*(n(4:end-1,4:end-1)-n(2:end-3,4:end-1));

    N_4=1/4*(N(4:end-1,2:end-3)+N(2:end-3,2:end-3)+N(3:end-2,3:end-2)+N(3:end-2,1:end-4))-alpha/2*(M(3:end-2,3:end-2).*N(3:end-2,3:end-2)./(H(3:end-2,3:end-2)+n(3:end-2,3:end-2))-M(3:end-2,1:end-4).*N(3:end-2,1:end-4)./(H(3:end-2,1:end-4)+n(3:end-2,1:end-4)))...
    -beta/2*(N(4:end-1,2:end-3).^2./(H(4:end-1,2:end-3)+n(4:end-1,2:end-3))-N(2:end-3,2:end-3).^2./(H(2:end-3,2:end-3)+n(2:end-3,2:end-3)))...
    -g*beta/2*(H(3:end-2,2:end-3)+n(3:end-2,2:end-3)).*(n(4:end-1,2:end-3)-n(2:end-3,2:end-3));
    
    nnew(3:end-2,3:end-2)=n(3:end-2,3:end-2)-alpha*(M_3-M_4)-beta*(N_1-N_2);
    Mnew(3:end-2,3:end-2)=M(3:end-2,3:end-2)-alpha*(M_3.^2./(H(3:end-2,4:end-1)+n_3)-M_4.^2./(H(3:end-2,2:end-3)+n_4))...
    -beta*(M_1.*N_1./(H(4:end-1,3:end-2)+n_1)-M_2.*N_2./(H(2:end-3,3:end-2)+n_2))...
    -g*alpha*(H(3:end-2,3:end-2)+(n_1+n_2+n_3+n_4)/4).*(n_3-n_4);
    Nnew(3:end-2,3:end-2)=N(3:end-2,3:end-2)-alpha*(M_3.*N_3./(H(3:end-2,4:end-1)+n_3)-M_4.*N_4./(H(3:end-2,2:end-3)+n_4))...
    -beta*(N_1.^2./(H(4:end-1,3:end-2)+n_1)-N_2.^2./(H(2:end-3,3:end-2)+n_2))...
    -g*beta*(H(3:end-2,3:end-2)+(n_1+n_2+n_3+n_4)/4).*(n_1-n_2);
    
    % Boundary condition

    if BC==0
            
        % Null displacement free velocity

        Mnew(2:end-1,end)=Mnew(2:end-1,end-2);
        Mnew(2:end-1,1)=Mnew(2:end-1,3);
        Mnew(end,2:end-1)=Mnew(end-2,2:end-1);
        Mnew(1,2:end-1)=Mnew(3,2:end-1);
        Mnew(1,1)=Mnew(3,3);
        Mnew(end,end)=Mnew(end-2,end-2);
        Mnew(1,end)=Mnew(3,end-2);
        Mnew(end,1)=Mnew(end-2,3);
        Nnew(2:end-1,end)=Nnew(2:end-1,end-2);
        Nnew(2:end-1,1)=Nnew(2:end-1,3);
        Nnew(end,2:end-1)=Nnew(end-2,2:end-1);
        Nnew(1,2:end-1)=Nnew(3,2:end-1);
        Nnew(1,1)=Nnew(3,3);
        Nnew(end,end)=Nnew(end-2,end-2);
        Nnew(1,end)=Nnew(3,end-2);
        Nnew(end,1)=Nnew(end-2,3);

    elseif BC==1
       
        % Free displacement null velocity

        nnew(2:end-1,end)=nnew(2:end-1,end-2);
        nnew(2:end-1,1)=nnew(2:end-1,3);
        nnew(end,2:end-1)=nnew(end-2,2:end-1);
        nnew(1,2:end-1)=nnew(3,2:end-1);
        nnew(1,1)=nnew(3,3);
        nnew(end,end)=nnew(end-2,end-2);
        nnew(1,end)=nnew(3,end-2);
        nnew(end,1)=nnew(end-2,3);
                
   
    end 

    n=nnew;
    M=Mnew;
    N=Nnew;
    
    if mod(k,snap)==0
        n_final(:,:,fix(k/snap))=n;
        M_final(:,:,fix(k/snap))=M;
        N_final(:,:,fix(k/snap))=N;
    end
end

close(f)

%% Ploting

% Gaussian source

figure(1);
imagesc(x*1.e-3,y*1.e-3,s);
title('Gaussian source');
xlabel('X (km)');
ylabel('Y (km)');


% Propagation - Vertical displacement

figure(2)
i=1;
for j=[1 190 250 360 440 520]
    subplot(3,2,i);
    imagesc(x*1.e-3,y*1.e-3,n_final(:,:,j));
    set(gca,'DataAspectRatio',[1 1 1]);
    colormap('jet');
    c=colorbar;
    c.Label.String = 'Vertical displacement (m)';
    set(c.Label,'Interpreter','latex');
    title('T = '+string(round(dt*snap*j))+' s','Interpreter','latex');
    xlabel('X Distance (km)','Interpreter','latex');
    ylabel('Y Distance (km)','Interpreter','latex');
    i=i+1;
end
saveas(gcf,strcat(root4,'n'),'png')

% Propagation - X Flux

figure(3)
i=1;
for j=[1 190 250 360 440 520]
    subplot(3,2,i);
    imagesc(x*1.e-3,y*1.e-3,M_final(:,:,j));
    set(gca,'DataAspectRatio',[1 1 1]);
    colormap('jet');
    c=colorbar;
    c.Label.String = 'X Flux ($m^{2}/s$)';
    set(c.Label,'Interpreter','latex');
    title('T = '+string(round(dt*snap*j))+' s','Interpreter','latex');
    xlabel('X Distance (km)','Interpreter','latex');
    ylabel('Y Distance (km)','Interpreter','latex');
    i=i+1;
end
% saveas(gcf,strcat(root4,'FX'),'png')

% Propagation - Y Flux

figure(4)
i=1;
for j=[1 190 250 360 440 520]
    subplot(3,2,i);
    imagesc(x*1.e-3,y*1.e-3,N_final(:,:,j));
    set(gca,'DataAspectRatio',[1 1 1]);
    colormap('jet');
    c=colorbar;
    c.Label.String = 'Y Flux ($m^{2}/s$)';
    set(c.Label,'Interpreter','latex');
    title('T = '+string(round(dt*snap*j))+' s','Interpreter','latex');
    xlabel('X Distance (km)','Interpreter','latex');
    ylabel('Y Distance (km)','Interpreter','latex');
    i=i+1;
end
% saveas(gcf,strcat(root4,'FY'),'png')

% Bathymetry 

figure(5)
imagesc(x*1.e-3,y*1.e-3,-H);
xlabel('X Distance (km)','Interpreter','latex')
ylabel('Y Distance (km)','Interpreter','latex')
title('Bathymetry','Interpreter','latex')
colormap('jet');
c4=colorbar;
c4.Label.String = 'Elevation (m)';
set(c4.Label,'Interpreter','latex');
set(gca,'DataAspectRatio',[1 1 1]);

% saveas(gcf,strcat(root4,'Bathymetry'),'png')

%% Initialize video

myVideo = VideoWriter(strcat(root2,'Nonlinear_2D_EXPLICIT'), 'MPEG-4'); %open video file
myVideo.FrameRate = 15;  
open(myVideo)
figure()
for j=1:length(n_final(1,1,:))
    subplot(2,2,1)
    imagesc(x*1.e-3,y*1.e-3,squeeze(n_final(:,:,j)));
    xlabel('X (km)')
    ylabel('Y (km)')
    title('2D Propagation')
    colormap('jet')
    c1=colorbar;
    c1.Label.String = 'Vertical displacement (m)';
    set(gca,'DataAspectRatio',[1 1 1]);
    text(950,100,string(round(dt*snap*j))+' s', 'Color', 'white','FontWeight','Bold','FontSize',15);
    subplot(2,2,2)
    imagesc(x*1.e-3,y*1.e-3,-H);
    xlabel('X (km)')
    ylabel('Y (km)')
    title('Bathymetry')
    colormap('jet')
    c2=colorbar;
    c2.Label.String = 'Bathymetry (m)';
    set(gca,'DataAspectRatio',[1 1 1]);
    subplot(2,2,3)
    imagesc(x*1.e-3,y*1.e-3,squeeze(M_final(:,:,j)));
    xlabel('X (km)')
    ylabel('Y (km)')
    title('2D Propagation')
    colormap('jet')
    c3=colorbar;
    c3.Label.String = 'X velocity (m/s)';
    set(gca,'DataAspectRatio',[1 1 1]);
    text(950,100,string(round(dt*snap*j))+' s', 'Color', 'white','FontWeight','Bold','FontSize',15);
    subplot(2,2,4)
    imagesc(x*1.e-3,y*1.e-3,squeeze(N_final(:,:,j)));
    xlabel('X (km)')
    ylabel('Y (km)')
    title('2D Propagation')
    colormap('jet')
    c4=colorbar;
    c4.Label.String = 'Y velocity (m/s)';
    set(gca,'DataAspectRatio',[1 1 1]);
    text(950,100,string(round(dt*snap*j))+' s', 'Color', 'white','FontWeight','Bold','FontSize',15);
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

