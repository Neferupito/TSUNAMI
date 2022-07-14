clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TSUNAMI MODEL ANALYTICAL SOLUTION LINEAR 1D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROOTS

%root='/Users/francoisdallasta/Documents/Numerical_project/txt_files/'; % Txt files folder

%% Parameters

% Space

dx=1250; % Spatial step (m)
L=1000000; % Model length (m)
x=0:dx:L; % Space vector
Nx=length(x); % Number of x samples

% Time

dt=10; % Time step (s)
T=4000; % Model duration (s)
t=0:dt:T; % Time vector

% Constants

g=9.81; % Gravity acceleration (N/kg)
H=1000; % (Constant!!!) Bathymetry (+ down) (m)
c=sqrt(g*H); % Propagation speed (m/s)

% Gaussian source
% 6 sigma = source width

width=100000; % Gaussian width (m)
sigma=width/6; % Gaussian standard deviation (m)
center=round((x(end)-x(1))/2); % Gaussian mean (m)
h=1; % Gaussian height (m)

%% Calculations

n=zeros(length(x),length(t)); % Vertical displacement storage vector

for j=1:length(t) % Iteration over time
    n(:,j)=h*1/2*(exp(-(x+c*t(j)-center).^2./(2*sigma^2))+exp(-(x-c*t(j)-center).^2/(2*sigma^2)));
end

%% Ploting

% Propagation - Vertical displacement

figure(1)
hold on;
for j=[1 16 61 141 261 381]
    plot(x*1.e-3,n(:,j),'DisplayName','t = '+string(t(j))+' s');
    title('Propagation');
    xlabel('Distance (km)');
    ylabel('Vertical displacement (m)');
    ylim([0 1*h]);
end
hold off;
legend;


%% Save matrix as txt files

%writematrix(x,strcat(root,'Linear_1D_x.txt'));
%writematrix(t,strcat(root,'Linear_1D_t.txt'));
%writematrix(n,strcat(root,'Linear_1D_ANALYTICAL.txt'));
