clear all;
close all;

%Read matrices
LW=readmatrix('Nonlinear_1D_LW.txt');
CN=readmatrix('Nonlinear_1D_CN.txt');
X=readmatrix('Linear_1D_x.txt');
T=readmatrix('Linear_1D_t.txt');


%Plot
figure(1)
for i=1:length(T)*10
    plot(X,CN(:,i),X,LW(:,i));
    xlabel('Distance (m)');
    ylabel('Amplitude');
    title('Nonlinear propagation');
    ylim([0 1]);
    pause(0.01);
end