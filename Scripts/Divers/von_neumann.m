clear all;
close all;



%VON NEUMANN STABILITY ANALYSIS EXPLICIT

%ROOT
%%%%%

root2 = '/Users/francoisdallasta/Documents/Numerical_project/Figures/Linear_1D/VON_NEUMANN/'; % Save figure


x=-5:0.001:5;
a=ones(1,length(x));
b=-2*x;
c=ones(1,length(x));
delta=b.^2-4*a.*c;
N1=[];
N2=[];
R1=[];
R2=[];

for j=1:length(x)
    if delta(j) < 0
      C1=(-b(j)-1i*sqrt(-delta(j)))/(2*a(j));
      C2=(-b(j)+1i*sqrt(-delta(j)))/(2*a(j));
      N1=[N1,sqrt(real(C1)^2+imag(C1)^2)];
      N2=[N2,sqrt(real(C2)^2+imag(C2)^2)];
      R1=[R1,NaN];
      R2=[R2,NaN];
    else
      R1=[R1,(-b(j)-sqrt(delta(j)))/(2*a(j))];
      R2=[R2,(-b(j)+sqrt(delta(j)))/(2*a(j))];
      N1=[N1,NaN];
      N2=[N2,NaN];
    end
end


%PLOTTING
figure(1)
hold on;
plot(x,N1,'black','LineWidth',2,'DisplayName','Norm Complex roots');
plot(x,R1,'red','LineWidth',2,'DisplayName','Real root 1');
plot(x,R2,'blue','LineWidth',2,'DisplayName','Real root 2');
ylim([-5,5]);
hold off;
xlabel('$\beta$','Interpreter','latex')
ylabel('$G$','Interpreter','latex')
title('$G^2-2 \beta G+1=0$','Interpreter','latex')
legend('Location','Southeast');
saveas(gcf,strcat(root2,'EXPLICIT'),'png')

%VON NEUMANN STABILITY ANALYSIS IMPLICIT

x=-5:0.001:5;
a=x;
b=-4+2*x;
c=x;
delta=b.^2-4*a.*c;
N1=[];
N2=[];
R1=[];
R2=[];

for j=1:length(x)
    if delta(j) < 0
      C1=(-b(j)-1i*sqrt(-delta(j)))/(2*a(j));
      C2=(-b(j)+1i*sqrt(-delta(j)))/(2*a(j));
      N1=[N1,sqrt(real(C1)^2+imag(C1)^2)];
      N2=[N2,sqrt(real(C2)^2+imag(C2)^2)];
      R1=[R1,NaN];
      R2=[R2,NaN];
    else
      R1=[R1,(-b(j)-sqrt(delta(j)))/(2*a(j))];
      R2=[R2,(-b(j)+sqrt(delta(j)))/(2*a(j))];
      N1=[N1,NaN];
      N2=[N2,NaN];
    end
end


%PLOTTING
figure(2)
hold on;
plot(x,N1,'black','LineWidth',2,'DisplayName','Norm Complex roots');
plot(x,R1,'red','LineWidth',2,'DisplayName','Real root 1');
plot(x,R2,'blue','LineWidth',2,'DisplayName','Real root 2');
ylim([-5,5]);
hold off;
xlabel('$\beta$','Interpreter','latex')
ylabel('$G$','Interpreter','latex')
title('$\beta G^2+(2\beta-4)G+\beta=0$','Interpreter','latex')

legend;
saveas(gcf,strcat(root2,'IMPLICIT'),'png')