clc; clear all; 

n=10;           %Number of Assets  
p_actual=1;           %Lag order
q_actual=3;            %Number of Factors
N=20;          %Number of Days

T=[2:N+1];

%Initialisations 
e_Cexact=[];
e_Cest=[];
Ynorm_actual=[];
Ynorm_Cexact=[];
Ynorm_Cest=[];

%% Number of Days

tic
for i=1:N
    [Ynorm_actual(:,i),Ynorm_Cexact(:,i),Ynorm_Cest(:,i), e_Cexact(:,i),e_Cest(:,i),C_exact(:,i),C0_FLM(:,i),e_Cmatrix(:,i)]=modelData_pest(T(i),n,q_actual,p_actual);
end
toc

%Ynorm_actual(1:5:N)
%Ynorm_Cexact(1:5:N)
%Ynorm_Cest(1:5:N)
figure();
plot([1:N],Ynorm_Cest,'r','Linewidth',1);
hold on;
plot([1:N],Ynorm_actual,'m','Linewidth',0.75);
plot([1:N],Ynorm_Cexact,'-.b','Linewidth',1);
title (['Norm. Returns for ',num2str(N),' days,' num2str(n),' assets,', 'q=', num2str(q_actual), 'p=', num2str(p_actual)],'Linewidth',10);
xlabel('Number of Days');
ylabel('Normalised Returns');
ax=gca;
ax.YAxis.Exponent = 0;
grid(gca,'minor');
grid on;
legend('Y_{Cest}','Y_{actual}','Y_{Cexact}','Location', 'Best', 'Fontsize',11,'NumColumns',3);
hold off;


%Save figure plot as a pdf
figures = gcf;
position = figures.PaperPosition;
figures.PaperSize = [position(3) position(4)];
%print(gcf, 'Testing_InitialReturns_BreakError','-dpdf','-fillpage');


figure();
plot([1:N],e_Cexact,'b','Linewidth',1.5);
hold on;
plot([1:N],e_Cest,'r', 'Linewidth', 1.0);

title (['Norm. Error for ',num2str(N),' days,' num2str(n),' assets,','q=', num2str(q_actual), 'p=', num2str(p_actual)],'Linewidth',10);
xlabel('Number of Days');
ylabel('Normalised Error');
ax=gca;
ax.YAxis.Exponent = 0;
grid(gca,'minor');
grid on;
legend('e_{Cexact}','e_{Cest}','Location', 'northwest','NumColumns',2, 'Fontsize',11)
hold off; 


%Save figure plot as a pdf
figures = gcf;
position = figures.PaperPosition;
figures.PaperSize = [position(3) position(4)];
%print(gcf, 'Testing_InitialReturnsError_BreakError','-dpdf','-fillpage');

figure();
plot([1:N],e_Cmatrix,'Linewidth',1);
title('Normalised Error between exact and estimated Factor Loading Matrix C','Fontsize',11);
xlabel('Number of Days');
ylabel('Normalised Error');
ax=gca;
ax.YAxis.Exponent = 0;
grid(gca,'minor');
grid on;
legend('e_{Cmatrix}','Location', 'Best')

%Save figure plot as a pdf
figures = gcf;
position = figures.PaperPosition;
figures.PaperSize = [position(3) position(4)];
%print(gcf, 'Testing_InitialNormCmatrix_NoiseFreeError','-dpdf','-fillpage');

