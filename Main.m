clc; clear all; 

n=10;           %Number of Assets  
p_actual=2;           %Lag order
q_actual=2;            %Number of Factors
N=30;          %Number of Days

T=[2:N+1];
%T=[1:q_actual];
%Initialisations 
e_Cexact=[];
e_Cest=[];
Ynorm_actual=[];
Yhat_Cexact=[];
Yhat_Cest=[];

%% Number of Days
q_test=4;
p_test=4;

tic
for i=1:N
    for q_est=1:q_test
        [Ynorm_actual(:,i),Ynorm_Cexact(:,i),Ynorm_Cest(:,i,q_est),e_Cexact(:,i),e_Cest(:,i,q_est),epsilon(:,i)]=modelData_pqtest(T(i),n,q_actual,p_actual,N,q_test,p_test,q_est); 
    end
end
toc

[p_estimate, q_estimate, p_est_min,diff_vec]=estimate_pq(e_Cest,N);

error_est=mean(abs(Ynorm_actual-Ynorm_Cest(p_estimate,:,q_estimate)))
%fprintf('The estimated factors are found to be q_est= %.0f\n',q_estimate);
%fprintf('The estimated lag order is found to be p_est= %.0f\n',p_estimate);

%{
q_estimate=2;
p_estimate=2;
for i=1:N
    [Ynorm_actual(:,i),Ynorm_Cexact(:,i),Ynorm_Cest(:,i), e_Cexact(:,i),e_Cest(:,i)]=modelData_pest(T(i),n,q_estimate,p_estimate);
end
%}

figure();
plot([1:N],Ynorm_Cest(p_estimate,:,q_estimate),'r','Linewidth',1);
hold on;
plot([1:N],Ynorm_actual,'m','Linewidth',0.75);
plot([1:N],Ynorm_Cexact,'-.b','Linewidth',1);
plot([1:N],Ynorm_Cest(p_est_min,:,q_estimate),'-.g','Linewidth',0.75);
title(['Norm. Returns for ',num2str(N),' days,' num2str(n),' assets,',' $$\hat{q}=$$', num2str(q_estimate), ' $$\hat{p}=$$', num2str(p_estimate), ' q=', num2str(q_actual), ' p=',num2str(p_actual) ],'Interpreter','Latex','Linewidth',10);
xlabel('Number of Days');
ylabel('Normalised Returns');
ax=gca;
ax.YAxis.Exponent = 0;
grid(gca,'minor');
grid on;
legend(strcat('$$Y_{Cest}, \hat{p}= $$',num2str(p_estimate)),'$$ Y_{actual} $$','$$ Y_{Cexact} $$',strcat('$$ Y_{Cest}, \hat{p}_{min}= $$',num2str(p_est_min)),'Location', 'Best','Interpreter','Latex','Fontsize',10, 'NumColumns', 3);

hold off;


%Save figure plot as a pdf
figures = gcf;
position = figures.PaperPosition;
figures.PaperSize = [position(3) position(4)];
print(gcf, 'ImprovedReturns_TimeHorizon50_HighError','-dpdf','-fillpage');


figure();
plot([1:N],e_Cexact,'b','Linewidth',1.5);
hold on;
plot([1:N],e_Cest(p_estimate,:,q_estimate),'r', 'Linewidth', 1.0);
plot([1:N],e_Cest(p_est_min,:,q_estimate),'-.g', 'Linewidth', 0.75);
title(['Norm. Error for ',num2str(N),' days,' num2str(n),' assets,',' $$\hat{q}=$$', num2str(q_estimate), ' $$\hat{p}=$$', num2str(p_estimate), ' q=', num2str(q_actual), ' p=',num2str(p_actual) ],'Interpreter','Latex','Linewidth',10);
xlabel('Number of Days');
ylabel('Normalised Error in Returns');
ax=gca;
ax.YAxis.Exponent = 0;
grid(gca,'minor');
grid on;
legend(' $$e_{Cexact}$$',strcat(' $$e_{Cest} \hat{p}=$$',num2str(p_estimate)),strcat(' $$e_{Cest} \hat{p}_{min}$$=', num2str(p_est_min)),'Interpreter','Latex','Location', 'Best')
hold off; 


%Save figure plot as a pdf
figures = gcf;
position = figures.PaperPosition;
figures.PaperSize = [position(3) position(4)];
print(gcf, 'ImprovedReturnsError_TimeHorizon50_HighError','-dpdf','-fillpage');
figure;
plot([1:N],epsilon)
