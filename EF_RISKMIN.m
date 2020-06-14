% Solving the risk minimization portfolio optimization problem 
%using QUADPROG - quadratic programming 
clear all;
close all;
clc;
set(groot, 'defaultAxesFontSize',  15);
set(groot, 'defaultLegendFontSize',  15);
set(groot, 'defaultFigurePosition',  [0, 0, 800, 400]);
set(groot, 'defaultLegendFontSizeMode',  'manual');
set(groot,'defaultAxesLooseInset',[0,0,0,0]);

%example works with moments for monthly total returns of a universe of 30 "blue-chip" stocks.
load BlueChipStockMoments
mrkt_return = MarketMean;
mrkt_risk = sqrt(MarketVar);
cash_ret = CashMean;
cash_rsk = sqrt(CashVar);

% n_assets risky assets
%get size of the assets list i.e. num of assets 30
[~,n_assets] = size(AssetList);

mean_asset = AssetMean;
cov_asset = AssetCovar;

% equal weights function for comparison purposes 
weight_f = 1/n_assets*ones(n_assets,1);
return_weight_f = mean_asset'*weight_f; 
risk_weight_f = weight_f'*cov_asset*weight_f;
% Objective quadprog function 0.5*x'Hx + f'x subject to:  A*x <= b  and Aeq*x constraints 
% X = quadprog(H,f,A,b,Aeq,beq,LB,UB) defines a set of lower and upper

H = 2*cov_asset;
f = zeros(n_assets,1);

% Constraints:
% weights must sum to one
A = ones(1,n_assets);
b = 1;

% set no short selling constraint 
lower_bound = zeros(n_assets,1);
upper_bound = ones(n_assets,1);
% return is set to pre-defined level
A = -mean_asset';
b = 0;
% setting return parameter to be equal to weight_f with equal weights
set_return_par = abs(mean_asset'*weight_f);

%Initialise parameters
iterations = 1000;
w_optimal = 0;
s_optimal = 0;
weightd_ret_v = zeros(iterations,1);%weighted returns vector initialisation
w_rsk_v = zeros(iterations,1);%weighted risk vector initialisation

for k = 1 : iterations
% solving optimization problem resulting in optimal portfolio weights vec
[w_optimal, weightedRisk] = quadprog(H,f,A,b,A,b,lower_bound,upper_bound);
weightd_ret_v(k) = w_optimal'*mean_asset; %weighted returns vector
w_rsk_v(k)= weightedRisk;

% Set paramter constraint
b =  -set_return_par*k/(0.5*iterations);
end
% plot efficient frontier Curve with the setparamers 
figure;

plot(w_rsk_v, weightd_ret_v, '','LineWidth',1.5);

legend('Efficient Frontier Curve','Location','best')
title('Efficient Frontier Curve for Risk Minimization','FontSize',16)
xlabel('\sigma - % Standard Deviation of Returns (Annualized)','Fontsize',14)
ylabel('\mu - % Mean of Returns (Annualized)', 'Fontsize',14)

grid(gca,'minor')
grid on
ticks = get(gca,'YTick');
set(gca, 'YMinorTick','on', 'YMinorGrid','on')
set(gca,'YTickLabel',num2str(ticks'*100))
ticks = get(gca,'XTick');
set(gca,'XTickLabel',num2str(ticks'*100))

%print plot as a pdf
figures = gcf;
position = figures.PaperPosition;
figures.PaperSize = [position(3) position(4)];

print(gcf,'EF_RiskMin','-dpdf','-fillpage')
