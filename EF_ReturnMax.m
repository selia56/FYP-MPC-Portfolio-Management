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
cash_return = CashMean;
cash_rsk = sqrt(CashVar);
% n_assets risky assets
%get size of the assets list i.e. num of assets 30
[~,n_assets] = size(AssetList);

mean_asset = AssetMean;
cov_asset = AssetCovar;
% assign a function with same weights
weight_f = 1/n_assets*ones(n_assets,1);
return_weight_f = mean_asset'*weight_f;
risk_weight_f = weight_f'*cov_asset*weight_f;

%Initialize parameters
iterations = 80; 
% set return equal to value of an equally weighted portfolio
set_mrkt_risk_par = abs(weight_f'*AssetCovar*weight_f);
w_risk_vctr = zeros(iterations,1);
w_return_vctr = zeros(iterations,1);
w_p = zeros(n_assets,1);

%initialise max risk level function
max_risk_lvl_fun=0; 
k_var=4.3; %variable for max risk level fun
k_const=0.0015; %constant for max risk level fun
for j = 1 : iterations
%setting max risk level function = k+ar
max_risk_lvl_fun = k_const + set_mrkt_risk_par*k_var*j/iterations;
epsilon = ones(n_assets,1);
cvx_begin sdp quiet
cvx_precision high
variable w_p(n_assets);
maximize(mean_asset'*w_p);

%adding constraints
subject to
epsilon'*w_p == 1;
for k=1:n_assets
    w_p(k)>=0;
end
[max_risk_lvl_fun, w_p'; w_p, inv(cov_asset)] >= 0;
cvx_end
w_risk_vctr(j)  = w_p'*cov_asset*w_p;
w_return_vctr(j) = w_p'*mean_asset;
epsilon = (eig([max_risk_lvl_fun, w_p'; w_p, inv(cov_asset)]));
end

% plot efficient frontier curve for return maximization
figure();
plot(w_risk_vctr, w_return_vctr, '','LineWidth',1.5)
legend('Efficient Frontier Curve','Location','best')
title('Efficient Frontier Curve for Returns Maxmization','FontSize',16)
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

print(gcf,'EF_ReturnMax','-dpdf','-fillpage')