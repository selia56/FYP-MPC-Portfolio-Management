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
mrkt_mean = MarketMean;
mrkt_risk = sqrt(MarketVar);
cash_ret = CashMean;
cash_risk = sqrt(CashVar);


% number of risky assets
[~,n_assets] = size(AssetList);

asset_mean = AssetMean;
asset_cov = AssetCovar;

weights_f = 1/n_assets*ones(n_assets,1);
weight_f_return = asset_mean'*weights_f;
weight_f_risk = weights_f'*asset_cov*weights_f;


H = 2*asset_cov;
f = -asset_mean';
% Defining constraints
A = ones(1,n_assets);
b = 1;

% no short-selling assumption-define upper and lower bounds
lower_bound = zeros(n_assets,1);
upper_bound = ones(n_assets,1);

%Initialise parameters
iterations = 100;
w_optimal = 0;
cov_optimal = 0;
w_return_vctr = zeros(iterations,1);
w_risk_vctr = zeros(iterations,1);

% initialise Arrow Pratt risk aversion index
lambda = 0;


%Quadprog optimziation solver
for i = 1 : iterations
%slving using quadratic programming
[w_optimal, output] = quadprog(lambda*H,f,[],[],A, b, lower_bound, upper_bound);
w_return_vctr(i) = w_optimal'*asset_mean;
w_risk_vctr(i)=w_optimal'*asset_cov*w_optimal;
lambda = i/15;
end

figure;
plot(w_risk_vctr, w_return_vctr, '','LineWidth',1.5)
legend('Efficient Frontier Curve','Location','best')
title('Efficient Frontier Curve for Risk Aversion','FontSize',16)
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

print(gcf,'EF_RiskAversion','-dpdf','-fillpage')