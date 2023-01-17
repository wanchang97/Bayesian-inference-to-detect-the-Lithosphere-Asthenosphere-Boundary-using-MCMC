% Illustration of modified TMCMC 
% understanding the literature Ching, J. and Chen, Y.C., 2007. Transitional Markov chain Monte Carlo method for Bayesian model updating, model class selection, and model averaging. Journal of engineering mechanics, (7), pp.816-832.
% Generate 50 random number
clear all; close all;clc
seed = 0; N = 50;
rng(seed);
x = rand(N,2);
X1 = x;
% %% Plot the samples when t = 0
% figure
% txt = 't = 0, p_0 = 0, N_0 = 50';
% scatter(x(:,1),x(:,2),'black','o');
% hold on
% text(0,1.2,txt,'FontSize', 15);
% xlim([-0.3,1.3]);
% ylim([-0.3,1.5]);
% set(gca,'visible','off')
%% Add weights
weights = rand(N,1); % assume that the weights are random, but it should be calculated from the likelihood function
weights_normalized = weights / sum(weights,1);
%% generate samples following the weights
L = cumsum(weights_normalized);
X_update = zeros(N,2);
index = zeros(N,1);
% assume the target pdf
alpha = 2.43;beta = 1;
delta = .05;
%pdf = @(x)gampdf(x,alpha,beta); % Target distribution
pdf = @(x) mvnpdf(x);
proppdf = @(x,y) unifpdf(y-x,-delta,delta);
proprnd = @(x) x + rand*2*delta; 
chosen_index = 0;
X_chosen = zeros(N,2);
for i = 1:N
    index(i) = find(rand <= L,1);
    x_chosen = x(index(i),:);
    x_MHsample = mhsample(x_chosen,1,'pdf',pdf,'proprnd',proprnd,'symmetric',1);
    x(index(i),:) = x_MHsample;
    X_chosen(i,:) = x_chosen;
    X_update(i,:) = x_MHsample;
end
%% Analyse the samples
figure('units','normalized','outerposition',[0 0 1 0.7])
p = panel();
p.pack(1,4);
p(1,1).select();
scatter(X1(:,1),X1(:,2),'black','o');
txt = '$t = 0, p_0 = 0, N_0 = 50$';
title(txt,'Interpreter','latex','FontSize',22);
axis([-0.3,1.3 -0.05 1.15]);
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[])
set(gca,'xtick',[]); set(gca,'ytick',[])
set(gca,'box','on')
%
p(1,2).select();
scatter(X1(:,1),X1(:,2),'black','o');
hold on
scatter(X_chosen(:,1),X_chosen(:,2),'red','s');
legend('$\{\textbf{m}_{0,j}:j = 1,\cdots,N_0\}$','chosen samples','Interpreter','latex','Location','Best');
txt = 'Choose samples using weights';
title(txt,'Interpreter','latex','FontSize',22);
axis([-0.3,1.3 -0.05 1.15]);
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[])
set(gca,'xtick',[]); set(gca,'ytick',[])
set(gca,'box','on')
%
p(1,3).select();
scatter(X_chosen(:,1),X_chosen(:,2),'red','s');
hold on
scatter(X_update(:,1),X_update(:,2),'b','p');
legend('chosen samples','$\{\textbf{m}_{1,j}:j = 1,\cdots,N_1\}$','Interpreter','latex','Location','Best');
txt = 'Generate MH sample on the chosen sample';
title(txt,'Interpreter','latex','FontSize',22);
axis([-0.3,1.3 -0.05 1.15]);
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[])
set(gca,'xtick',[]); set(gca,'ytick',[])
set(gca,'box','on')
%
p(1,4).select();
scatter(X_update(:,1),X_update(:,2),'b','p');
legend('$\{\textbf{m}_{1,j}:j = 1,\cdots,N_1\}$','Interpreter','latex','Location','Best');
txt = '$t = 1, p_1, N_1 = 50$';
title(txt,'Interpreter','latex','FontSize',22);
axis([-0.3,1.3 -0.05 1.15]);
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[])
set(gca,'xtick',[]); set(gca,'ytick',[])
set(gca,'box','on')
% p(1,3).select();
% scatter(X1(:,1),X1(:,2),'black','o');
% hold on
% scatter(X_chosen(:,1),X_chosen(:,2),'red','s');
% scatter(X_update(:,1),X_update(:,2),'b','p');
% legend('$\{\textbf{m}_{0,j}:j = 1,\cdots,N_0\}$','chosen samples','$\{\textbf{m}_{1,j}:j = 1,\cdots,N_1\}$','Interpreter','latex','Location','Best');
% txt = 't = 1, p_1 , N_1 = 50';
% title(txt,'Interpreter','latex','FontSize',22);
% axis([-0.3,1.3 -0.05 1.15]);
% set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[])
% set(gca,'xtick',[]); set(gca,'ytick',[])
% set(gca,'box','on')
p.fa.margin = 4;
p.fa.margintop = 12;
p.fontsize = 19;
saveas(gcf,'modified_TMCMC_illustration_panel_firsttwo.png');

