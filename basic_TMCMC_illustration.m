% Illustration of basic TMCMC
% understanding the literature Ching, J. and Chen, Y.C., 2007. Transitional Markov chain Monte Carlo method for Bayesian model updating, model class selection, and model averaging. Journal of engineering mechanics, (7), pp.816-832.
% Generate 50 random number
clear all; close all; clc
seed= 0;N = 50;
rng(seed);
x = rand(N,2);
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
weights = rand(N,1); % assume that the weights are random
weights_normalized = weights / sum(weights,1);
%% generate samples following the weights
L = cumsum(weights_normalized);
X_update = zeros(N,2);
index = zeros(N,1);
for i = 1:N
    index(i) = find(rand <= L,1);
    X_update(i,:) = x(index(i),:);
end
%% Plot the samples when t = 1
figure
subplot(1,2,1)
txt = 't = 0, p_0 = 0, N_0 = 50';
scatter(x(:,1),x(:,2),'black','o');
hold on
text(0,1.2,txt,'FontSize', 15);
xlim([-0.3,1.3]);
ylim([-0.3,1.5]);
set(gca,'visible','off')

subplot(1,2,2)
txt = 't = 1, p_1 , N_1 = 50';
scatter(x(:,1),x(:,2),'black','o');
hold on
scatter(X_update(:,1),X_update(:,2),'red','*');
legend('$\{\textbf{m}_{0,j}:j = 1,\cdots,N_0\}$','$\{\textbf{m}_{1,j}:j = 1,\cdots,N_1\}$','Interpreter','latex','Location','Best');
text(0,1.2,txt,'FontSize', 15);
xlim([-0.3,1.3]);
ylim([-0.3,1.5]);
set(gca,'visible','off')
%% plot the first two stages using panel
p = panel();
p.pack(1,2);
p(1,1).select();
scatter(x(:,1),x(:,2),'black','o');
txt = '$t = 0, p_0 = 0, N_0 = 50$';
title(txt,'Interpreter','latex','FontSize',22);
axis([-0.3,1.3 -0.05 1.15]);
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[])
set(gca,'xtick',[]); set(gca,'ytick',[])
set(gca,'box','on')
%
p(1,2).select();
scatter(x(:,1),x(:,2),'black','o');
hold on
scatter(X_update(:,1),X_update(:,2),'b','p');
legend('$\{\textbf{m}_{0,j}:j = 1,\cdots,N_0\}$','$\{\textbf{m}_{1,j}:j = 1,\cdots,N_1\}$','Interpreter','latex','Location','Best');
txt = '$t = 1, p_1 , N_1 = 50$';
title(txt,'Interpreter','latex','FontSize',22);
axis([-0.3,1.3 -0.05 1.15]);
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[])
set(gca,'xtick',[]); set(gca,'ytick',[])
set(gca,'box','on')
p.fa.margin = 2;
p.fa.margintop = 12;
p.fontsize = 20;
saveas(gcf,'basic_TMCMC_illustration_panel_firsttwo.png');

%% Code sequentially for t = 0,1,...,s
s = 5;
X = zeros(s+1,N,2);
X(1,:,:) = x;
Index = zeros(s+1,N);
Index(1,:) = linspace(1,N,N);
Numberofunique = zeros(s+1,1);
x_table = table(x(:,1),x(:,2));
[C,ia,ic]= unique(x_table);
Numberofunique(1) = length(ia);
for t = 2:s+1
    weights = rand(N,1); % assume that the weights are random
    weights_normalized = weights / sum(weights,1);
    L = cumsum(weights_normalized);
    X_update = zeros(N,2);
    index = zeros(N,1);
    for i = 1:N
        index(i) = find(rand <= L,1);
        X_update(i,:) = X(t-1,index(i),:);
    end
    X_update_table = table(X_update(:,1),X_update(:,2));
    [C,ia,ic]= unique(X_update_table,"stable");
    num = length(ia);
    Numberofunique(t) = num;
    X(t,:,:) = X_update;
    Index(t,:) = index;
end

%% Plot 5 resampling stages
% color
C = {'k','b','r','g','m','c','#7E2F8E','#D95319','y','c',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.
M = {'o','p','h','x','*','.','d','s','^','>'};
figure('units','normalized','outerposition',[0 0 1 0.7])
for t = 1:s+1
    subplot(1,s+1,t)
    if t > 1
        scatter(X(t-1,:,1),X(t-1,:,2),C{t-1},M{t-1});
        hold on
    end
    scatter(X(t,:,1),X(t,:,2),C{t},M{t});
    txt_t = ['t = ' num2str(t-1) '; $N_{diff}$ = ' num2str(Numberofunique(t))];
    txt_p = ['$p_' num2str(t-1) '$'];
    txt_N = ['$N_' num2str(t-1) '$'];
    text(-0.2,-0.2,txt_t,'FontSize', 15,'Interpreter','Latex');
    % text(0,1.4,txt_p,'FontSize', 15,'Interpreter','latex');
    % text(0,1.2,txt_N,'FontSize', 15,'Interpreter','latex');
    xlim([-0.3,1.3]);
    ylim([-0.3,1.5]);
    txt_mt1 = ['$\{\textbf{m}_{' num2str(t-1) ',j}:j = 1,\cdots,N_' num2str(t-1) '\}$'];
    txt_mt2 = ['$\{\textbf{m}_{' num2str(t) ',j}:j = 1,\cdots,N_' num2str(t) '\}$'];
    legend(txt_mt1,txt_mt2,'Interpreter','latex');
    set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[])
    set(gca,'xtick',[]); set(gca,'ytick',[])
    set(gca,'box','on')
end
sgtitle('basic TMCMC illustration')
saveas(gcf,'basic_TMCMC_illustration.png');


%% Plot use panel
figure('units','normalized','outerposition',[0 0 1 0.7])
p = panel();
p.pack(1,s+1);
% plot into each panel in turn
for t = 1:s+1
    p(1,t).select();
    scatter(X(t,:,1),X(t,:,2),C{t},M{t});
    hold on
    if t > 1
        %scatter(X(t-1,:,1),X(t-1,:,2),C{t-1},M{t-1});
        scatter(X(t-1,:,1),X(t-1,:,2),C{t-1},M{t-1});
    end
    
    
    txt_t = ['t = ' num2str(t-1) '; $N_{diff}$ = ' num2str(Numberofunique(t))];
    txt_p = ['$p_' num2str(t-1) '$'];
    txt_N = ['$N_' num2str(t-1) '$'];
    title(txt_t,'Interpreter','latex','FontSize',22)
    axis([-0.3,1.3 -0.05 1.15]);
    txt_mt1 = ['$\{\textbf{m}_{' num2str(t-2) ',j}:j = 1,\cdots,N_' num2str(t-2) '\}$'];
    txt_mt2 = ['$\{\textbf{m}_{' num2str(t-1) ',j}:j = 1,\cdots,N_' num2str(t-1) '\}$'];
    legend(txt_mt2,txt_mt1,'Interpreter','latex','FontSize',19);
    set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[])
    set(gca,'xtick',[]); set(gca,'ytick',[])
    set(gca,'box','on')
end
p.fa.margin = 2;
p.fa.margintop = 12;
p.fa.marginbottom= 12;
p.fontsize = 19;
saveas(gcf,'basic_TMCMC_illustration_panel.png');
