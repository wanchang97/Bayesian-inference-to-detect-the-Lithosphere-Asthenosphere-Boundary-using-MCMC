function TraceHistogramlnLPlot(lnL,figurepath,fileName)


%% Trace plot and the histogram plot

figure
subplot(1,2,1) % histogram plot
histogram(lnL,'Normalization','probability');
hold on

xlabel('steps')
%ylabel('steps')
title('histogram plot with burn-in period')
camroll(90)

subplot(1,2,2) % trace plot
plot(lnL);
hold on
xlabel('steps');
ylabel('log likelihood');
title('trace plot')
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);

end