function histplot_results_RW(M,figurepath,nRW,sigmaRW_ratio)
E = mean(M')'; % equivalent to mean(M,2)
s = std(M')';
nOfParam = size(M,1);
figure
for i = 1: nOfParam
    subplot(nOfParam,1,i);
    histogram(M(i,:),'Normalization','pdf','FaceAlpha',0.1,'EdgeAlpha',0.1);
    title("M(" + i + ")")
    hold on
    xline([E(i)-s(i) E(i) E(i)+s(i)],"-",{"E_{RW}-\sigma_{M}","E_{RW}","E_{RW}+\sigma_{M}"});
end
sgtitle("E_{RW} with nRW = " + nRW + " and \sigma_{RW} = " + sigmaRW_ratio + "m0")
fileName = ['Histgram plot of Random Walk Sampler with nRW = ' num2str(nRW) ' and sigmaRW_ratio = ' num2str(sigmaRW_ratio)  ' m0.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
