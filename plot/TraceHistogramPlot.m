function TraceHistogramPlot(M,true_m,delta,m0,figurepath,fwd,fileName)


%% Trace plot and the histogram plot
for i = 1:fwd.nm
    figure
    subplot(1,2,1) % histogram plot
    histogram(M(:,i),'Normalization','probability');
    hold on
%     plot(true_m(i),0,'*','MarkerSize',20)
%     plot(true_m(i)+delta(1,i),0,'|','MarkerSize',20)
%     plot(true_m(i)-delta(2,i),0,'|','MarkerSize',20)
    xline(true_m(i),'r','LineStyle','-')
    xline(true_m(i)+delta(1,i),'r','LineStyle','-.')
    xline(true_m(i)-delta(2,i),'r','LineStyle','-.')
    xlabel(['m_',num2str(i)])
    ylabel(['\pi(m_',num2str(i),')'])
    title('histogram plot with burn-in period')
    camroll(90)

    subplot(1,2,2) % trace plot
    plot(M(:,i));
    hold on
    plot(0,m0(i),'ro');
    yline(true_m(i),'r','LineStyle','-')
    yline(true_m(i)+delta(1,i),'r','LineStyle','-.')
    yline(true_m(i)-delta(2,i),'r','LineStyle','-.')
    xlabel('steps');
    ylabel(['m',num2str(i)]);
    title('trace plot')
    fn = fullfile(figurepath, ['m',num2str(i)],fileName);
    saveas(gcf,fn);
end

end