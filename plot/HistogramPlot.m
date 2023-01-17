function HistogramPlot(M,true_m,delta,error_m,nburn,n,m,acc,figurepath,fwd,stringinput,sigmam,sigmaObratio,fileName,priorInfo)
fig = figure;
k = fwd.k;
switch fwd.mtype
    case 0
        for i = 1:k
            subplot(1,k+1,i)
            histogram(M(nburn:end,i),'Normalization','probability');
            hold on
            title(['histogram of d',num2str(i)])
            xline(true_m(i),'r','LineStyle','-')
            xline(true_m(i)+delta(1,i),'r','LineStyle','-.')
            xline(true_m(i)-delta(2,i),'r','LineStyle','-.')
            txtd = ['Error_{d_' num2str(i) '} = ',num2str(error_m(i),'%.4f%%')];
            text(0.6,0.7,txtd,'Units','normalized')
            xlabel(stringinput(i));
        end
        % add relative errors
        subplot(1,k+1,k+1)
        text(0,1,['Acceptance ratio = ',num2str(acc/n)]);
        text(0,1-1/(2*k)*(k),['\sigma_{RW_d} = ',num2str(sigmam.d)]);
        text(0,0,['\sigma_{ob} = ',num2str(sigmaObratio) '*true_{ob}']);
        axis('off')
    case 1
        B = M(:,1:k-1);
        D = M(:,k:end);
        b = m(1:k-1);
        d = m(k:end);
        error_b = error_m(1:k-1);
        error_d = error_m(k:end);
        true_b = true_m(1:k-1);
        true_d = true_m(k:end);
        sigmab = sigmam.b;
        sigmad = sigmam.d;
        for i = 1:k-1
            subplot(2,k,i)
            histogram(B(nburn:end,i),'Normalization','probability');
            hold on
            title(['histogram of b',num2str(i)])
            xline(true_b(i),'r','LineStyle','-')
            xline(true_b(i)+delta(1,i),'r','LineStyle','-.')
            xline(true_b(i)-delta(2,i),'r','LineStyle','-.')
            xlabel(stringinput(i))
        end
        for i = 1:k
            subplot(2,k,k+i)
            histogram(D(nburn:end,i),'Normalization','probability');
            hold on
            title(['histogram of d',num2str(i)])
            xline(true_d(i),'r','LineStyle','-')
            xline(true_d(i)+delta(1,i+k-1),'r','LineStyle','-.')
            xline(true_d(i)-delta(2,i+k-1),'r','LineStyle','-.')
            xlabel(stringinput(i+k-1))
        end
        % add relative errors
        subplot(2,k,k)
        text(-0.3,1,'Relative error of b :')
        for i = 1: k-1
            txtb =['Error_{b_' num2str(i) '} = ',num2str(error_b(i),'%.4f%%')];
            text(-0.2,1-1/(2*k)*(i),txtb)
        end
        text(-0.3,1-1/(2*k)*(k),'Relative error of d :')
        for i = 1:k
            txtd = ['Error_{d_' num2str(i) '} = ',num2str(error_d(i),'%.4f%%')];
            text(-0.2,1-1/(2*k)*(i+k),txtd)
        end
        text(0.5,1,['Acceptance ratio = ',num2str(acc/n)]);
        text(0.5,0.8,['\sigma_{RW_b} = ',num2str(sigmab)]);
        text(0.5,1-1/(2*k)*(k),['\sigma_{RW_d} = ',num2str(sigmad)]);
        text(0.5,0,['\sigma_{ob} = ',num2str(sigmaObratio) '*true_{ob}']);
        axis('off')
end
han = axes(fig,'visible','off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
sgtitle('Normalized histogram plot of RW samples without burn-in period')
fig.Position = [100 100 1000 600];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
end