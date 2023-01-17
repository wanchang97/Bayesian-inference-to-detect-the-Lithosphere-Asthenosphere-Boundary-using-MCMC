function HistogramPlot1(M,true_m,delta,error_m,nburn,n,m,acc,figurepath,fwd,stringinput,sigmam,sigmaObratio,fileName,priorInfo)
fig = figure;
k = fwd.k;
switch fwd.mtype
    case 0
        for i = 1:k
            subplot(1,k+1,i)
            h = histfit(M(nburn:end,i));
            h(1).YData = h(1).YData/max(h(1).YData);
            h(2).YData = h(2).YData/max(h(2).YData);
            h(1).FaceColor = [0 0.4470 0.7410];
            h(2).Color = [0 1 1];
            %histogram(M(nburn:end,i),'Normalization','probability');
            hold on
            % Plot the prior pdf
            xmin = priorInfo.mdist(1).Par(1);
            xmax = priorInfo.mdist(1).Par(2);
            x = xmin:(xmax-xmin)/1000:xmax;
            y = priorInfo.mdist(1).pdf(x);
            plot(x,y,'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
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
            h = histfit(B(nburn:end,i));    
            h(1).YData = h(1).YData/max(h(1).YData);
            h(2).YData = h(2).YData/max(h(2).YData);
            h(1).FaceColor = [0 0.4470 0.7410];
            h(2).Color = [0 1 1];
            %histogram(B(nburn:end,i),'Normalization','probability');
            hold on
            % Plot the prior pdf
            xmin = priorInfo.mdist(1).Par(1);
            xmax = priorInfo.mdist(1).Par(2);
            x = xmin:(xmax-xmin)/1000:xmax;
            y = priorInfo.mdist(1).pdf(x);
            plot(x,y,'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
            title(['histogram of b',num2str(i)])
            %plot(true_b(i,1),0,'rx')
            xline(true_b(i),'r','LineStyle','-')
            xline(true_b(i)+delta(1,i),'r','LineStyle','-.')
            xline(true_b(i)-delta(2,i),'r','LineStyle','-.')
            %xlim([fwd.param.bmin,fwd.param.bmax])
            xlabel(stringinput(i))
        end
        for i = 1:k
            subplot(2,k,k+i)
            h = histfit(D(nburn:end,i));
            h(1).YData = h(1).YData/max(h(1).YData);
            h(2).YData = h(2).YData/max(h(2).YData);
            h(1).FaceColor = [0 0.4470 0.7410];
            h(2).Color = [0 1 1];
            %histogram(D(nburn:end,i),'Normalization','probability');
            hold on
            % Plot the prior pdf
            xmin = priorInfo.mdist(2).Par(1);
            xmax = priorInfo.mdist(2).Par(2);
            x = xmin:(xmax-xmin)/1000:xmax;
            y = priorInfo.mdist(2).pdf(x);
            plot(x,y,'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);

            title(['histogram of d',num2str(i)])
            xline(true_d(i),'r','LineStyle','-')
            xline(true_d(i)+delta(1,i+k-1),'r','LineStyle','-.')
            xline(true_d(i)-delta(2,i+k-1),'r','LineStyle','-.')
            %plot(true_d(i,1),0,'rx')
            %xlim([fwd.param.dmin,fwd.param.dmax])
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
%title(han,'Normalized histogram plot of RW samples','Position',[0.5, -0.1, 0])
sgtitle('Normalized histogram plot of RW samples without burn-in period')
fig.Position = [100 100 1000 600];
%fileName = ['normalized histogram plot of RW samples without burn-in period np = ' num2str(fwd.np) ' sigmaRatio = ' num2str(sigmaratio) ' sigmaObRatio = ' num2str(sigma_ob_ratio) ' n = ' num2str(n) '.png'];
%fileName = ['normalized histogram plot of RW samples nRW = ' num2str(nRW) ' sigmaratio = ' num2str(sigmaratio) ' np = ' num2str(fwd.np) ' and mtype = ' num2str(fwd.mtype) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
end