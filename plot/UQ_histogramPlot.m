function UQ_histogramPlot(M,true_m,nburn,n,m,acc,figurepath,fwd,stringinput,sigmam,fileName)
fig = figure;
k = fwd.k;
switch fwd.mtype
    case 0
        for i = 1:k
            subplot(1,k,i)
            uq_histogram(M(nburn:end,i));
            hold on
            title(['histogram of d',num2str(i)])
            plot(true_m(i,1),0,'rx')
            txtd = ['Error_{d_' num2str(i) '} = ',num2str(abs((m(i)-true_m(i)) ./ true_m(i))*100,'%.4f%%')];
            text(0.6,0.7,txtd,'Units','normalized')
            %xlim([fwd.param.dmin,fwd.param.dmax])
            xlabel(stringinput(i));
        end
    case 1
        B = M(:,1:k-1);
        D = M(:,k:end);
        b = m(1:k-1);
        d = m(k:end);
        true_b = true_m(1:k-1);
        true_d = true_m(k:end);
        sigmab = sigmam.b;
        sigmad = sigmam.d;
        for i = 1:k-1
            subplot(2,k,i)
            uq_histogram(B(nburn:end,i));
            hold on
            title(['histogram of b',num2str(i)])
            plot(true_b(i,1),0,'rx')
            %xlim([fwd.param.bmin,fwd.param.bmax])
            xlabel(stringinput(i))
        end
        for i = 1:k
            subplot(2,k,k+i)
            uq_histogram(D(nburn:end,i));
            hold on
            title(['histogram of d',num2str(i)])
            plot(true_d(i,1),0,'rx')
            %xlim([fwd.param.dmin,fwd.param.dmax])
            xlabel(stringinput(i+k-1))
        end
        % add relative errors
        subplot(2,k,k)
        text(-0.3,1,'Relative error of b :')
        for i = 1: k-1
            txtb =['Error_{b_' num2str(i) '} = ',num2str(abs((b(i)-true_b(i)) ./ true_b(i))*100,'%.4f%%')];
            text(-0.2,1-1/(2*k)*(i),txtb)
        end
        text(-0.3,1-1/(2*k)*(k),'Relative error of d :')
        for i = 1:k
            txtd = ['Error_{d_' num2str(i) '} = ',num2str(abs((d(i)-true_d(i)) ./ true_d(i))*100,'%.4f%%')];
            text(-0.2,1-1/(2*k)*(i+k),txtd)
        end
        text(0.5,1,['Acceptance ratio = ',num2str(acc/n)]);
        text(0.5,0.8,['\sigma_b = ',num2str(sigmab)]);
        text(0.5,1-1/(2*k)*(k),['\sigma_d = ',num2str(sigmad)]);
        axis('off')
end
han = axes(fig,'visible','off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
%title(han,'Normalized histogram plot of RW samples','Position',[0.5, -0.1, 0])
sgtitle('Normalized UQhistogram plot of RW samples without burn-in period')
fig.Position = [100 100 1000 600];
fileName = ['normalized UQhistogram plot of RW samples without burn-in period np = ' num2str(fwd.np) ' sigmaRatio = ' num2str(sigmaratio) ' sigmaObRatio = ' num2str(sigma_ob_ratio) ' n = ' num2str(n) '.png'];
%fileName = ['normalized histogram plot of RW samples nRW = ' num2str(nRW) ' sigmaratio = ' num2str(sigmaratio) ' np = ' num2str(fwd.np) ' and mtype = ' num2str(fwd.mtype) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
end