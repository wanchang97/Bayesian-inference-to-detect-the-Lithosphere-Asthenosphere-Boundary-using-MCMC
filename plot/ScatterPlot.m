function ScatterPlot(M,nburn,t,m,figurepath,fwd,stringinput,plotmin,plotmax,fileName)
% ScatterPlot(M,true_m,nburn,n,m,acc,figurepath,fwd,stringinput,plotmin,plotmax,sigmam,sigma_ob_ratio,fileName)
figscatter = figure;
k = fwd.k;
for i = 1:fwd.nm
    for j = 1:i
        num = fwd.nm*(i-1)+j;
        subplot(fwd.nm,fwd.nm,num)
        if j == i
            histogram(M(nburn:t:end,j),'Normalization','probability');
            hold on
            xline(m(j),'r','Fontsize',17);
            xlim([plotmin(i),plotmax(i)])
            xlabel(stringinput(j))
        else
            scatter(M(nburn:t:end,i),M(nburn:t:end,j));
            hold on
            plot(m(i),m(j),'r+')
            xlim([plotmin(i),plotmax(i)])
            ylim([plotmin(j),plotmax(j)])
            xlabel(stringinput(j))
            ylabel(stringinput(i))
        end
    end
%     % add relative errors
%     
%     switch fwd.mtype
%         case 0
%             sigmad = sigmam.d;
%             text(-0.3,1,'Relative error of d :');
%             for i = 1:k
%                 num = rem(i,2)+1;
%                 subplot(fwd.nm,fwd.nm,num);
%                 txtd = ['Error_{d_' num2str(i) '} = ',num2str(abs((m(i)-true_m(i)) ./ true_m(i))*100,'%.4f%%')];
%                 text(-0.2,1-1/(2*k)*(i),txtd)
%             end
%             subplot(fwd.nm,fwd.nm,fwd.nm);
%             text(0.5,1,['Acceptance ratio = ',num2str(acc/n)]);
%             text(0.5,0.8,['\sigma_d = ',num2str(sigmad)]);
%             text(0.5,0,['\sigma_{ob} = ',num2str(sigma_ob_ratio) '*true_{ob}']);
%         case 1
%             B = M(:,1:k-1);
%             D = M(:,k:end);
%             b = m(1:k-1);
%             d = m(k:end);
%             true_b = true_m(1:k-1);
%             true_d = true_m(k:end);
%             sigmab = sigmam.b;
%             sigmad = sigmam.d;
%             subplot(fwd.nm,fwd.nm,2)
%             text(-0.3,1,'Relative error of b :');
%             for i = 1: k-1
%                 num = i;
%                 subplot(fwd.nm,fwd.nm,num);
%                 txtb =['Error_{b_' num2str(i) '} = ',num2str(abs((b(i)-true_b(i)) ./ true_b(i))*100,'%.4f%%')];
%                 text(-0.2,1-1/(2*k)*(i),txtb);
%             end
%             text(-0.3,1-1/(2*k)*(k),'Relative error of d :');
%             for i = 1:k
%                 num = rem(i,2)+1;
%                 subplot(fwd.nm,fwd.nm,num);
%                 txtd = ['Error_{d_' num2str(i) '} = ',num2str(abs((d(i)-true_d(i)) ./ true_d(i))*100,'%.4f%%')];
%                 text(-0.2,1-1/(2*k)*(i+k),txtd);
%             end
%             subplot(fwd.nm,fwd.nm,fwd.nm);
%             text(0.5,1,['Acceptance ratio = ',num2str(acc/n)]);
%             text(0.5,0.8,['\sigma_b = ',num2str(sigmab)]);
%             text(0.5,1-1/(2*k)*(k),['\sigma_d = ',num2str(sigmad)]);
%             text(0.5,0,['\sigma_{ob} = ',num2str(sigma_ob_ratio),'*true_{ob}']);
%     end
    %axis('off')
end
sgtitle('Scatter plot of steady MCMC samples (without burn-in periods)')
figscatter.Position = [100 100 1000 800];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);