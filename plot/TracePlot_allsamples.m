function TracePlot_allsamples(M,M_all,true_m,delta,m0,n,m,acc,figurepath,fwd,sigmam,fileName)
fig = figure;
k = fwd.k;
switch fwd.mtype
    case 0
        sigmad = sigmam.d;
        for i = 1:k
            subplot(1,k,i)
            plot(M(:,i));
            hold on
            plot(M_all(:,i));
            plot(0,m0(i),'ro');
            txtd = ['Error_{d_' num2str(i) '} = ',num2str(abs((m(i)-true_m(i)) ./ true_m(i))*100,'%.4f%%')];
            text(0.3,0.7,txtd,'Units','normalized')
            xlabel('steps')
            ylabel(['d',num2str(i)]);
            yline(true_m(i),'r-');
            yline(true_m(i)+delta(1,i), 'r-.')
            yline(true_m(i)-delta(2,i), 'r-.')
        end
    case 1
        B = M(:,1:k-1);
        D = M(:,k:end);
        b = m(1:k-1);
        d = m(k:end);
        b0 = m0(1:k-1);
        d0 = m0(k:end);
        true_b = true_m(1:k-1);
        true_d = true_m(k:end);
        sigmab = sigmam.b;
        sigmad = sigmam.d;
        for i = 1:k-1
            subplot(2,k,i)
            plot(B(:,i));
            hold on
            plot(M_all(:,i));
            plot(0,b0(i),'ro');
            xlabel('steps')
            ylabel(['b',num2str(i)]);
            yline(true_b(i),'r-');
            yline(true_b(i)+delta(1,i),'r-.')
            yline(true_b(i)-delta(2,i),'r-.')
        end

        for i = 1:k
            subplot(2,k,i+k)
            plot(D(:,i));
            hold on
            plot(M_all(:,i+k));
            plot(0,d0(i),'ro');
            xlabel('steps')
            ylabel(['d',num2str(i)]);
            yline(true_d(i),'r-');
            yline(true_d(i)+delta(1,i+k-1),'r-.')
            yline(true_d(i)-delta(2,i+k-1), 'r-.')
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
% Common xlabel
han = axes(fig,'visible','off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
%xlabel(han,'number of MCMC states')
fig.Position = [100 100 1000 600];
sgtitle('Markov Chain states plot')
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
end