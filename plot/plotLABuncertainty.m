function plotLABuncertainty(E,delta,fwd)
%{
Input: 
E      : nm x 1  or (nm-1)x 1
Cov    : nm x nm or (nm-1)x (nm-1)
Output:
Figure
%}
alpha = 0.05;% 95% confidence interval
nm = fwd.nm;
k = fwd.k;
LABmesh = setParameterization(E,fwd);
switch fwd.mtype
    case 0
        %d_uncer = E_uncertainty;
        for I= 1:k
            %delta_d = d_uncer(I)*norminv(1-alpha/2);
            % Plot the LAB depth
            xx = [LABmesh.X(I) LABmesh.X(I+1)];
            zz = LABmesh.d(I) * [1 1];
            zz_lower = (LABmesh.d(I)-delta(I)) * [1 1];
            zz_upper = (LABmesh.d(I)+delta(I)) * [1 1];
            plot(xx,zz_lower,'g','LineWidth',2);% plot the lab Depth
            hold on
            plot(xx,zz_upper,'g','LineWidth',2);% plot the lab Depth
            x2 = [xx, fliplr(xx)];
            inBetween = [zz_upper, fliplr(zz_lower)];
            fill(x2, inBetween, 'cyan','FaceAlpha',0.2);
            plot(xx,zz,'k','LineWidth',3,'LineStyle','-.');% plot the lab Depth
        end
    case 1
        %b_uncer = [E_uncertainty(1:k-1)]; % standard deviation
        %d_uncer = E_uncertainty(k:end);
        for I= 1:k
            %delta_d = d_uncer(I)*sqrt(norminv(1-alpha));
            % Plot the LAB depth
            xx = [LABmesh.X(I) LABmesh.X(I+1)];
            zz = LABmesh.d(I) * [1 1];
            zz_lower = (LABmesh.d(I)-delta(I)) * [1 1];
            zz_upper = (LABmesh.d(I)+delta(I)) * [1 1];
            plot(xx,zz_lower,'g','LineWidth',2);% plot the lab Depth
            hold on
            plot(xx,zz_upper,'g','LineWidth',2);% plot the lab Depth
            x2 = [xx, fliplr(xx)];
            inBetween = [zz_upper, fliplr(zz_lower)];
            fill(x2, inBetween, 'cyan','FaceAlpha',0.2);
            plot(xx,zz,'k','LineWidth',3,'LineStyle','-.');% plot the lab Depth
        end
        for I = 1:k-1
            %delta_b = b_uncer(I)*sqrt(norminv(1-alpha));
            xx = LABmesh.X(I+1) * [1 1];
            xx_left = (LABmesh.X(I+1)-delta(I+k-1))*[1 1];
            xx_right = (LABmesh.X(I+1)+delta(I+k-1))*[1 1];
            zz1 = [LABmesh.d(I) LABmesh.d(I+1)];
            plot(xx_left,zz1,'r','LineWidth',2);
            hold on
            plot(xx_right,zz1,'r','LineWidth',2);
            z1 = [zz1, fliplr(zz1)];
            inBetween1 = [xx_left, fliplr(xx_right)];
            fill(inBetween1,z1, 'red','FaceAlpha',0.2);
            plot(xx,zz1,'k','LineWidth',3,'LineStyle','-.');% plot the lab Depth
        end
end

