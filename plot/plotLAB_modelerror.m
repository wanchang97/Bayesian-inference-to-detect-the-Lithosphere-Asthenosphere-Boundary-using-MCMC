function plotLAB_modelerror(fwd,labMesh,delta,denotemodelerror)
%{
---------------------------------------------------------------------------
Plot LAB and the FEM model error
---------------------------------------------------------------------------
Input: 
    fwd        : fwd a struct containing fwd.mtype, fwd.k
    labMesh    : labMesh a struct containing labMesh.X, labMesh.d labMesh.T
    delta      : containing all the upper and lower limit of each input
    denotemodelerror: 1: plot 0 : not plot
Output:
---------------------------------------------------------------------------           
%}
X = labMesh.X;
d = labMesh.d;
k = fwd.k;
nm = fwd.nm;
%plotSites = getValueFromVarargin('plotSites',varargin,0);
switch fwd.mtype
    case 0
        dplus = delta(1,:);
        dminus = delta(2,:);
        % plot the depths with its FEM model error
        for I = 1:k
            xx = [X(I) X(I+1)];
            zz = d(I) * [1 1];
            zz_minus = (d(I)-dminus(I)) * [1 1];
            zz_plus = (d(I)+dplus(I)) * [1 1];
            plot(xx,zz_minus,'b-.','LineWidth',1);
            hold on
            plot(xx,zz_plus,'b-.','LineWidth',1);
            x2 = [xx, fliplr(xx)];
            inBetween = [zz_minus, fliplr(zz_plus)];
            fill(x2, inBetween, 'blue','FaceAlpha',0.2);
            plot(xx,zz,'b','LineWidth',1);
            if denotemodelerror
                txtd = ['\delta_{d_' num2str(I) '}: (-' num2str(dminus(I),'%.1i') ', +' num2str(dplus(I),'%.1i') ')'];
                text(fwd.geometry.lx/k*(I-1),d(I)*0.9,txtd);
            else 
                continue
            end
        end

    case 1
        dplus = delta(1,:);
        dminus = delta(2,:);
        % plot the depths with its FEM model error
        for I = 1:k
            xx = [X(I) X(I+1)];
            zz = d(I) * [1 1];
            zz_minus = (d(I)-dminus(I+k-1)) * [1 1];
            zz_plus = (d(I)+dplus(I+k-1)) * [1 1];
            plot(xx,zz_minus,'b-.','LineWidth',1);
            hold on
            plot(xx,zz_plus,'b-.','LineWidth',1);
            x2 = [xx, fliplr(xx)];
            inBetween = [zz_minus, fliplr(zz_plus)];
            fill(x2, inBetween, 'blue','FaceAlpha',0.2);
            plot(xx,zz,'b','LineWidth',1);
            if denotemodelerror
                txtd = ['\delta_{d_' num2str(I) '}: (-' num2str(dminus(I+k-1),'%.1i') ', +' num2str(dplus(I+k-1),'%.1i') ')'];
                text(X(I)*1.2,d(I)*0.9,txtd);
            else 
                continue
            end
        end

        % plot the boundary associated with model error
        for I = 1:k-1
            xx = X(I+1) * [1 1];
            xx_left = (X(I+1)-dminus(I)) * [1 1];
            xx_right = (X(I+1)+dplus(I)) * [1 1];
            zz = [d(I) d(I+1)];
            plot(xx_left,zz,'b-.','LineWidth',1);
            plot(xx_right,zz,'b-.','LineWidth',1);
            z2 = [zz, fliplr(zz)];
            inBetween = [xx_left, fliplr(xx_right)];
            fill(inBetween,z2, 'blue','FaceAlpha',0.2);
            plot(xx,zz,'b','LineWidth',1);
            if denotemodelerror
               txtb = ['\delta_{b_' num2str(I) '}: (-' num2str(dminus(I),'%.1i') ', +' num2str(dplus(I),'%.1i') ')'];
               text(X(I+1)*0.9,d(I)*0.5,txtb);
            else
                continue
            end
        end
        %text(x,z,'The true LAB considered FEM model error')
end
end

