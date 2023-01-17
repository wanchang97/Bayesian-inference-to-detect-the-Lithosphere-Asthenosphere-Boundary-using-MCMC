classdef Dist
    %{
This is modified from the ERADist
https://www.cee.ed.tum.de/era/software/eradist/ 
The modification is to fit our special problem case with two types of input
parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PriorPdf definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    name:            name of the distribution 'normal' or 'uniform'
    opt :            'MOM' or 'PAR'
    val :            the parameter used to define the distribution
    M:               samples to be evaluated | dim: np x nsamples
    mu:              prior mean | dim: np x 1
    cov:             prior covmatrix | dim :  np x np

Outputs: 
    dist:            Priorpdf evaluated at the points M              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The following distribution types are available:
opt = "PAR", specification of the distribution by its parameters:
Normal:                     Obj = ERADist('normal','PAR',[mean,std])
Uniform:                    Obj = ERADist('uniform','PAR',[lower,upper])
Log-normal:                 Obj = ERADist('lognormal','PAR',[mu_lnx,sig_lnx])
opt = "MOM", specification of the distribution by its moments:
Normal:                     Obj = ERADist('normal','MOM',[mean,std])
Uniform:                    Obj = ERADist('uniform','MOM',[mean,std])
Log-normal:                 Obj = ERADist('lognormal','MOM',[mean,std])
    %}
    properties
        Name % Name of the distribution
        Par  % parameter of the distribution
    end
    methods
        function Obj = Dist(name,opt,val)
            Obj.Name = lower(name);
            switch upper(opt)

                %------------------------------------------------------------------------------------------------------------------
                % If PAR is chosen, the validity of the given parameters is
                % checked.

                case 'PAR'
                    Obj.Par = val;
                    switch lower(Obj.Name)
                        case 'normal'
                            if Obj.Par(2) > 0
                            else
                                error('The Normal distribution is not defined for your parameters');
                            end
                        case 'uniform'
                            if val(2) <= val(1)
                                error('The upper bound b must be larger than the lower bound a.');
                            end
                        case 'lognormal'
                            if Obj.Par(2) > 0
                            else
                                error('The Lognormal distribution is not defined for your parameters');
                            end
                    end
                case 'MOM'
                    if length(val) > 1 && val(2) < 0
                        error('The standard deviation must be non-negative.');
                    else
                        switch lower(Obj.Name)
                            case {'normal','gaussian'}
                                if val(2) > 0
                                    Obj.Par  = val;
                                else
                                    error('The Normal distribution is not defined for your parameters');
                                end
                            case 'uniform'
                                % compute parameters
                                Obj.Par(1) = val(1) - sqrt(12)*val(2)/2;
                                Obj.Par(2) = val(1) + sqrt(12)*val(2)/2;
                            case 'lognormal'
                                % Solve two equations for the parameters of the distribution
                                Obj.Par(1) = log(val(1)) - log(sqrt(1+(val(2)/val(1))^2));  % mean normal
                                Obj.Par(2) = sqrt(log(1+(val(2)/val(1))^2));   % sigma normal
                        end
                    end
            end

        end
        function Mean = mean(Obj) % return the mean of the distribution
            switch lower(Obj.Name)
                case {'normal','gaussian'}
                    Mean = Obj.Par(1);
                case 'uniform'
                    Mean = (Obj.Par(2)+Obj.Par(1))/2;
                case 'lognormal'
                    Mean = exp(Obj.Par(1)+(Obj.Par(2)^2/2));
            end
        end
        function Standarddeviation = std(Obj) % return the standard deviation of the distribution
            switch lower(Obj.Name)
                case {'normal','gaussian'}
                    Standarddeviation = Obj.Par(2);
                case 'uniform'
                    Standarddeviation = sqrt(((Obj.Par(2)-Obj.Par(1))^2)/12);
                case 'lognormal'
                    Standarddeviation = sqrt((exp(Obj.Par(2)^2)-1)*exp(2*Obj.Par(1)+Obj.Par(2)^2));
            end
        end
        function CDF = cdf(Obj,x) % return the CDF value at x
            switch lower(Obj.Name)
                case {'normal','gaussian'}
                    CDF = normcdf(x,Obj.Par(1),Obj.Par(2));
                case 'uniform'
                    CDF = unifcdf(x,Obj.Par(1),Obj.Par(2));
                case 'lognormal'
                    CDF = logncdf(x,Obj.Par(1),Obj.Par(2));
            end
        end
        function PDF = pdf(Obj,x) % return the PDF value at x
            switch lower(Obj.Name)
                case {'normal','gaussian'}
                    PDF = normpdf(x,Obj.Par(1),Obj.Par(2));
                case 'uniform'
                    PDF = unifpdf(x,Obj.Par(1),Obj.Par(2));
                case 'lognormal'
                    PDF = lognpdf(x,Obj.Par(1),Obj.Par(2));
            end
        end
        function logPDF = logpdf(Obj,x) % return the logPDF value at x
            switch lower(Obj.Name)
                case {'normal','gaussian'}
                    logPDF = logGauss(x,Obj.Par(1),Obj.Par(2));
                case 'uniform'
                    logPDF = log(unifpdf(x,Obj.Par(1),Obj.Par(2)));
                case 'lognormal'
                    logPDF = log(lognpdf(x,Obj.Par(1),Obj.Par(2)));
            end
        end
        function Random = random(Obj,m,n) % Generate random samples according to the distribution of the object
        if nargin == 2
            switch lower(Obj.Name)
                case {'normal','gaussian','standardnormal','standardgaussian'}
                Random = normrnd(Obj.Par(1),Obj.Par(2),m);
                case 'uniform'
                Random = random('uniform',Obj.Par(1),Obj.Par(2),m);
                case 'lognormal'
                Random = lognrnd(Obj.Par(1),Obj.Par(2),m);
            end
        end
        if nargin == 3
            switch lower(Obj.Name)
                case {'normal','gaussian','standardnormal','standardgaussian'}
                Random = normrnd(Obj.Par(1),Obj.Par(2),m,n);
                case 'uniform'
                Random = random('uniform',Obj.Par(1),Obj.Par(2),m,n);
                case 'lognormal'
                Random = lognrnd(Obj.Par(1),Obj.Par(2),m,n);
            end
        end
        
        end

    end
end

% addpath("Classes");
% switch type
%     case 'normal'
%         pi = mvnpdf(M,mu,cov);
%     case 'uniform'
%         dist_uniform = ERADist('uniform','MOM',[mu,Sigma]);
%         y = dist_uniform.pdf;
%         % ? Fixme: only for 1p case
%         %Sigma = diag(cov);
%         %a = mu - sqrt(3*Sigma);
%         %b = mu + sqrt(3*Sigma);
%         %pi = unipdf(M,a,b);  %  np x nMC
%     end
% end