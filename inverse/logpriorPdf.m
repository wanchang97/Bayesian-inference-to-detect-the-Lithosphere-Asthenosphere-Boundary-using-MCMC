function logPi = logpriorPdf(m,m_dist,fwd)
% It only works for the independent input parameter case
% m : c,d size np x nsamples = 2k x nsamples
%     c dim n x nsamples
%     d dim n x nsamples
% m_dist : c_dist ,d_dist
% logPi
k = fwd.k;
nsamples = size(m,2);
switch fwd.mtype
    case 0
        d = m;
        d_dist = m_dist;
        logPi = zeros(nsamples,1);
        for j = 1:nsamples
            logpid = 0;
            for i = 1:k-1
                logpid = logpid + d_dist.logpdf(d(i,j));
            end
            logpid = logpid + d_dist.logpdf(d(i+1,j));
            logPi(j,1) = logpid;
        end
    case 1   
        b = m(1:k-1,:);
        d = m(k:end,:);
        b_dist = m_dist(1);
        d_dist = m_dist(2);
        logPi = zeros(nsamples,1);
        for j = 1:nsamples
            logpib = 0; % initialize
            logpid = 0;
            for i = 1:k-1
                logpib = logpib + b_dist.logpdf(b(i,j));
                logpid = logpid + d_dist.logpdf(d(i,j));
            end
            logpid = logpid + d_dist.logpdf(d(i+1,j));
            logpi = logpib + logpid;
            logPi(j,1) = logpi;
        end
end
end