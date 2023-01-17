function Pi = priorPdf(m,m_dist,fwd)
% It only works for the independent input parameter case
% m : c,d size np x nsamples = 2k x nsamples
%     c dim n x nsamples
%     d dim n x nsamples
% m_dist : c_dist ,d_dist
% Pi
k = fwd.k;
nsamples = size(m,2);
switch fwd.mtype
    case 0
        d = m;
        d_dist = m_dist;
        Pi = zeros(nsamples,1);
        for j = 1:nsamples
            pid = 1;
            for i = 1:k-1
                pid = pid * d_dist.pdf(d(i,j));
            end
            pid = pid * d_dist.pdf(d(i+1,j));
            Pi(j,1) = pid;
        end
    case 1
        b = m(1:k-1,:);
        d = m(k:end,:);
        b_dist = m_dist(1);
        d_dist = m_dist(2);
        Pi = zeros(nsamples,1);

        for j = 1:nsamples
            pib = 1; % initialize
            pid = 1;
            for i = 1:k-1
                pib = pib * b_dist.pdf(b(i,j));
                pid = pid * d_dist.pdf(d(i,j));
            end
            pid = pid * d_dist.pdf(d(i+1,j));
            pi = pib * pid;
            Pi(j,1) = pi;
        end
end
end
