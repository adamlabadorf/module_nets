function [samples] = randsamp(v,p,n)
    % take n samples with replacement from v in proportion to p
    dist = cumsum(p)/sum(p);
    for ii = 1:n
        r = rand;
        samples(ii) = v(dist > r)(1);
    end
end
