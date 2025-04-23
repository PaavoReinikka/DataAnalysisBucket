function [p_BY, h] = correcting_with_correlations(dataX, dataY, alpha)
    if nargin < 3
        alpha = 0.05;
    end
    n_tests = size(dataX,1);
    c = harmonic(n_tests);

    [h,p,ci] = ttest2(dataX,dataY,'Vartype','unequal','Dim',2);
    [p_BY, inds] = sort(p,1);
    p_BY = p_BY(inds);

    correction =  1:length(p) * alpha  / (n_tests*c);
    h = (p_BY <= correction);

end
    