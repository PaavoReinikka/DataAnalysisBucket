function [D, variances, indices] = var_filter(data, th)
    all_variances = var(data, 0, 2);
    bools     = all_variances >= th;
    rnge      = 1:size(data,1);
    indices   = transpose(rnge(bools));
    variances = all_variances(bools);
    D         = data(bools,:);
end