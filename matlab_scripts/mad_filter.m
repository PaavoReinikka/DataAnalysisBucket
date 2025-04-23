function [D, mads, indices] = mad_filter(data, th)
    all_mads  = mad(data, 0, 2);
    bools     = all_mads >= th;
    rnge      = 1:size(data,1);
    indices   = transpose(rnge(bools));
    mads      = all_mads(bools);
    D         = data(bools,:);
end