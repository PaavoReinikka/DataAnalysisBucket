function [D, values, sparse_inds] = pre_filter(data, th, method)

    if(method=='var')
        [D,values,sparse_inds] = var_filter(data, th);
    elseif(method=='mad')
        [D,values,sparse_inds] = mad_filter(data, th);
    else
        error("Pass var or mad as method (3) argument!");
    end

end