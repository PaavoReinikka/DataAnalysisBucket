function y = kth_minmargin(dataX, dataY, k, relative)
    % Rough approximation of quantile margin.
    % Can also be used to calculate minmargin by setting k=0.
    % Can also take a relative measure wrt. fold-change.

    if nargin < 3
        k=0;
    end
    if nargin < 4
        relative = false;
    end
    
    sortedX = sort(dataX,2);
    sortedY = sort(dataY,2);

    if relative
        multiplier = abs(1./fchange(dataX, dataY));
    else
        multiplier = 1;
    end    
    % drop the extremum points, take the margin, and scale.
    y = multiplier .* minmargin(sortedX(:,1+k:end-k),sortedY(:,1+k:end-k));
end