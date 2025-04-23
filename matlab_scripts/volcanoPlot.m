function y = volcanoPlot(p, fc, spread, th, tle)

    if nargin == 3
        tle = spread;
        spread = ones(size(p));
        th = 0;
    end
    
    h = spread >= th;
    plot(fc(h), -log10(p(h)), '.', 'MarkerSize',1)
    hold on

    plot(fc(~h), -log10(p(~h)), 'k.', 'MarkerSize',1)
    yline(-log10(0.05), 'r--')
    xline([-2,2],'r--')
    xlabel('log-fc (log2 ratio)')
    ylabel('-log10(p)')
    title(tle)
    hold off

end