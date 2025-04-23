function y = dump_single(featureInfo, p_benjamini, p_uncorrected, ...
    fc, median_fc, mmg, q90_mmg, spread, patterns, TAG, spread_method, fc_th)
    %%DUMP save all information to csv file
    % featureInfo holds all features you want to save as table.
    % p_benjamini/_holm/_uncorrected hold corresponding p values.
    % fc is log fold change (distance of mean(log2(area)) between tested
    % groups.
    % median fc is same except medians of groups.
    % mmg is minimum margin.
    % q90_mmg is minimum margin between 90th quantiles.
    % spread is variance or mad.
    % patterns holds string representations of patterns.
    % TAG is to identify experiment and possibly set save path.
    % keep_all is boolean, whether to save all (significant or not).
    % spread_method tells whether var or mad.
    if nargin < 12
        fc_th = 0;
    end
    
    fname_out = sprintf("%s_dump.csv",TAG);
    n_tests = length(p_benjamini);
    n_features = size(featureInfo,1);
    n_patterns = length(patterns);
    p_bonf = min(n_tests*p_uncorrected,1);
    h = p_benjamini<=0.05;
    hh = abs(fc) >= fc_th;
    sz = [n_features,n_patterns];

    [useless, inds] = sort(p_benjamini);
    ind = @(i) inds(i);

    % ASSERT1
    assert(n_tests == n_patterns*n_features, "Assert1");
        
    header = sprintf("Name;masstime;masstime_HP" + ...
        ";TEST;P_ORIG;P_FDR;P_BONF;LOG-FC;median LOG-FC;MMARGIN;90th Q MARGIN;%s\n",spread_method);

    fid = fopen(fname_out, 'w');
    fprintf(fid, header);
    for i=1:n_tests
        if(h(ind(i)) && hh(ind(i)))
            [feat_ind, pattern_ind] = ind2sub(sz, ind(i));
            fprintf(fid,'%s;',featureInfo.VarName1(feat_ind));

            %fprintf(fid,'%s;',featureInfo.Var2(feat_ind));
            %fprintf(fid,'%s;',featureInfo.Var3(feat_ind));
            %fprintf(fid,'%s;',featureInfo.Var4(feat_ind));
            fprintf(fid,'%.3f@',featureInfo.Var5(feat_ind));%mass
            fprintf(fid,'%.1f;',featureInfo.Var6(feat_ind));%time
            fprintf(fid,'%s@',featureInfo.Var5(feat_ind));%mass
            fprintf(fid,'%s;',featureInfo.Var6(feat_ind));%time
            %fprintf(fid,'%s;',featureInfo.Var7(feat_ind));
            %fprintf(fid,'%s;',featureInfo.Var8(feat_ind));
            
            fprintf(fid,'%s;%.4e;%.4e;%.4e;%.3f;%.3f;%.3f;%.3f;%.3f\n', ...
                patterns(pattern_ind), p_uncorrected(ind(i)), p_benjamini(ind(i)),p_bonf(ind(i)),fc(ind(i)), ...
                median_fc(ind(i)),mmg(ind(i)), q90_mmg(ind(i)), spread(ind(i)));
        end
    end
    fclose(fid);

end