

%% INITIALISE THE ENVIRONMENT, LOAD AND COMPILE DATA, FILTER, ETC:
clear
addpath("univariate_utils\")

% BEGIN BY DEFINING THE EXPERIMENT:
% SINGLE_MODE = true refers to reading data from original csv files (mode per
% experiment). NOTE: now values start at field 'I4'
% SINGLE_MODE = false refers to reading from a compiled 
% csv file, e.g., all modes/medium/cells, with one added column indicating
% the sourcefile. NOTE: now values start at field 'J2'


%SINGLE_MODE = false;
%fname='../celldata/clean/AllModes.csv';experiment="ALL_VF0";FILTER = true;filter_method='var';filter_th = 0;fc_th=1.5;
%fname='../celldata/clean/cellsAllModes.csv';experiment="cells_VF0";FILTER = true;filter_method='var';filter_th =0;fc_th=1.5;
%fname='../celldata/clean/mediumAllModes.csv';experiment="medium_VF0";FILTER = true;filter_method='var';filter_th =0;fc_th=1.5;

SINGLE_MODE = true;
%fname='../celldata/clean/cellslipidomicspos_log2.csv';experiment="lipidpos_VF025FC15";FILTER= true;filter_method='var';filter_th =0.25;fc_th=1.5;
%fname='../celldata/clean/cellslipidomicsneg_log2.csv';experiment="lipidneg_VF025FC15";FILTER = true;filter_method='var';filter_th =0.25;fc_th=1.5;

%fname='../celldata/clean/cellshilicneg_log2.csv';experiment="cellshilicneg_VF025FC15";FILTER = true;filter_method='var';filter_th =0.25;fc_th=1.5;
%fname='../celldata/clean/cellshilicpos_log2.csv';experiment="cellshilicpos_VF025FC15";FILTER = true;filter_method='var';filter_th =0.25;fc_th=1.5;
%fname='../celldata/clean/cellsRPneg_log2.csv';experiment="cellsRPneg_VF025FC15";FILTER = true;filter_method='var';filter_th =0.25;fc_th=1.5;
%fname='../celldata/clean/cellsRPpos_log2.csv';experiment="cellsRPpos_VF025FC15";FILTER = true;filter_method='var';filter_th =0.25;fc_th=1.5;

%fname='../celldata/clean/mediumHILICneg_log2.csv';experiment="mediumhilicneg_VF025FC15";FILTER = true;filter_method='var';filter_th =0.25;fc_th=1.5;
%fname='../celldata/clean/mediumHILICpos_log2.csv';experiment="mediumhilicpos_VF025FC15";FILTER = true;filter_method='var';filter_th =0.25;fc_th=1.5;
%fname='../celldata/clean/mediumRPneg_log2.csv';experiment="mediumRPneg_VF025FC15";FILTER = true;filter_method='var';filter_th =0.25;fc_th=1.5;
fname='../celldata/clean/mediumRPpos_log2.csv';experiment="mediumRPpos_VF025FC15";FILTER = true;filter_method='var';filter_th =0.25;fc_th=1.5;

%minmargin threshold
mmg_th=0;

PATH_DUMPS = "../cellresults/ttest/dumps/";
PATH_SUMMARY = "../cellresults/ttest/";
PATH_PLOT = "../cellresults/ttest/volcanoplots/";
PATH_NEW_DATA = "../celldata/clean/ttested/";
PATH_NEW_IDS =  "../cellresults/ttest/significant_IDs/";

% SAVE THE RESULTS AND COMPILE NEW DATA WITH SIGNF FEATURES ID's
SAVE_SUMMARIES = false;
SAVE_DUMPS = true;
SAVE_VOLCANO = false;
PLOT_VOLCANO = false;
GENERATE_NEW_IDS = false;

if(~SINGLE_MODE)
    D        = readmatrix(fname, 'Range','J2');
    features = readFeatureInfo(fname,9,2); % feature info (other than area values)
    fullData = importfile(fname);
else
    D        = readmatrix(fname, 'Range','I4');
    features = readFeatureInfo(fname,8,4); % feature info (other than area values)
end

n_tests_bfr = size(D,1);

% Apply filter:
% replace the initial datamatrix with new filtered one, and update features
if FILTER
    % spread refers to either variance or mad, depending on method
    [D, spread, sparse_inds] = pre_filter(D, filter_th, filter_method);
    features = features(sparse_inds,:);
else
    filter_th = nan;
end

% Compile groups
aSYN = D(:,1:10);
comb = D(:,11:20);
INFg = D(:,21:30);
UT   = D(:,31:40);

% Number of patterns.
patterns = ["aSYN--UT", "comb.--UT", "INFg--UT", "aSYN--INFg","comb.--INFg.", "aSYN--comb."];
n_patterns = length(patterns);
n_tests_bfr = n_tests_bfr * n_patterns;

% expands the spread/variance to match the testsize
fullSpread = repmat(spread, n_patterns,1);

% Compile patterns
dataX = cat(1, aSYN, comb, INFg, aSYN, comb, aSYN);
dataY = cat(1, UT,   UT,   UT,   INFg, INFg, comb);
dataXY = cat(2, dataX, dataY);

n_tests = size(dataX,1);

%TODO: add some asserts here and in results
assert(n_tests == size(dataY,1));

% Calculate fold-change and margin
fc = fchange(dataX, dataY); % fold-change
mfc = median_fchange(dataX,dataY); % median fold change
mmg = minmargin(dataX, dataY); % minmargin
q_mmg = kth_minmargin(dataX,dataY,1); % quantile margin


% RUN TESTS IN ONE BATCH

[h_,p_uncorrected,ci_] = ttest2(dataX,dataY,'Vartype','unequal','Dim',2);

% CORRECTIONS

% Holm-Bonferroni
[p_holm, c_alpha, h_holm] = fwer_holmbonf(p_uncorrected, 0.05);
p_holm = transpose(p_holm);
h_holm = transpose(h_holm);

% B-Y correction
[h_benjamini, crit_p, adj_ci_cvrg, p_benjamini]=fdr_bh(p_uncorrected,0.05,'dep');
% Cap p-values to 1
p_benjamini = min(p_benjamini,1);

% RESULTS

h = h_benjamini;
% get linear indices of result
linear_inds = bool2ind(h);

% map linear indices to results -- [row=features, col=patterns]
sz = [size(D,1),n_patterns];
[feat_inds, pattern_inds] = ind2sub(sz, linear_inds);

% FC and margin for significant features
fc_BYsignificant       = fc(linear_inds);
mfc_BYsignificant      = mfc(linear_inds);
mmg_BYsignificant      = mmg(linear_inds);
q_mmg_BYsignificant    = q_mmg(linear_inds);



% number of significant (uncorrected/corrected)
count_uncorrected = sum(p_uncorrected <= 0.05);
count_bonf = sum(length(p_uncorrected)*p_uncorrected <= 0.05); %signif after Bonferroni
count_holm = sum(h_holm); %significant after Holm-Bonferroni correction
count_benjamini = sum(h);
count_unique_benjamini = length(unique(feat_inds)); %same features across all comparisons.

% Check
assert(count_benjamini == sum(p_benjamini <= 0.05), 'Problem with alpha!');
assert(n_tests == length(h));

% strong effects
numFC = sum(abs(fc)>=fc_th); %how often |FC|>= th
numMFC = sum(abs(mfc)>=fc_th); %how often |meadianFC|>= th
numMMG = sum(abs(mmg)>mmg_th); %how often |msize|> th
numQMMG = sum(abs(q_mmg)>mmg_th);


% strong effects among significant
numFC_BYsignificant = sum(abs(fc_BYsignificant)>=fc_th);
numMFC_BYsignificant = sum(abs(mfc_BYsignificant)>=fc_th);
numMMG_BYsignificant = sum(abs(mmg_BYsignificant)>mmg_th);
num_q_mmg_BYsignificant = sum(abs(q_mmg_BYsignificant)>mmg_th);
numFC_MMG_BYsignificant = sum(abs(fc_BYsignificant)>=fc_th & abs(mmg_BYsignificant)>mmg_th);
numMFC_Q_MMG_BYsignificant = sum(abs(mfc_BYsignificant)>=fc_th & abs(q_mmg_BYsignificant)>mmg_th);

% collect significant features (B-Y)
%features_BYsignificant = features(feat_inds); %-- [row=feat_inds, col=pattern_inds]

%% WRITE DETAILED RESULTS TO FILE

if(SAVE_DUMPS)
    
    if(SINGLE_MODE)
        % DUMP IF SINGLE MODES
        tag2 = sprintf("%s%s_SGNF", PATH_DUMPS, experiment);
        dump_single(features, p_benjamini, p_uncorrected, fc, mfc, mmg, q_mmg, fullSpread, patterns, tag2, filter_method, fc_th);
    else
        % DUMP IF ALLMODES
        tag2 = sprintf("%s%s_SGNF", PATH_DUMPS, experiment);
        dump(features, p_benjamini, p_uncorrected, fc, mfc, mmg, q_mmg, fullSpread, patterns, tag2, filter_method, 0);
    end
end


%% PRINT 
fprintf("Number of tests: %d\n", n_tests);
fprintf("Filtered (%s) proportion: %d\n", filter_method, 1 - n_tests / n_tests_bfr);
fprintf("Number of discoveries: %d\n", count_benjamini);

fprintf('Bonferroni corrected p <= 0.05: %d\n',count_bonf);
fprintf('B-Y corrected p <= 0.05: %d\n',count_benjamini);
fprintf('Unique and significant: %d\n',count_unique_benjamini);
    
fprintf('|FC| >= %.3f: %d\n', fc_th, numFC);
fprintf('|medianFC| >= %.3f: %d\n', fc_th, numMFC);
fprintf('|minmargin| > %.3f: %d\n', mmg_th,numMMG);
fprintf('|q90_margin| > %.3f: %d\n', mmg_th,numQMMG);
fprintf('Out of B-Y significant, |FC| >= %.3f: %d\n', fc_th, numFC_BYsignificant);
fprintf('Out of B-Y significant, |medianFC| >= %.3f: %d\n', fc_th, numMFC_BYsignificant);
fprintf('Out of B-Y significant, |minmargin| > %.3f: %d\n', mmg_th, numMMG_BYsignificant);
fprintf('Out of B-Y significant, |q90_margin| > %.3f: %d\n', mmg_th, num_q_mmg_BYsignificant);
fprintf('Out of B-Y significant, |FC| >= %.3f and |minmargin| > %.3f: %d\n',...
        fc_th, mmg_th, numFC_MMG_BYsignificant);
fprintf('Out of B-Y significant, |medianFC| >= %.3f and |q90_margin| > %.3f: %d\n',...
        fc_th, mmg_th,numMFC_Q_MMG_BYsignificant);

for i=1:n_patterns
    count = sum(pattern_inds==i);
    fprintf("Pattern %s: %d\n",patterns(i),count)
end


%% PLOT FOR SPREAD (VAR OR MAD) AGAINST LOG-P VALUES
if(FILTER & filter_th==0)
    plot(fullSpread, -log(p_benjamini),'r.');
    hold on
    xlim([0,inf])
    yline(-log(0.05));
    xline(filter_th);
    xlabel(sprintf('%s',filter_method));
    ylabel('log p-value');
    title(sprintf('Filtering th: %d.',filter_th));
    hold off
end

%% SAVE BASIC RESULTS

if(SAVE_SUMMARIES)
    fname_summary = sprintf("%s%s_summary.txt", PATH_SUMMARY,experiment);
    fid_summary = fopen(fname_summary,"w");
    fprintf(fid_summary,'Experiment: %s\n',experiment);
    fprintf(fid_summary,'Filter (bool): %d\nFilter_th: %f\n',FILTER, filter_th);
    fprintf(fid_summary,'Filter method: %s\n',filter_method);
    fprintf(fid_summary,'Number of tests: %d\n',n_tests);
    fprintf(fid_summary,'Number of patterns: %d\n',n_patterns);
    
    fprintf(fid_summary,'Patterns: (');
    for i=1:n_patterns
        fprintf(fid_summary,'%s ',patterns(i));
    end
    fprintf(fid_summary,')\n');

    fprintf(fid_summary,'Uncorrected p <= 0.05: %d\n',count_uncorrected);
    fprintf(fid_summary,'Bonferroni corrected p <= 0.05: %d\n',count_bonf);
    fprintf(fid_summary,'Holm-Bonferroni corrected p <= 0.05: %d\n',count_holm);
    fprintf(fid_summary,'B-Y corrected p <= 0.05: %d\n',count_benjamini);
    fprintf(fid_summary,'Unique and significant: %d\n',count_unique_benjamini);

    fprintf(fid_summary,'|FC| >= %.3f: %d\n', fc_th, numFC);
    fprintf(fid_summary,'|medianFC| >= %.3f: %d\n', fc_th, numMFC);
    fprintf(fid_summary,'|minmargin| > %.3f: %d\n', mmg_th,numMMG);
    fprintf(fid_summary,'|q90_margin| > %.3f: %d\n', mmg_th,numQMMG);
    fprintf(fid_summary,'Out of B-Y significant, |FC| >= %.3f: %d\n', fc_th, numFC_BYsignificant);
    fprintf(fid_summary,'Out of B-Y significant, |medianFC| >= %.3f: %d\n', fc_th, numMFC_BYsignificant);
    fprintf(fid_summary,'Out of B-Y significant, |minmargin| > %.3f: %d\n', mmg_th, numMMG_BYsignificant);
    fprintf(fid_summary,'Out of B-Y significant, |q90_margin| > %.3f: %d\n', mmg_th, num_q_mmg_BYsignificant);
    fprintf(fid_summary,'Out of B-Y significant, |FC| >= %.3f and |minmargin| > %.3f: %d\n',...
            fc_th, mmg_th, numFC_MMG_BYsignificant);
    fprintf(fid_summary,'Out of B-Y significant, |medianFC| >= %.3f and |q90_margin| > %.3f: %d\n',...
            fc_th, mmg_th,numMFC_Q_MMG_BYsignificant);
    
    fprintf(fid_summary,"------------------------\n");
    for i=1:n_patterns
        count = sum(pattern_inds==i);
        fprintf(fid_summary,"Pattern %s: %d\n",patterns(i),count);
    end

    fclose(fid_summary);
end



%%

if(GENERATE_NEW_IDS)
    % select only significant (h) or add other conditions
    hh = h & (abs(fc)>fc_th);
    if(SINGLE_MODE)
        linear_inds = bool2ind(hh);
        % map linear indices to results -- [row=features, col=patterns]
        sz = [size(D,1),n_patterns];
        [feat_inds2, pattern_inds] = ind2sub(sz, linear_inds);
    
        tag4 = sprintf("%s%s.csv",PATH_NEW_IDS, experiment);
        T = table(unique(feat_inds2));
        writetable(T, tag4, 'Delimiter',';','WriteVariableNames',0);
    else
        fullData=repmat(fullData(sparse_inds,:), n_patterns,1);
        tag4 = sprintf("%s%s.",PATH_NEW_IDS, experiment);
        write_SGNF_ids(fullData,hh, tag4);
    end
end


%% PLOT RESULTS
if(SAVE_VOLCANO || PLOT_VOLCANO)
    figure()
    tiledlayout(2,2)
    nexttile
    plot(fc, -log10(p_uncorrected), '.', 'MarkerSize',1)
    yline(-log10(0.05), 'r--')
    xline([-2,2],'r--')
    xlabel('log-fc (log2 ratio)')
    ylabel('-log10(p)')
    title('Uncorrected')
    
    nexttile
    p_bonf = min(p_uncorrected.*n_tests,1);
    plot(fc, -log10(p_bonf), '.', 'MarkerSize',1)
    ylim([0,inf])
    yline(-log10(0.05), 'r--')
    xline([-2,2],'r--')
    xlabel('log-fc (log2 ratio)')
    ylabel('-log10(p)')
    title('Bonferroni')
    
    
    nexttile
    plot(fc, -log10(p_holm), '.', 'MarkerSize',1)
    yline(-log10(0.05), 'r--')
    xline([-2,2],'r--')
    xlabel('log-fc (log2 ratio)')
    ylabel('-log10(p)')
    title('Holm-Bonferroni')
    
    if(filter_th>0)
        nexttile
        plot(fc, -log10(p_benjamini), '.', 'MarkerSize',1)
        ylim([0,inf])
        yline(-log10(0.05), 'r--')
        xline([-2,2],'r--')
        xlabel('log-fc (log2 ratio)')
        ylabel('-log10(p)')
        title('B-Y')   
    else
        nexttile
        th = 1;
        volcanoPlot(p_benjamini, fc, fullSpread, th, sprintf("B-Y   (%s < %.2f)", filter_method, th));
    end
end

if(SAVE_VOLCANO)
    fname_volcano = sprintf("%s%s%s", PATH_PLOT, experiment, "_TILES.png");
    saveas(gcf,fname_volcano);
    
    figure()
    volcanoPlot(p_benjamini, fc,'B-Y')
    fname_volcano = sprintf("%s%s%s", PATH_PLOT, experiment, "_MAIN.png");
    saveas(gcf,fname_volcano);
end
