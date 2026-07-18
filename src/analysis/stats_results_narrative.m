% ========================================================================
% stats_results_narrative.m
%
% PURPOSE
%   Prints the summary numbers used in the Results narrative (the Fig. 4 box
%   plots): the CAM16-UCS observer-metamerism magnitude (dE') as a function of
%   illuminant class, primary bandwidth (SD = sigma of the Gaussian display
%   primaries, in nm), and degree of chromatic adaptation D. Also prints the
%   per-observer percentage reduction relative to the unadapted state (D = 0)
%   and the daylight D = 0 vs D = 1 Wilcoxon signed-rank tests.
%
% INPUT (.mat)
%   results/compare_nl_adapt_ucs.mat
%     .Res      struct with one field per illuminant class (e.g. daylight,...);
%               each class R holds:
%                 R.sd     row of primary bandwidths sigma (nm)
%                 R.Dlist  row of adaptation degrees D
%                 R.pD     1 x nD cell; pD{iD} is (nSD x nObsIll) matrix of the
%                          CAM16-UCS dE' metamerism magnitude at degree D(iD),
%                          one column per pooled (observer x illuminant) sample
%     .classes  cell/char list of illuminant class names
%
% PRINTS (nothing saved to disk)
%   - median dE' pooled over observers x illuminants, per SD row and D column
%   - per-observer % reduction of dE' vs the D = 0 baseline, per SD
%   - daylight D = 0 vs D = 1 paired signrank tests (per SD and pooled), with
%     Cohen's dz effect size and the signed-rank z statistic
%
% METHOD
%   dE' is the between-observer spread of the per-observer metameric display
%   matches, read out in CAM16-UCS. Reduction is 100*(1 - dE'(D)/dE'(D=0))
%   computed per sample (so a positive value means adaptation shrank
%   metamerism). Paired significance uses MATLAB's signrank (Wilcoxon).
% ========================================================================

% Numbers + stats for the rewritten Results narrative (Fig 4 box plots).
clearvars; clc;
projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
% Load the precomputed metamerism magnitudes (already reduced to dE' values).
S = load(fullfile(projectRoot,'results','compare_nl_adapt_ucs.mat'));
% Res: results grouped by illuminant class; classes: list of class names.
Res = S.Res; classes = cellstr(S.classes);
% Loop over illuminant classes and print one summary block each.
for ci=1:numel(classes)
    % Skip any class name that has no field in Res.
    cls=classes{ci}; if ~isfield(Res,cls), continue; end
    % Pull this class's bandwidth list (sd), adaptation degrees (D), and count.
    R=Res.(cls); sd=R.sd(:)'; D=R.Dlist(:)'; nD=numel(D);
    fprintf('\n================= %s =================\n', cls);
    % --- Table of median dE' : rows = bandwidth SD, columns = adaptation D ---
    fprintf('median dE'' (pooled over observers x illuminants):\n');
    fprintf('  SD |');
    % Column header, one per adaptation degree D.
    for iD=1:nD, fprintf('  D=%-4g |', D(iD)); end
    fprintf('\n');
    % One row per primary bandwidth sigma.
    for si=1:numel(sd)
        fprintf(' %3g |', sd(si));
        % Median dE' pooled over all (observer x illuminant) samples at each D.
        for iD=1:nD, fprintf(' %6.3f |', median(R.pD{iD}(si,:),'omitnan')); end
        fprintf('\n');
    end
    % per-observer reduction medians vs D=0
    % Locate the unadapted baseline column (D = 0).
    iD0=find(D==0,1);
    fprintf('per-unit reduction vs D=0 [%%] (median):\n');
    % Report the % shrinkage of dE' at each adapted D relative to D = 0.
    for iD=1:nD
        % Nothing to report for the baseline against itself.
        if D(iD)==0, continue; end
        fprintf('  D=%g: ', D(iD));
        vals=[];
        for si=1:numel(sd)
            % Per-sample reduction: 100*(1 - dE'(D)/dE'(D=0)); positive = smaller.
            red=100*(1 - R.pD{iD}(si,:)./R.pD{iD0}(si,:));
            % Store the per-bandwidth median so we can print the SD-range after.
            vals(end+1)=median(red,'omitnan');
            fprintf('SD%g=%.1f  ', sd(si), median(red,'omitnan'));
        end
        % Range of the median reduction across the three bandwidths.
        fprintf(' [range %.1f-%.1f]\n', min(vals), max(vals));
    end
end

% ---- daylight D=0 vs D=1 paired test (per SD) ----
% Wilcoxon signed-rank tests that adaptation (D=1) lowers dE' below D=0, on the
% daylight class only. a = dE' unadapted, b = dE' fully adapted (paired samples).
fprintf('\n=========== daylight D=0 vs D=1 paired tests ===========\n');
R=Res.daylight; sd=R.sd(:)'; D=R.Dlist(:)';
% Column indices of the unadapted (D=0) and fully adapted (D=1) states.
i0=find(D==0,1); i1=find(D==1,1);
% One paired test per primary bandwidth sigma.
for si=1:numel(sd)
    % Paired D=0 and D=1 dE' vectors for this bandwidth.
    a=R.pD{i0}(si,:); b=R.pD{i1}(si,:);
    % Keep only pairs finite in both states.
    ok=isfinite(a)&isfinite(b); a=a(ok); b=b(ok);
    % Wilcoxon signed-rank p-value and stats (st.zval = normal-approx z).
    [p,~,st]=signrank(a,b);
    dz=mean(a-b)/std(a-b);              % Cohen's dz
    rbc=st.zval/sqrt(numel(a));         % rank-biserial-ish (z/sqrt n)
    fprintf('SD=%2g: n=%d, median D0=%.2f D1=%.2f, signrank p=%.2e, z=%.1f, dz=%.2f\n', ...
        sd(si), numel(a), median(a), median(b), p, st.zval, dz);
end
% pooled across SD (stack)
% Stack all three bandwidths into one paired sample and test once more.
a=[]; b=[];
for si=1:numel(sd), a=[a R.pD{i0}(si,:)]; b=[b R.pD{i1}(si,:)]; end
% Drop non-finite pairs from the pooled sample.
ok=isfinite(a)&isfinite(b); a=a(ok); b=b(ok);
% Pooled signed-rank test and Cohen's dz effect size.
[p,~,st]=signrank(a,b); dz=mean(a-b)/std(a-b);
fprintf('POOLED SD: n=%d, signrank p=%.2e, z=%.1f, dz=%.2f\n', numel(a), p, st.zval, dz);
