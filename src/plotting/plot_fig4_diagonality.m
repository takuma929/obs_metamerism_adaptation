% ========================================================================
% MANUSCRIPT FIGURE 6 (Optics Express, "Adaptation reduces inter-observer
%   metamerism with narrow-band primaries").  IMPORTANT: the script name says
%   "fig4" but the figure was RENUMBERED to Figure 6 in the manuscript; the
%   script number does NOT track the manuscript number.
%
% INPUT  .mat : results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat
%               (only D_set and params are loaded - see note near the load).
% OUTPUTS     : figs/fig4/fig4_diagonality.pdf (vector) and .png (600 dpi, for
%               the Word manuscript), plus a between-observer spread table to
%               the Command Window.
% ------------------------------------------------------------------------
% ========================================================================
% plot_fig4_diagonality.m   —   Figure 4
%
% The cross-observer 3x3 map D and its spread ACROSS observers.
%
% For a display with primaries P, an observer's cone excitations are
% C_obs = T_obs*P and the standard observer's are C_std = T_std*P, so
%     D = C_std * inv(C_obs)
% maps an individual observer's LMS to the standard observer's LMS for the same
% display stimulus.  Both C_std and C_obs depend on the primaries (hence on
% sigma), but for the standard observer itself C_obs = C_std, so
%     D_std = C_std * inv(C_std) = I
% EXACTLY, for any primary set — the sigma dependence cancels.  The identity is
% therefore a sigma-invariant reference (diagonal = 1, off-diagonal = 0), and a
% single reference line serves both sigma = 5 and sigma = 20: sigma acts only on
% where the 10,000 individual observers fall relative to it.
%
% Averaging D over observers hides the very thing we care about — the
% between-observer variability cancels out.  This figure therefore shows the FULL
% distribution: one histogram per matrix cell (3x3 panels), with the narrow
% (sigma = 5 nm) and broad (sigma = 20 nm) primaries overlaid, and a vertical
% line marking the standard observer (D = I) so the spread away from it is
% directly visible.
%
% D is illuminant-independent (it depends only on the display primaries and the
% observer cone fundamentals), so only the small D_set/params fields are loaded —
% not the large per-observer chromaticity arrays.
%
% Run on its own:  >> plot_fig4_diagonality
% Output: figs/fig4/fig4_diagonality.pdf   (vector, two-column width)
% ========================================================================
clearvars; close all; clc;

% Resolve project root from this script's own path (.../src/plotting/<file>):
% three fileparts() strip the file, the plotting/ and the src/ folders.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist(fullfile(projectRoot,'src'),'dir')
    projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
end
addpath(genpath(fullfile(projectRoot,'src')));
figDir = fullfile(projectRoot,'figs','fig4'); if ~exist(figDir,'dir'); mkdir(figDir); end
T = load(fullfile(projectRoot,'data','fig_parameters.mat')); figp = T.figp;
inMat = fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
assert(exist(inMat,'file')==2,'Missing results: %s',inMat);

% ---- config ------------------------------------------------------------
CFG.SDpair  = [5 20];                       % narrow / broad primaries to overlay
CFG.nBins   = 55;                           % bins per panel (shared by both SDs)
CFG.clipQ   = [0.002 0.998];                % percentile clip for the x-range
sdCol = [0.298 0.447 0.690; 0.769 0.306 0.322];  % sigma = 5 (steel blue) / sigma = 20 (muted red)
sdAlp = 0.55;
coneNm = {'L','M','S'};

% ---- load small fields only -------------------------------------------
% D_set{si} holds the precomputed 3x3xnObs cross-observer maps D for the
% si-th bandwidth; params.sdList lists the available sigmas. (D is illuminant-
% independent, so the large per-illuminant chromaticity arrays are skipped.)
fprintf('Loading D_set...\n');
S = load(inMat, 'D_set','params');
sdList = S.params.sdList(:)';
% Map the requested narrow/broad sigmas to their column indices in sdList.
iSD = arrayfun(@(s) find(sdList==s,1), CFG.SDpair, 'UniformOutput', false);
assert(all(~cellfun(@isempty,iSD)), 'sdList (%s) lacks requested SDs %s', ...
       mat2str(sdList), mat2str(CFG.SDpair));
iSD = cell2mat(iSD);
% Dm{1} = narrow (sigma=5), Dm{2} = broad (sigma=20); each 3 x 3 x nObs.
Dm  = {S.D_set{iSD(1)}, S.D_set{iSD(2)}};   % each 3 x 3 x nObs
nObs = size(Dm{1},3);
fprintf('N = %d observers;  sigma = %g and %g nm.\n', nObs, CFG.SDpair);

% ---- pass 1: x ranges --------------------------------------------------
% The 3 diagonal cells share one range, so their spreads are directly
% comparable.  The 6 off-diagonal cells keep individual ranges: they span two
% orders of magnitude (D_LM ~ +/-0.5 vs D_SL ~ +/-0.005), so a shared range
% would collapse the S row onto the D = 0 line and hide its structure.
LO = nan(3); HI = nan(3);
for i = 1:3
    for j = 1:3
        % Pool both sigmas for this matrix cell, then clip to the 0.2/99.8
        % percentiles so a few outliers do not stretch the axis.
        both = [squeeze(Dm{1}(i,j,:)); squeeze(Dm{2}(i,j,:))];
        both = both(isfinite(both));
        LO(i,j) = quantile(both,CFG.clipQ(1));
        HI(i,j) = quantile(both,CFG.clipQ(2));
    end
end
% Diagonal cells share one range (centred on the reference D = 1) so their
% spreads are directly comparable; pad so D = 1 stays inside the axis.
dg = logical(eye(3));
[dLo,dHi] = pad_range(min(LO(dg)), max(HI(dg)), 1);   % diagonal: shared, reference D = 1
fprintf('shared diagonal x range: [%.4g %.4g]\n',dLo,dHi);

% ---- figure: 3x3 histograms -------------------------------------------
fig = figure('Color','w','ToolBar','none','MenuBar','none');
tl  = tiledlayout(fig,3,3,'Padding','compact','TileSpacing','compact');
hLeg = gobjects(1,3);

% D maps individual-observer LMS (input, col j) to standard-observer LMS
% (output, row i); the 3x3 grid of panels mirrors that matrix layout.
for i = 1:3           % row = standard-observer cone (output)
    for j = 1:3       % col = individual-observer cone (input)
        ax = nexttile(tl); hold(ax,'on');
        % The nObs values of entry D(i,j) for the narrow and broad primaries.
        v5  = squeeze(Dm{1}(i,j,:));
        v20 = squeeze(Dm{2}(i,j,:));
        % Standard observer has D = I, so this cell's reference is 1 on the
        % diagonal and 0 off it.
        stdVal = double(i==j);                      % standard observer: D = I

        if i==j                                   % diagonal: shared range
            lo = dLo; hi = dHi;
        else                                      % off-diagonal: its own range
            [lo,hi] = pad_range(LO(i,j), HI(i,j), 0);
        end
        edges = linspace(lo,hi,CFG.nBins+1);

        % both sigmas are scaled by the SAME per-panel maximum, so only the taller
        % distribution reaches 1 and their relative peak heights are preserved
        % (a sharply peaked sigma = 5 stays visibly taller than a spread-out
        % sigma = 20).  Normalising each to its own peak would destroy that.
        % Raw counts, then divide BOTH by the shared per-panel maximum (see
        % header comment) so relative peak heights between sigmas survive.
        c5  = histcounts(v5, edges);
        c20 = histcounts(v20,edges);
        mx  = max([c5 c20]);   c5 = c5/mx;   c20 = c20/mx;
        % Overlay narrow (blue) and broad (red) distributions as translucent bars.
        h1 = histogram(ax,'BinEdges',edges,'BinCounts',c5, ...
             'FaceColor',sdCol(1,:),'FaceAlpha',sdAlp,'EdgeColor','none');
        h2 = histogram(ax,'BinEdges',edges,'BinCounts',c20, ...
             'FaceColor',sdCol(2,:),'FaceAlpha',sdAlp,'EdgeColor','none');

        % Dashed vertical reference line at the standard observer (D = I),
        % raised above the bars so it stays visible.
        ylim(ax,[0 1.05]); yticks(ax,0:0.5:1); ytickformat(ax,'%.1f');   % 0.0 / 0.5 / 1.0
        h3 = plot(ax,[stdVal stdVal],[0 1.05],'--','Color',[0.10 0.10 0.10],'LineWidth',0.9);
        uistack(h3,'top');

        xlim(ax,[lo hi]);
        % Panel title D_{ij}, e.g. D_{LM}, via TeX subscripts.
        title(ax,sprintf('D_{%s%s}',coneNm{i},coneNm{j}), ...
              'FontWeight','normal','FontSize',figp.fontsize,'Interpreter','tex');
        % Axis labels only on the outer edge (left column / bottom row).
        if j==1, ylabel(ax,'Norm. frequency','FontWeight','Bold'); end   % peak-normalised
        if i==3, xlabel(ax,'Matrix entry','FontWeight','Bold'); end
        fig_axis_style(ax, figp);
        % per-cell scales differ by orders of magnitude (S-cone rows are ~1e-3),
        % so let each panel pick its own tick format; just keep them unrotated.
        ax.XTickLabelRotation = 0; ax.XAxis.Exponent = 0;
        % ~3 round tick values spanning this panel, guaranteed to include the
        % reference; format decimals so all labels in the panel align.
        tv = nice_ticks(lo,hi,stdVal);
        xticks(ax, tv);
        % same number of decimals for every label within a panel (panels differ
        % by ~3 orders of magnitude, so the count is derived per panel)
        xtickformat(ax, sprintf('%%.%df', decimals_for(tv)));
        ticklengthcm(ax,0.08);
        % Capture the first panel's handles once to seed the shared legend.
        if i==1 && j==1, hLeg = [h1 h2 h3]; end
    end
end

% Single shared legend (two sigmas + the standard-observer line) placed in the
% layout's north strip above all nine panels.
lg = legend(hLeg, {sprintf('\\sigma = %g nm',CFG.SDpair(1)), ...
                   sprintf('\\sigma = %g nm',CFG.SDpair(2)), ...
                   'standard observer'}, ...
            'Box','off','FontSize',figp.fontsize,'Orientation','horizontal');
lg.Layout.Tile = 'north';

% Two-column width; export a vector PDF plus a 600-dpi PNG for the Word draft.
fig.Units='centimeters'; fig.InvertHardcopy='off';
fig.Position = [6,6,figp.twocolumn, figp.twocolumn*0.80];
exportgraphics(fig, fullfile(figDir,'fig4_diagonality.pdf'),'ContentType','vector');
exportgraphics(fig, fullfile(figDir,'fig4_diagonality.png'),'Resolution',600);  % for the .docx
fprintf('Saved %s (+ .png)\n', fullfile(figDir,'fig4_diagonality.pdf'));

% ---- report the spread that the mean would have hidden ------------------
% Print mean +/- std of each D entry across observers - the between-observer
% variability that averaging D would have collapsed (the point of the figure).
fprintf('\nBetween-observer spread of D (std across %d observers):\n', nObs);
for s = 1:2
    fprintf('  sigma = %2g nm\n', CFG.SDpair(s));
    for i = 1:3
        fprintf('    ');
        for j = 1:3
            v = squeeze(Dm{s}(i,j,:));
            fprintf('D_%s%s: %6.3f +/- %.3f   ', coneNm{i},coneNm{j}, mean(v), std(v));
        end
        fprintf('\n');
    end
end

% ===================== local helpers =====================
function [lo,hi] = pad_range(lo,hi,ref)
% shared range for a group of cells, padded so the reference (D = I) stays visible.
% Degenerate/non-finite input -> a tiny symmetric window around the reference.
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo, lo = ref-1e-3; hi = ref+1e-3; return; end
% Add 6% padding and force the reference value inside the range.
pad = 0.06*(hi-lo);
lo  = min(lo, ref-pad);
hi  = max(hi, ref+pad);
end

function nd = decimals_for(tv)
% fewest decimals that still render every tick exactly (so labels line up).
% Round to 10 dp first to kill floating-point fuzz, then find the smallest d
% for which rounding the ticks to d decimals is a no-op.
tv = round(tv(:),10);
for d = 0:6
    if all(abs(tv - round(tv,d)) <= 1e-9*max(1,max(abs(tv)))), nd = d; return; end
end
nd = 6;
end

function t = nice_ticks(lo,hi,ref)
% ~3 round ticks spanning [lo,hi]; falls back to the reference value.
r = hi - lo;
if ~isfinite(r) || r <= 0, t = ref; return; end
% Target step ~ range/3, snapped to a "nice" 1/2/2.5/5/10 multiple of the
% appropriate power of ten (mag), then generate ticks on that grid.
raw = r/3;  mag = 10^floor(log10(raw));  mult = [1 2 2.5 5 10];
[~,ix] = min(abs(mult - raw/mag));  step = mult(ix)*mag;
t = (ceil(lo/step)*step) : step : (floor(hi/step)*step);
% Guarantee at least two ticks (fall back to the endpoints).
if numel(t) < 2, t = [lo hi]; end
end
