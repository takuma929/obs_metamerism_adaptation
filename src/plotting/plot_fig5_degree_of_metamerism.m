% ========================================================================
% plot_fig5_degree_of_metamerism.m   —   Figure 4: observer metamerism vs degree of adaptation (all panels)
%
% NOTE ON FIGURE NUMBERING: this SCRIPT is named "fig5", but it produces the
% MANUSCRIPT's Figure 4 (script and manuscript numbers differ because the
% figures were renumbered during revision). Output files are written to
% figs/fig5/ (the script-number folder), not figs/fig4/.
%
% Observer metamerism as a function of the degree of chromatic adaptation D,
% scored in CAM16-UCS. At each D every observer AND the standard-observer
% reference are placed at the same adaptation state (own reproduced white for the
% individuals, reference white for the standard observer), so the global
% adaptation-state cast cancels and only the between-observer metamerism remains.
% Three levels are shown:
%   D = 0    : no chromatic adaptation (white not discounted)
%   D = 0.5  : partial adaptation
%   D = 1    : complete adaptation to the adopted white
% All three use the identical full CIECAM16 -> CAM16-UCS read-out; D controls only
% the CAT16 white-discounting gain.
%   (a)  one compact box per D at each primary SD, per class.
%   (a2) per-observer metamerism reduction 100*(1 - d_D/d_0) at D = 0.5 and D = 1.
%
% Variants of the adaptation's functional form (post-compression gain, a
% Nayatani-type power law, an opponent-stage second site) are NOT plotted: a gain
% commutes with the near-power-law compression, so moving it across the
% compression is the same operation up to a constant. Reported in the text as an
% algebraic identity. See compare_nonlinear_adapt_ucs.m for the numbers.
%
% Input : results/compare_nl_adapt_ucs.mat (from compare_nonlinear_adapt_ucs.m)
% Output: figs/fig5/fig5{a,a_reduction}_*.pdf (vector)
% ========================================================================
clearvars; close all; clc;

% --- Resolve the project root from this file's location -------------------
% mfilename('fullpath') is .../src/plotting/plot_fig5_...; three fileparts()
% strip filename -> plotting -> src, leaving the repository root.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
% Fallback to an absolute path if the script is run in a way where the
% inferred root has no src/ folder (e.g. copy-pasted into the console).
if ~exist(fullfile(projectRoot,'src'),'dir')
    projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
end
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
% Vector PDFs for this figure go here; created on first run.
figDir  = fullfile(projectRoot,'figs','fig5');
if ~exist(figDir,'dir'); mkdir(figDir); end
% Precomputed metamerism scores (produced by compare_nonlinear_adapt_ucs.m).
inMat = fullfile(projectRoot,'results','compare_nl_adapt_ucs.mat');
assert(exist(inMat,'file')==2, 'Missing %s (run compare_nonlinear_adapt_ucs first).', inMat);
% Res: per-class results struct; classes: cellstr/string of class names.
S = load(inMat);  Res = S.Res;  classes = S.classes;
% figp: shared figure-style parameters (fonts, column widths, etc.).
T = load(fullfile(dataDir,'fig_parameters.mat')); figp = T.figp;

% Human-readable titles keyed by internal class name (map is available for
% annotation; the panels themselves are labelled in the caption).
clsTitle = containers.Map({'daylight','smooth','peaky'}, ...
    {'daylight','smooth daylight','peaky daylight'});

%% ---- (a) CAM16-UCS metamerism vs degree of adaptation, one figure per class
% One compact box per degree of adaptation D (0, 0.5, 1) at each primary SD, side
% by side. Each illuminant class is an independent figure. Sequential single-hue
% ramp: light (D=0, no adaptation) -> dark (D=1, complete adaptation).
% Degree-of-adaptation levels D (same list for every class); one box per D.
Dlist        = Res.(char(classes(1))).Dlist(:)';
% Legend labels "D = 0", "D = 0.5", "D = 1".
methodNames  = arrayfun(@(d) sprintf('D = %g',d), Dlist, 'uni', 0);
methodCols   = [0.63 0.09 0.11; 0.90 0.42 0.13; 0.98 0.73 0.48];   % sequential red->orange (D 0->1)
for ci = 1:numel(classes)
    cls = char(classes(ci)); if ~isfield(Res,cls), continue; end
    % R.sd: primary spectral bandwidths (SD, nm); R.pD{mi}: nSD x N array of
    % per-observer CAM16-UCS dE' at the mi-th D level (N = observers).
    R = Res.(cls);  sdList = R.sd(:)';  nM = numel(Dlist);

    fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');
    % Lay out the nM boxes of one SD as a tight cluster centred on that SD tick.
    % scale ties box spacing to the SD-axis gap so it looks the same at any set
    % of SDs; off = symmetric per-box offsets; bw = box width (data units).
    % fixed inter-box spacing, centred on each SD -> geometry independent of nM
    scale = min(diff(sdList))/5;  spacing = 1.0;
    off = ((1:nM) - (nM+1)/2) * spacing;
    bw  = spacing*scale*0.72;
    for mi = 1:nM
        A = R.pD{mi};                                     % nSD x N
        % One compact box per SD showing the spread of dE' over all observers.
        for si = 1:numel(sdList)
            draw_cbox(ax, sdList(si)+off(mi)*scale, A(si,:), bw, methodCols(mi,:));
        end
    end
    ylabel(ax,'\DeltaE''','FontWeight','Bold');   % CAM16-UCS distance, as defined in §2.3
    xlim(ax,[min(sdList)-3, max(sdList)+3]); ylim(ax,[0 12]);
    xticks(ax,sdList); ax.XTickLabel=string(sdList); ax.XTickLabelRotation=0;
    yticks(ax,0:2:12); ytickformat(ax,'%.0f');
    fig_axis_style(ax, figp);
    % Shared x-axis label placed manually below the axes (normalized units,
    % clipping off so it can sit outside the axes box).
    text(ax,0.5,-0.12,'Primary SD [nm]','Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontWeight','bold','FontName','Arial','FontSize',figp.fontsize,'Clipping','off');
    % --- Absolute sizing so every panel is byte-identical across classes -----
    % InvertHardcopy off preserves the white background on export.
    fig.Units='centimeters'; fig.InvertHardcopy='off';
    % One third of the journal two-column width per figure; inner = drawable area.
    figW  = figp.twocolumn/3;  inner = figW - 1.10;
    figH  = 0.90 + inner*0.75 + 0.17;
    axL = 0.90;  axB = 0.63;  axW = inner*1.08;  axH = inner*0.825;
    fig.Position=[10,10,figW,figH];
    ax.Units='centimeters'; ax.Position=[axL,axB,axW,axH]; ticklengthcm(ax,0.1);
    % Let the renderer settle before exporting (avoids clipped/partial PDFs).
    pause(0.5)
    exportgraphics(fig, fullfile(figDir,sprintf('fig5a_metamerism_%s.pdf',cls)),'ContentType','vector');
    close(fig);
end

%% ---- (a2) metamerism reduction [%], one figure per class -----------------
% Per-observer percentage reduction of CAM16-UCS metamerism relative to the
% unadapted state (D=0): 100*(1 - d_D/d_0), at D = 0.5 and D = 1.
% Reduction panel omits D=0 (its reduction vs itself is 0 by definition).
redIdx        = find(Dlist>0);                     % non-zero D levels (0.5, 1)
rMethodNames  = methodNames(redIdx);
rMethodCols   = methodCols(redIdx,:);
for ci = 1:numel(classes)
    cls = char(classes(ci)); if ~isfield(Res,cls), continue; end
    R = Res.(cls);  sdList = R.sd(:)';  nM = numel(redIdx);
    % Index of the unadapted baseline (D=0) used as the reference denominator.
    iD0 = find(Dlist==0,1);

    fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');
    % Zero reference line: no change relative to the unadapted state.
    yline(ax,0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
    scale = min(diff(sdList))/5;  spacing = 1.0;      % same geometry as panel (a)
    off = ((1:nM) - (nM+1)/2) * spacing;
    bw  = spacing*scale*0.72;
    for mi = 1:nM
        % Per-observer % reduction of metamerism vs D=0: element-wise so every
        % observer is compared with its OWN unadapted dE' (kept as a matched pair).
        redPct = 100*(1 - R.pD{redIdx(mi)} ./ R.pD{iD0});   % nSD x N per-observer reduction
        for si = 1:numel(sdList)
            draw_cbox(ax, sdList(si)+off(mi)*scale, redPct(si,:), bw, rMethodCols(mi,:));
        end
    end
    ylabel(ax,'Metamerism reduction [%]','FontWeight','Bold');
    xlim(ax,[min(sdList)-3, max(sdList)+3]); ylim(ax,[0 100]);
    xticks(ax,sdList); ax.XTickLabel=string(sdList); ax.XTickLabelRotation=0;
    yticks(ax,0:20:100); ytickformat(ax,'%.0f');
    fig_axis_style(ax, figp);
    % Shared manual x-axis label (see panel (a)).
    text(ax,0.5,-0.12,'Primary SD [nm]','Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontWeight','bold','FontName','Arial','FontSize',figp.fontsize,'Clipping','off');
    % Identical absolute sizing to panel (a) so the two panels are matched.
    fig.Units='centimeters'; fig.InvertHardcopy='off';
    figW  = figp.twocolumn/3;  inner = figW - 1.10;
    figH  = 0.90 + inner*0.75 + 0.17;
    axL = 0.90;  axB = 0.63;  axW = inner*1.08;  axH = inner*0.825;
    fig.Position=[10,10,figW,figH];
    ax.Units='centimeters'; ax.Position=[axL,axB,axW,axH]; ticklengthcm(ax,0.1);
    pause(0.5)
    exportgraphics(fig, fullfile(figDir,sprintf('fig5a_reduction_%s.pdf',cls)),'ContentType','vector');
    close(fig);
end
%% ---- median values per panel, for annotating the figure ------------------
% Panel (a): median dE' at each D (= the white line of each box).
% Panel (a2): the box shows the PER-OBSERVER reduction vs D=0, so its median is
% not identical to 100*(1 - median/median); both are printed.
for ci = 1:numel(classes)
    cls = char(classes(ci)); if ~isfield(Res,cls), continue; end
    R = Res.(cls);  sdList = R.sd(:)';  iD0 = find(Dlist==0,1);
    fprintf('\n================ [%s] ================\n', cls);
    % Panel (a) numbers: box median = median dE' across observers, per SD and D.
    fprintf('(a) median dE'' at each D (= box median)\n');
    fprintf('   SD |');  fprintf(' %10s |', methodNames{:});  fprintf('\n');
    for si = 1:numel(sdList)
        fprintf('  %3g |', sdList(si));
        for mi = 1:numel(Dlist)
            % 'omitnan' ignores observers with undefined matches for this SD/D.
            fprintf(' %10.3f |', median(R.pD{mi}(si,:),'omitnan'));
        end
        fprintf('\n');
    end
    fprintf('(a2) metamerism reduction vs D=0 [%%]\n');
    fprintf('   SD |');
    for mi = redIdx, fprintf(' %-9s per-obs | med/med |', methodNames{mi}); end
    fprintf('\n');
    for si = 1:numel(sdList)
        fprintf('  %3g |', sdList(si));
        for mi = redIdx
            % 'per-obs' = median of the per-observer reduction (what the box shows);
            red = 100*(1 - R.pD{mi}(si,:) ./ R.pD{iD0}(si,:));
            % 'med/med' = reduction of the two medians. They differ because the
            % median of a ratio is not the ratio of medians; both are printed.
            fprintf(' %17.1f | %7.1f |', median(red,'omitnan'), ...
                100*(1 - median(R.pD{mi}(si,:),'omitnan')/median(R.pD{iD0}(si,:),'omitnan')));
        end
        fprintf('\n');
    end
end

fprintf('\nFig 5 done (CAM16-UCS metamerism vs degree of adaptation D = 0, 0.5, 1).\n');

% ------------------------------------------------------------------------
function draw_cbox(ax, xc, d, w, col)
% Compact filled box for one group: IQR box, white median line, 1*IQR whiskers.
% Inputs: ax = target axes; xc = box centre on the x-axis; d = data vector
% (observers); w = box width in data units; col = fill colour.
% Drop NaNs (observers without a valid value); nothing to draw if all missing.
d = d(~isnan(d));
if isempty(d), return; end
% Quartiles Q1/median/Q3 and the interquartile range.
q  = quantile(d,[0.25 0.5 0.75]);  iqr = q(3)-q(1);
% Whisker ends = most extreme data still within 1*IQR of the box edges.
lo = min(d(d >= q(1)-iqr));  hi = max(d(d <= q(3)+iqr));
% Vertical whisker line spanning [lo hi].
line(ax,[xc xc],[lo hi],'Color',col,'LineWidth',0.8);
% Filled rectangle spanning the IQR (Q1..Q3), centred at xc.
patch(ax,'XData',xc+[-1 1 1 -1]*w/2,'YData',[q(1) q(1) q(3) q(3)], ...
      'FaceColor',col,'EdgeColor','none');
% White median tick across the box.
line(ax,xc+[-1 1]*w/2,[q(2) q(2)],'Color','w','LineWidth',0.8);
end
