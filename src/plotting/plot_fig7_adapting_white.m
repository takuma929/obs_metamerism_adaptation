% ========================================================================
% plot_fig7_adapting_white.m   —   Figure 7: sensitivity to the adapting-white definition
%
% Median metamerism reduction at complete adaptation (D = 1), relative to the
% unadapted state, for three definitions of the adapting white, each estimated
% per observer from that observer's own reproduced stimuli:
%   brightest   : highest-luminance surface            (paper default)
%   mean LMS    : mean cone excitation over the surfaces (grey world,
%                 luminance-weighted mean chromaticity)
%   mean chrom. : mean xy chromaticity over the surfaces at the mean luminance
%                 (every surface weighted equally)
% The numerically optimised adapting white (optimal_white_upperlimit.m) is
% drawn as a thin dotted line over each primary-SD group: the upper bound that any
% white-anchored diagonal correction can reach.
%
% One panel per illuminant class, same absolute panel geometry as Figure 4 so
% the three panels tile to the two-column width. Bars carry their median value
% written inside the bar (white, rotated, one decimal)
% (median of the per-(observer, illuminant) reduction pooled over the 12
% illuminants of the class). The spread across observers is the subject of
% Figure 4 and is not repeated here; with 120,000 samples per bar the 95 %%
% confidence interval of each median is a fraction of a point (printed below).
%
% Inputs : results/adapting_white_variants.mat   (adapting_white_variants.m)
%          results/optimal_white_upperlimit.mat  (optimal_white_upperlimit.m)
%          results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat
%          (only the 'illuminants' variable, for the class of each illuminant)
% Output : figs/fig7/fig7_reduction_{daylight,smooth,peaky}.pdf   (one panel each)
%          figs/fig7/fig7_composite.pdf / .png                   (a-c + key)
% ========================================================================
clearvars; close all; clc;

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist(fullfile(projectRoot,'src'),'dir')
    projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
end
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
figDir  = fullfile(projectRoot,'figs','fig7');
if ~exist(figDir,'dir'); mkdir(figDir); end
T = load(fullfile(dataDir,'fig_parameters.mat')); figp = T.figp;

% ---- data ---------------------------------------------------------------
V = load(fullfile(projectRoot,'results','adapting_white_variants.mat'));
O = load(fullfile(projectRoot,'results','optimal_white_upperlimit.mat'));
M = load(fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat'),'illuminants');
catVec = string({M.illuminants.category});           % class of each illuminant (1..36)
sdList = V.sdList(:)'; nSD = numel(sdList);
[~, ~, nIll, ~] = size(V.d1);
assert(numel(catVec) == nIll, 'illuminant count mismatch');
classes = ["daylight","smooth","peaky"];
% variants to draw as bars, in fixed order (index into V.variants)
barVars  = ["brightest","meanLMS","meanChrom"];
barNames = {'brightest surface','mean LMS','mean chromaticity'};
barIdx   = arrayfun(@(v) find(V.variants == v,1), barVars);
nB = numel(barIdx);
% optimised-white reductions were pooled in illuminant order, nSub observers each
nSub = O.nSub;
assert(all(cellfun(@numel, O.redOpt) == nIll*nSub), 'optimal_white_upperlimit pooling order unexpected');

% ---- colours (identity, fixed order; grey = paper default) ---------------
barCols = [0.45 0.45 0.45;      % brightest surface (default)
           0.00 0.45 0.70;      % mean LMS
           0.80 0.47 0.65];     % mean chromaticity
optCol  = [0.84 0.37 0.00];     % optimised white (upper bound), distinct from the bars

% ---- geometry shared by every panel (identical to Figure 4) ---------------
figW  = figp.twocolumn/3;  inner = figW - 1.10;
figH  = 0.90 + inner*0.75 + 0.17;
axL = 1.15;  axB = 0.63;  axW = figW - axL - 0.22;  axH = inner*0.825;   % wider left margin: two-row y label
% categorical x: one group per SD at 1,2,3; bars spaced 0.27 within a group
xg = 1:nSD;  bsp = 0.27;  bw = 0.22;
off = ((1:nB) - (nB+1)/2) * bsp;
lblFS = figp.fontsize - 1;      % value labels one point smaller than tick labels

% ---- per-class statistics -------------------------------------------------
Stat = struct();
for ci = 1:numel(classes)
    sel = find(catVec == classes(ci));
    med = nan(nSD,nB); q1 = med; q3 = med; ciLo = med; ciHi = med; medOpt = nan(nSD,1);
    for si = 1:nSD
        a0 = reshape(V.d0(si,:,sel),[],1);
        for b = 1:nB
            a1 = reshape(V.d1(si,:,sel,barIdx(b)),[],1);
            r  = 100*(1 - a1./a0);
            qq = quantile(r,[0.25 0.5 0.75]);
            q1(si,b) = qq(1); med(si,b) = qq(2); q3(si,b) = qq(3);
            % distribution-free 95 %% CI of the median from order statistics
            rs = sort(r(~isnan(r))); n = numel(rs);
            lo = max(1,floor(n/2 - 1.96*sqrt(n)/2)); hi = min(n,ceil(n/2 + 1.96*sqrt(n)/2));
            ciLo(si,b) = rs(lo); ciHi(si,b) = rs(hi);
        end
        rOpt = reshape(O.redOpt{si},nSub,nIll);          % nSub x nIll (illuminant order)
        medOpt(si) = median(reshape(rOpt(:,sel),[],1),'omitnan');
    end
    Stat.(classes(ci)) = struct('med',med,'q1',q1,'q3',q3,'ciLo',ciLo,'ciHi',ciHi,'medOpt',medOpt);
end

% ---- one panel per class ---------------------------------------------------
for ci = 1:numel(classes)
    cls = classes(ci); S = Stat.(cls);
    fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');
    draw_panel(ax, S, xg, off, bw, bsp, barCols, optCol, lblFS, figp);
    xlabel_manual(ax, figp);
    fig.Units='centimeters'; fig.InvertHardcopy='off';
    fig.Position=[10,10,figW,figH];
    ax.Units='centimeters'; ax.Position=[axL,axB,axW,axH]; ticklengthcm(ax,0.1);
    pause(0.3)
    exportgraphics(fig, fullfile(figDir,sprintf('fig7_reduction_%s.pdf',cls)),'ContentType','vector');
    close(fig);
end

% ---- composite: (a)-(c) side by side with a key above ---------------------
keyH = 0.45; keyGap = 0.50;                    % key height and gap above the panel letters [cm]
fig = figure('Color','w'); fig.Units='centimeters'; fig.InvertHardcopy='off';
fig.Position=[5,5,figp.twocolumn,figH+keyGap+keyH];
letters = {'(a) daylight','(b) smooth artificial','(c) peaky artificial'};
for ci = 1:numel(classes)
    cls = classes(ci); S = Stat.(cls);
    ax = axes(fig); hold(ax,'on');
    draw_panel(ax, S, xg, off, bw, bsp, barCols, optCol, lblFS, figp);
    xlabel_manual(ax, figp);
    ax.Units='centimeters'; ax.Position=[(ci-1)*figW+axL, axB, axW, axH]; ticklengthcm(ax,0.1);
    text(ax,-0.16,1.12,letters{ci},'Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','bottom','FontWeight','bold','FontName','Arial','FontSize',figp.fontsize,'Clipping','off');
    if ci > 1, ylabel(ax,''); end
end
% key: three filled swatches + the dashed upper-bound line, one row, centred
kax = axes(fig,'Units','centimeters','Position',[0, figH+keyGap, figp.twocolumn, keyH], ...
    'Visible','off','XLim',[0 figp.twocolumn],'YLim',[0 1]); hold(kax,'on');
entries = [barNames, {'optimized white (upper bound)'}];
swW = 0.32; gap = 0.18; pad = 0.45;
% estimate widths from character counts (Arial ~0.11 cm per char at 7 pt)
wTxt = cellfun(@(s) 0.11*numel(s), entries);
wTot = sum(swW + gap + wTxt) + pad*(numel(entries)-1);
x0 = (figp.twocolumn - wTot)/2;
for e = 1:numel(entries)
    if e <= nB
        patch(kax,'XData',x0+[0 swW swW 0],'YData',[0.3 0.3 0.7 0.7],'FaceColor',barCols(e,:),'EdgeColor','none');
    else
        line(kax,x0+[0 swW],[0.5 0.5],'Color',optCol,'LineStyle',':','LineWidth',0.5);
    end
    text(kax,x0+swW+gap,0.5,entries{e},'FontName','Arial','FontSize',figp.fontsize, ...
        'VerticalAlignment','middle','HorizontalAlignment','left');
    x0 = x0 + swW + gap + wTxt(e) + pad;
end
pause(0.3)
exportgraphics(fig, fullfile(figDir,'fig7_composite.pdf'),'ContentType','vector');
exportgraphics(fig, fullfile(figDir,'fig7_composite.png'),'Resolution',400);
close(fig);

% ---- numbers for the caption / text ---------------------------------------
fprintf('\nMedian reduction at D = 1 [%%] (95%% CI of the median; IQR across observer x illuminant), per class and SD\n');
for ci = 1:numel(classes)
    S = Stat.(classes(ci));
    fprintf('[%s]\n   SD |', classes(ci));
    for b = 1:nB, fprintf(' %-34s |', barNames{b}); end
    fprintf(' optimised\n');
    for si = 1:nSD
        fprintf('  %3g |', sdList(si));
        for b = 1:nB, fprintf(' %5.1f [%5.1f-%5.1f] (%4.1f-%4.1f) |', S.med(si,b), S.ciLo(si,b), S.ciHi(si,b), S.q1(si,b), S.q3(si,b)); end
        fprintf(' %5.1f\n', S.medOpt(si));
    end
end
fprintf('\nFig 7 done (adapting-white definitions).\n');

% ------------------------------------------------------------------------
function draw_panel(ax, S, xg, off, bw, bsp, barCols, optCol, lblFS, figp)
% Grouped bars (median reduction) with IQR error bars and value labels, plus the
% optimised-white upper bound as a dashed line across each SD group.
nSD = numel(xg); nB = size(S.med,2);
for si = 1:nSD
    for b = 1:nB
        xc = xg(si) + off(b); m = S.med(si,b);
        patch(ax,'XData',xc+[-1 1 1 -1]*bw/2,'YData',[0 0 m m], ...
            'FaceColor',barCols(b,:),'EdgeColor','none');
        % value written inside the bar, rotated 90 deg, ending just below the top
        text(ax,xc,m-1.5,sprintf('%.1f',m),'Rotation',90,'HorizontalAlignment','right', ...
            'VerticalAlignment','middle','FontName','Arial','FontSize',lblFS,'Color','w');
    end
    % upper bound: dotted line spanning the group, value centred above it
    xl = xg(si)+off(1)-bw/2; xr = xg(si)+off(end)+bw/2;
    line(ax,[xl xr],[1 1]*S.medOpt(si),'Color',optCol,'LineStyle',':','LineWidth',0.5);
    text(ax,xg(si),S.medOpt(si)+1,sprintf('%.1f',S.medOpt(si)),'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontName','Arial','FontSize',lblFS,'Color',optCol);
end
ylabel(ax,{'Metamerism','reduction [%]'},'FontWeight','Bold');
xlim(ax,[xg(1)-0.55, xg(end)+0.55]); ylim(ax,[0 100]);
xticks(ax,xg); ax.XTickLabel = string([5 10 20]); ax.XTickLabelRotation = 0;
yticks(ax,0:20:100); ytickformat(ax,'%.0f');
fig_axis_style(ax, figp);
end

function xlabel_manual(ax, figp)
text(ax,0.5,-0.12,'Primary SD [nm]','Units','normalized','HorizontalAlignment','center', ...
    'VerticalAlignment','top','FontWeight','bold','FontName','Arial','FontSize',figp.fontsize,'Clipping','off');
end
