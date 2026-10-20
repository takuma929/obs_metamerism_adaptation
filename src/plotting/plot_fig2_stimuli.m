% ========================================================================
% plot_fig2_stimuli.m  -  Figure 2: reflectances (a); illuminant classes as
%   three independent panels daylight (b), smooth (c), peaky (d) each 25% of the
%   two-column width; primaries spectra e=SD5 nm and f=SD20 nm (one-column wide,
%   half height, stacked vertically); xy gamut (g). Also saves axis-free,
%   transparent-background icon versions (vector PDF) of the b/c/d representative
%   illuminant spectra and the e/f primaries, for use as small figure symbols.
% Self-contained; run on its own:  >> plot_fig2_stimuli
% Loads only small fields (no 9 GB byIll) so it is fast.
%
% MANUSCRIPT FIGURE
%   Optics Express: "Adaptation reduces inter-observer metamerism with
%   narrow-band primaries." Figure 2 documents the stimuli: the reflectance
%   set (CIE 99 CFI), the three illuminant classes (daylight / smooth / peaky),
%   the Gaussian virtual-display primaries (peaks 450/540/610 nm) at two
%   bandwidths, and the display gamut versus the retained reflectance
%   chromaticities in the CIE 1931 xy plane.
%
% INPUT .mat
%   results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat
%   Only small fields are loaded (wls, params, display, keepIdx_global,
%   illuminants, cone_fund); the multi-gigabyte per-illuminant field is skipped.
%
% OUTPUT FILES (written into figs/fig2/)
%   fig2a_reflectances.png                       - (a) retained reflectances,
%       each coloured by its own sRGB rendering under D65.
%   fig2b/c/d_illuminant_<class>.png             - (b,c,d) each illuminant class,
%       members max-normalised, ~6500 K member bold. Labelled panels.
%   fig2b/c/d_illuminant_icon_<class>.pdf        - axis-free transparent icons
%       of the bold representative SPD of each class.
%   fig2e/f_primaries_spectrum_SD<n>.png         - (e,f) Gaussian RGB primaries
%       at sigma = 5 and 20 nm. Labelled panels.
%   fig2e/f_primaries_icon_SD<n>.pdf             - axis-free transparent icons
%       of the same primaries.
%   fig2g_gamut.png                              - (g) CIE 1931 xy chromaticity
%       diagram with the reflectance cloud and the per-bandwidth display gamut
%       triangles.
%
% NOTE: local helper functions live at the bottom of this file. Several
%   (xyz2lab_wref, deltaE2000, renderFilledSpectrumTile, makeFilledGraphOnly,
%   make_srgb_mosaic) are shared colorimetry/plotting utilities kept alongside
%   the figure scripts; this figure actively uses xyz_to_srgb_d65_wref (panel a)
%   and resample_T_to_wls (setup). Each helper is documented at its definition.
% ========================================================================
clearvars; close all; clc;

% Resolve project root from this script's own location so it can be run
% from any working directory (e.g. directly from src/plotting/).
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist(fullfile(projectRoot,'src'),'dir')               % section run: mfilename points at a temp file
    projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
end
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
[~,thisName] = fileparts(mfilename('fullpath'));            % e.g. plot_fig1_observers
subFig  = regexp(thisName,'^plot_(fig[^_]*)','tokens','once'); % -> 'fig1' ... 'figS'
if isempty(subFig), subFig = {'fig2'}; end                 % fallback for section runs
figDir  = fullfile(projectRoot,'figs',subFig{1});     % per-figure subfolder
if ~exist(figDir,'dir'); mkdir(figDir); end
inMat=fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
assert(exist(inMat,'file')==2,'Missing results: %s.',inMat);
% ---- load results (small fields only; no 9 GB byIll) ----
S = load(inMat,'wls','params','display','keepIdx_global','illuminants','cone_fund');
% wls: common wavelength sampling grid (column) shared by every spectrum here.
wls       = S.wls(:);
% sdList: Gaussian primary bandwidths sigma (nm) analysed, as a row vector.
sdList    = S.params.sdList(:)';
nW        = numel(wls);
% Pset: stored Gaussian display primaries (spectral power) per bandwidth.
Pset      = S.display.Pset;
% Mdisp_set: per-bandwidth 3x3 matrix of primary tristimulus values (XYZ cols).
Mdisp_set = S.display.Mdisp_set;
% keepIdx: indices of the reflectance samples retained after screening.
keepIdx   = S.keepIdx_global(:);
% cone_fund: [Nobs x nW x 3] individual-observer L/M/S cone fundamentals.
cone_fund = S.cone_fund;

% ---- figure style parameters ----
T = load(fullfile(dataDir,'fig_parameters.mat'));
figp = T.figp;

% ---- reflectance spectra (99 CIE test samples, resampled to wls) ----
% CSV column 1 = wavelength (1 nm), columns 2..end = the 99 CIE CFI reflectances.
Rmat     = readmatrix(fullfile(dataDir,'CIE_srf_cfi_1nm.csv'));
wl_ref   = Rmat(:,1);
R_99_1nm = Rmat(:,2:end);
% Resample onto wls (pchip) and clamp reflectance to the physical range [0,1].
R_wlsx99 = min(max(interp1(wl_ref,R_99_1nm,wls,'pchip','extrap'),0),1);
% Keep only the screened subset used throughout the analysis.
R_keep   = R_wlsx99(:,keepIdx);

% ---- CIE 1931 2-deg colour-matching functions ----
% Reference (standard) observer used to compute XYZ tristimulus values.
Axyz  = readmatrix(fullfile(dataDir,'CIE1931XYZ.txt'));
% Interpolate the xbar/ybar/zbar CMFs onto wls; T_xyz is 3 x nW (rows X,Y,Z).
T_xyz = interp1(Axyz(:,1),Axyz(:,2:4),wls,'pchip',0)';

% ---- representative illuminants (daylight D65 and peaky) ----
% Each illuminant carries a category label and an isRep flag marking the
% representative member of its class; pick the daylight (D65) representative.
illCats   = string({S.illuminants.category});
isRep     = [S.illuminants.isRep];
idxD65    = find(isRep & illCats=="daylight",1,'first');
d65       = S.illuminants(idxD65).spd(:);
% Also fetch the representative "peaky" illuminant if present (unused here).
idxHouser = find(isRep & illCats=="peaky",1,'first');
if ~isempty(idxHouser), spdHouser = S.illuminants(idxHouser).spd(:); else, spdHouser = []; end

% ---- world colours under D65 and the adopted (brightest) white ----
% Colour signal = illuminant SPD times reflectance (element-wise per wavelength).
nSamples      = size(R_keep,2);
SPD_keep_D65  = d65 .* R_keep;
% Tristimulus of every retained sample under D65: 3 x nSamples.
XYZ_world_D65 = T_xyz * SPD_keep_D65;
% Adopted white = the brightest sample (max luminance Y); used as the
% reference white for sRGB rendering of sample colours in panel (a).
[~,idxWhite]  = max(XYZ_world_D65(2,:));
XYZ_white_D65 = XYZ_world_D65(:,idxWhite);
% Guard against non-positive components before any division/normalisation.
XYZ_white_D65(XYZ_white_D65<=0) = eps;

% ---- standard-observer cone fundamentals (Stockman-Sharpe 2-deg) ----
% Reference LMS fundamentals, resampled from their native grid onto wls.
Lstd = load(fullfile(dataDir,'T_cones_ss2.mat'));
T_cones_std = resample_T_to_wls(Lstd.S_cones_ss2,Lstd.T_cones_ss2,wls);

%% ============================================================
% (1) Figure: Reflectance spectra (kept only; GLOBAL keep-set)
%% ============================================================
fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

% Colour signal and tristimulus of each retained reflectance under D65.
SPD_keep_D65 = d65 .* R_keep;                 % nW x nSamples
XYZ_keep     = T_xyz * SPD_keep_D65;          % 3 x nSamples

% Render each sample to sRGB relative to the adopted white so every
% reflectance curve is drawn in its own perceived colour.
XYZ_white = XYZ_white_D65(:);                 % 3 x 1
rgbLines  = xyz_to_srgb_d65_wref(XYZ_keep, XYZ_white);  % nSamples x 3

% Draw each reflectance spectrum; *0.9 slightly darkens for on-white legibility.
for i = 1:nSamples
    plot(ax, wls, R_keep(:,i), 'LineWidth', 0.35, 'Color', rgbLines(i,:)*0.9);
end

xlabel(ax,'Wavelength (nm)','FontWeight','Bold');
ylabel(ax,'Reflectance','FontWeight','Bold');

fig.Units          = 'centimeters';
fig.Color          = 'w';
fig.InvertHardcopy = 'off';

ax.FontName  = 'Arial';
ax.FontSize  = figp.fontsize;
ax.XColor    = 'k';
ax.YColor    = 'k';
ax.LineWidth = 0.5;
ax.Units     = 'centimeters';

% panel_geom: helper returning figure width, axes position, figure height for
% a panel occupying 1/4 of the two-column width (so a-d tile horizontally).
[fw,ap,fh]   = panel_geom(4,figp,0.9);   % 25% of two-column width (row panels a-d)
fig.Position = [10,10,fw,fh];
ax.Position  = ap;

xlim([390-10 780+10]); ylim([0 1]);
xticks([400,500 600 700 780]);
yticks([0 0.5 1.0]);
ax.XTickLabel = {'400','500','600','700','780'};
ax.YTickLabel = {'0.0','0.5','1.0'};

ticklengthcm(ax,0.05);
grid(ax,'off'); box(ax,'off');
% Light-grey axes background to keep the coloured curves legible.
ax.Color = ones(1,3)*0.97;

exportgraphics(fig, fullfile(figDir,'fig2a_reflectances.png'), 'Resolution',300);

%% ============================================================
% (b,c,d) Three CCT-matched illuminant classes, each an INDEPENDENT panel
%   (max-normalised SPDs; the member closest to 6500 K is drawn bold).
%   Each panel is 25% of the two-column width so a-d tile horizontally.
%% ============================================================
% Class name, plot colour and panel letter for the three illuminant classes.
clsN   = {'daylight','smooth','peaky'};
clsC   = {[.10 .35 .85],[.10 .55 .20],[.85 .10 .10]};
clsLet = {'b','c','d'};

for cc = 1:3
    fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

    % gather this class's member SPDs, max-normalise each
    % (peak of every SPD scaled to 1 so shape, not level, is compared).
    memb = find(illCats==clsN{cc});
    spdM = zeros(nW,numel(memb));
    for kk = 1:numel(memb), spdM(:,kk) = S.illuminants(memb(kk)).spd(:); end
    spdM = spdM ./ max(spdM,[],1);

    % highlight the member closest to 6500 K
    % (CCT-matched representative; ir indexes the nearest-to-6500 K SPD).
    cctv = [S.illuminants(memb).CCT];
    [~,ir] = min(abs(cctv-6500));
    % Faint members (alpha 0.28) form a band; the representative is drawn bold.
    for kk = 1:numel(memb)
        plot(ax,wls,spdM(:,kk),'-','Color',[clsC{cc} 0.28],'LineWidth',0.4);
    end
    plot(ax,wls,spdM(:,ir),'-','Color',clsC{cc},'LineWidth',1.1);

    % axes styling
    xlim(ax,[390 730]); ylim(ax,[0 1.05]);
    xticks(ax,400:100:700); ytickformat(ax,'%.1f');
    ax.FontName = 'Arial'; ax.Color = ones(3,1)*0.97; ax.FontSize = figp.fontsize;
    ax.XColor = 'k'; ax.YColor = 'k'; ax.LineWidth = 0.5;
    box(ax,'off'); grid(ax,'on'); ax.GridAlpha = 0.1; ax.XGrid = 'off';

    % labels (y-label on every panel so a-d export to identical width; no title)
    xlabel(ax,'Wavelength [nm]','FontWeight','bold');
    ylabel(ax,'Rel. power','FontWeight','bold');

    % geometry: 25% of two-column width (matches panels a-d)
    fig.Units = 'centimeters'; fig.InvertHardcopy = 'off';
    [fw,ap,fh] = panel_geom(4,figp,0.9);
    fig.Position = [10,10,fw,fh];
    ax.Units = 'centimeters'; ax.Position = ap; ticklengthcm(ax,0.08);

    % One PNG per class: fig2b_illuminant_daylight.png, etc.
    outName = sprintf('fig2%s_illuminant_%s.png', clsLet{cc}, clsN{cc});
    exportgraphics(fig, fullfile(figDir,outName),'Resolution',300);
    close(fig);
end

% (b,c,d icons) The bold representative SPD of each class (member closest to
%   6500 K), as an axis-free, transparent-background symbol filled with a
%   strongly desaturated tint of the class colour (same hue as the line).
%   Vector PDF, matching the primary-icon style/size.
% iconW/iconH: symbol size in cm. mixClass blends the class colour toward
% white for a low-saturation fill (fill = 15% class colour + 85% white).
iconW = 1.8;  iconH = 0.9;  mixClass = 0.85;   % fill = 15% class colour + 85% white
for cc = 1:3
    % Rebuild this class's max-normalised member SPDs and its representative.
    memb = find(illCats==clsN{cc});
    spdM = zeros(nW,numel(memb));
    for kk = 1:numel(memb), spdM(:,kk) = S.illuminants(memb(kk)).spd(:); end
    spdM = spdM ./ max(spdM,[],1);
    cctv = [S.illuminants(memb).CCT];
    [~,ir] = min(abs(cctv-6500));               % the bold representative member
    yRep = spdM(:,ir);

    % Outline in the class colour; fill in a desaturated tint of the same hue.
    lineCol = clsC{cc};
    fillCol = (1-mixClass)*lineCol + mixClass*[1 1 1];   % low-saturation same-hue fill

    % Transparent figure/axes so the exported glyph has no background.
    fig = figure('Color','none'); ax = axes(fig); hold(ax,'on');
    % Closed polygon fills the area under the representative SPD.
    fill(ax, [wls; flipud(wls)], [yRep; zeros(size(yRep))], fillCol, 'EdgeColor','none','FaceAlpha',1.0);
    plot(ax, wls, yRep, '-', 'LineWidth', 0.3, 'Color', lineCol);

    xlim(ax,[390 730]); ylim(ax,[0.01 1.02]);
    axis(ax,'off'); ax.Color = 'none';

    fig.Units='centimeters'; fig.InvertHardcopy='off';
    fig.Position = [10,10,iconW,iconH];
    ax.Units='centimeters'; ax.Position = [0 0 iconW iconH];

    % Vector PDF with no background: a clean, scalable class symbol.
    outName = sprintf('fig2%s_illuminant_icon_%s.pdf', clsLet{cc}, clsN{cc});
    exportgraphics(fig, fullfile(figDir,outName), 'ContentType','vector','BackgroundColor','none');
    close(fig);
end

% (e,f) Display RGB primaries spectra as TWO independent labelled panels
%   e = SD 5 nm, f = SD 20 nm. Each is one-column wide and half that height,
%   so e can be stacked vertically above f (bottom-left of the final layout).
%% ============================================================
% Primary channel colours: B blue, G green, R red (peaks at 450/540/610 nm).
colB = [0.10 0.30 0.85];
colG = [0.10 0.60 0.10];
colR = [0.85 0.10 0.10];

% Desaturated same-hue fills (15% channel colour + 85% white) for the areas.
mixToWhite = 0.85;
fillB = (1-mixToWhite)*colB + mixToWhite*[1 1 1];
fillG = (1-mixToWhite)*colG + mixToWhite*[1 1 1];
fillR = (1-mixToWhite)*colR + mixToWhite*[1 1 1];

primSD  = [5 20];               % e -> SD 5 nm, f -> SD 20 nm
primLet = {'e','f'};

for k = 1:numel(primSD)
    sdTarget = primSD(k);
    % Confirm the requested bandwidth exists in the analysed set.
    si = find(sdList == sdTarget, 1, 'first');
    assert(~isempty(si), 'sdTarget=%g not found in sdList.', sdTarget);

    % Rebuild the Gaussian primaries on a 1 nm grid so the narrowband curves
    % render smoothly (the stored Pset lives on the coarser analysis grid).
    wlF   = (wls(1):1:wls(end))';                    % 1 nm wavelength axis
    peaks = S.params.peakBase(:)';                   % [B G R] peak wavelengths
    % Gaussian per channel: exp(-0.5*((lambda-peak)/sigma)^2), sigma = sdTarget.
    G     = exp(-0.5*((wlF - peaks)./sdTarget).^2);  % nWf x 3
    G     = G ./ max(G,[],1);                         % unit-max per channel
    yB = G(:,1); yG = G(:,2); yR = G(:,3);

    fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

    % Filled area under each primary (closed polygon curve-to-baseline).
    fill(ax, [wlF; flipud(wlF)], [yB; zeros(size(yB))], fillB, 'EdgeColor','none', 'FaceAlpha', 1.0);
    fill(ax, [wlF; flipud(wlF)], [yG; zeros(size(yG))], fillG, 'EdgeColor','none', 'FaceAlpha', 1.0);
    fill(ax, [wlF; flipud(wlF)], [yR; zeros(size(yR))], fillR, 'EdgeColor','none', 'FaceAlpha', 1.0);

    % Outline each primary on top of its fill.
    plot(ax, wlF, yB, '-', 'LineWidth', 0.5, 'Color', colB);
    plot(ax, wlF, yG, '-', 'LineWidth', 0.5, 'Color', colG);
    plot(ax, wlF, yR, '-', 'LineWidth', 0.5, 'Color', colR);

    xlabel(ax,'Wavelength (nm)','FontWeight','Bold');
    ylabel(ax,'Rel. power','FontWeight','Bold');

    xlim([390-10 780+10]); ylim([0 1.00]);
    xticks([400 500 600 700 780]);
    yticks([0 0.5 1.0]);
    ax.XTickLabel = {'400','500','600','700','780'};
    ax.YTickLabel = {'0.0','0.5','1.0'};

    ax.FontName='Arial'; ax.FontSize=figp.fontsize;
    ax.XColor='k'; ax.YColor='k'; ax.LineWidth=0.5;
    ticklengthcm(ax,0.05);
    grid(ax,'off'); box(ax,'off'); ax.Color = ones(1,3)*0.97;

    % geometry: one-column width, half that height (for vertical e/f stacking)
    fig.Units='centimeters'; fig.InvertHardcopy='off';
    fig.Position = [10,10,figp.onecolumn, figp.onecolumn/2];
    ax.Units='centimeters';
    ax.Position  = [1.0 0.8 figp.onecolumn-1.2 figp.onecolumn/2-1.0];

    % One PNG per bandwidth: fig2e_primaries_spectrum_SD5.png, ...SD20.png.
    outName = sprintf('fig2%s_primaries_spectrum_SD%g.png', primLet{k}, sdTarget);
    exportgraphics(fig, fullfile(figDir, outName), 'Resolution',300);
    close(fig);
end

%% ============================================================
% (e,f icons) The SAME SD 5 nm and SD 20 nm RGB primaries, saved as small
%   axis-free, transparent-background symbols (no ticks/labels/frame), so each
%   can be dropped in as a glyph — e.g. next to a legend entry or an axis tick
%   — to denote the primary bandwidth. Transparent PNG at high resolution.
%% ============================================================
iconW = 1.8;  iconH = 0.9;      % symbol size in cm (landscape tile)
for k = 1:numel(primSD)
    sdTarget = primSD(k);

    % rebuild the Gaussian primaries on a 1 nm grid (smooth narrowband curves)
    % (identical construction to the e/f panels above, just no axes).
    wlF   = (wls(1):1:wls(end))';
    peaks = S.params.peakBase(:)';
    G     = exp(-0.5*((wlF - peaks)./sdTarget).^2);
    G     = G ./ max(G,[],1);
    yB = G(:,1); yG = G(:,2); yR = G(:,3);

    % Transparent figure/axes so the exported glyph has no background.
    fig = figure('Color','none'); ax = axes(fig); hold(ax,'on');

    % Filled areas plus thin outlines for the three primaries.
    fill(ax, [wlF; flipud(wlF)], [yB; zeros(size(yB))], fillB, 'EdgeColor','none', 'FaceAlpha', 1.0);
    fill(ax, [wlF; flipud(wlF)], [yG; zeros(size(yG))], fillG, 'EdgeColor','none', 'FaceAlpha', 1.0);
    fill(ax, [wlF; flipud(wlF)], [yR; zeros(size(yR))], fillR, 'EdgeColor','none', 'FaceAlpha', 1.0);
    plot(ax, wlF, yB, '-', 'LineWidth', 0.3, 'Color', colB);
    plot(ax, wlF, yG, '-', 'LineWidth', 0.3, 'Color', colG);
    plot(ax, wlF, yR, '-', 'LineWidth', 0.3, 'Color', colR);

    xlim(ax, [390-10 780+10]); ylim(ax, [0.01 1.02]);
    axis(ax,'off'); ax.Color = 'none';               % no axes, transparent axes

    fig.Units = 'centimeters'; fig.InvertHardcopy = 'off';
    fig.Position = [10,10,iconW,iconH];
    ax.Units = 'centimeters'; ax.Position = [0 0 iconW iconH];   % curves fill the tile

    outName = sprintf('fig2%s_primaries_icon_SD%g.pdf', primLet{k}, sdTarget);
    % vector PDF: reliably transparent AND scalable -> the best form of a symbol.
    exportgraphics(fig, fullfile(figDir,outName), 'ContentType','vector','BackgroundColor','none');
    close(fig);
end

%% ============================================================
% (4) Figure: xy chromaticity (Psychtoolbox background + reflectance gamut polyshape + SD triangles)
%   Upgrade based on your "Figure 8" style:
%     - use plotxyChromaticity (PsychToolbox) for the background
%     - convex hull + polyshape for reflectance gamut (outline)
%     - scatter reflectance chromaticities with nicer styling
%     - keep your SD-dependent display gamut triangles
% ============================================================

fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

% --- Use PsychToolbox chromaticity background (same as your nice script)
% NOTE: plotxyChromaticity creates/uses current axes; keep hold on.
% Draws the spectral-locus horseshoe background; capture that axes as ax.
plotxyChromaticity; hold on;
ax = gca;

% --- Reflectance samples chromaticity under D65 (already computed in your script)
% CIE 1931 xy: x = X/(X+Y+Z), y = Y/(X+Y+Z); eps guard avoids /0 for black.
XYZ_world = XYZ_world_D65;  % 3 x nSamples
sumXYZ = sum(XYZ_world,1); sumXYZ(sumXYZ<=0) = eps;
x_samp = XYZ_world(1,:) ./ sumXYZ;
y_samp = XYZ_world(2,:) ./ sumXYZ;

% --- Convex hull + polyshape gamut (like your example)
% Convex hull of the sample chromaticities = boundary of the reflectance gamut
% (built as a polyshape; the outline plot itself is commented out below).
k_gamut = convhull(x_samp(:), y_samp(:));
C_gamut = polyshape(x_samp(k_gamut), y_samp(k_gamut));

% Draw reflectance gamut outline
% gamutColor = [0.2 0.2 0.2];
% plot(C_gamut, ...
%     'FaceColor','none', ...
%     'EdgeColor', gamutColor, ...
%     'LineWidth', 0.6, ...
%     'LineStyle','-', ...
%     'HandleVisibility','off');

% --- Scatter reflectance points (nicer style)
% Semi-transparent black dots (face alpha 0.20, no edge) show the reflectance
% chromaticity cloud without swamping the background.
sc = scatter(x_samp, y_samp, ...
    14, 'filled', ...
    'MarkerFaceColor',[0 0 0], ...
    'MarkerEdgeColor',[0.2 0.2 0.2], ...
    'LineWidth',0.1, ...
    'HandleVisibility','off');
sc.MarkerFaceAlpha = 0.20;
sc.MarkerEdgeAlpha = 0.00;

% --- D65 white point: plain black cross at the standard D65 chromaticity
% (0.3127, 0.3290); the same black '+' is used for the key.
scatter(0.3127, 0.3290, 42, [0 0 0], '+', 'LineWidth', 1.3, 'HandleVisibility','off');

% --- Display gamut triangles for the analysed bandwidths only (σ = 5, 10, 20 nm)
% One triangle per bandwidth; brewermap 'Reds' ramp shades them by sigma.
sdShow = intersect(sdList, [5 10 20], 'stable');
col = brewermap(max(5, numel(sdShow)+2), 'Reds');
lsList = repmat({'-'}, 1, max(1,numel(sdShow)));

for si = 1:numel(sdShow)
    % Locate this bandwidth in the full sdList to index its primary matrix.
    siG = find(sdList == sdShow(si), 1, 'first');
    Mdisp = Mdisp_set{siG}; % 3 x 3, XYZ of primaries as columns (B,G,R)
    % Convert each primary's XYZ to its xy chromaticity (triangle vertices).
    Xp = Mdisp(1,:); Yp = Mdisp(2,:); Zp = Mdisp(3,:);
    denp = (Xp + Yp + Zp); denp(denp<=0) = eps;

    xP = Xp ./ denp;
    yP = Yp ./ denp;

    % connect B-G-R-B
    % Close the triangle by repeating the first (blue) vertex.
    xTri = [xP(1) xP(2) xP(3) xP(1)];
    yTri = [yP(1) yP(2) yP(3) yP(1)];

    % Draw the display gamut triangle (colour deepens with bandwidth index).
    plot(ax, xTri, yTri, ...
        'LineStyle', lsList{si}, ...
        'LineWidth', 0.8, ...
        'Color', col(min(si+1, size(col,1)),:), ...
        'HandleVisibility','off');
end

% --- Axes limits/ticks (match your nice example)
ax.XLim = [0 0.9];
ax.YLim = [0 0.9];
ax.XTick = [0 0.45 0.9];
ax.YTick = [0 0.45 0.9];
ax.XTickLabel = char('0.00','0.45','0.90');
ax.YTickLabel = char('0.00','0.45','0.90');

xlabel('CIE1931 x', 'FontSize', figp.fontsize, 'FontWeight','normal');
ylabel('CIE1931 y', 'FontSize', figp.fontsize, 'FontWeight','normal');

% --- Figure styling (closer to your published-style script)
fig.Units          = 'centimeters';
fig.Color          = 'w';
fig.InvertHardcopy = 'off';

ax.FontName  = 'Arial';
ax.FontSize  = figp.fontsize;
ax.Color     = ones(3,1);   % white background
ax.TickDir   = 'out';
ax.LineWidth = 0.3;

% Square axes so the xy plane is undistorted (equal x/y scaling).
axis(ax,'square');
box(ax,'off');
grid(ax,'off');

% Size/position similar to your previous two-column/2 style
fig.Position = [10,10,figp.twocolumn/2,figp.twocolumn/2];
ax.Units = 'centimeters';
ax.Position = [1 1 5.2 5.2];   % same "nice" sizing feel

% Export at 600 dpi (the chromaticity background benefits from high resolution).
exportgraphics(fig, fullfile(figDir,'fig2g_gamut.png'), ...
    'ContentType','image','Resolution',600);

fprintf('Fig 2 done.\n');

%% ============================================================
% Helper: XYZ -> Lab (CIE 1976) with reference white
%   Converts tristimulus values to CIELAB (L*a*b*) relative to a given white.
%   XYZ: 3x1 or 3xN
%   XYZ_white: 3x1
%   returns Lab: 3xN
%% ============================================================
function Lab = xyz2lab_wref(XYZ, XYZ_white)

XYZ = double(XYZ);
% Accept Nx3 input by transposing to the 3xN convention.
if size(XYZ,1) ~= 3
    XYZ = XYZ'; % allow Nx3
end
XYZ_white = double(XYZ_white(:));
% Normalise each channel by the white point (X/Xn, Y/Yn, Z/Zn).
XYZn = XYZ ./ max(XYZ_white, eps);

% f(t) per CIE Lab: cube-root above the linear-segment threshold, linear below.
% delta = 6/29 is the break point; kappa is the linear-segment slope.
delta = 6/29;
thr   = delta^3;
kappa = 1/(3*delta^2);   % = (29/6)^2 / 3

f = zeros(size(XYZn));
mask = XYZn > thr;
f(mask)  = XYZn(mask).^(1/3);
f(~mask) = kappa*XYZn(~mask) + 4/29;

fx = f(1,:); fy = f(2,:); fz = f(3,:);

% Standard CIELAB combination of the f() values.
L = 116*fy - 16;
a = 500*(fx - fy);
b = 200*(fy - fz);

Lab = [L; a; b];
end


%% ============================================================
% Helper: CIEDE2000 ΔE00 between two Lab colors
%   Full CIEDE2000 colour-difference formula (Sharma et al. 2005
%   implementation), vectorised over columns.
%   Lab1, Lab2 are 3x1 (or 3xN matched)
%   returns dE00 as scalar (or 1xN)
%% ============================================================
function dE00 = deltaE2000(Lab1, Lab2)

Lab1 = double(Lab1); Lab2 = double(Lab2);
% Enforce the 3xN column convention for both inputs.
if size(Lab1,1) ~= 3, Lab1 = Lab1'; end
if size(Lab2,1) ~= 3, Lab2 = Lab2'; end
assert(all(size(Lab1) == size(Lab2)), 'Lab1 and Lab2 must be same size.');

L1 = Lab1(1,:); a1 = Lab1(2,:); b1 = Lab1(3,:);
L2 = Lab2(1,:); a2 = Lab2(2,:); b2 = Lab2(3,:);

% Chroma of each colour and their mean.
C1 = hypot(a1,b1);
C2 = hypot(a2,b2);
Cbar = (C1 + C2) / 2;

% G scales a* to compensate the formula's known bias in the blue/neutral region.
Cbar7 = Cbar.^7;
G = 0.5 * (1 - sqrt(Cbar7 ./ (Cbar7 + 25^7)));

% Adjusted a' and the resulting chroma C' and hue-angle h' (wrapped to [0,2pi)).
a1p = (1+G).*a1;
a2p = (1+G).*a2;

C1p = hypot(a1p,b1);
C2p = hypot(a2p,b2);

h1p = atan2(b1, a1p);
h2p = atan2(b2, a2p);
h1p(h1p<0) = h1p(h1p<0) + 2*pi;
h2p(h2p<0) = h2p(h2p<0) + 2*pi;

% Lightness and chroma differences.
dLp = L2 - L1;
dCp = C2p - C1p;

% Hue-angle difference, wrapped to (-pi,pi]; zero when either chroma is zero.
dhp = h2p - h1p;
dhp(C1p.*C2p==0) = 0;
dhp(dhp >  pi) = dhp(dhp >  pi) - 2*pi;
dhp(dhp < -pi) = dhp(dhp < -pi) + 2*pi;

% Hue difference expressed in the a*b* plane.
dHp = 2*sqrt(C1p.*C2p) .* sin(dhp/2);

% Mean lightness, chroma and hue used by the weighting functions.
Lbp = (L1 + L2)/2;
Cbp = (C1p + C2p)/2;

hbp = h1p + h2p;
idx = abs(h1p - h2p) > pi;
hbp(idx) = hbp(idx) + 2*pi;
hbp = hbp/2;
hbp(C1p.*C2p==0) = h1p(C1p.*C2p==0) + h2p(C1p.*C2p==0); % doesn't matter; set to sum

% T: hue-dependent weighting term.
T = 1 ...
    - 0.17*cos(hbp - deg2rad(30)) ...
    + 0.24*cos(2*hbp) ...
    + 0.32*cos(3*hbp + deg2rad(6)) ...
    - 0.20*cos(4*hbp - deg2rad(63));

% Rotation term (RT = -sin(2*dTheta)*RC) that couples chroma and hue in blue.
dTheta = deg2rad(30) * exp(-((rad2deg(hbp) - 275)/25).^2);
RC = 2*sqrt((Cbp.^7) ./ (Cbp.^7 + 25^7));

% Weighting (compensation) functions for lightness, chroma and hue.
SL = 1 + (0.015*(Lbp - 50).^2) ./ sqrt(20 + (Lbp - 50).^2);
SC = 1 + 0.045*Cbp;
SH = 1 + 0.015*Cbp.*T;

RT = -sin(2*dTheta) .* RC;

% Parametric factors (unity = reference viewing conditions).
kL = 1; kC = 1; kH = 1;

% Combine weighted lightness/chroma/hue differences plus the rotation term.
dE00 = sqrt( ...
    (dLp./(kL*SL)).^2 + ...
    (dCp./(kC*SC)).^2 + ...
    (dHp./(kH*SH)).^2 + ...
    RT .* (dCp./(kC*SC)) .* (dHp./(kH*SH)) );
end


%% ============================================================
% Helper: render a filled spectrum into a square tile (no axes/labels)
%   Rasterises a single filled spectral curve to a tilePx-by-tilePx RGB image,
%   e.g. for use as a small inline glyph. Returns the pixel array (uint8).
%   - wls: nW x 1
%   - y  : nW x 1 (normalized 0..1)
%   - rgbFill: 1 x 3 in [0,1]
%   - tilePx: scalar
%% ============================================================
function tile = renderFilledSpectrumTile(wls, y, rgbFill, tilePx)

wls = wls(:);
y   = y(:);
% Clamp curve to [0,1] so the fill stays inside the tile.
y   = max(min(y,1),0);

% create a hidden figure sized to tilePx x tilePx pixels
fig = figure('Visible','off','Color','w','Units','pixels','Position',[10 10 tilePx tilePx]);
ax  = axes(fig); hold(ax,'on');

% fill under curve (graph-only)
% Closed polygon [wls; reversed wls] over [0; reversed y] fills the area
% between the curve and the baseline in the solid fill colour.
fill(ax, [wls; flipud(wls)], [zeros(size(y)); flipud(y)], ...
    rgbFill, 'EdgeColor','none', 'FaceAlpha', 1.0);

% subtle outline for definition (same idea as your reflected spectrum plot)
plot(ax, wls, y, 'k-', 'LineWidth', 0.8);

xlim(ax, [min(wls) max(wls)]);
ylim(ax, [0.01 1]);

% Axis-free, filling the whole tile so the curve occupies the full image.
axis(ax,'off');
ax.Position = [0 0 1 1];  % fill entire tile

% capture to image
% getframe/frame2im grabs the rendered pixels as an RGB array.
frame = getframe(fig);
tile  = frame2im(frame);

% ensure exact size (sometimes off by 1 px depending on renderer)
tile = imresize(tile, [tilePx tilePx]);

close(fig);
end



%% ============================================================
% Local helper: filled graph-only plot in the "illuminant style"
%   Draws a single filled spectrum with no axes and exports it as a
%   transparent-background vector PDF icon (outFile) into figDir.
%% ============================================================
function makeFilledGraphOnly(wls, y, rgbFill, figp, figDir, outFile)
nW = numel(wls);
y = y(:);
y = max(y,0);

fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

% Filled area under the curve plus a thin black outline (same idiom as above).
fill(ax, [wls(:); flipud(wls(:))], [zeros(nW,1); flipud(y)], ...
    rgbFill, 'EdgeColor','none', 'FaceAlpha', 1.0);
plot(ax, wls, y, 'k-', 'LineWidth', 0.5);

xlim(ax, [390-10 780+10]);
ylim(ax, [0.1 1.0]);

fig.Units          = 'centimeters';
fig.Color          = 'w';
fig.InvertHardcopy = 'off';

ax.FontName  = 'Arial';
ax.FontSize  = figp.fontsize;
ax.XColor    = 'k';
ax.YColor    = 'k';
ax.LineWidth = 0.5;
ax.Units     = 'centimeters';

fig.Position = [10, 10, figp.twocolumn/3, figp.twocolumn/3];
ax.Position  = [0.2 0.2, figp.twocolumn/3*0.9, figp.twocolumn/3*0.9];
axis(ax,'off');

% BackgroundColor 'none' -> transparent PDF suitable as a figure symbol.
exportgraphics(fig, fullfile(figDir, outFile), ...
    'Resolution',300,'BackgroundColor','none');
end


%% ============================================================
% Helper: build mosaic image for sRGB patches
%   Tiles N sRGB colours (rgb: N x 3) into an nRows-by-nCols grid of colour
%   patches separated by gapPx-wide white gutters; returns a uint8 RGB image.
%% ============================================================
function img = make_srgb_mosaic(rgb, nCols, nRows, tilePx, gapPx)
N = size(rgb,1);
% Total canvas size = tiles plus the gaps around and between them.
W = nCols*tilePx + (nCols+1)*gapPx;
H = nRows*tilePx + (nRows+1)*gapPx;

img = ones(H, W, 3); % white background

% Fill the grid row-major; each cell's pixel block is painted its colour.
k = 1;
for r = 1:nRows
    for c = 1:nCols
        % Pixel bounds of tile (r,c), offset by the gutters that precede it.
        y0 = (r-1)*tilePx + r*gapPx + 1;
        y1 = y0 + tilePx - 1;
        x0 = (c-1)*tilePx + c*gapPx + 1;
        x1 = x0 + tilePx - 1;

        % Paint only while colours remain (trailing cells stay white).
        if k <= N
            col = reshape(rgb(k,:), 1,1,3);
            img(y0:y1, x0:x1, :) = repmat(col, tilePx, tilePx, 1);
        end
        k = k + 1;
    end
end

img = im2uint8(img);
end

%% ============================================================
% Helper: XYZ -> sRGB (D65) using a reference white (your pipeline style)
%   Converts tristimulus values to display sRGB for on-screen colour rendering
%   of the reflectance samples: normalise by the white, apply the sRGB matrix,
%   clip to gamut, then apply the sRGB transfer function (gamma).
%% ============================================================
function rgb = xyz_to_srgb_d65_wref(XYZ, XYZ_white)
% XYZ: 3 x N (absolute scale arbitrary)
% XYZ_white: 3 x 1 reference white under same scaling
% rgb: N x 3 in [0,1]

% Scale so the reference white maps to relative luminance 1 (Y = 1).
Yw = max(XYZ_white(2), eps);
XYZn = XYZ ./ Yw;

% Standard linear-sRGB (D65) XYZ->RGB matrix.
M = [ 3.2406 -1.5372 -0.4986;
     -0.9689  1.8758  0.0415;
      0.0557 -0.2040  1.0570];

% Linear RGB, then clip out-of-gamut / negative values into [0,1].
rgb_lin = (M * XYZn)';         % N x 3
rgb_lin = max(min(rgb_lin,1),0);

% sRGB opto-electronic transfer function: linear near black, ~gamma 2.4 above.
a = 0.055;
thr = 0.0031308;

rgb = zeros(size(rgb_lin));
low = rgb_lin <= thr;
rgb(low)  = 12.92 * rgb_lin(low);
rgb(~low) = (1+a) * (rgb_lin(~low).^(1/2.4)) - a;

% Final safety clip to displayable [0,1].
rgb = max(min(rgb,1),0);
end

%% ============================================================
% Helper: resample (S, T) to your wls vector (for standard observer cones)
%   Reinterpolates a spectral table stored in Psychtoolbox "S" format onto wls.
%   S: [wl0 step n]  (Psychtoolbox style)
%   T: 3 x n OR n x 3
% ============================================================
function T_out = resample_T_to_wls(S, T_in, wls)
% Accept n x 3 tables by transposing to the 3 x n convention.
if size(T_in,1) ~= 3 && size(T_in,2) == 3
    T_in = T_in';
end
% Reconstruct the source wavelength axis from [start, step, count].
wl0 = S(1); step = S(2); n = S(3);
wl = (wl0:step:(wl0 + step*(n-1)))';
% Interpolate each row onto wls (pchip, zero outside the source range).
T_out = interp1(wl, T_in', wls(:), 'pchip', 0)'; % 3 x numel(wls)
end