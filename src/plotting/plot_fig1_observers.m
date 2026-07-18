% ========================================================================
% plot_fig1_observers.m  -  Figure 1: individual observer cone fundamentals
% Self-contained; run on its own:  >> plot_fig1_observers
% Loads only small fields (no 9 GB byIll) so it is fast.
%
% MANUSCRIPT FIGURE
%   Optics Express: "Adaptation reduces inter-observer metamerism with
%   narrow-band primaries." Figure 1 shows the observer model: the family
%   of individual-observer cone fundamentals produced by the Asano model on
%   the CIE 2006 basis (10,000 simulated observers, each with discrete
%   gene-determined lambda-max). Curves are peak-normalised L/M/S sensitivities.
%
% INPUT .mat
%   results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat
%   Only the small fields are loaded (wls, params, display, keepIdx_global,
%   illuminants, cone_fund); the multi-gigabyte per-illuminant field is skipped.
%
% OUTPUT FILES (written into figs/fig1/)
%   fig1a_cone_fundamentals_family.png  - subsampled family (every 10th observer)
%                                         plus bold median L/M/S, WITH axes/ticks.
%   fig1b_cone_fundamentals_all.png     - ALL 10,000 observers, graph-only
%                                         (no axes) plus bold median, high-res PNG.
%
% PANELS
%   (a) Cone-fundamental family with labelled axes: thin L/M/S curves for a
%       1-in-10 subsample of observers overlaid with the population median.
%   (b) Same family drawn for all 10,000 observers as an axis-free tile, for
%       use as a compact illustration of inter-observer spread.
%
% NOTE: local helper functions live at the bottom of this file. Several
%   (xyz2lab_wref, deltaE2000, renderFilledSpectrumTile, makeFilledGraphOnly,
%   make_srgb_mosaic, xyz_to_srgb_d65_wref) are shared colorimetry/plotting
%   utilities kept alongside the figure scripts; only resample_T_to_wls and
%   xyz_to_srgb_d65_wref are exercised by this particular figure. Each helper
%   is documented at its definition below.
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
if isempty(subFig), subFig = {'fig1'}; end                 % fallback for section runs
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
% Mdisp_set: per-bandwidth 3x3 matrix of primary tristimulus values (XYZ).
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
% Adopted white = the brightest sample (max luminance Y); used later as the
% reference white for sRGB rendering of sample chromaticities.
[~,idxWhite]  = max(XYZ_world_D65(2,:));
XYZ_white_D65 = XYZ_world_D65(:,idxWhite);
% Guard against non-positive components before any division/normalisation.
XYZ_white_D65(XYZ_white_D65<=0) = eps;

% ---- standard-observer cone fundamentals (Stockman-Sharpe 2-deg) ----
% Reference LMS fundamentals, resampled from their native grid onto wls.
Lstd = load(fullfile(dataDir,'T_cones_ss2.mat'));
T_cones_std = resample_T_to_wls(Lstd.S_cones_ss2,Lstd.T_cones_ss2,wls);

%% ============================================================
% (5) Figure: Cone fundamentals across observers (Asano family)
%% ============================================================
assert(isfield(S,'cone_fund'), 'S.cone_fund not found in MAT file.');
cone_fund = S.cone_fund;
assert(size(cone_fund,2) == numel(wls), 'cone_fund wavelength dim != numel(wls).');
assert(size(cone_fund,3) == 3, 'cone_fund third dim must be 3 (L,M,S).');

Nobs = size(cone_fund,1);

% Peak-normalise each observer's each cone to a maximum of 1 (max over the
% wavelength dimension, dim 2), so curves are compared by shape not absolute
% sensitivity. eps guard avoids divide-by-zero for any all-zero curve.
cf = cone_fund;
mx = max(cf, [], 2);
mx(mx<=0) = eps;
cfN = cf ./ mx;

% Population median curve across observers (dim 1) -> [nW x 3] for L,M,S.
cfMed = squeeze(median(cfN, 1));      % [nW x 3]

% Cone colour convention: L red, M green, S blue.
colL = [0.85 0.10 0.10];
colM = [0.10 0.60 0.10];
colS = [0.10 0.30 0.85];

fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

% Draw a 1-in-10 subsample of observers as thin curves (transpose so each
% column-per-observer is plotted as its own line against wls).
step = 10; idx = 1:step:Nobs;
plot(ax, wls, squeeze(cfN(idx,:,1))', 'LineWidth', 0.15, 'Color', colL);
plot(ax, wls, squeeze(cfN(idx,:,2))', 'LineWidth', 0.15, 'Color', colM);
plot(ax, wls, squeeze(cfN(idx,:,3))', 'LineWidth', 0.15, 'Color', colS);

% Overlay the population median for each cone as a thicker line.
plot(ax, wls, cfMed(:,1), '-', 'LineWidth', 0.8, 'Color', colL);
plot(ax, wls, cfMed(:,2), '-', 'LineWidth', 0.8, 'Color', colM);
plot(ax, wls, cfMed(:,3), '-', 'LineWidth', 0.8, 'Color', colS);

xlabel(ax,'Wavelength (nm)','FontWeight','Bold');
ylabel(ax,'Normalized sensitivity','FontWeight','Bold');

% Pad the visible-range x-limits by 10 nm each side; y just above 1 (peak-norm).
xlim([390-10 780+10]); ylim([0 1.05]);

% InvertHardcopy off keeps the white figure background on export (no auto-white).
fig.Units          = 'centimeters';
fig.Color          = 'w';
fig.InvertHardcopy = 'off';

ax.FontName  = 'Arial';
ax.FontSize  = figp.fontsize;
ax.XColor    = 'k';
ax.YColor    = 'k';
ax.LineWidth = 0.5;
ax.Units     = 'centimeters';

% Physical figure/axes size in cm, keyed to the journal two-column width.
fig.Position = [10,10,figp.twocolumn/3*1.2,figp.twocolumn/3*1.2];
ax.Position  = [0.80 0.74, figp.twocolumn/3*0.97, figp.twocolumn/3*0.97];

xticks([400 500 600 700 780]);
yticks([0 0.5 1.0]);
ax.XTickLabel = {'400','500','600','700','780'};
ax.YTickLabel = {'0.0','0.5','1.0'};

% ticklengthcm: helper that sets tick length in cm (renderer-independent).
ticklengthcm(ax,0.05);
grid(ax,'off'); box(ax,'off');
% Light-grey axes background so the many thin curves stay legible.
ax.Color = ones(1,3)*0.97;

% Export panel (a) as a 300-dpi raster PNG.
exportgraphics(fig, fullfile(figDir,'fig1a_cone_fundamentals_family.png'), ...
    'ContentType','image','Resolution',300);

%% ============================================================
% (5b) Cone fundamentals (ALL 10,000 observers) graph-only + bold median
%   - no axis / ticks / labels (like illuminants)
%   - WARNING: vector PDF with 10k curves can be huge; PNG is safer.
%% ============================================================
fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

% cfN is [Nobs x nW x 3], already normalized peak=1 per observer & cone
Nobs = size(cfN,1);
assert(Nobs == 10000, 'Expected 10000 observers, got %d.', Nobs);

% --- draw ALL observers (thin) ---
% One line per observer (columns of the transposed slice); very thin so the
% 10,000-curve overlay reads as a density band rather than solid fill.
plot(ax, wls, squeeze(cfN(:,:,1))', 'LineWidth', 0.05, 'Color', colL); % L
plot(ax, wls, squeeze(cfN(:,:,2))', 'LineWidth', 0.05, 'Color', colM); % M
plot(ax, wls, squeeze(cfN(:,:,3))', 'LineWidth', 0.05, 'Color', colS); % S

% --- bold median on top ---
plot(ax, wls, cfMed(:,1), '-', 'LineWidth', 1.2, 'Color', colL);
plot(ax, wls, cfMed(:,2), '-', 'LineWidth', 1.2, 'Color', colM);
plot(ax, wls, cfMed(:,3), '-', 'LineWidth', 1.2, 'Color', colS);

xlim(ax, [390-10 780+10]);
ylim(ax, [0 1.05]);

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
% Graph-only tile: strip axes, ticks and labels (curves only).
axis(ax,'off');

% --- EXPORT: prefer raster to avoid massive vector files ---
% 600-dpi PNG: a vector export of 10k curves would be enormous/slow to open.
exportgraphics(fig, fullfile(figDir,'fig1b_cone_fundamentals_all.png'), ...
    'Resolution', 600);

% If you *really* want vector (can be enormous / slow):
% exportgraphics(fig, fullfile(figDir,'fig1b_cone_fundamentals_all.png'), ...
%     'Resolution',300,'BackgroundColor','none');



fprintf('Fig 1 done.\n');

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