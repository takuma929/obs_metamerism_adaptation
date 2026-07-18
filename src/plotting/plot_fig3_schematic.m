% ========================================================================
% MANUSCRIPT FIGURE 3 (Optics Express, "Adaptation reduces inter-observer
%   metamerism with narrow-band primaries").  NOTE: the script "fig3" label
%   here happens to match manuscript Figure 3, but in general the figure
%   SCRIPT numbers do NOT track the renumbered manuscript figures.
%
% WHAT THIS FIGURE SHOWS - a schematic of the metamerism-quantification
%   pipeline for ONE fixed example surface (origIdx = 83):
%     (a) the "world" spectrum = D65 illuminant x surface reflectance;
%     (b) the surface reflectance on its own;
%     (c) the D65 illuminant on its own;
%     (d) for a few representative individual observers PLUS the standard
%         observer, the reproduced display spectrum (sum of 3 Gaussian
%         RGB primaries).  The RGB weights are solved so the reproduced
%         stimulus matches the observer's own cone LMS; the standard
%         observer then reads that stimulus and we report its LMS, XYZ and
%         the CAM16-UCS delta-E' vs the target.
%   Tristimulus values are also printed to the Command Window for the paper.
%
% INPUT  .mat : results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat
%               (only small fields loaded: wls, params, display, keepIdx,
%                illuminants, cone_fund - never the ~9 GB per-illuminant array)
% OTHER data  : data/fig_parameters.mat (figure style), CIE_srf_cfi_1nm.csv
%               (99 CIE test reflectances), CIE1931XYZ.txt (2-deg CMFs),
%               T_cones_ss2.mat (Stockman-Sharpe 2-deg standard cones)
% OUTPUTS     : vector PDFs into figs/fig3/  (fig3a_..., fig3b_..., fig3c_...,
%               fig3_rgbRepro_...) plus Command-Window tables.
%
% Reproducibility: the example surface and observer indices are PINNED, with a
% seeded rng(0) fallback, so the schematic is identical across runs.
% LOCAL HELPERS at the bottom (xyz2lab_wref, deltaE2000, cam16_s3, the sRGB /
% tile / mosaic renderers, resample_T_to_wls) are commented in place.
% ------------------------------------------------------------------------
% ========================================================================
% plot_fig3_schematic.m  -  Figure 3: metamerism-quantification schematic
%   (a) reflected-light spectrum, (b) surface reflectance, (c) D65 illuminant
%   spectrum, and the per-observer reproduced RGB spectra (2:1 panels), origIdx=83.
% Self-contained; run on its own:  >> plot_fig3_schematic
% Loads only small fields (no 9 GB byIll) so it is fast.
% ========================================================================
clearvars; close all; clc;

% Resolve project root from this script's own location so it can be run
% from any working directory (e.g. directly from src/plotting/).
% mfilename('fullpath') is .../src/plotting/plot_fig3_schematic, so three
% fileparts() calls strip the file + the plotting/ + the src/ folder.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist(fullfile(projectRoot,'src'),'dir')               % section run: mfilename points at a temp file
    projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
end
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
% Derive the output subfolder name ("fig3") from this script's own name so a
% renamed copy writes to the matching figs/<name>/ directory automatically.
[~,thisName] = fileparts(mfilename('fullpath'));            % e.g. plot_fig1_observers
subFig  = regexp(thisName,'^plot_(fig[^_]*)','tokens','once'); % -> 'fig1' ... 'figS'
if isempty(subFig), subFig = {'fig3'}; end                 % fallback for section runs
figDir  = fullfile(projectRoot,'figs',subFig{1});     % per-figure subfolder
if ~exist(figDir,'dir'); mkdir(figDir); end
inMat=fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
assert(exist(inMat,'file')==2,'Missing results: %s.',inMat);
% ---- load results (small fields only; no 9 GB byIll) ----
S = load(inMat,'wls','params','display','keepIdx_global','illuminants','cone_fund');
% wls: sampling wavelengths (column); sdList: primary bandwidths sigma (nm).
wls       = S.wls(:);
sdList    = S.params.sdList(:)';
nW        = numel(wls);
% Pset{si}: nW x 3 Gaussian primary SPDs for the si-th bandwidth; columns [B G R].
% Mdisp_set: matching display-to-cone matrices (loaded but not used in this figure).
Pset      = S.display.Pset;
Mdisp_set = S.display.Mdisp_set;
% keepIdx: original 1..99 indices of the in-gamut CIE test surfaces kept.
keepIdx   = S.keepIdx_global(:);
% cone_fund: [Nobs x nW x 3] individual Asano CIE 2006 cone fundamentals, [L M S].
cone_fund = S.cone_fund;

% ---- figure style parameters ----
T = load(fullfile(dataDir,'fig_parameters.mat'));
figp = T.figp;

% ---- reflectance spectra (99 CIE test samples, resampled to wls) ----
% First CSV column is wavelength; remaining 99 columns are reflectances at 1 nm.
Rmat     = readmatrix(fullfile(dataDir,'CIE_srf_cfi_1nm.csv'));
wl_ref   = Rmat(:,1);
R_99_1nm = Rmat(:,2:end);
% Resample all 99 reflectances onto wls and clamp to the physical [0,1] range.
R_wlsx99 = min(max(interp1(wl_ref,R_99_1nm,wls,'pchip','extrap'),0),1);
% Restrict to the in-gamut surfaces kept for this observer/display run.
R_keep   = R_wlsx99(:,keepIdx);

% ---- CIE 1931 2-deg colour-matching functions ----
% Reference (standard) observer's CMFs; transpose to 3 x nW so T_xyz*spd -> XYZ.
Axyz  = readmatrix(fullfile(dataDir,'CIE1931XYZ.txt'));
T_xyz = interp1(Axyz(:,1),Axyz(:,2:4),wls,'pchip',0)';

% ---- representative illuminants (daylight D65 and peaky) ----
% Pick the flagged representative daylight (D65) and peaky illuminant SPDs.
illCats   = string({S.illuminants.category});
isRep     = [S.illuminants.isRep];
idxD65    = find(isRep & illCats=="daylight",1,'first');
d65       = S.illuminants(idxD65).spd(:);
idxHouser = find(isRep & illCats=="peaky",1,'first');
% (peaky SPD grabbed for completeness; this schematic only uses D65 below)
if ~isempty(idxHouser), spdHouser = S.illuminants(idxHouser).spd(:); else, spdHouser = []; end

% ---- world colours under D65 and the adopted (brightest) white ----
% World spectra = illuminant x reflectance (column-wise); XYZ via the CMFs.
nSamples      = size(R_keep,2);
SPD_keep_D65  = d65 .* R_keep;
XYZ_world_D65 = T_xyz * SPD_keep_D65;
% Adopt the highest-luminance (max Y) surface as the scene white for the
% sRGB rendering and the CAM16-UCS adopting white below.
[~,idxWhite]  = max(XYZ_world_D65(2,:));
XYZ_white_D65 = XYZ_world_D65(:,idxWhite);
XYZ_white_D65(XYZ_white_D65<=0) = eps;

% ---- standard-observer cone fundamentals (Stockman-Sharpe 2-deg) ----
% Standard observer's LMS cones (used to read reproduced spectra); resampled
% from Psychtoolbox [wl0 step n] format onto wls (3 x nW).
Lstd = load(fullfile(dataDir,'T_cones_ss2.mat'));
T_cones_std = resample_T_to_wls(Lstd.S_cones_ss2,Lstd.T_cones_ss2,wls);

%% ============================================================
% (6b) Reflected-light spectrum of a RANDOM kept surface under D65
%   - graph-only, fill with sRGB render color
%   - print standard observer LMS (normalized max=1) to Command Window
%% ============================================================
% Fixed target surface to reproduce the original figure (origIdx = 83).
% Pinning the surface avoids the random pick shifting whenever the kept
% set changes (a random draw over a different-sized kept set lands elsewhere).
TARGET_ORIGIDX = 83;
j = find(keepIdx == TARGET_ORIGIDX, 1);
if isempty(j)                    % fallback if that surface is out-of-gamut this run
    rng(0); j = randi(nSamples);
end
origIdx = keepIdx(j);            % original index in 1..99

% World spectrum for this example surface = D65 x reflectance; the "_n"
% version is peak-normalised to 1 for plotting (shape only).
spd_rand = d65 .* R_keep(:,j);  % nW x 1 (relative radiance)
spd_rand_n = spd_rand ./ max(spd_rand + eps);

% Fill color: sRGB rendering of this surface under D65 (using your reference white)
XYZ_rand = T_xyz * spd_rand;                     % 3 x 1
rgbFill  = xyz_to_srgb_d65_wref(XYZ_rand, XYZ_white_D65); % 1 x 3

% Standard observer LMS (normalize max=1)
LMS = T_cones_std * spd_rand;        % 3 x 1
LMSn = LMS ./ max(LMS + eps);

fprintf('\n[RANDOM SURFACE under D65]\n');
fprintf('  kept-set index j = %d / %d\n', j, nSamples);
fprintf('  original reflectance index (1..99) = %d\n', origIdx);
fprintf('  standard-observer LMS (max-normalized): [L M S] = [%.4f  %.4f  %.4f]\n', ...
    LMSn(1), LMSn(2), LMSn(3));

fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

% Filled area under the normalised spectrum, coloured by the surface's own
% sRGB render; thin black outline for definition. The [wls; flipud(wls)] /
% [0; flipud(y)] pattern builds a closed polygon for fill().
fill(ax, [wls; flipud(wls)], [zeros(nW,1); flipud(spd_rand_n)], ...
    rgbFill, 'EdgeColor','none', 'FaceAlpha', 1.0);
plot(ax, wls, spd_rand_n, 'k-', 'LineWidth', 0.8); % subtle outline

% Visible-range window (with a small margin) and full 0..1 height.
xlim(ax, [390-10 780+10]);
ylim(ax, [0 1.0]);

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

% Panel (a): world spectrum PDF, filename carries the example surface index.
% Transparent-background vector export keeps the axis-free tile crisp.
outName = sprintf('fig3a_target_spectrum_idx%02d.pdf', origIdx);
exportgraphics(fig, fullfile(figDir, outName), ...
    'ContentType','vector','BackgroundColor','none');

%% ============================================================
% (6c) ALSO plot SURFACE REFLECTANCE of the chosen random sample
%   - graph-only (no axis/labels)
%   - filled under curve using SAME sRGB color as reflected light (rgbFill)
%% ============================================================
% Panel (b): the surface reflectance alone (illuminant removed), same square
% tile style. Filled with the SAME sRGB colour as (a) so the two panels read
% as the same physical surface.
R_rand = R_keep(:,j);                 % reflectance (0..1)
R_rand = min(max(R_rand,0),1);

fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

fill(ax, [wls(:); flipud(wls(:))], [zeros(nW,1); flipud(R_rand(:))], ...
    rgbFill, 'EdgeColor','none', 'FaceAlpha', 1.0);
plot(ax, wls, R_rand, 'k-', 'LineWidth', 0.8);  % subtle outline

xlim(ax, [390-10 780+10]);
ylim(ax, [0 1.0]);

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

outName = sprintf('fig3b_target_reflectance_idx%02d.pdf', origIdx);
exportgraphics(fig, fullfile(figDir, outName), ...
    'ContentType','vector','BackgroundColor','none');


%% ============================================================
% (6c-2) Illuminant spectrum (D65) that lights the surface
%   - graph-only (no axis/labels), same square style as (6b)/(6c)
%   - filled with a pale-yellow "light source" colour
%   - illuminant is surface-independent, so no origIdx in the name
%% ============================================================
% Panel (c): the D65 illuminant SPD alone, peak-normalised. Surface-independent,
% so the filename has no origIdx. Filled pale grey to read as a light source.
y_ill    = d65 ./ max(d65 + eps);      % normalised SPD (shape only)
illFill  = [0.9 0.9 0.9];              % light grey = light source

fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

fill(ax, [wls(:); flipud(wls(:))], [zeros(nW,1); flipud(y_ill(:))], ...
    illFill, 'EdgeColor','none', 'FaceAlpha', 1.0);
plot(ax, wls, y_ill, 'k-', 'LineWidth', 0.8);  % subtle outline

xlim(ax, [390-10 780+10]);
ylim(ax, [0 1.0]);

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

exportgraphics(fig, fullfile(figDir, 'fig3c_illuminant_D65.pdf'), ...
    'ContentType','vector','BackgroundColor','none');


%% ============================================================
% (6d) RGB reproduction (SD targets) for:
%        - 3 requested observers  +  standard observer
%      - filled spectra (same style as reflected light)
%      - outputs SUM only (R+G+B)
%      - prints:
%           * LMS(original) and LMS(reproduced) for each observer
%           * XYZ1931(target) and XYZ1931(reproduced) (STANDARD OBSERVER)
%           * ΔE00 (STANDARD OBSERVER Lab, relative to XYZ_white_D65)
%
% REQUIREMENTS already in your script:
%   - spd_rand (reflected-light spectrum, nW x 1) and origIdx exist
%   - Pset, sdList, wls, nW, figp, figDir exist
%   - T_xyz (3 x nW), XYZ_white_D65 (3 x 1) exist
%   - xyz_to_srgb_d65_wref, makeFilledGraphOnly exist
%   - standard observer cones already loaded as:
%       T_cones_std  (3 x nW)
%% ============================================================

sdTargets = [5 10 20];   % primary bandwidths shown in the worked example

% CAM16-UCS parameters (match the paper's metric, §2.3) for the ΔE′ column.
%   LA = adapting luminance, Yb = background luminance factor, F/c/Nc = the
%   CIECAM16 "average" surround constants. M16cam is the XYZ->RGB (CAT16)
%   sharpened-cone matrix used inside the appearance model.
camCFG = struct('LA',20,'Yb',20,'F',1.0,'c',0.69,'Nc',1.0);
M16cam = [ 0.401288 0.650173 -0.051461; -0.250268 1.204414 0.045854; -0.002079 0.048952 0.953127 ];
% Scale the adopting white to Y = 100 (paper convention) and reuse that same
% scale for every stimulus so the CAM16 lightness reference is consistent.
wsc    = 100 / max(XYZ_white_D65(2), eps);
Wref_c = XYZ_white_D65(:) * wsc;

% ---- observer cone fundamentals (raw) ----
assert(isfield(S,'cone_fund'), 'S.cone_fund missing.');
cone_fund = S.cone_fund;                 % [Nobs x nW x 3], order [L M S]
Nobs = size(cone_fund,1);
assert(size(cone_fund,2)==nW && size(cone_fund,3)==3, 'cone_fund size mismatch.');

% ---- observers shown in the figure (three requested individuals) ----
% A greedy max-min sweep in cone space originally picked five observers; the
% paper figure keeps just these three observer indices.
sel = [8255; 2438; 119];
assert(all(sel>=1 & sel<=Nobs), 'requested observer index out of range (Nobs=%d).', Nobs);

fprintf('\n[RGB REPRODUCTION for SD targets]\n');
fprintf('Selected observer indices (%d): %s\n', numel(sel), mat2str(sel'));

spd_target = spd_rand(:);  % nW x 1 (reflected-light spectrum)

% ============================================================
% STANDARD OBSERVER reference metrics for the TARGET (world) spectrum.
%   These are the "ground truth" values the reproduced stimuli are compared
%   against: XYZ under CIE 1931, CIELAB, and the CAM16-UCS [J';a';b'] vector.
% ============================================================
XYZ_tgt_std = T_xyz * spd_target;                         % 3x1
Lab_tgt_std = xyz2lab_wref(XYZ_tgt_std(:), XYZ_white_D65(:));
% D = 1 -> full chromatic adaptation to the scene white (von Kries-style).
s3_tgt_std  = cam16_s3(XYZ_tgt_std(:)*wsc, Wref_c, 1, M16cam, camCFG);  % [J';a';b'] CAM16-UCS

fprintf('\n[Standard-observer target metrics]\n');
fprintf('  XYZ1931(target) = [%.4f  %.4f  %.4f]\n', XYZ_tgt_std(1), XYZ_tgt_std(2), XYZ_tgt_std(3));
fprintf('  Lab*(target)    = [%.4f  %.4f  %.4f]  (ref white = XYZ_white_D65)\n', Lab_tgt_std(1), Lab_tgt_std(2), Lab_tgt_std(3));

% ============================================================
% Build a list of "observers to test":
%   1..5: selected individuals
%   6:    standard observer (T_cones_std)
% ============================================================
nSel = numel(sel);

% collect per-observer metrics, then print a transposed table (observers as
% columns, metrics as rows) to match the Fig 3c layout.
store = struct();

% Loop over the selected individuals (kk = 1..nSel) plus one extra pass
% (kk = nSel+1) for the standard observer.
for kk = 1:(nSel+1)

    if kk <= nSel
        o = sel(kk);
        % To: the observer's colour-matching space used to SOLVE the RGB mix.
        % For an individual that is their own CIE 2006 cone fundamentals.
        To = squeeze(cone_fund(o,:,:))';     % 3 x nW  (Asano CIE 2006 individual cones)
        To_lms = To;                         % show this individual's own cone LMS
        obsLabel = sprintf('observer%05d', o);
    else
        % Standard observer: solve the mix in CIE 1931 XYZ so the reproduced
        % stimulus is an exact XYZ match -> its ΔE00/ΔE' vs target is ~0.
        To = T_xyz;                          % solve & ΔE00 reference = CIE 1931 (keeps ΔE00 = 0)
        To_lms = T_cones_std;                % show its LMS via the CIE 2006 / SS 2-deg standard cones
        obsLabel = 'standardObserver';
    end

    % LMS_tgt drives the RGB solve (in the observer's matching space);
    % LMS_ind is the cone LMS shown in the table (Asano individual / standard cones, never XYZ).
    LMS_tgt   = To * spd_target;             % 3 x 1  (for the RGB solve)
    LMS_ind   = To_lms * spd_target;         % 3 x 1  (cone LMS for display)
    LMS_ind_n = LMS_ind ./ max(LMS_ind + eps);

    % One reproduction per primary bandwidth sigma.
    for ii = 1:numel(sdTargets)
        sd = sdTargets(ii);

        si = find(sdList == sd, 1, 'first');
        assert(~isempty(si), 'SD=%g not found in sdList.', sd);

        P = Pset{si};   % nW x 3, columns assumed [B G R]

        % Metameric match: find RGB primary weights whose mixture excites THIS
        % observer's cones exactly like the target. A = To*P is the 3x3
        % display-to-observer matrix; solving A\LMS_tgt inverts it.
        % Solve (To*P)*rgb = LMS_tgt  (3x3 system)
        A = To * P;                 % 3 x 3
        rgb_w = A \ LMS_tgt;        % 3 x 1 weights for [B G R]

        % Reproduced display spectrum = weighted sum of the 3 primaries.
        % Reproduced spectrum (sum only)
        spd_sum = P * rgb_w;        % nW x 1
        spd_sum = max(spd_sum, 0);

        % Now read that reproduced stimulus with the STANDARD observer: this is
        % the mismatch an individual's metamer creates for everyone else.
        % Standard-observer LMS and XYZ1931 of the reproduced spectrum
        LMS_std     = T_cones_std * spd_sum;                   % 3x1
        LMS_std_n   = LMS_std ./ max(LMS_std + eps);
        % For the standard observer, individual and Std LMS differ only by the
        % (tiny) metamer residual; align them to the individual LMS to avoid confusion.
        if kk > nSel, LMS_std_n = LMS_ind_n; end
        XYZ_rep_std = T_xyz * spd_sum;                         % 3x1
        % CAM16-UCS ΔE′ between the reproduced stimulus and the target, both read
        % by the standard observer and adapted to the scene white (D = 1).
        % Euclidean distance in [J';a';b'] IS the CAM16-UCS colour difference.
        s3_rep_std  = cam16_s3(XYZ_rep_std(:)*wsc, Wref_c, 1, M16cam, camCFG);
        dEp_std     = norm(s3_rep_std(:) - s3_tgt_std(:));

        % stash metrics for the transposed table. rgb_w is in [B G R] order
        % (matching Pset columns), so flip to R,G,B for display.
        store.(obsLabel).LMS_ind(:,ii) = LMS_ind_n(:);
        store.(obsLabel).rgb(:,ii)     = [rgb_w(3); rgb_w(2); rgb_w(1)];
        store.(obsLabel).LMS_std(:,ii) = LMS_std_n(:);
        store.(obsLabel).XYZ(:,ii)     = XYZ_rep_std(:);
        store.(obsLabel).dE(ii)        = dEp_std;

        % ---- fill color for spectrum = its sRGB rendering (same style) ----
        rgbFill_sum = xyz_to_srgb_d65_wref(XYZ_rep_std, XYZ_white_D65); % 1x3

        % ---- normalize for plotting (shape only) ----
        ySum = spd_sum ./ max(spd_sum + eps);

        % ---- OUTPUT individual plot (PDF): one 2:1 panel per observer x sigma
        outPdf = sprintf('fig3_rgbRepro_SD%02d_%s_idx%02d.pdf', sd, obsLabel, origIdx);
        makeFilledGraphOnly(wls, ySum, rgbFill_sum, figp, figDir, outPdf);
    end
end

% ---- transposed table (columns = observers, rows = metrics; matches Fig 3c) ----
% Column order pins the standard observer first, then the three individuals.
colKeys  = {'standardObserver','observer00119','observer02438','observer08255'};
colNames = {'std observer','00119','02438','08255'};
% Anonymous formatter: render a 3-vector with a run-time-chosen decimal count.
fmtV = @(v,dec) sprintf(['[%.' dec 'f %.' dec 'f %.' dec 'f]'], v(1), v(2), v(3));
%   label                    field(in store)  decimals
metrics = { 'Ind LMS (L, M, S)','LMS_ind','2' ; 'RGB weights (R, G, B)','rgb','2' ; ...
            'Std LMS (L, M, S)','LMS_std','2' ; 'XYZ1931 (X, Y, Z)','XYZ','1' };
% One printed sub-table per primary bandwidth sigma.
for ii = 1:numel(sdTargets)
    fprintf('\nSD = %d nm  (columns = observers, rows = metrics)\n', sdTargets(ii));
    % Header row: blank corner cell then the observer column names.
    fprintf('%-22s', '');
    for c = 1:numel(colKeys), fprintf(' | %-17s', colNames{c}); end
    fprintf('\n');
    % One row per metric (LMS / RGB / std LMS / XYZ), each column an observer.
    for m = 1:size(metrics,1)
        fprintf('%-22s', metrics{m,1});
        for c = 1:numel(colKeys)
            v = store.(colKeys{c}).(metrics{m,2})(:,ii);
            fprintf(' | %-17s', fmtV(v, metrics{m,3}));
        end
        fprintf('\n');
    end
    % Final row: the CAM16-UCS colour difference for each observer.
    fprintf('%-22s', 'ΔE'' (CAM16-UCS)');
    for c = 1:numel(colKeys), fprintf(' | %-17.2f', store.(colKeys{c}).dE(ii)); end
    fprintf('\n');
end

fprintf('\nFig 3 done.\n');

%% ============================================================
% Helper: XYZ -> Lab (CIE 1976) with reference white
%   XYZ: 3x1 or 3xN
%   XYZ_white: 3x1
%   returns Lab: 3xN
%% ============================================================
function Lab = xyz2lab_wref(XYZ, XYZ_white)
% CIE 1976 L*a*b* from tristimulus XYZ relative to a reference white.

XYZ = double(XYZ);
if size(XYZ,1) ~= 3
    XYZ = XYZ'; % allow Nx3
end
XYZ_white = double(XYZ_white(:));
% Normalise each channel by the white point (X/Xn, Y/Yn, Z/Zn).
XYZn = XYZ ./ max(XYZ_white, eps);

% Nonlinear CIELAB compression f(t): cube root above the linear-toe threshold,
% linear below it. delta = 6/29 is the knee; kappa is the toe slope.
% f(t) per CIE Lab
delta = 6/29;
thr   = delta^3;
kappa = 1/(3*delta^2);   % = (29/6)^2 / 3

f = zeros(size(XYZn));
mask = XYZn > thr;
f(mask)  = XYZn(mask).^(1/3);
f(~mask) = kappa*XYZn(~mask) + 4/29;

fx = f(1,:); fy = f(2,:); fz = f(3,:);

% Standard L*a*b* combinations of the compressed channels.
L = 116*fy - 16;
a = 500*(fx - fy);
b = 200*(fy - fz);

Lab = [L; a; b];
end


%% ============================================================
% Helper: CIEDE2000 ΔE00 between two Lab colors
%   Lab1, Lab2 are 3x1 (or 3xN matched)
%   returns dE00 as scalar (or 1xN)
%% ============================================================
function dE00 = deltaE2000(Lab1, Lab2)
% CIEDE2000 colour difference between two CIELAB colours. (Standard reference
% implementation, kept as a helper; the paper's reported metric is CAM16-UCS
% via cam16_s3, so this may be unused in the current figure flow.)

Lab1 = double(Lab1); Lab2 = double(Lab2);
if size(Lab1,1) ~= 3, Lab1 = Lab1'; end
if size(Lab2,1) ~= 3, Lab2 = Lab2'; end
assert(all(size(Lab1) == size(Lab2)), 'Lab1 and Lab2 must be same size.');

L1 = Lab1(1,:); a1 = Lab1(2,:); b1 = Lab1(3,:);
L2 = Lab2(1,:); a2 = Lab2(2,:); b2 = Lab2(3,:);

% Chroma of each colour and the mean chroma (drives the a* re-scaling G).
C1 = hypot(a1,b1);
C2 = hypot(a2,b2);
Cbar = (C1 + C2) / 2;

% G boosts a* for low-chroma colours (CIEDE2000 grey-axis correction).
Cbar7 = Cbar.^7;
G = 0.5 * (1 - sqrt(Cbar7 ./ (Cbar7 + 25^7)));

a1p = (1+G).*a1;
a2p = (1+G).*a2;

% Re-scaled chroma and hue angles (radians, wrapped to [0,2*pi)).
C1p = hypot(a1p,b1);
C2p = hypot(a2p,b2);

h1p = atan2(b1, a1p);
h2p = atan2(b2, a2p);
h1p(h1p<0) = h1p(h1p<0) + 2*pi;
h2p(h2p<0) = h2p(h2p<0) + 2*pi;

% Lightness / chroma / hue differences (hue wrapped to the short way round).
dLp = L2 - L1;
dCp = C2p - C1p;

dhp = h2p - h1p;
dhp(C1p.*C2p==0) = 0;
dhp(dhp >  pi) = dhp(dhp >  pi) - 2*pi;
dhp(dhp < -pi) = dhp(dhp < -pi) + 2*pi;

dHp = 2*sqrt(C1p.*C2p) .* sin(dhp/2);

% Mean lightness / chroma / hue used by the weighting functions below.
Lbp = (L1 + L2)/2;
Cbp = (C1p + C2p)/2;

hbp = h1p + h2p;
idx = abs(h1p - h2p) > pi;
hbp(idx) = hbp(idx) + 2*pi;
hbp = hbp/2;
hbp(C1p.*C2p==0) = h1p(C1p.*C2p==0) + h2p(C1p.*C2p==0); % doesn't matter; set to sum

% Hue-dependent weighting T (four cosine terms per CIEDE2000).
T = 1 ...
    - 0.17*cos(hbp - deg2rad(30)) ...
    + 0.24*cos(2*hbp) ...
    + 0.32*cos(3*hbp + deg2rad(6)) ...
    - 0.20*cos(4*hbp - deg2rad(63));

% Rotation term for the blue region (~275 deg) that couples chroma and hue.
dTheta = deg2rad(30) * exp(-((rad2deg(hbp) - 275)/25).^2);
RC = 2*sqrt((Cbp.^7) ./ (Cbp.^7 + 25^7));

% Lightness / chroma / hue scaling factors (SL, SC, SH) and the RT interaction.
SL = 1 + (0.015*(Lbp - 50).^2) ./ sqrt(20 + (Lbp - 50).^2);
SC = 1 + 0.045*Cbp;
SH = 1 + 0.015*Cbp.*T;

RT = -sin(2*dTheta) .* RC;

% Parametric weights (unity = reference viewing conditions).
kL = 1; kC = 1; kH = 1;

% Combined CIEDE2000 distance (weighted L/C/H plus the RT cross term).
dE00 = sqrt( ...
    (dLp./(kL*SL)).^2 + ...
    (dCp./(kC*SC)).^2 + ...
    (dHp./(kH*SH)).^2 + ...
    RT .* (dCp./(kC*SC)) .* (dHp./(kH*SH)) );
end


%% ============================================================
% Helper: render a filled spectrum into a square tile (no axes/labels)
%   - wls: nW x 1
%   - y  : nW x 1 (normalized 0..1)
%   - rgbFill: 1 x 3 in [0,1]
%   - tilePx: scalar
%% ============================================================
function tile = renderFilledSpectrumTile(wls, y, rgbFill, tilePx)
% Rasterise a single filled-spectrum square as a tilePx x tilePx RGB image.
% (Utility helper for building image mosaics; not the vector-PDF export path.)

wls = wls(:);
y   = y(:);
y   = max(min(y,1),0);

% create a hidden figure sized to tilePx x tilePx pixels
fig = figure('Visible','off','Color','w','Units','pixels','Position',[10 10 tilePx tilePx]);
ax  = axes(fig); hold(ax,'on');

% fill under curve (graph-only)
fill(ax, [wls; flipud(wls)], [zeros(size(y)); flipud(y)], ...
    rgbFill, 'EdgeColor','none', 'FaceAlpha', 1.0);

% subtle outline for definition (same idea as your reflected spectrum plot)
plot(ax, wls, y, 'k-', 'LineWidth', 0.8);

xlim(ax, [min(wls) max(wls)]);
ylim(ax, [0.01 1]);

axis(ax,'off');
ax.Position = [0 0 1 1];  % fill entire tile

% Grab the rendered figure as a raster frame and convert to an image array.
% capture to image
frame = getframe(fig);
tile  = frame2im(frame);

% ensure exact size (sometimes off by 1 px depending on renderer)
tile = imresize(tile, [tilePx tilePx]);

close(fig);
end



%% ============================================================
% Local helper: filled graph-only plot in the "illuminant style"
%% ============================================================
function makeFilledGraphOnly(wls, y, rgbFill, figp, figDir, outFile)
% Draw one reproduced-spectrum panel (filled, axis-free, 2:1 aspect) and export
% it as a transparent vector PDF. Used for every observer x sigma reproduction.
nW = numel(wls);
y = y(:);
y = max(y,0);

fig = figure('Color','w'); ax = axes(fig); hold(ax,'on');

% Closed polygon fill under the curve + thin outline (same style as panels a-c).
fill(ax, [wls(:); flipud(wls(:))], [zeros(nW,1); flipud(y)], ...
    rgbFill, 'EdgeColor','none', 'FaceAlpha', 1.0);
plot(ax, wls, y, 'k-', 'LineWidth', 0.5);

% narrow wavelength window framing the primaries (peaks 450/540/610 nm), and a
% wide 2:1 (width:height) panel for the reproduced spectra.
xlim(ax, [410 670]);
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

W = figp.twocolumn/3;  H = W/2;                 % 2:1 width:height
fig.Position = [10, 10, W, H];
ax.Position  = [0 0, W, H];                     % axes fills the whole panel
axis(ax,'off');

% Transparent vector PDF (the tight xlim keeps the drawn spectrum close to the
% 2:1 panel proportions after exportgraphics trims the transparent margin).
exportgraphics(fig, fullfile(figDir, outFile), ...
    'ContentType','vector','BackgroundColor','none');
end


%% ============================================================
% Helper: build mosaic image for sRGB patches
%% ============================================================
function img = make_srgb_mosaic(rgb, nCols, nRows, tilePx, gapPx)
% Tile a list of solid sRGB colours into a gap-separated grid image (utility
% helper for colour-patch mosaics).
N = size(rgb,1);
% Overall canvas size = tiles + inter-tile and border gaps.
W = nCols*tilePx + (nCols+1)*gapPx;
H = nRows*tilePx + (nRows+1)*gapPx;

img = ones(H, W, 3); % white background

% Fill row-major; leave any trailing tiles white when k exceeds N.
k = 1;
for r = 1:nRows
    for c = 1:nCols
        y0 = (r-1)*tilePx + r*gapPx + 1;
        y1 = y0 + tilePx - 1;
        x0 = (c-1)*tilePx + c*gapPx + 1;
        x1 = x0 + tilePx - 1;

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
%% ============================================================
function rgb = xyz_to_srgb_d65_wref(XYZ, XYZ_white)
% Convert XYZ to gamma-encoded sRGB for on-screen swatch colours only (a
% rendering aid, not part of the colorimetric measurement).
% XYZ: 3 x N (absolute scale arbitrary)
% XYZ_white: 3 x 1 reference white under same scaling
% rgb: N x 3 in [0,1]

% Scale so the white's luminance maps to 1 (Y-relative normalisation).
Yw = max(XYZ_white(2), eps);
XYZn = XYZ ./ Yw;

% Standard sRGB (D65) XYZ->linear-RGB matrix.
M = [ 3.2406 -1.5372 -0.4986;
     -0.9689  1.8758  0.0415;
      0.0557 -0.2040  1.0570];

% Linear RGB, clipped to the displayable [0,1] cube.
rgb_lin = (M * XYZn)';         % N x 3
rgb_lin = max(min(rgb_lin,1),0);

% sRGB transfer function: linear 12.92x toe, then the (1+a)*t^(1/2.4)-a curve.
a = 0.055;
thr = 0.0031308;

rgb = zeros(size(rgb_lin));
low = rgb_lin <= thr;
rgb(low)  = 12.92 * rgb_lin(low);
rgb(~low) = (1+a) * (rgb_lin(~low).^(1/2.4)) - a;

rgb = max(min(rgb,1),0);
end

%% ============================================================
% Helper: resample (S, T) to your wls vector (for standard observer cones)
%   S: [wl0 step n]  (Psychtoolbox style)
%   T: 3 x n OR n x 3
% ============================================================
function T_out = resample_T_to_wls(S, T_in, wls)
% Reinterpolate a Psychtoolbox-style spectral table onto the wls grid.
% Orient to 3 x n if it came in as n x 3.
if size(T_in,1) ~= 3 && size(T_in,2) == 3
    T_in = T_in';
end
% Reconstruct the source wavelength axis from S = [start step nSamples].
wl0 = S(1); step = S(2); n = S(3);
wl = (wl0:step:(wl0 + step*(n-1)))';
% pchip interpolation, zero outside the original support.
T_out = interp1(wl, T_in', wls(:), 'pchip', 0)'; % 3 x numel(wls)
end

%% ============================================================
% CIECAM16 forward -> CAM16-UCS [J'; a'; b'] (same as the paper's metric).
% XYZ, W are 3 x N (here 3 x 1); D in [0,1].
% ============================================================
function s3 = cam16_s3(XYZ, W, D, M16, cam)
% CIECAM16 forward model reduced to the CAM16-UCS coordinates [J'; a'; b'].
% Euclidean distance between two such vectors IS the CAM16-UCS colour
% difference used as the paper's metamerism metric (§2.3).
%   XYZ : 3 x N stimulus tristimulus (Y ~ 100 scale)
%   W   : 3 x N adopting white (same scale); D in [0,1] adaptation degree.
% Unpack surround/viewing constants; Yw = white luminance factor.
LA=cam.LA; Yb=cam.Yb; c=cam.c; Nc=cam.Nc; Yw=W(2,1,:);
% Luminance-adaptation factor FL and its intermediate k (partial adaptation).
k=1/(5*LA+1); FL=0.2*k^4*(5*LA)+0.1*(1-k^4)^2*(5*LA)^(1/3);
% Background induction factors: n = relative background, z = base exponent,
% Nbb/Ncb = brightness/chromatic background constants.
n=max(Yb./Yw,1e-6); z=1.48+sqrt(n); Nbb=0.725.*(1./n).^0.2; Ncb=Nbb;
% White through CAT16 sharpened cones (M16), then the von Kries-style degree-of-
% adaptation gains Dr; RGBcw = adapted white cone signals.
RGBw=pagemtimes(M16,W); Dr=D.*Yw./RGBw+(1-D); RGBcw=Dr.*RGBw;
% Hunt-Pointer-Estevez nonlinear response compression of the adapted white,
% and the achromatic white response Aw (reference for lightness J).
tw=(FL.*abs(RGBcw)/100).^0.42; RGBaw=sign(RGBcw).*400.*tw./(27.13+tw)+0.1;
Aw=(2*RGBaw(1,1,:)+RGBaw(2,1,:)+RGBaw(3,1,:)/20-0.305).*Nbb;
% Same pipeline for the stimulus: sharpen, adapt, then compress.
RGB=pagemtimes(M16,XYZ); RGBc=Dr.*RGB;
t=(FL.*abs(RGBc)/100).^0.42; RGBa=sign(RGBc).*400.*t./(27.13+t)+0.1;
R1=RGBa(1,:,:); R2=RGBa(2,:,:); R3=RGBa(3,:,:);
% Opponent signals a,b (red-green / yellow-blue) and achromatic response A.
a=R1-12*R2/11+R3/11; b=(R1+R2-2*R3)/9; A=(2*R1+R2+R3/20-0.305).*Nbb;
% Lightness J (stimulus-to-white achromatic ratio raised to c*z).
J=100.*sign(A).*(abs(A)./Aw).^(c.*z);
% Hue angle h and its eccentricity weighting et.
h=atan2d(b,a); et=0.25*(cos(deg2rad(h)+2)+3.8);
% Chroma C and colourfulness M from the opponent magnitude and lightness.
denC=R1+R2+(21/20)*R3; tC=(50000/13*Nc.*Ncb).*et.*hypot(a,b)./max(denC,1e-9);
Cc=(max(tC,0)).^0.9.*sqrt(max(J,0)/100).*(1.64-0.29.^n).^0.73; Mc=Cc.*FL^0.25;
% CAM16-UCS compression: J' and M' plus the Cartesian a'/b' from M' and hue.
Jp=(1+100*0.007).*J./(1+0.007.*J); Mp=(1/0.0228)*log(1+0.0228.*Mc);
s3=cat(1,Jp,Mp.*cosd(h),Mp.*sind(h));
end