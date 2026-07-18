% ########################################################################
% run_metamerism_analysis.m
% ------------------------------------------------------------------------
% MASTER SIMULATION SCRIPT for the paper
%   "Adaptation reduces inter-observer metamerism with narrow-band primaries"
%   (Optics Express).
%
% PURPOSE
%   Quantify how much individual observers DISAGREE about colour when a
%   virtual narrow-band RGB display is driven to a per-observer metameric
%   match, and how chromatic adaptation reduces that disagreement.
%
% ONE-LINE METHOD
%   For every (observer, surface, illuminant, primary bandwidth):
%     1. Compute the observer's LMS response to the "world" spectrum
%        (illuminant x surface reflectance).
%     2. Solve the 3 display drives (RGB) so the 3-primary display reproduces
%        THAT observer's LMS exactly (a per-observer metameric match).
%     3. Read every observer's reproduced stimulus in a COMMON colour space.
%     The spread BETWEEN observers of those readouts is the "observer
%     metamerism"; it is scored as a colour difference and optionally after
%     chromatic adaptation at several degrees of adaptation D.
%
% KEY INGREDIENTS
%   Observers : 10,000 Asano Individual Colorimetric Observers on the
%               CIE 2006 basis (age + gene-determined discrete lambda_max
%               shifts). Reference = CIE 1931 2-deg standard observer.
%   Surfaces  : CIE 99 colour-fidelity reflectance samples (CIE_srf_cfi_1nm).
%   Illuminants: three classes x 12 CCTs each (daylight, smooth-artificial,
%               peaky-artificial), CCT-matched, from illuminant_sets.mat.
%   Display   : 3 Gaussian primaries (peaks 450/540/610 nm), bandwidth
%               sigma in {5,10,20} nm.
%   Scoring   : CIEDE2000 in CIELAB, RMS cone-contrast error, u'v' and
%               MacLeod-Boynton scatter, plus a full CIECAM16/CAT16
%               multi-stage decomposition (adaptation degrees D in Dlist).
%
% COLORIMETRIC PIPELINE (per observer o, surface s)
%   C_obs   = T_obs * P            (3x3 : cone response to each primary)
%   LMS_o   = T_obs * world_s      (observer LMS to the world spectrum)
%   RGB_o   = inv(C_obs) * LMS_o   (metameric drives that reproduce LMS_o)
%   XYZ_o   = Mdisp * RGB_o        (what the STANDARD XYZ observer sees on
%                                   the display driven by RGB_o)
%   The between-observer spread of XYZ_o (and derived quantities) is the
%   signal of interest. Adaptation (von Kries / CAT16 / CIECAM16) is applied
%   with the observer stimulus AND the reference placed at the SAME white so
%   the global adaptation cast cancels and only residual disagreement remains.
%
% INPUTS (all under ./data, resolved relative to pwd = repo root)
%   CIE2006/Asano tables : cie2006_RelativeMacularDensity.txt,
%                          cie2006_docul.txt, cie2006_LMSAbsorbance.txt
%   CIE_srf_cfi_1nm.csv  : [wl, 99 reflectances]
%   T_cones_ss2.mat      : Stockman-Sharpe 2-deg standard-observer cones
%   CIE1931XYZ.txt       : [wl X Y Z] CIE 1931 2-deg CMFs
%   illuminant_sets.mat  : IllSet (daylight/smooth/peaky SPDs + metadata),
%                          built by build_illuminant_sets
%   Helper functions on path: generate_coneFund_family_AsanoCIE2006,
%                          get_daylightspd, build_gaussian_primaries, etc.
%
% OUTPUT
%   results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat  (-v7.3)
%   One struct with params, observer family, display family, per-illuminant
%   results (byIll), gamut bookkeeping, and all adaptation-variant metrics.
%
% MAJOR STAGES (section banners below)
%   1  Paths
%   2  Key parameters (observers, display, gamut, adaptation config)
%   3  Generate the Asano CIE2006 observer family
%   4  Load reflectances
%   5  Load standard-observer cones + XYZ CMFs
%   6  Build the illuminant list (3 classes, CCT-matched, scaled to Y=100)
%   7  Build the display primaries per bandwidth
%   8  Precompute inv(C_obs) and D matrices (illuminant-independent)
%   8a Strict per-illuminant gamut filtering + common in-gamut surface set
%   9  Main loop over illuminants (SD sweep, adaptation variants, metrics)
%   10 Save everything into one .mat
%   + Local helper functions (colorimetry, adaptation, CIEDE2000, CIECAM16).
%
% HOW TO RUN
%   From the repository root, after adding the source tree to the path:
%     >> addpath(genpath('src'));            % (plus wherever data/ lives)
%     >> run_metamerism_analysis
%   Set environment variable VERIFY_CAM=1 to run only the paged-vs-scalar
%   CIECAM16 self-test and return.
%
% REPRODUCIBILITY / PERFORMANCE
%   rng/seed for the observer ensemble is ASANO.seed = 0 (passed to the
%   generator). Uses parfor for the per-observer optimal-diagonal search when
%   the Parallel Computing Toolbox is available; the parfor body is
%   deterministic (no random draws inside), so results are pool-independent.
% ########################################################################
%
% ========================================================================
% run_AsanoCIE2006_metamerism_calc_multiIlluminant_RMSconeContrast_uv_MB.m
%
% (Original development header; retained for provenance.)
%
% UPDATE (requested):
%   - Illuminants: two CCT-matched reference classes from illuminant_sets.mat
%       (1) DAYLIGHT   : 12 CIE daylight phases (4000-8000 K)
%       (2) ARTIFICIAL : 12 real high-fidelity/high-gamut sources (TM-30
%                        Rf/Rg screened from the Houser 401 set, CCT-matched)
%   - Compute RMS cone-contrast error (SD-dependent, observer-dependent)
%   - For each surface, save u'v' for all 10000 observers (per SD, per illuminant)
%     before and after adaptation:
%       uv_indiv_no{si} : nSamples x 2 x nObs
%       uv_indiv_ad{si} : nSamples x 2 x nObs
%       uv_std_no{si}   : nSamples x 2
%       uv_std_ad{si}   : nSamples x 2
%   - Compute MacLeod–Boynton chromaticity (r,b) and save coordinates
%     for ellipse drawing (per SD, per illuminant), before and after adaptation:
%       mb_indiv_no{si} : nSamples x 2 x nObs   (r,b)
%       mb_indiv_ad{si} : nSamples x 2 x nObs   (r,b)
%       mb_std_no{si}   : nSamples x 2          (r,b)
%       mb_std_ad{si}   : nSamples x 2          (r,b)
%   - Save ALL results into ONE .mat file (-v7.3)
%
% Notes:
%   - Gamut filtering is strict across observers (per illuminant, across SDs)
%   - u'v' for individuals is computed from observer-specific XYZ (via pRGB and Mdisp per SD)
%   - Adaptation is applied in XYZ (von Kries on white patch; Y preserved)
%   - RMS cone-contrast error definition (meaningful SD dependence):
%       For each observer and SD:
%         RGB_indiv      = invC_obs * LMS_indiv_world
%         LMS_std_seen   = C_std * RGB_indiv
%         CC_seen        = (LMS_std_seen - mean_surfaces)/mean_surfaces  (per obs)
%         CC_std_world   = standard observer on TRUE surfaces, same contrast rule
%         RMS across surfaces of ||CC_seen - CC_std_world||_2
% ========================================================================

% Fresh workspace so a re-run cannot inherit stale variables; start timer.
clearvars; close all; clc;
tic;

%% ============================================================
% 1) Paths
%% ============================================================
% All inputs live under ./data and outputs under ./results, both relative to
% pwd (the repo root when launched as instructed).
dataDir = fullfile(pwd,'data');
outDir  = fullfile(pwd,'results');
if ~exist(outDir,'dir'); mkdir(outDir); end

% Single consolidated output file for the whole study.
outMat = fullfile(outDir,'results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
fprintf('  -> %s\n', outMat);

% optional self-test: paged CIECAM16 == per-observer scalar CIECAM16
% (enabled by environment variable VERIFY_CAM=1, which survives the clearvars above)
% Environment variables are unaffected by clearvars, so this switch works even
% though the workspace was just cleared; it runs the self-test and returns
% without doing the full simulation.
if ~isempty(getenv('VERIFY_CAM'))
    verify_ciecam16_paged();
    return;
end

%% ============================================================
% 2) Key parameters
%% ============================================================
% Gaussian-primary bandwidth sweep (sigma in nm): narrower primaries -> more
% inter-observer metamerism. sigma=40 nm is deliberately excluded so that a
% large surface set stays in-gamut for EVERY observer at EVERY bandwidth.
sdList   = [5 10 20];             % primary SD sweep (nm); sigma=40 dropped so a
                                  % large common in-gamut surface set survives
% Fixed peak wavelengths of the three Gaussian primaries (Blue, Green, Red).
peakBase = [450 540 610];         % fixed peaks for B,G,R (nm)

% Observer ensemble size. The Asano generator splits the population into a
% 2-deg field-size group (N2) and a 10-deg group (N10); here both halves are
% ~5000 to build a 10,000-strong ensemble.
TOTAL_OBS = 10000;
N2  = floor(TOTAL_OBS/2);
N10 = TOTAL_OBS - N2;

% Asano model options passed straight to the generator.
ASANO = struct();
ASANO.seed      = 0;              % fixes the random physiological parameter draws
ASANO.age_range = [18 70];        % sampled observer ages (years)
ASANO.normalize = 'max';          % normalise each cone fundamental to unit peak

% Representative observer to store per-SD example spectra for (clamped to nObs below)
REP_OBS_IDX = 993;

% Gamut precheck configuration (see stage 8a / gamut_precheck_all_SDs_fromPset):
%   chunkObs : observers processed per block to bound peak memory
%   tolNeg   : a drive more negative than -tolNeg counts as out-of-gamut
CFG_GAMUT = struct();
CFG_GAMUT.chunkObs = 2000;
CFG_GAMUT.tolNeg   = 1e-10;

% ------------------------------------------------------------------
% NEW adaptation validations (in addition to white-patch von Kries)
%   (1) OPT_DIAG : per-observer optimal diagonal (XYZ) that MINIMISES
%                  that observer's median DeltaE00 across surfaces.
%                  This is the ceiling of a diagonal-in-XYZ correction.
%   (2) CAT16    : CIECAM16-style chromatic adaptation (CAT16). Diagonal
%                  gain applied in the sharpened M16 cone-like space with
%                  a degree-of-adaptation D. Run at D=1 (complete, directly
%                  comparable to von Kries) and at the CIECAM16 default D.
% ------------------------------------------------------------------
% Master on/off switches for the adaptation-model variants evaluated per SD.
CFG_ADAPT = struct();
CFG_ADAPT.runOptDiag    = true;  % (1) optimal diagonal in XYZ
CFG_ADAPT.runOptM16     = true;  % (1b) optimal diagonal in the sharpened M16 basis
CFG_ADAPT.runOptPowM16  = true;  % (1c) optimal diagonal + shared power exponent in M16 (minimal nonlinear)
CFG_ADAPT.runCAT16      = true;  % (2) CIECAM16-style CAT16 (linear, CIELAB readout)
CFG_ADAPT.runCAT16stages= true;  % (3) full CIECAM16 multi-stage decomposition

% CIECAM16 viewing parameters (for the multi-stage decomposition). These fix
% the appearance-model constants (luminance level, background, surround) that
% drive the nonlinear compression and lightness/chroma exponents.
CFG_ADAPT.cam.LA = 20;           % adapting-field luminance (cd/m^2)
CFG_ADAPT.cam.Yb = 20;           % relative background luminance factor
CFG_ADAPT.cam.F  = 1.0;          % surround: average
CFG_ADAPT.cam.c  = 0.69;         % surround: average
CFG_ADAPT.cam.Nc = 1.0;          % surround: average
% Degrees of adaptation swept in the CIECAM16 stage decomposition: D=0 (no
% adaptation), the CIECAM16 default (~0.858 here), and D=1 (complete).
CFG_ADAPT.cam.Dlist = [0 0.858 1.0];  % degrees of adaptation to sweep (D=0 baseline .. D=1 complete)

% Optimiser settings for the per-observer optimal diagonal (fminsearch). The
% search is warm-started (see optimize_diag_xyz) so a tight budget suffices.
CFG_ADAPT.optObjective = 'median';        % objective on DeltaE00: 'median' or 'mean'
CFG_ADAPT.optMaxFunEvals = 200;           % fminsearch budget (warm-started, converges fast)
CFG_ADAPT.optTolX = 1e-4;
CFG_ADAPT.optTolFun = 1e-4;

% CIECAM16 degree-of-adaptation for the CAT16 variant "CIE-default"
% D = F*(1 - (1/3.6)*exp((-LA-42)/92)), clamped to [0,1]
% Two operating points are compared: complete adaptation (D_full=1, directly
% comparable to a white-patch von Kries correction) and the CIECAM16 default
% partial adaptation D_cie computed below from F and LA.
CFG_ADAPT.cat16_F  = 1.0;                 % average surround
CFG_ADAPT.cat16_LA = 20;                  % adapting-field luminance (cd/m^2); ~Lwhite/5 for a 100 cd/m^2 white
CFG_ADAPT.cat16_D_full = 1.0;             % complete adaptation (comparable to von Kries)

% CAT16 sharpened cone matrix (CIECAM16) and its inverse (XYZ with Y~100).
% M16 maps XYZ into a "sharpened" cone-like RGB space where the diagonal
% von-Kries-style adaptation gains are applied; M16inv maps back to XYZ.
M16 = [ 0.401288  0.650173 -0.051461; ...
       -0.250268  1.204414  0.045854; ...
       -0.002079  0.048952  0.953127 ];
M16inv = inv(M16);
CFG_ADAPT.M16    = M16;
CFG_ADAPT.M16inv = M16inv;

% CIECAM16 default degree of adaptation (standard CIECAM16 formula), clamped
% to the valid [0,1] range.
CFG_ADAPT.cat16_D_cie = min(1, max(0, ...
    CFG_ADAPT.cat16_F * (1 - (1/3.6)*exp((-CFG_ADAPT.cat16_LA - 42)/92)) ));
fprintf('CAT16 degrees of adaptation: D_full=%.3f, D_CIEdefault=%.3f (F=%.2f, LA=%g)\n', ...
    CFG_ADAPT.cat16_D_full, CFG_ADAPT.cat16_D_cie, CFG_ADAPT.cat16_F, CFG_ADAPT.cat16_LA);

%% ============================================================
% 3) Generate observer family (Asano CIE2006)
%% ============================================================
% Build the 10,000-observer ensemble of cone fundamentals. Each observer is a
% CIE 2006 template perturbed by age and gene-determined lambda_max shifts
% (Asano Individual Colorimetric Observer Model).
fprintf('Generating Asano CIE2006 observer family...\n');
assert(exist('generate_coneFund_family_AsanoCIE2006','file')==2, ...
    'Missing generate_coneFund_family_AsanoCIE2006.m on path.');

% Data tables the generator needs; searched for across several candidate
% subdirectories of ./data so the layout is flexible.
needFiles = {'cie2006_RelativeMacularDensity.txt','cie2006_docul.txt','cie2006_LMSAbsorbance.txt'};
candDirs  = { ...
    dataDir, ...
    fullfile(dataDir,'cie2006'), ...
    fullfile(dataDir,'CIE2006'), ...
    fullfile(dataDir,'asano'), ...
    fullfile(dataDir,'asano_cie2006'), ...
    fullfile(dataDir,'tables'), ...
    fullfile(dataDir,'cie_tables') ...
};

% Pick the first candidate directory that contains ALL required tables.
asanoDir = '';
for k = 1:numel(candDirs)
    d = candDirs{k};
    if ~exist(d,'dir'), continue; end
    ok = true;
    for f = 1:numel(needFiles)
        if exist(fullfile(d,needFiles{f}),'file')~=2
            ok = false; break;
        end
    end
    if ok
        asanoDir = d;
        break;
    end
end

assert(~isempty(asanoDir), ...
    'Could not find CIE2006/Asano tables. Put required files under ./data (or ./data/cie2006 etc.).');
fprintf('  Using CIE2006 tables from: %s\n', asanoDir);

coneFam = generate_coneFund_family_AsanoCIE2006( ...
    'data_dir',   asanoDir, ...
    'N2',         N2, ...
    'N10',        N10, ...
    'seed',       ASANO.seed, ...
    'age_range',  ASANO.age_range, ...
    'normalize',  ASANO.normalize ...
);

% The generator's output wavelength grid is the SINGLE source of truth: every
% reflectance, illuminant and CMF below is resampled onto wls.
wls      = coneFam.wl_out(:);                         % SOURCE OF TRUTH
nW       = numel(wls);
% Reorder the cone tensor to 3 (LMS) x nW (wavelengths) x nObs so that per
% observer T_indiv(:,:,o) is a 3 x nW matrix ready for matrix products.
T_indiv  = permute(coneFam.fund, [3 2 1]);            % 3 x nW x nObs
obs_meta = coneFam.meta;

nObs = size(T_indiv,3);
fprintf('Using nObs=%d (N2=%d, N10=%d), nW=%d spectral channels\n', nObs, N2, N10, nW);

% Clamp representative-observer index to the available observer count
repObsIdx = min(REP_OBS_IDX, nObs);
if repObsIdx ~= REP_OBS_IDX
    fprintf('  REP_OBS_IDX=%d exceeds nObs; using repObsIdx=%d instead.\n', REP_OBS_IDX, repObsIdx);
end

%% ============================================================
% 4) Load reflectances (99) -> interpolate to wls
%% ============================================================
fprintf('Loading reflectances: CIE_srf_cfi_1nm.csv ...\n');
refFile = fullfile(dataDir,'CIE_srf_cfi_1nm.csv');
assert(exist(refFile,'file')==2, 'Missing reflectance file: %s', refFile);

% CSV layout: column 1 = wavelength (nm), columns 2..100 = the 99 CIE
% colour-fidelity reflectance samples at 1 nm.
Rmat     = readmatrix(refFile);
wl_ref   = Rmat(:,1);
R_99_1nm = Rmat(:,2:end);
assert(size(R_99_1nm,2)==99, 'Expected 99 reflectance samples.');

% Resample reflectances onto the observer wavelength grid; clamp to physical
% [0,1] to remove any interpolation overshoot.
R_wlsx99 = interp1(wl_ref, R_99_1nm, wls, 'pchip', 'extrap');  % nW x 99
R_wlsx99 = min(max(R_wlsx99, 0), 1);

%% ============================================================
% 5) Load standard observer cones and XYZ CMFs (interpolate to wls)
%% ============================================================
fprintf('Loading standard observer: T_cones_ss2.mat ...\n');
stdFile = fullfile(dataDir,'T_cones_ss2.mat');
assert(exist(stdFile,'file')==2, 'Missing standard observer file: %s', stdFile);
% S_cones_ss2 is a PTB-style [start step nWl] wavelength descriptor; the helper
% expands it and resamples the 2-deg cone fundamentals onto wls.
load(stdFile,'S_cones_ss2','T_cones_ss2');
T_cones_std = resample_cmf_to_wls(S_cones_ss2, T_cones_ss2, wls);  % 3 x nW

fprintf('Loading CIE1931XYZ.txt ...\n');
xyzFile = fullfile(dataDir,'CIE1931XYZ.txt');
assert(exist(xyzFile,'file')==2, 'Missing CMF file: %s', xyzFile);
Axyz = readmatrix(xyzFile);
assert(size(Axyz,2) >= 4, 'CIE1931XYZ.txt must be [wl X Y Z] in cols 1..4');

% CIE 1931 2-deg CMFs: column 1 = wavelength, cols 2..4 = X,Y,Z bars. Resample
% to wls (0 outside the tabulated range) and transpose to 3 x nW so that
% XYZ = T_xyz * spectrum for any nW x 1 spectral vector.
wl_xyz     = Axyz(:,1);
T_xyz_1nm  = Axyz(:,2:4);
T_xyz      = interp1(wl_xyz, T_xyz_1nm, wls, 'pchip', 0)';        % 3 x nW

%% ============================================================
% 6) Build illuminant list: daylight class + artificial class (CCT-matched)
%% ============================================================
fprintf('Preparing illuminant classes (daylight + artificial) from illuminant_sets.mat...\n');

% --- load the two CCT-matched reference classes (built by build_illuminant_sets) ---
illSetFile = fullfile(dataDir,'illuminant_sets.mat');
assert(exist(illSetFile,'file')==2, ...
    'Missing %s. Run build_illuminant_sets first.', illSetFile);
Sset   = load(illSetFile,'IllSet');
IllSet = Sset.IllSet;
wl_set = IllSet.wl(:);

% --- standalone D65 (daylight model), for surface sRGB rendering only ---
%     (illuminant-independent reference used to colour the u'v' scatter points)
CCT_D65 = 6504;
assert(exist('get_daylightspd','file')==2, 'get_daylightspd.m not found on path.');
spd_D65 = max(get_daylightspd(CCT_D65, wls), 0);
spd_D65 = spd_D65(:);

% --- assemble the illuminant list: daylight + smooth + peaky (CCT-matched) ---
% Flatten the three illuminant classes (12 CCTs each) into one struct array
% Ill(1..nIll). Each entry carries its resampled SPD plus metadata (CCT, Duv,
% TM-30 Rf/Rg, spectral peakiness) for later grouping and reporting.
Ill  = struct([]);
ii   = 0;
cats = {'daylight','smooth','peaky'};
for c = 1:numel(cats)
    catName = cats{c};
    spdC  = IllSet.(catName).spd;      % 81 x M
    metaC = IllSet.(catName).meta;     % [CCT Duv Rf Rg peakiness]
    Mcat  = size(spdC,2);
    for m = 1:Mcat
        % Resample this SPD onto the common wavelength grid and clip negatives.
        spd_m = interp1(wl_set, spdC(:,m), wls, 'linear', 'extrap');
        spd_m = max(spd_m(:), 0);
        ii = ii + 1;
        Ill(ii).name      = sprintf('%s%02d', catName(1:3), m);   % day01.. smo01.. pea01..
        Ill(ii).category  = string(catName);
        Ill(ii).type      = string(catName);
        Ill(ii).source    = "illuminant_sets.mat";
        Ill(ii).CCT       = metaC(m,1);
        Ill(ii).Duv       = metaC(m,2);
        Ill(ii).Rf        = metaC(m,3);
        Ill(ii).Rg        = metaC(m,4);
        Ill(ii).peakiness = metaC(m,5);
        Ill(ii).label     = sprintf('%s #%d (CCT=%.0fK, Rf=%.0f, pk=%.3f)', ...
                                    catName, m, metaC(m,1), metaC(m,3), metaC(m,5));
        Ill(ii).spd       = spd_m;
        Ill(ii).isRep     = false;
    end
end

% --- mark one representative illuminant per class (nearest 6500 K) for the ---
%     illustrative per-illuminant figures (u'v' scatter, spectra) ---
catVec = string({Ill.category});
cctVec = [Ill.CCT];
for c = 1:numel(cats)
    inC = find(catVec == string(cats{c}));
    [~, jj] = min(abs(cctVec(inC) - 6500));
    Ill(inC(jj)).isRep = true;
end

nIll = numel(Ill);

% Normalise every illuminant to a common luminance level so that a perfect
% (unit) diffuser yields Y=100 under each. T_xyz(2,:) is the Ybar row, so
% T_xyz(2,:)*spd is the diffuser luminance; the scale factor sets it to 100.
% This makes CIELAB / CIECAM readouts and colour differences comparable across
% illuminants of different absolute power.
for k = 1:nIll
    Yk = T_xyz(2,:) * Ill(k).spd;
    sY = 100 / max(Yk, eps);
    Ill(k).scaleToY100 = sY;
    Ill(k).spd = Ill(k).spd * sY;

    % Cache the illuminant's own white point (XYZ and chromaticity xy).
    XYZk = T_xyz * Ill(k).spd;
    Ill(k).XYZ = XYZk;
    Ill(k).xy  = xyz_to_xy(XYZk);
end

% cache scaled D65 SPD for ellipse colors (illuminant-independent reference,
% scaled to Y(perfect diffuser)=100 to match the illuminant scaling above)
spd_D65_scaled = spd_D65 * (100 / max(T_xyz(2,:) * spd_D65, eps));

%% ============================================================
% 7) Define display primaries (outside illuminant loop) and scale white to Y=100
%% ============================================================
nSD = numel(sdList);

% Per-bandwidth caches, built once and reused across all illuminants:
%   Pset{si}      : nW x 3 primary SPD matrix (columns = B,G,R primaries)
%   Mdisp_set{si} : 3 x 3 display->XYZ matrix (XYZ = Mdisp * RGB drives)
%   scale_set     : luminance scale applied to the primaries
%   Ywhite0_set   : unscaled full-white luminance (diagnostic)
Pset        = cell(nSD,1);   % each: nW x 3
Mdisp_set   = cell(nSD,1);   % each: 3 x 3
scale_set   = zeros(nSD,1);
Ywhite0_set = zeros(nSD,1);

for si = 1:nSD
    sd = sdList(si);

    % Unit-peak Gaussian primaries at this bandwidth.
    P0 = build_gaussian_primaries(wls, peakBase, sd);   % nW x 3 (unit max)
    % Full-white SPD = sum of all three primaries at max drive; its luminance
    % sets the display scale so that RGB=[1 1 1] produces Y=100.
    SPD_white0 = P0 * ones(3,1);                        % nW x 1
    Ywhite0 = T_xyz(2,:) * SPD_white0;                  % scalar
    Ywhite0_set(si) = Ywhite0;

    s = 100 / max(Ywhite0, eps);
    scale_set(si) = s;

    % Scaled primaries and the resulting display primary matrix in XYZ.
    P = P0 * s;
    Pset{si}      = P;
    Mdisp_set{si} = T_xyz * P;                          % 3 x 3
end

%% ============================================================
% 8) Precompute invC and D (illuminant-independent)
%% ============================================================
% The cone-to-primary matrix C_obs = T_obs*P and its inverse depend only on the
% observer and the display primaries, NOT on the illuminant, so they are built
% once per bandwidth and reused inside the illuminant loop.
%   C_obs      : each observer's LMS response to the 3 primaries (3x3)
%   invC       : maps a target LMS to the primary drives (RGB = invC*LMS)
%   C_std      : the SAME for the standard-observer cones (used for cone
%                contrast and MacLeod-Boynton readouts)
%   D = C_std*invC : the 3x3 map from an observer's drives to the standard
%                cone excitations they produce; how far D departs from a pure
%                diagonal (off-diagonal energy) indexes the metameric mismatch.
fprintf('\nPrecomputing invC and D (independent of illuminant)...\n');

saveD = true;

D_set          = cell(nSD,1);          % each: 3 x 3 x nObs (per SD)
offDiag_perObs = nan(nSD, nObs);       % nSD x nObs off-diagonal ratio
invC_set       = cell(nSD,1);          % each: 3 x 3 x nObs (per SD)
Cstd_set       = cell(nSD,1);          % each: 3 x 3 (per SD), used for cone-contrast error & MB

for si = 1:nSD
    P = Pset{si};

    % pagemtimes multiplies each 3 x nW observer slice by the nW x 3 primaries,
    % giving a 3 x 3 x nObs stack of per-observer cone-to-primary matrices.
    C_obs = pagemtimes(T_indiv, P);      % 3 x 3 x nObs
    invC  = pageinv_or_loop(C_obs);      % 3 x 3 x nObs
    invC_set{si} = invC;

    C_std = T_cones_std * P;             % 3 x 3
    Cstd_set{si} = C_std;

    if saveD
        % D_o = C_std * invC_o for every observer (batched via pagemtimes).
        D = pagemtimes(repmat(C_std,[1 1 nObs]), invC);   % 3 x 3 x nObs
        D_set{si} = D;
        offDiag_perObs(si,:) = offdiag_energy_ratio(D);   % 1 x nObs
    end
end

%% ============================================================
% 8a) PER-ILLUMINANT strict gamut precheck (across observers & all SDs)
%     Each illuminant keeps the surfaces that are in-gamut (reproducible on
%     the display for EVERY observer, at EVERY SD) under THAT illuminant.
%     A global intersection across a broad, diverse illuminant set would
%     collapse to near-empty, so surfaces are screened per illuminant; the
%     between-observer metamerism metrics are medians over surfaces and are
%     therefore comparable across illuminants despite differing surface sets.
%% ============================================================
fprintf('\nPer-illuminant gamut filtering (strict across observers & SDs)...\n');

% A surface is IN-GAMUT only if its per-observer metameric drives RGB=invC*LMS
% are non-negative for EVERY observer at EVERY bandwidth under this illuminant.
% Negative drives cannot be physically displayed, so such surfaces are dropped.
allIdx = (1:size(R_wlsx99,2))';   % 1..99
gamutByIll = cell(nIll,1);        % per-illuminant precheck summaries
keepByIll  = cell(nIll,1);        % per-illuminant kept-surface indices
badUnionAll = [];                 % union of bad across illuminants (for reference)

for kIll = 1:nIll
    % World spectra (illuminant x each reflectance) fed to the gamut checker.
    SPD_world_99 = Ill(kIll).spd .* R_wlsx99;  % nW x 99 (apply this illuminant)
    gamutByIll{kIll} = gamut_precheck_all_SDs_fromPset(T_indiv, SPD_world_99, Pset, CFG_GAMUT);

    % A surface out-of-gamut at ANY bandwidth is dropped for this illuminant.
    badThisIll = unique(vertcat(gamutByIll{kIll}.badReflectanceIdx{:})); % bad for ANY SD
    keepByIll{kIll} = setdiff(allIdx, badThisIll);
    badUnionAll = union(badUnionAll, badThisIll);

    if isempty(keepByIll{kIll})
        error('No reflectances remain in-gamut for illuminant %s (strict).', ...
            string(Ill(kIll).name));
    end

    fprintf('  %s: keep=%d/%d (strict across observers; union across SDs)\n', ...
        string(Ill(kIll).name), numel(keepByIll{kIll}), numel(allIdx));
end

% --- enforce a COMMON surface set across EVERY illuminant and SD ---
% Take the intersection of the per-illuminant in-gamut sets so that every
% illuminant class and every primary bandwidth is evaluated on the SAME
% surfaces. This turns the class/bandwidth comparison into a within-surface
% design and removes the gamut-composition confound (different illuminants
% otherwise retain different, and differently saturated, surface sets).
% Intersect the per-illuminant in-gamut sets down to the surfaces common to
% ALL illuminants and bandwidths, then overwrite each illuminant's keep-list
% with that common set (within-surface design; see banner above).
commonKeep = keepByIll{1};
for kIll = 2:nIll, commonKeep = intersect(commonKeep, keepByIll{kIll}); end
if isempty(commonKeep)
    error('No surfaces are in-gamut in common across all illuminants and SDs.');
end
fprintf('COMMON surface set across all %d illuminants x %d SDs: %d/%d surfaces.\n', ...
    nIll, numel(sdList), numel(commonKeep), numel(allIdx));
for kIll = 1:nIll, keepByIll{kIll} = commonKeep; end

% keepIdx_global here means "kept for the representative daylight illuminant";
% used only by the illustrative spectra/surface-rendering figures.
repIll = find([Ill.isRep] & (string({Ill.category})=="daylight"), 1);
if isempty(repIll), repIll = 1; end
keepIdx_global = keepByIll{repIll};
fprintf('Representative daylight (%s) keeps %d/%d surfaces (used for illustrative figures).\n', ...
    string(Ill(repIll).name), numel(keepIdx_global), numel(allIdx));

%% ============================================================
% 9) Main loop over illuminants
%% ============================================================
fprintf('\nRunning all illuminant conditions (%d total)...\n', nIll);

% start a parallel pool for the per-observer optimal-diagonal search (if the
% Parallel Computing Toolbox is available; parfor falls back to serial otherwise)
if license('test','Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
    try
        parpool('local', feature('numcores'));
    catch ME
        fprintf('  (parpool unavailable: %s; running serially)\n', ME.message);
    end
end

% One entry of IllResults per illuminant condition (populated at loop end).
IllResults = struct([]);
for kIll = 1:nIll
    fprintf('\n==============================\n');
    fprintf('Illuminant %d/%d: %s\n', kIll, nIll, Ill(kIll).label);
    fprintf('xy = [%.6f, %.6f]\n', Ill(kIll).xy(1), Ill(kIll).xy(2));

    % Apply illuminant to reflectances
    SPD_world_99 = Ill(kIll).spd .* R_wlsx99;     % nW x 99

    % --- use this illuminant's own in-gamut surface set ---
    % (After stage 8a this equals the common set, so nSamples is identical
    %  across illuminants.)
    keepIdx = keepByIll{kIll};
    SPD_world = Ill(kIll).spd .* R_wlsx99(:, keepIdx);   % nW x nSamples
    nSamples  = size(SPD_world,2);

    % (optional) keep per-illuminant gamut summary for reporting/debug
    gamut = gamutByIll{kIll};
    badThisIll = unique(vertcat(gamut.badReflectanceIdx{:}));
    badUnion   = badUnionAll;

    %% --- compute sRGB surface colors under D65 (for ellipse coloring) ---
    % Render each kept surface under a fixed D65 (illuminant-independent) purely
    % to colour the scatter/ellipse markers in later figures; not used in any
    % metric. srgb_D65_cell just replicates the colours for every bandwidth.
    SPD_D65_keep = spd_D65_scaled .* R_wlsx99(:, keepIdx);   % nW x nSamples
    XYZ_D65_keep = T_xyz * SPD_D65_keep;                     % 3 x nSamples (Y white ~100)
    srgbSurf_D65 = xyz_to_srgb_d65(XYZ_D65_keep);            % nSamples x 3 in [0,1]
    srgb_D65_cell = repmat({srgbSurf_D65}, nSD, 1);

    %% --- World XYZ/Lab reference + brightest-surface white ---
    % Standard-observer XYZ of the true world spectra = the ground-truth target
    % every observer's reproduction is compared against. The brightest surface
    % (max Y) serves as the adopted white for CIELAB and for the von Kries /
    % CAT16 adaptation anchors.
    XYZ_world = T_xyz * SPD_world;              % 3 x nSamples
    [~, idxWhite] = max(XYZ_world(2,:));
    XYZ_white = XYZ_world(:, idxWhite);
    XYZ_white(XYZ_white<=0) = eps;

    Lab_ref = xyz_to_lab(XYZ_world, XYZ_white); % 3 x nSamples

    %% --- Precompute world LMS for all observers ---
    % Each observer's true LMS response to the world spectra (batched slice mult
    % over the 3 x nW x nObs cone tensor): 3 x nSamples x nObs.
    LMS_indiv_world = pagemtimes(T_indiv, SPD_world);   % 3 x nSamples x nObs

    %% --- Standard observer cone-contrast reference on TRUE surfaces ---
    % Reference cone contrast: standard-observer LMS relative to its mean across
    % surfaces. The per-observer reproduced contrast is later compared to this.
    LMS_std_world = T_cones_std * SPD_world;            % 3 x nSamples
    mu_std_world  = mean(LMS_std_world, 2, 'omitnan');  % 3 x 1
    mu_std_world  = max(mu_std_world, 1e-12);
    CC_std_world  = (LMS_std_world - mu_std_world) ./ mu_std_world;  % 3 x nSamples

    %% --- SD sweep outputs ---
    % Preallocate every per-(bandwidth x observer) metric collected in the SD
    % loop below. medDE_no/ad = median CIEDE2000 vs the world reference, before
    % ("no") and after ("ad") white-patch von Kries adaptation.
    medDE_no   = nan(nSD, nObs);
    medDE_ad   = nan(nSD, nObs);

    % NEW adaptation variants (median DeltaE00 per SD x observer)
    medDE_opt        = nan(nSD, nObs);   % (1) optimal diagonal in XYZ
    meanDE_opt       = nan(nSD, nObs);   % companion mean DeltaE00 for the optimal diagonal
    medDE_opt_m16    = nan(nSD, nObs);   % (1b) optimal diagonal in the M16 sharpened basis
    medDE_optpow_m16 = nan(nSD, nObs);   % (1c) optimal diagonal + shared power exponent in M16 (nonlinear)
    medDE_cat16full  = nan(nSD, nObs);   % (2a) CAT16, D=1 (complete)
    medDE_cat16cie   = nan(nSD, nObs);   % (2b) CAT16, D=CIECAM16 default
    gains_opt        = nan(3, nObs, nSD);% stored optimal diagonal gains (gX,gY,gZ)

    % u'v' for new variants (for optional ellipse plots)
    uv_indiv_opt       = cell(nSD,1);    % each: nSamples x 2 x nObs
    uv_indiv_cat16full = cell(nSD,1);
    uv_std_opt         = cell(nSD,1);    % each: nSamples x 2
    uv_std_cat16full   = cell(nSD,1);

    % (3) CIECAM16 multi-stage decomposition:
    %   median between-observer distance at each CIECAM16 stage, for each D
    %   stage 1 = linear post-CAT16 (M16 adapted signals)
    %   stage 2 = post-nonlinear compressed cone responses (R'a G'a B'a)
    %   stage 3 = CAM16-UCS [J' a' b']  (perceptual color difference)
    nD_cam   = numel(CFG_ADAPT.cam.Dlist);
    cam_d1   = nan(nSD, nObs, nD_cam);
    cam_d2   = nan(nSD, nObs, nD_cam);
    cam_d3   = nan(nSD, nObs, nD_cam);
    cam_d3c  = nan(nSD, nObs, nD_cam);   % control: observer forced to common (reference) white
    cam_d3lig= nan(nSD, nObs, nD_cam);   % own-white stage-3 lightness-only |ΔJ'|
    cam_d3chr= nan(nSD, nObs, nD_cam);   % own-white stage-3 chromatic-only hypot(Δa',Δb')

    rmsConeContrastErr = nan(nSD, nObs);   % per SD x per observer

    % u'v'
    uv_indiv_no = cell(nSD,1);  % each: nSamples x 2 x nObs
    uv_indiv_ad = cell(nSD,1);  % each: nSamples x 2 x nObs
    uv_std_no   = cell(nSD,1);  % each: nSamples x 2
    uv_std_ad   = cell(nSD,1);  % each: nSamples x 2

    % MacLeod–Boynton (r,b)
    mb_indiv_no = cell(nSD,1);  % each: nSamples x 2 x nObs
    mb_indiv_ad = cell(nSD,1);  % each: nSamples x 2 x nObs
    mb_std_no   = cell(nSD,1);  % each: nSamples x 2
    mb_std_ad   = cell(nSD,1);  % each: nSamples x 2

    meanOffDiagRatio = nan(nSD,1);
    semOffDiagRatio  = nan(nSD,1);

    fprintf('  Running SD sweep on nSamples=%d...\n', nSamples);

    for si = 1:nSD
        P     = Pset{si};          % nW x 3
        Mdisp = Mdisp_set{si};     % 3 x 3
        invC  = invC_set{si};      % 3 x 3 x nObs
        C_std = Cstd_set{si};      % 3 x 3

        % --- Per-observer metameric match ---
        % Solve the display drives that reproduce each observer's world LMS:
        %   RGB_o = invC_o * LMS_o. This is the heart of the method: every
        %   observer gets a display stimulus that is a metameric match to the
        %   world FOR THEM.
        pRGB = pagemtimes(invC, LMS_indiv_world);                   % 3 x nSamples x nObs

        % --- Common-space readout ---
        % Push each observer's drives through the display->XYZ matrix to get
        % what the STANDARD XYZ observer sees on that stimulus. The spread of
        % these XYZ across observers IS the inter-observer metamerism.
        XYZ  = pagemtimes(repmat(Mdisp,[1 1 nObs]), pRGB);          % 3 x nSamples x nObs
        uv_indiv_no{si} = xyz_pages_to_upvp(XYZ);                   % nSamples x 2 x nObs

        % Standard observer reference u'v' on TRUE surfaces
        uv_std_no{si} = xyz_to_upvp(XYZ_world);                     % nSamples x 2

        % ------------------------------------------------------------
        % NEW: MacLeod–Boynton (r,b) BEFORE adaptation
        % Use standard-cone excitations "seen" when displaying RGB_indiv:
        %   LMS_std_seen = C_std * RGB_indiv
        % ------------------------------------------------------------
        % LMS_std_seen = standard-observer cone excitations produced by each
        % observer's drives (C_std * RGB_o); the between-observer scatter of
        % these in MB space is another view of the metamerism.
        LMS_std_seen = pagemtimes(repmat(C_std,[1 1 nObs]), pRGB);   % 3 x nSamples x nObs
        mb_indiv_no{si} = lms_pages_to_macleodboynton_rb(LMS_std_seen); % nSamples x 2 x nObs
        mb_std_no{si}   = lms_to_macleodboynton_rb(LMS_std_world);      % nSamples x 2

        % Median CIEDE2000 per observer BEFORE adaptation: convert each
        % observer's reproduced XYZ to CIELAB (adopted white XYZ_white) and
        % compare to the world reference Lab_ref.
        Lab = xyz_to_lab_pages(XYZ, XYZ_white);
        medDE_no(si,:) = median_de00_per_observer(Lab, Lab_ref);

        % White-patch von Kries in XYZ: scale each observer's XYZ so their
        % reproduced white patch lands on the reference white (global cast
        % removed); Y of the white is preserved by construction.
        XYZ_ad = vonkries_xyz_whitePatch(XYZ, XYZ_white, idxWhite);
        uv_indiv_ad{si} = xyz_pages_to_upvp(XYZ_ad);                % nSamples x 2 x nObs

        % The reference is adapted with the SAME rule so any residual is genuine
        % inter-observer disagreement, not a shared adaptation shift.
        XYZ_world_ad = vonkries_xyz_whitePatch_std(XYZ_world, XYZ_white, idxWhite);
        uv_std_ad{si} = xyz_to_upvp(XYZ_world_ad);                  % nSamples x 2
        
        % repObsIdx is defined once above (clamped to nObs)
        % Save only rep observer per SD
        % Allocate the representative-observer caches on the first bandwidth,
        % then store that single observer's drives and XYZ (pre/post adapt) for
        % illustrative figures.
        if si == 1
            rep_XYZ_no = cell(nSD,1);
            rep_XYZ_ad = cell(nSD,1);
            rep_RGB    = cell(nSD,1);
        end

        rep_RGB{si}    = squeeze(pRGB(:,:,repObsIdx));      % 3 x nSamples
        rep_XYZ_no{si} = squeeze(XYZ(:,:,repObsIdx));       % 3 x nSamples
        rep_XYZ_ad{si} = squeeze(XYZ_ad(:,:,repObsIdx));    % 3 x nSamples

        % ------------------------------------------------------------
        % NEW: MacLeod–Boynton (r,b) AFTER adaptation
        % Apply VK gains in LMS space using the same white-patch index.
        % ------------------------------------------------------------
        % Same white-patch von Kries, but applied directly in LMS (gains that
        % send the white-patch LMS to unity) for the MacLeod-Boynton readout.
        LMS_std_seen_ad = vonkries_lms_whitePatch(LMS_std_seen, idxWhite);       % 3 x nSamples x nObs
        mb_indiv_ad{si} = lms_pages_to_macleodboynton_rb(LMS_std_seen_ad);       % nSamples x 2 x nObs

        LMS_std_world_ad = vonkries_lms_whitePatch_std(LMS_std_world, idxWhite); % 3 x nSamples
        mb_std_ad{si}    = lms_to_macleodboynton_rb(LMS_std_world_ad);            % nSamples x 2

        % Median CIEDE2000 per observer AFTER white-patch von Kries adaptation.
        Lab_ad = xyz_to_lab_pages(XYZ_ad, XYZ_white);
        medDE_ad(si,:) = median_de00_per_observer(Lab_ad, Lab_ref);

        % ============================================================
        % NEW (1): OPTIMAL DIAGONAL in XYZ (per observer)
        %   Find g=(gX,gY,gZ) minimising that observer's median (or mean)
        %   DeltaE00 across surfaces. Warm-started at the white-patch gains
        %   and a least-squares diagonal, then polished with fminsearch.
        %   This is the ceiling of any diagonal-in-XYZ correction.
        % ============================================================
        if CFG_ADAPT.runOptDiag
            optSet = optimset('Display','off', ...
                'MaxFunEvals', CFG_ADAPT.optMaxFunEvals, ...
                'MaxIter',     CFG_ADAPT.optMaxFunEvals, ...
                'TolX',        CFG_ADAPT.optTolX, ...
                'TolFun',      CFG_ADAPT.optTolFun);
            % per-observer optimal-diagonal search (parallel over observers).
            % Results collected in per-observer sliced arrays, then stored.
            % Each column iObs is written independently, so the loop is a valid
            % parfor with sliced outputs; scalars pulled out of the CFG struct
            % (runM16, objT, M16_, ...) become parfor-broadcast variables and
            % keep the parfor transparent (no struct field access inside).
            XYZ_opt   = zeros(size(XYZ));   % 3 x nSamples x nObs
            g_col     = zeros(3, nObs);
            med_col   = zeros(1, nObs);
            mean_col  = zeros(1, nObs);
            m16_col   = nan(1, nObs);
            pow_col   = nan(1, nObs);
            runM16    = CFG_ADAPT.runOptM16;
            runPow    = CFG_ADAPT.runOptPowM16;
            M16_ = CFG_ADAPT.M16; M16inv_ = CFG_ADAPT.M16inv; objT = CFG_ADAPT.optObjective;
            parfor iObs = 1:nObs
                % Slice this observer's reproduced XYZ and find the diagonal
                % gain g that minimises their median/mean CIEDE2000 vs the world
                % reference (ceiling of any diagonal-in-XYZ correction).
                XYZo = XYZ(:,:,iObs);                              % 3 x nSamples (sliced input)
                [g, medv, meanv] = optimize_diag_xyz( ...
                    XYZo, Lab_ref, XYZ_world, XYZ_white, idxWhite, objT, optSet);
                g_col(:,iObs)     = g;
                med_col(iObs)     = medv;
                mean_col(iObs)    = meanv;
                XYZ_opt(:,:,iObs) = g .* XYZo;                     % broadcast 3x1 .* 3xN

                % (1b) optimal diagonal in the sharpened M16 basis
                if runM16
                    m16_col(iObs) = optimize_diag_m16( ...
                        XYZo, Lab_ref, XYZ_world, XYZ_white, idxWhite, ...
                        M16_, M16inv_, objT, optSet);
                end

                % (1c) optimal diagonal + shared power exponent in M16 (nonlinear)
                if runPow
                    pow_col(iObs) = optimize_diag_power_m16( ...
                        XYZo, Lab_ref, XYZ_white, idxWhite, ...
                        M16_, M16inv_, objT, optSet);
                end

            end
            % Gather the parfor sliced results into the per-SD stores.
            gains_opt(:,:,si)   = g_col;
            medDE_opt(si,:)     = med_col;
            meanDE_opt(si,:)    = mean_col;
            medDE_opt_m16(si,:) = m16_col;
            medDE_optpow_m16(si,:) = pow_col;
            uv_indiv_opt{si} = xyz_pages_to_upvp(XYZ_opt);        % nSamples x 2 x nObs
            uv_std_opt{si}   = xyz_to_upvp(XYZ_world);            % std maps to itself
        end

        % ============================================================
        % NEW (2): CIECAM16-style CAT16 adaptation
        %   Diagonal gain in the sharpened M16 cone-like space, mapping each
        %   observer's reproduced white to the reference white, with degree
        %   of adaptation D. D=1 is complete (comparable to von Kries);
        %   D=CIECAM16-default models realistic partial adaptation.
        % ============================================================
        if CFG_ADAPT.runCAT16
            % Complete adaptation (D=1): CAT16 gains in M16 space map each
            % observer's reproduced white exactly to the reference white; score
            % the residual median CIEDE2000.
            [XYZ_c_full, ~] = cat16_adapt_pages(XYZ, XYZ_white, idxWhite, ...
                CFG_ADAPT.M16, CFG_ADAPT.M16inv, CFG_ADAPT.cat16_D_full);
            Lab_c_full = xyz_to_lab_pages(XYZ_c_full, XYZ_white);
            medDE_cat16full(si,:) = median_de00_per_observer(Lab_c_full, Lab_ref);
            uv_indiv_cat16full{si} = xyz_pages_to_upvp(XYZ_c_full);
            uv_std_cat16full{si}   = xyz_to_upvp(XYZ_world);      % std maps to itself

            % Partial adaptation at the CIECAM16 default D (realistic incomplete
            % adaptation); same scoring.
            [XYZ_c_cie, ~] = cat16_adapt_pages(XYZ, XYZ_white, idxWhite, ...
                CFG_ADAPT.M16, CFG_ADAPT.M16inv, CFG_ADAPT.cat16_D_cie);
            Lab_c_cie = xyz_to_lab_pages(XYZ_c_cie, XYZ_white);
            medDE_cat16cie(si,:) = median_de00_per_observer(Lab_c_cie, Lab_ref);
        end

        % ============================================================
        % NEW (3): CIECAM16 MULTI-STAGE decomposition
        %   Run each observer's reproduced stimuli (and the standard-observer
        %   world reference) through the full CIECAM16 forward model and read
        %   out the between-observer distance at three successive stages, for
        %   each degree of adaptation D in CFG_ADAPT.cam.Dlist. The change from
        %   D=0 to D=1 isolates the chromatic-adaptation stage; the change of
        %   the ratio across stages isolates the nonlinear/opponent/UCS stages.
        % ============================================================
        if CFG_ADAPT.runCAT16stages
            % Sweep degree of adaptation; for each D read the between-observer
            % distance at CIECAM16 stages 1/2/3 (plus common-white control and
            % lightness/chromatic splits of stage 3).
            for di = 1:numel(CFG_ADAPT.cam.Dlist)
                Dc = CFG_ADAPT.cam.Dlist(di);
                DEs = ciecam16_stage_metamerism(XYZ, XYZ_world, XYZ_white, ...
                    idxWhite, Dc, CFG_ADAPT.M16, CFG_ADAPT.M16inv, CFG_ADAPT.cam);
                cam_d1(si,:,di) = DEs.d1;
                cam_d2(si,:,di) = DEs.d2;
                cam_d3(si,:,di) = DEs.d3;
                cam_d3c(si,:,di)   = DEs.d3c;
                cam_d3lig(si,:,di) = DEs.d3lig;
                cam_d3chr(si,:,di) = DEs.d3chr;
            end
        end

        % ============================================================
        % RMS cone-contrast error (SD-dependent)
        %   CC_seen computed per observer relative to mean across surfaces
        %   Compare vs CC_std_world
        % ============================================================
        % Per-observer cone contrast of the standard-observer excitations
        % actually seen (LMS_std_seen), relative to that observer's own mean
        % across surfaces.
        mu_seen = mean(LMS_std_seen, 2, 'omitnan');                 % 3 x 1 x nObs
        mu_seen = max(mu_seen, 1e-12);

        CC_seen = (LMS_std_seen - mu_seen) ./ mu_seen;              % 3 x nSamples x nObs

        % Error vs the standard-observer world cone contrast, summed over LMS,
        % then RMS across surfaces -> one SD-dependent scalar per observer.
        err = CC_seen - repmat(CC_std_world, [1 1 nObs]);           % 3 x nSamples x nObs
        err2 = squeeze(sum(err.^2, 1));                             % nSamples x nObs
        rmsConeContrastErr(si,:) = sqrt(mean(err2, 1, 'omitnan'));   % 1 x nObs

        % Off-diagonal energy ratio of D = C_std*invC: how non-diagonal the
        % observer->standard cone map is (0 = pure gain change, correctable by
        % adaptation; larger = irreducible metameric distortion). Reuse the
        % precomputed values when saveD, else recompute here.
        if saveD
            offEach = offDiag_perObs(si,:);
        else
            D = pagemtimes(repmat(C_std,[1 1 nObs]), invC);
            offEach = offdiag_energy_ratio(D);
        end

        meanOffDiagRatio(si) = mean(offEach,'omitnan');
        nEff = sum(isfinite(offEach));
        semOffDiagRatio(si)  = std(offEach,0,'omitnan') ./ max(sqrt(nEff),1);

        fprintf('    SD=%g nm | mean RMS cone-contrast err = %.6f\n', sdList(si), ...
            mean(rmsConeContrastErr(si,:), 'omitnan'));
    end

    % Post/pre adaptation ratio per (SD x observer): <1 means adaptation
    % reduced that observer's metamerism relative to no-adapt.
    ratio_ad_over_no = medDE_ad ./ (medDE_no + eps);

    % Ratios for the new adaptation variants (post/pre, relative to No-Adapt)
    ratio_opt_over_no       = medDE_opt       ./ (medDE_no + eps);
    ratio_optm16_over_no    = medDE_opt_m16   ./ (medDE_no + eps);
    ratio_optpowm16_over_no = medDE_optpow_m16./ (medDE_no + eps);
    ratio_cat16full_over_no = medDE_cat16full ./ (medDE_no + eps);
    ratio_cat16cie_over_no  = medDE_cat16cie  ./ (medDE_no + eps);

    %% --- store per-illuminant results ---
    % Pack every metric and diagnostic computed for this illuminant into one
    % IllResults element. (ok flags a completed condition.)
    IllResults(kIll).ok = true;

    IllResults(kIll).illuminant = Ill(kIll);
    IllResults(kIll).gamutSummary = gamut;
    IllResults(kIll).keepIdx      = keepIdx;
    IllResults(kIll).badUnion     = badUnion;
    IllResults(kIll).nSamples     = nSamples;

    IllResults(kIll).idxWhite   = idxWhite;
    IllResults(kIll).XYZ_white  = XYZ_white;

    IllResults(kIll).XYZ_world  = XYZ_world;
    IllResults(kIll).Lab_ref    = Lab_ref;

    IllResults(kIll).medDE_no         = medDE_no;
    IllResults(kIll).medDE_ad         = medDE_ad;
    IllResults(kIll).ratio_ad_over_no = ratio_ad_over_no;

    % NEW adaptation variants
    IllResults(kIll).medDE_opt              = medDE_opt;
    IllResults(kIll).meanDE_opt             = meanDE_opt;
    IllResults(kIll).medDE_opt_m16          = medDE_opt_m16;
    IllResults(kIll).ratio_optm16_over_no   = ratio_optm16_over_no;
    IllResults(kIll).medDE_optpow_m16       = medDE_optpow_m16;
    IllResults(kIll).ratio_optpowm16_over_no= ratio_optpowm16_over_no;
    IllResults(kIll).medDE_cat16full        = medDE_cat16full;
    IllResults(kIll).medDE_cat16cie         = medDE_cat16cie;
    IllResults(kIll).ratio_opt_over_no      = ratio_opt_over_no;
    IllResults(kIll).ratio_cat16full_over_no= ratio_cat16full_over_no;
    IllResults(kIll).ratio_cat16cie_over_no = ratio_cat16cie_over_no;
    IllResults(kIll).gains_opt              = gains_opt;      % 3 x nObs x nSD
    IllResults(kIll).cat16_D_full           = CFG_ADAPT.cat16_D_full;
    IllResults(kIll).cat16_D_cie            = CFG_ADAPT.cat16_D_cie;

    IllResults(kIll).uv_indiv_opt        = uv_indiv_opt;
    IllResults(kIll).uv_indiv_cat16full  = uv_indiv_cat16full;
    IllResults(kIll).uv_std_opt          = uv_std_opt;
    IllResults(kIll).uv_std_cat16full    = uv_std_cat16full;

    % (3) CIECAM16 multi-stage decomposition (median distance per SD x obs x D)
    IllResults(kIll).cam_Dlist  = CFG_ADAPT.cam.Dlist;
    IllResults(kIll).cam_params = CFG_ADAPT.cam;
    IllResults(kIll).cam_d1     = cam_d1;   % stage 1: linear post-CAT16
    IllResults(kIll).cam_d2     = cam_d2;   % stage 2: post-nonlinear compressed cone responses
    IllResults(kIll).cam_d3     = cam_d3;   % stage 3: CAM16-UCS [J' a' b']
    IllResults(kIll).cam_d3c    = cam_d3c;   % control: observer at common (reference) white
    IllResults(kIll).cam_d3lig  = cam_d3lig; % stage-3 lightness-only |dJ'|
    IllResults(kIll).cam_d3chr  = cam_d3chr; % stage-3 chromatic-only hypot(da',db')

    IllResults(kIll).meanOffDiagRatio = meanOffDiagRatio;
    IllResults(kIll).semOffDiagRatio  = semOffDiagRatio;

    % NEW fields
    IllResults(kIll).rmsConeContrastErr = rmsConeContrastErr;  % nSD x nObs

    IllResults(kIll).uv_indiv_no = uv_indiv_no;
    IllResults(kIll).uv_indiv_ad = uv_indiv_ad;
    IllResults(kIll).uv_std_no   = uv_std_no;
    IllResults(kIll).uv_std_ad   = uv_std_ad;

    IllResults(kIll).mb_indiv_no = mb_indiv_no;
    IllResults(kIll).mb_indiv_ad = mb_indiv_ad;
    IllResults(kIll).mb_std_no   = mb_std_no;
    IllResults(kIll).mb_std_ad   = mb_std_ad;

    IllResults(kIll).srgb_D65    = srgb_D65_cell;
    
    IllResults(kIll).repObsIdx  = repObsIdx;
    IllResults(kIll).rep_RGB    = rep_RGB;
    IllResults(kIll).rep_XYZ_no = rep_XYZ_no;
    IllResults(kIll).rep_XYZ_ad = rep_XYZ_ad;

    IllResults(kIll).keepIdx      = keepIdx;          % this illuminant's kept surfaces
    IllResults(kIll).badUnion     = badUnionAll;
    IllResults(kIll).badThisIll   = badThisIll;      % optional
    IllResults(kIll).gamutSummary = gamut;           % optional

end

%% ============================================================
% 10) Save EVERYTHING (single MAT)
%% ============================================================
% Assemble one results struct: run metadata, parameters, the wavelength grid,
% the observer family and its metadata, the display family, all illuminants,
% the per-illuminant results (byIll), and gamut bookkeeping.
results = struct();

results.timestamp = string(datetime("now"));
results.dataDir   = dataDir;
results.outDir    = outDir;

results.params.sdList      = sdList(:);
results.params.peakBase    = peakBase(:);
results.params.ASANO       = ASANO;
results.params.N2          = N2;
results.params.N10         = N10;
results.params.CFG_GAMUT   = CFG_GAMUT;

results.wls = wls(:);
results.nW  = nW;

results.nObs     = nObs;
results.obs_meta = obs_meta;
results.asanoDir = asanoDir;

results.cone_fund = coneFam.fund;     % nObs x nW x 3

results.display.Pset       = Pset;
results.display.Mdisp_set  = Mdisp_set;
results.display.scale_set  = scale_set;
results.display.Ywhite0    = Ywhite0_set;

results.illuminants = Ill;
results.byIll       = IllResults;

if saveD
    results.D_set = D_set;
    results.offDiag_perObs = offDiag_perObs;
end

results.keepIdx_global = keepIdx_global;
results.badUnionAll    = badUnionAll;
results.gamutByIll     = gamutByIll;   % optional (can be big-ish but usually fine)


% '-struct' expands each top-level field of results into its own variable in
% the file; '-v7.3' (HDF5) is required because the data exceeds the 2 GB limit
% of the older MAT formats.
save(outMat, '-struct', 'results', '-v7.3');
fprintf('\nSaved results:\n  %s\n', outMat);

fprintf('\nALL DONE (CALC).\n');
toc;

%% ============================================================
% Helper functions
%% ============================================================

% Resample cone/CMF curves given as a PTB-style [start step nWl] descriptor S
% and a matrix T_in onto the target grid wls; returns 3 x numel(wls).
function T = resample_cmf_to_wls(S, T_in, wls)
if size(T_in,1) ~= 3 && size(T_in,2) == 3, T_in = T_in'; end
wl0 = S(1); step = S(2); n = S(3);
wl = (wl0:step:(wl0 + step*(n-1)))';
T = interp1(wl, T_in', wls, 'pchip', 0)'; % 3 x numel(wls)
end

% Build three unit-peak Gaussian display primaries (columns) at the given peak
% wavelengths and common standard deviation sd (nm), sampled on wls.
function P = build_gaussian_primaries(wls, peaks, sd)
P = zeros(numel(wls), 3);
for i = 1:3
    g = exp(-0.5*((wls - peaks(i))/sd).^2);
    if max(g) > 0, g = g./max(g); end
    P(:,i) = g;
end
end

% Batched 3x3 inverse over the third dimension: use pageinv when available,
% otherwise fall back to a per-page inv loop (same result).
function invA = pageinv_or_loop(A)
nObs = size(A,3);
if exist('pageinv','file')==2
    invA = pageinv(A);
else
    invA = zeros(3,3,nObs);
    for i=1:nObs
        invA(:,:,i) = inv(A(:,:,i));
    end
end
end

% White-patch von Kries adaptation in XYZ, per observer page.
% For each observer, scale X,Y,Z by the ratio (target white / that observer's
% reproduced white patch) so their white lands on XYZ_target_white; a simple
% diagonal (channel-independent) correction.
% XYZ_in: 3 x nSamples x nObs.
function XYZ_corr = vonkries_xyz_whitePatch(XYZ_in, XYZ_target_white, idxWhite)
nObs = size(XYZ_in,3);
XYZ_corr = XYZ_in;
FLOOR = 1e-6;
for iObs = 1:nObs
    w_obs = XYZ_in(:, idxWhite, iObs);
    gX = XYZ_target_white(1) / max(w_obs(1), FLOOR);
    gY = XYZ_target_white(2) / max(w_obs(2), FLOOR);
    gZ = XYZ_target_white(3) / max(w_obs(3), FLOOR);
    XYZ_corr(1,:,iObs) = gX * XYZ_in(1,:,iObs);
    XYZ_corr(2,:,iObs) = gY * XYZ_in(2,:,iObs);
    XYZ_corr(3,:,iObs) = gZ * XYZ_in(3,:,iObs);
end
end

% Same white-patch von Kries as above but for the single standard-observer
% reference (no observer dimension). Applied so the reference is adapted with
% the identical rule and the shared cast cancels.
function XYZ_out = vonkries_xyz_whitePatch_std(XYZ_in, XYZ_target_white, idxWhite)
% XYZ_in: 3 x nSamples
XYZ_out = XYZ_in;
FLOOR = 1e-6;
w_obs = XYZ_in(:, idxWhite);
gX = XYZ_target_white(1) / max(w_obs(1), FLOOR);
gY = XYZ_target_white(2) / max(w_obs(2), FLOOR);
gZ = XYZ_target_white(3) / max(w_obs(3), FLOOR);
XYZ_out(1,:) = gX * XYZ_in(1,:);
XYZ_out(2,:) = gY * XYZ_in(2,:);
XYZ_out(3,:) = gZ * XYZ_in(3,:);
end

% Standard CIE XYZ -> CIELAB for a 3 x nSamples input, using XYZ_white as the
% reference white. f is the CIELAB compression with the linear toe below
% (6/29)^3. Returns 3 x nSamples [L; a; b].
function Lab = xyz_to_lab(XYZ, XYZ_white)
X = XYZ(1,:) ./ XYZ_white(1);
Y = XYZ(2,:) ./ XYZ_white(2);
Z = XYZ(3,:) ./ XYZ_white(3);
f = @(t) ((t > (6/29)^3) .* t.^(1/3)) + ((t <= (6/29)^3) .* (t*(29/6)^2/3 + 4/29));
fx = f(X); fy = f(Y); fz = f(Z);
L = 116*fy - 16;
a = 500*(fx - fy);
b = 200*(fy - fz);
Lab = [L; a; b];
end

% Paged version of xyz_to_lab: XYZ is 3 x nSamples x nObs, returns 3 x nSamples
% x nObs. Identical formula; cat(1,...) keeps the observer dimension.
function Lab = xyz_to_lab_pages(XYZ, XYZ_white)
X = XYZ(1,:,:) ./ XYZ_white(1);
Y = XYZ(2,:,:) ./ XYZ_white(2);
Z = XYZ(3,:,:) ./ XYZ_white(3);
f = @(t) ((t > (6/29)^3) .* t.^(1/3)) + ((t <= (6/29)^3) .* (t*(29/6)^2/3 + 4/29));
fx = f(X); fy = f(Y); fz = f(Z);
L = 116*fy - 16;
a = 500*(fx - fy);
b = 200*(fy - fz);
Lab = cat(1, L, a, b);
end

% Scalar CIEDE2000 colour difference between two CIELAB triplets (the CIE 2000
% total colour-difference formula with kL=kC=kH=1). Lab1,Lab2 are 3-vectors.
function dE = ciede2000(Lab1, Lab2)
L1=Lab1(1); a1=Lab1(2); b1=Lab1(3);
L2=Lab2(1); a2=Lab2(2); b2=Lab2(3);
kL=1; kC=1; kH=1;

C1 = hypot(a1,b1);
C2 = hypot(a2,b2);
Cm = (C1 + C2)/2;

G = 0.5 * (1 - sqrt((Cm^7)/(Cm^7 + 25^7)));
a1p = (1+G)*a1;
a2p = (1+G)*a2;

C1p = hypot(a1p,b1);
C2p = hypot(a2p,b2);

h1p = atan2d(b1,a1p); if h1p<0, h1p=h1p+360; end
h2p = atan2d(b2,a2p); if h2p<0, h2p=h2p+360; end

dLp = L2 - L1;
dCp = C2p - C1p;

dhp = h2p - h1p;
if C1p*C2p == 0
    dhp = 0;
elseif dhp > 180
    dhp = dhp - 360;
elseif dhp < -180
    dhp = dhp + 360;
end

dHp = 2*sqrt(C1p*C2p) * sind(dhp/2);

Lpm = (L1 + L2)/2;
Cpm = (C1p + C2p)/2;

if C1p*C2p == 0
    hpm = h1p + h2p;
else
    if abs(h1p-h2p) > 180
        hpm = (h1p + h2p + 360)/2;
    else
        hpm = (h1p + h2p)/2;
    end
end

T = 1 - 0.17*cosd(hpm-30) + 0.24*cosd(2*hpm) + 0.32*cosd(3*hpm+6) - 0.20*cosd(4*hpm-63);

dRo = 30*exp(-((hpm-275)/25)^2);
Rc  = 2*sqrt((Cpm^7)/(Cpm^7 + 25^7));
Rt  = -Rc*sind(2*dRo);

Sl = 1 + (0.015*(Lpm-50)^2)/sqrt(20 + (Lpm-50)^2);
Sc = 1 + 0.045*Cpm;
Sh = 1 + 0.015*Cpm*T;

dE = sqrt( (dLp/(kL*Sl))^2 + (dCp/(kC*Sc))^2 + (dHp/(kH*Sh))^2 + Rt*(dCp/(kC*Sc))*(dHp/(kH*Sh)) );
end

% For each observer, compute CIEDE2000 vs the reference at every surface and
% return the median across surfaces. Lab_pages: 3 x nSamples x nObs;
% Lab_ref: 3 x nSamples. Output: 1 x nObs.
function medPerObs = median_de00_per_observer(Lab_pages, Lab_ref)
nObs     = size(Lab_pages,3);
nSamples = size(Lab_pages,2);
medPerObs = nan(1,nObs);
for iObs = 1:nObs
    dE = zeros(nSamples,1);
    for j = 1:nSamples
        dE(j) = ciede2000(Lab_pages(:,j,iObs), Lab_ref(:,j));
    end
    medPerObs(iObs) = median(dE,'omitnan');
end
end

% Off-diagonal energy ratio of a stack of 3x3 matrices M(:,:,o): fraction of
% total squared entries that lies OFF the main diagonal. 0 = purely diagonal
% (a gain change), toward 1 = strongly cross-coupling. Returns 1 x nObs.
function r = offdiag_energy_ratio(M)
sumAll  = squeeze(sum(sum(M.^2, 1), 2));
sumDiag = squeeze(M(1,1,:).^2 + M(2,2,:).^2 + M(3,3,:).^2);
r = (sumAll - sumDiag) ./ (sumAll + eps);
r = r(:)';
end

% Strict gamut precheck across ALL observers and ALL bandwidths for a fixed set
% of world spectra. A reflectance is "bad" at a bandwidth if ANY observer's
% metameric drives RGB = invC*LMS go below -tolNeg (out of the physical
% display gamut). Observers are processed in chunks to bound memory.
% Returns per-SD counts, the bad reflectance indices, and the min drive seen.
function out = gamut_precheck_all_SDs_fromPset(T_indiv, SPD_world, Pset, CFG)
if ~isfield(CFG,'chunkObs'), CFG.chunkObs = 2000; end
if ~isfield(CFG,'tolNeg'),   CFG.tolNeg   = 1e-10; end

nObs = size(T_indiv,3);
nSD  = numel(Pset);
nRef = size(SPD_world,2);

out.nBadReflectances  = zeros(nSD,1);
out.badReflectanceIdx = cell(nSD,1);
out.nBadPairs         = zeros(nSD,1);
out.minRGB            = nan(nSD,3);

usePageInv = (exist('pageinv','file')==2);

for si = 1:nSD
    P = Pset{si};

    badRefMask = false(nRef,1);
    badPairs   = 0;
    minRGB     = [inf inf inf];

    % Process observers in chunks of CFG.chunkObs to bound peak memory.
    for i1 = 1:CFG.chunkObs:nObs
        i2  = min(i1 + CFG.chunkObs - 1, nObs);
        idx = i1:i2;

        % For this observer chunk: LMS to every world spectrum, and the
        % cone-to-primary matrices whose inverse gives the drives.
        Tchunk   = T_indiv(:,:,idx);                 % 3 x nW x Nc
        LMSchunk = pagemtimes(Tchunk, SPD_world);    % 3 x nRef x Nc
        Cobs     = pagemtimes(Tchunk, P);            % 3 x 3 x Nc

        if usePageInv
            invC = pageinv(Cobs);
        else
            Nc = size(Cobs,3);
            invC = zeros(3,3,Nc);
            for k = 1:Nc
                invC(:,:,k) = inv(Cobs(:,:,k));
            end
        end

        % Metameric drives for every (surface, observer) in the chunk.
        RGB = pagemtimes(invC, LMSchunk);            % 3 x nRef x Nc

        % Track the most-negative drive per channel (diagnostic).
        minRGB(1) = min(minRGB(1), min(RGB(1,:,:), [], 'all'));
        minRGB(2) = min(minRGB(2), min(RGB(2,:,:), [], 'all'));
        minRGB(3) = min(minRGB(3), min(RGB(3,:,:), [], 'all'));

        % A (surface,observer) pair is out-of-gamut if any channel is negative;
        % a surface is bad if it is out-of-gamut for any observer.
        negMask = (RGB < -CFG.tolNeg);
        badSO   = squeeze(any(negMask,1));           % nRef x Nc

        badPairs   = badPairs + nnz(badSO);
        badRefMask = badRefMask | any(badSO,2);
    end

    out.nBadReflectances(si)  = nnz(badRefMask);
    out.badReflectanceIdx{si} = find(badRefMask);
    out.nBadPairs(si)         = badPairs;
    out.minRGB(si,:)          = minRGB;
end
end

% CIE xy chromaticity of a single XYZ tristimulus vector.
function xy = xyz_to_xy(XYZ)
den = max(eps, sum(XYZ));
xy  = [XYZ(1)/den, XYZ(2)/den];
end

% CIE 1976 u'v' chromaticity for a 3 x nSamples XYZ input -> nSamples x 2.
function uv = xyz_to_upvp(XYZ)
% XYZ: 3 x nSamples -> uv: nSamples x 2
X = XYZ(1,:); Y = XYZ(2,:); Z = XYZ(3,:);
den = (X + 15*Y + 3*Z);
den = max(den, eps);
up = 4*X ./ den;
vp = 9*Y ./ den;
uv = [up(:), vp(:)];
end

% Paged CIE u'v': XYZ is 3 x nSamples x nObs -> uv_pages nSamples x 2 x nObs.
function uv_pages = xyz_pages_to_upvp(XYZ)
% XYZ: 3 x nSamples x nObs -> uv_pages: nSamples x 2 x nObs
X = squeeze(XYZ(1,:,:));   % nSamples x nObs
Y = squeeze(XYZ(2,:,:));
Z = squeeze(XYZ(3,:,:));

den = (X + 15*Y + 3*Z);
den = max(den, eps);

up = 4*X ./ den;
vp = 9*Y ./ den;

nSamples = size(X,1);
nObs = size(X,2);

uv_pages = zeros(nSamples, 2, nObs);
uv_pages(:,1,:) = reshape(up, [nSamples 1 nObs]);
uv_pages(:,2,:) = reshape(vp, [nSamples 1 nObs]);
end

% XYZ (3 x n, white ~ Y=100) -> gamma-encoded sRGB (n x 3), assuming a D65
% white. Used only to colour figure markers, not in any metric.
function srgb = xyz_to_srgb_d65(XYZ)
%XYZ_TO_SRGB_D65 Convert XYZ (3 x n) to sRGB (n x 3), assuming D65.
% XYZ assumed roughly Y~100 for white (we normalize by /100).

assert(size(XYZ,1)==3, 'XYZ must be 3 x n.');

X = XYZ(1,:)/100;
Y = XYZ(2,:)/100;
Z = XYZ(3,:)/100;

M = [ 3.2406  -1.5372  -0.4986;
     -0.9689   1.8758   0.0415;
      0.0557  -0.2040   1.0570];

RGB_lin = M * [X; Y; Z];
RGB_lin = max(RGB_lin, 0);

a = 0.055;
thr = 0.0031308;

RGB = zeros(size(RGB_lin));
mask = (RGB_lin <= thr);
RGB(mask)  = 12.92 * RGB_lin(mask);
RGB(~mask) = (1+a) * (RGB_lin(~mask).^(1/2.4)) - a;

RGB = min(max(RGB,0),1);
srgb = RGB.';
end

% ============================
% NEW: MacLeod–Boynton helpers
% ============================

% MacLeod-Boynton chromaticity (r,b) from LMS (3 x nSamples -> nSamples x 2).
% Lw, Mw, Sw are the standard MB scaling weights that put L+M on a luminance
% footing; r = L/(L+M) (red-green) and b = S/(L+M) (yellow-blue).
function rb = lms_to_macleodboynton_rb(LMS)
% LMS: 3 x nSamples
% rb : nSamples x 2  (r,b)
Lw = 0.689903;
Mw = 0.348322;
Sw = 0.0371597/0.0192;

L = Lw*LMS(1,:); M = Mw*LMS(2,:); S = Sw*LMS(3,:);

den = max(L + M, 1e-12);
r = L ./ den;
b = S ./ den;
rb = [r(:), b(:)];
end

% Paged MacLeod-Boynton: LMS 3 x nSamples x nObs -> rb_pages nSamples x 2 x nObs.
% Same weights and formulas as lms_to_macleodboynton_rb.
function rb_pages = lms_pages_to_macleodboynton_rb(LMS)
% LMS: 3 x nSamples x nObs
% rb_pages: nSamples x 2 x nObs

Lw = 0.689903;
Mw = 0.348322;
Sw = 0.0371597/0.0192;

L = Lw*squeeze(LMS(1,:,:));   % nSamples x nObs
M = Mw*squeeze(LMS(2,:,:));
S = Sw*squeeze(LMS(3,:,:));

den = max(L + M, 1e-12);
r = L ./ den;
b = S ./ den;

nSamples = size(L,1);
nObs     = size(L,2);

rb_pages = zeros(nSamples, 2, nObs);
rb_pages(:,1,:) = reshape(r, [nSamples 1 nObs]);
rb_pages(:,2,:) = reshape(b, [nSamples 1 nObs]);
end

% White-patch von Kries directly in LMS (paged, per observer). Gains send each
% observer's white-patch LMS to [1 1 1]; used for the MB after-adaptation plot.
function LMS_corr = vonkries_lms_whitePatch(LMS_in, idxWhite)
% LMS_in: 3 x nSamples x nObs
% Apply per-observer von Kries gains in LMS so that white-patch LMS -> [1 1 1]
FLOOR = 1e-12;
nObs  = size(LMS_in,3);

LMS_corr = LMS_in;
for iObs = 1:nObs
    w = LMS_in(:, idxWhite, iObs);     % 3 x 1
    g = 1 ./ max(w, FLOOR);           % 3 x 1
    LMS_corr(1,:,iObs) = g(1) * LMS_in(1,:,iObs);
    LMS_corr(2,:,iObs) = g(2) * LMS_in(2,:,iObs);
    LMS_corr(3,:,iObs) = g(3) * LMS_in(3,:,iObs);
end
end

% Standard-observer (single, non-paged) white-patch von Kries in LMS.
function LMS_out = vonkries_lms_whitePatch_std(LMS_in, idxWhite)
% LMS_in: 3 x nSamples
FLOOR = 1e-12;
w = LMS_in(:, idxWhite);
g = 1 ./ max(w, FLOOR);
LMS_out = LMS_in;
LMS_out(1,:) = g(1) * LMS_in(1,:);
LMS_out(2,:) = g(2) * LMS_in(2,:);
LMS_out(3,:) = g(3) * LMS_in(3,:);
end

% ============================================================
% NEW helpers: optimal diagonal (XYZ) + CAT16 + vectorized CIEDE2000
% ============================================================

function [g, medv, meanv] = optimize_diag_xyz(XYZo, Lab_ref, XYZ_world, XYZ_white, idxWhite, objType, optSet)
% Find the diagonal gain g=(gX,gY,gZ) applied to an observer's reproduced
% XYZ that MINIMISES that observer's DeltaE00 (median or mean) against the
% standard-observer world reference. Warm-started at the white-patch gains
% and a least-squares diagonal, then polished with fminsearch in log-space.
FLOOR = 1e-6;

% warm-start candidate 1: white-patch von Kries gains
w_o  = XYZo(:, idxWhite);
g_wp = XYZ_white(:) ./ max(w_o, FLOOR);

% warm-start candidate 2: least-squares diagonal in XYZ (min ||diag(g)XYZo - XYZ_world||^2)
g_lsq = sum(XYZo .* XYZ_world, 2) ./ max(sum(XYZo.^2, 2), FLOOR);
g_lsq = max(g_lsq, FLOOR);

objfun = @(gv) diag_objective(gv, XYZo, Lab_ref, XYZ_white, objType);

f_wp  = objfun(g_wp);
f_lsq = objfun(g_lsq);
if f_lsq < f_wp
    g0 = g_lsq;  f0 = f_lsq;
else
    g0 = g_wp;   f0 = f_wp;
end

% optimise in log-space (keeps gains positive, scale-stable)
l0   = log(max(g0, FLOOR));
lopt = fminsearch(@(l) objfun(exp(l)), l0, optSet);
g_opt = exp(lopt(:));
fopt  = objfun(g_opt);

% guard: never accept a solution worse than the warm start
if fopt <= f0
    g = g_opt;
else
    g = g0;
end

[medv, meanv] = diag_both(g, XYZo, Lab_ref, XYZ_white);
end

% Objective for the XYZ optimal diagonal: apply gains gv to XYZo, convert to
% Lab, and return the median (or mean) CIEDE2000 vs the reference.
function f = diag_objective(gv, XYZo, Lab_ref, XYZ_white, objType)
gv   = gv(:);
XYZa = gv .* XYZo;                       % 3 x N (implicit expansion)
Lab  = xyz_to_lab(XYZa, XYZ_white);      % 3 x N
dE   = ciede2000_vec(Lab, Lab_ref);      % 1 x N
if strcmpi(objType, 'mean')
    f = mean(dE, 'omitnan');
else
    f = median(dE, 'omitnan');
end
end

% Return BOTH the median and mean CIEDE2000 for a given XYZ diagonal gain
% (used to report the companion statistic alongside the optimised objective).
function [medv, meanv] = diag_both(gv, XYZo, Lab_ref, XYZ_white)
gv   = gv(:);
XYZa = gv .* XYZo;
Lab  = xyz_to_lab(XYZa, XYZ_white);
dE   = ciede2000_vec(Lab, Lab_ref);
medv  = median(dE, 'omitnan');
meanv = mean(dE, 'omitnan');
end

function medv = optimize_diag_m16(XYZo, Lab_ref, XYZ_world, XYZ_white, idxWhite, M16, M16inv, objType, optSet)
% Optimal diagonal in the sharpened M16 basis: minimise that observer's
% DeltaE00 over a per-channel gain applied in M16 space,
%   XYZa = M16inv * diag(g) * M16 * XYZo.
% Warm-started at the CAT16 white-anchored gains (Rw_ref./Rw_o) and an M16
% least-squares diagonal, then polished with fminsearch in log-space.
% Returns the resulting median (objType='median') or mean DeltaE00.
FLOOR = 1e-6;

RGBo   = M16 * XYZo;                 % 3 x N (observer stimuli in M16)
RGBref = M16 * XYZ_world;            % 3 x N (standard-observer world in M16)

% warm start 1: CAT16 white-anchored gains (map observer white -> reference white)
Rw_o   = M16 * XYZo(:, idxWhite);
Rw_ref = M16 * XYZ_white(:);
g_wp   = Rw_ref ./ sign_guard(Rw_o, FLOOR);
% warm start 2: least-squares diagonal in M16
g_lsq  = sum(RGBo .* RGBref, 2) ./ max(sum(RGBo.^2, 2), FLOOR);

objfun = @(gv) diag_objective_m16(gv, RGBo, M16inv, Lab_ref, XYZ_white, objType);

f_wp = objfun(g_wp); f_lsq = objfun(g_lsq);
if f_lsq < f_wp, g0 = g_wp*0 + g_lsq; f0 = f_lsq; else, g0 = g_wp; f0 = f_wp; end

l0    = log(max(abs(g0), FLOOR));    % optimise |g| in log-space; sign fixed from g0
sgn   = sign(g0); sgn(sgn==0) = 1;
lopt  = fminsearch(@(l) objfun(sgn.*exp(l)), l0, optSet);
g_opt = sgn .* exp(lopt(:));
fopt  = objfun(g_opt);
if fopt <= f0, gsel = g_opt; else, gsel = g0; end

medv = diag_objective_m16(gsel, RGBo, M16inv, Lab_ref, XYZ_white, 'median');
if strcmpi(objType,'mean')
    medv = diag_objective_m16(gsel, RGBo, M16inv, Lab_ref, XYZ_white, 'mean');
end
end

% Objective for the M16 optimal diagonal: apply per-channel gains in the
% sharpened M16 space, map back to XYZ, and score median/mean CIEDE2000.
function f = diag_objective_m16(gv, RGBo, M16inv, Lab_ref, XYZ_white, objType)
gv   = gv(:);
XYZa = M16inv * (gv .* RGBo);            % apply diagonal in M16, back to XYZ
Lab  = xyz_to_lab(XYZa, XYZ_white);
dE   = ciede2000_vec(Lab, Lab_ref);
if strcmpi(objType,'mean'), f = mean(dE,'omitnan'); else, f = median(dE,'omitnan'); end
end

function medv = optimize_diag_power_m16(XYZo, Lab_ref, XYZ_white, idxWhite, M16, M16inv, objType, optSet)
% Optimal diagonal + single shared power exponent in the sharpened M16 basis.
% Minimal nonlinear generalisation of the M16 diagonal: after white-normalising
% each M16 cone-like response, apply a per-channel gain and one shared exponent,
%   r_i   = (M16*XYZ)_i / (M16*XYZ_white)_i         (reference white ~ 1)
%   r_c,i = g_i .* sign(r_i).*|r_i|.^p
%   XYZa  = M16inv * ( r_c .* (M16*XYZ_white) )
% At p=1 this reduces EXACTLY to a linear diagonal in M16, so the improvement
% over Opt. diag (M16) isolates the nonlinear (response-compression) gain.
% Warm-started at the CAT16 white-anchored gains with p=1, then polished with
% fminsearch over [log|g|(3); log p]. Returns the median (or mean) DeltaE00.
FLOOR = 1e-6;
RGBo   = M16 * XYZo;                          % 3 x N  (observer stimuli in M16)
Rw_ref = M16 * XYZ_white(:);                  % 3 x 1  (reference white in M16)
rn     = RGBo ./ sign_guard(Rw_ref, FLOOR);   % 3 x N  (white-normalised responses)

% warm start: white-anchored gains at p=1 (map observer white -> reference white)
Rw_o = M16 * XYZo(:, idxWhite);
g0   = Rw_ref ./ sign_guard(Rw_o, FLOOR);     % 3 x 1
sgn  = sign(g0); sgn(sgn==0) = 1;

objfun = @(th) diag_power_objective_m16(sgn.*exp(th(1:3)), th(4), rn, Rw_ref, M16inv, Lab_ref, XYZ_white, objType);
th0    = [log(max(abs(g0),FLOOR)); 0];        % p = 1 -> log p = 0
f0     = objfun(th0);
thop   = fminsearch(objfun, th0, optSet);
fop    = objfun(thop);
thsel  = th0; if fop <= f0, thsel = thop; end

gsel = sgn.*exp(thsel(1:3)); lpsel = thsel(4);
medv = diag_power_objective_m16(gsel, lpsel, rn, Rw_ref, M16inv, Lab_ref, XYZ_white, 'median');
if strcmpi(objType,'mean')
    medv = diag_power_objective_m16(gsel, lpsel, rn, Rw_ref, M16inv, Lab_ref, XYZ_white, 'mean');
end
end

% Objective for the M16 diagonal-plus-power variant: apply per-channel gains g
% and one shared exponent p to the white-normalised M16 responses rn, restore
% the reference white, map to XYZ, and score median/mean CIEDE2000.
function f = diag_power_objective_m16(g, lp, rn, Rw_ref, M16inv, Lab_ref, XYZ_white, objType)
g   = g(:);  p = exp(lp);                     % exponent kept positive
rc  = g .* (sign(rn).*abs(rn).^p);            % 3 x N  (gain + shared power on white-normalised responses)
XYZa= M16inv * (rc .* Rw_ref);                % de-normalise, back to XYZ
Lab = xyz_to_lab(XYZa, XYZ_white);
dE  = ciede2000_vec(Lab, Lab_ref);
if strcmpi(objType,'mean'), f = mean(dE,'omitnan'); else, f = median(dE,'omitnan'); end
end

function y = sign_guard(x, floorMag)
% keep sign of x but floor its magnitude away from zero
s = sign(x); s(s==0) = 1;
y = s .* max(abs(x), floorMag);
end

function [XYZ_c, factor_all] = cat16_adapt_pages(XYZ, XYZ_white, idxWhite, M16, M16inv, D)
% CIECAM16 CAT16 chromatic adaptation.
% Maps each observer's reproduced white -> the reference white in the
% sharpened M16 cone-like space, with degree of adaptation D in [0,1].
%   D=1 : complete adaptation (observer white -> reference white exactly).
%   D=0 : no adaptation (identity).
% XYZ:       3 x nSamples x nObs (observer reproduced stimuli, before adapt)
% XYZ_white: 3 x 1 reference white (standard observer world white)
nObs  = size(XYZ,3);
Yw    = 100;             % adopting-white luminance convention (CIECAM16)
FLOOR = 1e-9;

Rw_ref = M16 * XYZ_white(:);                       % 3 x 1
den    = D*Yw ./ max(Rw_ref, FLOOR) + (1 - D);     % 3 x 1

W_o  = reshape(XYZ(:, idxWhite, :), 3, nObs);      % 3 x nObs (observer whites)
Rw_o = M16 * W_o;                                  % 3 x nObs
num  = D*Yw ./ max(Rw_o, FLOOR) + (1 - D);         % 3 x nObs
factor_all = num ./ den;                           % 3 x nObs (implicit expansion)

% Transform all stimuli into M16, apply the per-observer diagonal gains
% (broadcast over surfaces via reshape), and map back to XYZ. pagemtimes
% applies the shared 3x3 M16 / M16inv to every observer page at once.
RGB   = pagemtimes(M16, XYZ);                       % 3 x nSamples x nObs
RGB_c = RGB .* reshape(factor_all, [3, 1, nObs]);
XYZ_c = pagemtimes(M16inv, RGB_c);                  % 3 x nSamples x nObs
end

function dE = ciede2000_vec(Lab1, Lab2)
% Vectorized CIEDE2000. Lab1, Lab2: 3 x N -> dE: 1 x N
L1=Lab1(1,:); a1=Lab1(2,:); b1=Lab1(3,:);
L2=Lab2(1,:); a2=Lab2(2,:); b2=Lab2(3,:);
kL=1; kC=1; kH=1;

C1 = hypot(a1,b1);
C2 = hypot(a2,b2);
Cm = (C1 + C2)/2;

G = 0.5 * (1 - sqrt((Cm.^7)./(Cm.^7 + 25^7)));
a1p = (1+G).*a1;
a2p = (1+G).*a2;

C1p = hypot(a1p,b1);
C2p = hypot(a2p,b2);
prodC = C1p.*C2p;

h1p = atan2d(b1,a1p); h1p(h1p<0) = h1p(h1p<0)+360;
h2p = atan2d(b2,a2p); h2p(h2p<0) = h2p(h2p<0)+360;

dLp = L2 - L1;
dCp = C2p - C1p;

dhp = h2p - h1p;
dhp(dhp >  180) = dhp(dhp >  180) - 360;
dhp(dhp < -180) = dhp(dhp < -180) + 360;
dhp(prodC == 0) = 0;

dHp = 2*sqrt(C1p.*C2p).*sind(dhp/2);

Lpm = (L1 + L2)/2;
Cpm = (C1p + C2p)/2;

% mean hue (matches scalar ciede2000 convention in this file)
hpm = (h1p + h2p)/2;
mBig = (prodC ~= 0) & (abs(h1p - h2p) > 180);
hpm(mBig) = (h1p(mBig) + h2p(mBig) + 360)/2;
mZero = (prodC == 0);
hpm(mZero) = h1p(mZero) + h2p(mZero);

T = 1 - 0.17*cosd(hpm-30) + 0.24*cosd(2*hpm) + 0.32*cosd(3*hpm+6) - 0.20*cosd(4*hpm-63);

dRo = 30*exp(-((hpm-275)/25).^2);
Rc  = 2*sqrt((Cpm.^7)./(Cpm.^7 + 25^7));
Rt  = -Rc.*sind(2*dRo);

Sl = 1 + (0.015*(Lpm-50).^2)./sqrt(20 + (Lpm-50).^2);
Sc = 1 + 0.045*Cpm;
Sh = 1 + 0.015*Cpm.*T;

dE = sqrt( (dLp./(kL*Sl)).^2 + (dCp./(kC*Sc)).^2 + (dHp./(kH*Sh)).^2 ...
    + Rt.*(dCp./(kC*Sc)).*(dHp./(kH*Sh)) );
end

% ============================================================
% NEW helpers: CIECAM16 multi-stage decomposition
% ============================================================

function DEs = ciecam16_stage_metamerism(XYZ_obs_pages, XYZ_world, XYZ_white, idxWhite, D, M16, M16inv, cam)
% Between-observer metamerism read out at three CIECAM16 stages.
% For a given degree of adaptation D:
%   - reference  = standard observer on the world spectra, adopting white = XYZ_white
%   - observer o = reproduced stimuli, adopting white = that observer's reproduced white
% Distance at each stage is Euclidean; median across surfaces per observer.
%   d1 : linear post-CAT16 adapted signals (M16 space)
%   d2 : post-nonlinear compressed cone responses (R'a G'a B'a)
%   d3 : CAM16-UCS [J' a' b']  (perceptual color difference)
%
% Controls that separate front-end chromatic adaptation from the appearance
% model's own white normalization:
%   d3c   : stage-3 distance when the OBSERVER is forced to adopt the COMMON
%           reference white (so no per-observer adaptation to the display white);
%           the achromatic Aw normalization is still present for both parties.
%   d3lig : lightness-only part of the own-white stage-3 distance, median |ΔJ'|.
%   d3chr : chromatic-only part of the own-white stage-3 distance, median hypot(Δa',Δb').
O = size(XYZ_obs_pages,3);

ref = ciecam16_stages(XYZ_world, XYZ_white, D, M16, M16inv, cam);  % s1,s2,s3 : 3 x nSamples

% own-white adaptation: each observer adopts that observer's reproduced white
Wown = XYZ_obs_pages(:, idxWhite, :);                             % 3 x 1 x O
st   = ciecam16_stages_paged(XYZ_obs_pages, Wown, D, M16, cam);   % s1,s2,s3 : 3 x N x O
d1   = squeeze(median(sqrt(sum((st.s1 - ref.s1).^2, 1)), 2, 'omitnan'));
d2   = squeeze(median(sqrt(sum((st.s2 - ref.s2).^2, 1)), 2, 'omitnan'));
ds3  = st.s3 - ref.s3;                                            % 3 x N x O
d3   = squeeze(median(sqrt(sum(ds3.^2, 1)), 2, 'omitnan'));
d3lig= squeeze(median(abs(ds3(1,:,:)), 2, 'omitnan'));
d3chr= squeeze(median(hypot(ds3(2,:,:), ds3(3,:,:)), 2, 'omitnan'));

% common-white control: all observers forced to adopt the reference white
Wcom = repmat(XYZ_white(:), 1, 1, O);                             % 3 x 1 x O
stc  = ciecam16_stages_paged(XYZ_obs_pages, Wcom, D, M16, cam);
d3c  = squeeze(median(sqrt(sum((stc.s3 - ref.s3).^2, 1)), 2, 'omitnan'));

DEs.d1 = d1(:)'; DEs.d2 = d2(:)'; DEs.d3 = d3(:)';
DEs.d3c = d3c(:)'; DEs.d3lig = d3lig(:)'; DEs.d3chr = d3chr(:)';
end

function out = ciecam16_stages_paged(XYZ, XYZ_w, D, M16, cam)
% Vectorized CIECAM16 forward model over observer pages (see ciecam16_stages).
% XYZ:   3 x N x O stimuli;  XYZ_w: 3 x 1 x O adopting whites.
% Returns out.s1 (linear post-CAT16), out.s2 (compressed), out.s3 (CAM16-UCS),
% each 3 x N x O. Numerically identical to a per-observer ciecam16_stages call.
LA = cam.LA; Yb = cam.Yb; c = cam.c; Nc = cam.Nc;
Yw = XYZ_w(2,1,:);                              % 1 x 1 x O

k   = 1/(5*LA + 1);
FL  = 0.2*k^4*(5*LA) + 0.1*(1 - k^4)^2*(5*LA)^(1/3);   % scalar
n   = max(Yb ./ Yw, 1e-6);                     % 1 x 1 x O
z   = 1.48 + sqrt(n);                          % 1 x 1 x O
Nbb = 0.725.*(1./n).^0.2;                       % 1 x 1 x O
Ncb = Nbb;

RGBw = pagemtimes(M16, XYZ_w);                  % 3 x 1 x O
Dr   = D.*Yw ./ RGBw + (1 - D);                % 3 x 1 x O
RGBcw= Dr .* RGBw;                              % 3 x 1 x O
tw   = (FL .* abs(RGBcw) / 100).^0.42;
RGBaw= sign(RGBcw) .* 400 .* tw ./ (27.13 + tw) + 0.1;   % 3 x 1 x O
Aw   = (2*RGBaw(1,1,:) + RGBaw(2,1,:) + RGBaw(3,1,:)/20 - 0.305) .* Nbb;   % 1 x 1 x O

RGB  = pagemtimes(M16, XYZ);                    % 3 x N x O
RGBc = Dr .* RGB;                               % 3 x N x O
s1   = RGBc;

t    = (FL .* abs(RGBc) / 100).^0.42;
RGBa = sign(RGBc) .* 400 .* t ./ (27.13 + t) + 0.1;      % 3 x N x O
s2   = RGBa;

R1 = RGBa(1,:,:); R2 = RGBa(2,:,:); R3 = RGBa(3,:,:);    % 1 x N x O each
a = R1 - 12*R2/11 + R3/11;
b = (R1 + R2 - 2*R3)/9;
A = (2*R1 + R2 + R3/20 - 0.305) .* Nbb;                  % 1 x N x O

J   = 100 .* sign(A) .* (abs(A)./Aw).^(c.*z);            % 1 x N x O
h   = atan2d(b, a);
et  = 0.25*(cos(deg2rad(h) + 2) + 3.8);
denC= R1 + R2 + (21/20)*R3;
tC  = (50000/13 * Nc .* Ncb) .* et .* hypot(a, b) ./ max(denC, 1e-9);
Cc  = (max(tC,0)).^0.9 .* sqrt(max(J,0)/100) .* (1.64 - 0.29.^n).^0.73;
M   = Cc * FL^0.25;

Jp  = (1 + 100*0.007) .* J ./ (1 + 0.007.*J);
Mp  = (1/0.0228) * log(1 + 0.0228.*M);
ap  = Mp .* cosd(h);
bp  = Mp .* sind(h);
s3  = cat(1, Jp, ap, bp);                                % 3 x N x O

out.s1 = s1; out.s2 = s2; out.s3 = s3;
end

function out = ciecam16_stages(XYZ, XYZ_w, D, M16, M16inv, cam)
% CIECAM16 forward model returning intermediate stage representations.
% XYZ:   3 x N stimuli (Y on a 0..100 scale, white ~ Y=100)
% XYZ_w: 3 x 1 adopting white
% D:     degree of adaptation in [0,1]
% Returns out.s1 (linear post-CAT16), out.s2 ([J a b]), out.s3 ([J' a' b']).
LA = cam.LA; Yb = cam.Yb; c = cam.c; Nc = cam.Nc;
Yw = XYZ_w(2);

% viewing-condition constants
k   = 1/(5*LA + 1);
FL  = 0.2*k^4*(5*LA) + 0.1*(1 - k^4)^2*(5*LA)^(1/3);
n   = max(Yb / Yw, 1e-6);
z   = 1.48 + sqrt(n);
Nbb = 0.725*(1/n)^0.2;
Ncb = Nbb;

% white in sharpened space + adaptation factors
RGBw = M16 * XYZ_w(:);                         % 3 x 1
Dr   = D*Yw ./ RGBw + (1 - D);                 % 3 x 1 (per-channel gain)

% adapted white -> post-nonlinear -> achromatic Aw
RGBcw = Dr .* RGBw;                            % 3 x 1
tw    = (FL .* abs(RGBcw) / 100).^0.42;
RGBaw = sign(RGBcw) .* 400 .* tw ./ (27.13 + tw) + 0.1;
Aw    = (2*RGBaw(1) + RGBaw(2) + RGBaw(3)/20 - 0.305) * Nbb;

% stimuli: CAT16 adaptation (stage 1 = linear adapted signals)
RGB  = M16 * XYZ;                              % 3 x N
RGBc = Dr .* RGB;                              % 3 x N (implicit expansion)
s1   = RGBc;

% post-adaptation nonlinear compression (stage 2 = compressed cone responses)
t    = (FL .* abs(RGBc) / 100).^0.42;
RGBa = sign(RGBc) .* 400 .* t ./ (27.13 + t) + 0.1;   % 3 x N
s2   = RGBa;

% opponent + achromatic
a = RGBa(1,:) - 12*RGBa(2,:)/11 + RGBa(3,:)/11;
b = (RGBa(1,:) + RGBa(2,:) - 2*RGBa(3,:))/9;
A = (2*RGBa(1,:) + RGBa(2,:) + RGBa(3,:)/20 - 0.305) * Nbb;

% lightness, hue, chroma, colorfulness
J   = 100 * sign(A) .* (abs(A)/Aw).^(c*z);
h   = atan2d(b, a);
et  = 0.25*(cos(deg2rad(h) + 2) + 3.8);
denC = RGBa(1,:) + RGBa(2,:) + (21/20)*RGBa(3,:);
tC  = (50000/13 * Nc * Ncb) .* et .* hypot(a, b) ./ max(denC, 1e-9);
Cc  = (max(tC,0)).^0.9 .* sqrt(max(J,0)/100) .* (1.64 - 0.29^n)^0.73;
M   = Cc * FL^0.25;

% CAM16-UCS
Jp  = (1 + 100*0.007) .* J ./ (1 + 0.007*J);
Mp  = (1/0.0228) * log(1 + 0.0228*M);
ap  = Mp .* cosd(h);
bp  = Mp .* sind(h);
s3  = [Jp; ap; bp];                            % CAM16-UCS

out.s1 = s1; out.s2 = s2; out.s3 = s3;
end

function verify_ciecam16_paged()
% Assert that the vectorized ciecam16_stages_paged reproduces the per-observer
% scalar ciecam16_stages to numerical precision, across degrees of adaptation.
rng(0);
N = 40; O = 25;
XYZ = 100*rand(3,N,O);
W   = 80 + 40*rand(3,1,O);
cam = struct('LA',20,'Yb',20,'F',1,'c',0.69,'Nc',1);
M16 = [ 0.401288  0.650173 -0.051461; ...
       -0.250268  1.204414  0.045854; ...
       -0.002079  0.048952  0.953127 ];
M16inv = inv(M16);
for D = [0 0.858 1.0]
    P = ciecam16_stages_paged(XYZ, W, D, M16, cam);
    e1=0; e2=0; e3=0;
    for o = 1:O
        Sc = ciecam16_stages(XYZ(:,:,o), W(:,1,o), D, M16, M16inv, cam);
        e1 = max(e1, max(abs(P.s1(:,:,o)-Sc.s1),[],'all'));
        e2 = max(e2, max(abs(P.s2(:,:,o)-Sc.s2),[],'all'));
        e3 = max(e3, max(abs(P.s3(:,:,o)-Sc.s3),[],'all'));
    end
    fprintf('  D=%.3f: max|ds1|=%.2e  max|ds2|=%.2e  max|ds3|=%.2e\n', D, e1, e2, e3);
    assert(e1<1e-9 && e2<1e-9 && e3<1e-9, 'paged CIECAM16 mismatch at D=%.3f', D);
end
fprintf('VERIFY_CAM PASSED: paged CIECAM16 == scalar CIECAM16\n');
end
