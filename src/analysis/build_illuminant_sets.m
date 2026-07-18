function build_illuminant_sets(cfg)
% build_illuminant_sets  Construct three CCT-matched illuminant classes used
% by the observer-metamerism analysis, so the classes differ in spectral
% *peakiness* (smoothness), not white point:
%
%   (1) DAYLIGHT         : M CIE daylight phases (smooth by construction),
%                          CCTs spaced evenly in mireds over [cctMin,cctMax].
%   (2) SMOOTH ARTIFICIAL: M real, spectrally SMOOTH engineered sources.
%   (3) PEAKY  ARTIFICIAL: M real, spectrally PEAKY engineered sources.
%
% The two artificial classes are drawn from the Houser/TM-30 401-illuminant
% set. Peakiness is the RMS spectral gradient of the max-normalised SPD,
%   pk = sqrt(mean(diff(spd/max(spd)).^2)),
% which is ~0 for smooth daylight/Planckian spectra and large for spiky
% multi-LED sources. Candidates are pooled within the common CCT band, split
% at the median peakiness into a smooth half and a peaky half, and one source
% is chosen per daylight CCT target from each half (CCT-stratified), taking the
% smoothest (resp. peakiest) candidate near each target. No colour-rendition
% (Rf/Rg) quality screen is applied -- the point is to span smooth to peaky --
% except that the CIE daylight/Planckian *reference* illuminants embedded in the
% 401 set (Rf >= refRf, i.e. essentially non-manufactured references identical
% to the daylight class) are excluded from the artificial pool.
%
% All SPDs are returned on the 380:5:780 nm grid; the analysis pipeline
% interpolates them onto its working grid.
%
% INPUTS (all optional; supplied as fields of the struct 'cfg')
%   cfg.M       number of illuminants PER class (default 12)
%   cfg.cctMin  low  end of the shared CCT band, in kelvin (default 4000)
%   cfg.cctMax  high end of the shared CCT band, in kelvin (default 7000)
%   cfg.refRf   Rf threshold at/above which a 401-set source is treated as an
%               embedded reference (daylight/Planckian) and excluded from the
%               artificial candidate pool (default 99)
%
% OUTPUT
%   data/illuminant_sets.mat  containing the struct IllSet with fields:
%     .wl                 wavelength grid (81x1, 380:5:780 nm)
%     .cfg                the resolved configuration used
%     .cctTargets         Mx1 CCT targets (K), shared by all three classes
%     .metaCols           column labels for the per-source .meta matrices
%     .daylight/.smooth/.peaky  each a struct with:
%         .spd   81 x M matrix of SPDs
%         .meta  M x 5 matrix [CCT Duv Rf Rg peakiness]
%         .srcIdx (artificial classes only) index into the 401 set
%
% HOW TO RUN
%   >> build_illuminant_sets                 % defaults
%   >> build_illuminant_sets(struct('M',20)) % override any field(s)
%
% Output: data/illuminant_sets.mat  (struct IllSet)

% ---- resolve configuration: start from defaults, then overlay any fields
% the caller actually provided (so a partial cfg still gets all defaults). ----
if nargin < 1, cfg = struct(); end
def = struct('M',12, 'cctMin',4000, 'cctMax',7000, 'refRf',99);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(cfg, fn{i}), cfg.(fn{i}) = def.(fn{i}); end
end

% ---- locate the project root (three levels up from src/analysis/) and put
% the source tree + data folder on the path so helper functions and data
% files resolve regardless of the caller's working directory. ----
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(projectRoot,'src')));
addpath(fullfile(projectRoot,'data'));
dataDir = fullfile(projectRoot,'data');

% Load the TM-30 machinery (CMFs, daylight/Planckian basis, reference data)
% and take its wavelength grid as the master grid for everything below.
D  = load_tm30_data();
wl = D.wl;                                   % 81x1, 380:5:780

% CCT targets are spaced EVENLY IN MIREDS (reciprocal megakelvin), not in
% kelvin, because perceptually uniform CCT steps are uniform in mireds. The
% same targets are reused by all three classes so the classes are CCT-matched.
% Common CCT targets (even in mireds), shared by all three classes
miredTgt = linspace(1e6/cfg.cctMin, 1e6/cfg.cctMax, cfg.M);   % warm -> cool
cctTgt   = 1e6 ./ miredTgt;

%% ---------------- (1) DAYLIGHT phases at the target CCTs ----------------
% Class 1 is synthesised directly: for each CCT target reconstruct the CIE
% daylight SPD from the daylight eigenvector basis (S0,S1,S2). These are smooth
% by construction and serve as the "reference" smooth-spectrum class.
fprintf('Building %d CIE daylight phases at CCT targets %g-%g K...\n', ...
    cfg.M, cfg.cctMin, cfg.cctMax);
dCols = zeros(numel(wl), cfg.M);
dMeta = zeros(cfg.M, 5);                      % [CCT Duv Rf Rg peakiness]
for i = 1:cfg.M
    % reconstruct the daylight phase at this CCT, clamp any tiny negatives to 0
    v = max(cie_daylight_spd(cctTgt(i), wl, D.dayBasis), 0);
    % colorimetry / rendition metrics for this SPD (CCT, Duv, TM-30 Rf, Rg)
    r = compute_tm30(v, wl, D);
    dCols(:,i) = v;
    % store metrics plus the spectral peakiness measure defined below
    dMeta(i,:) = [r.CCT r.Duv r.Rf r.Rg spd_peakiness(v)];
end
fprintf('  daylight CCT [%.0f..%.0f] K, peakiness [%.3f..%.3f]\n', ...
    min(dMeta(:,1)), max(dMeta(:,1)), min(dMeta(:,5)), max(dMeta(:,5)));

%% ---------------- characterise the 401 set ----------------
% Load Houser's 401 real-world illuminant SPDs (column 1 is wavelength, the
% rest are the sources) and verify they sit on the same 380:5:780 grid.
S401 = load(fullfile(dataDir,'houser401.mat'),'i401');
i401 = S401.i401;
assert(max(abs(i401(:,1)-wl))<1e-9, '401 set must be on the 380:5:780 grid.');
spd401 = i401(:,2:end);                        % 81 x 401
nSrc = size(spd401,2);
% Compute CCT, Rf, Rg and peakiness for every source so we can later filter
% and split them. Preallocate as NaN so any all-zero/invalid source stays NaN.
mCCT = nan(nSrc,1); mRf = nan(nSrc,1); mRg = nan(nSrc,1); mPk = nan(nSrc,1);
for j = 1:nSrc
    % skip degenerate columns (no positive power)
    v = max(spd401(:,j),0); if ~any(v>0), continue; end
    r = compute_tm30(v, wl, D);
    mCCT(j)=r.CCT; mRf(j)=r.Rf; mRg(j)=r.Rg; mPk(j)=spd_peakiness(v);
end

% Candidate pool = sources with a valid CCT inside the shared band AND with
% Rf below refRf. The Rf<refRf test drops the CIE daylight/Planckian reference
% illuminants embedded in the 401 set (they would duplicate the daylight class).
% candidate pool: in CCT band, valid, and NOT an embedded reference illuminant
inBand = isfinite(mCCT) & mCCT>=cfg.cctMin & mCCT<=cfg.cctMax & mRf < cfg.refRf;
poolIdx = find(inBand);
fprintf('Artificial candidate pool in %g-%g K (refs Rf>=%g excluded): %d sources\n', ...
    cfg.cctMin, cfg.cctMax, cfg.refRf, numel(poolIdx));

% Split the pool by its OWN median peakiness: the smoother half becomes the
% "smooth artificial" class, the peakier half the "peaky artificial" class.
% Using the median guarantees a balanced two-way split of the real sources.
% split the pool at median peakiness -> smooth half / peaky half
medPk = median(mPk(poolIdx));
smoothPool = poolIdx(mPk(poolIdx) <= medPk);
peakyPool  = poolIdx(mPk(poolIdx) >  medPk);
fprintf('  median peakiness = %.3f; smooth pool n=%d, peaky pool n=%d\n', ...
    medPk, numel(smoothPool), numel(peakyPool));

%% ---------------- CCT-stratified selection from each half ----------------
% From each half pick one source per CCT target so the artificial classes end
% up CCT-matched to the daylight class. Within the smooth half take the
% SMOOTHEST candidate near each target; within the peaky half the PEAKIEST.
% This maximises the separation in peakiness between the two artificial classes
% while holding CCT fixed.
% smooth artificial: nearest-CCT to each target, pick the SMOOTHEST in-window
selSmooth = select_stratified(smoothPool, mCCT, mPk, miredTgt, 'min');
% peaky artificial: nearest-CCT to each target, pick the PEAKIEST in-window
selPeaky  = select_stratified(peakyPool,  mCCT, mPk, miredTgt, 'max');

% Gather the selected SPDs and their metrics. Duv (column 2) is left NaN here
% and filled in the loop below; the others come straight from the 401-set scan.
sCols = spd401(:,selSmooth); sMeta = [mCCT(selSmooth) nan(cfg.M,1) mRf(selSmooth) mRg(selSmooth) mPk(selSmooth)];
pCols = spd401(:,selPeaky);  pMeta = [mCCT(selPeaky)  nan(cfg.M,1) mRf(selPeaky)  mRg(selPeaky)  mPk(selPeaky)];
% Recompute TM-30 on the chosen sources only to obtain their Duv values.
% fill Duv
for i=1:cfg.M
    rs=compute_tm30(sCols(:,i),wl,D); sMeta(i,2)=rs.Duv;
    rp=compute_tm30(pCols(:,i),wl,D); pMeta(i,2)=rp.Duv;
end
fprintf('  SMOOTH artificial: CCT [%.0f..%.0f] K, Rf [%.0f..%.0f], peakiness [%.3f..%.3f]\n', ...
    min(sMeta(:,1)),max(sMeta(:,1)),min(sMeta(:,3)),max(sMeta(:,3)),min(sMeta(:,5)),max(sMeta(:,5)));
fprintf('  PEAKY  artificial: CCT [%.0f..%.0f] K, Rf [%.0f..%.0f], peakiness [%.3f..%.3f]\n', ...
    min(pMeta(:,1)),max(pMeta(:,1)),min(pMeta(:,3)),max(pMeta(:,3)),min(pMeta(:,5)),max(pMeta(:,5)));

%% ---------------- assemble + save ----------------
% Pack the three CCT-matched classes plus shared metadata into one struct.
% metaCols documents the column order of every .meta matrix; srcIdx records
% which 401-set columns the artificial sources came from (for traceability).
IllSet = struct();
IllSet.wl         = wl;
IllSet.cfg        = cfg;
IllSet.cctTargets = cctTgt(:);
IllSet.metaCols   = {'CCT','Duv','Rf','Rg','peakiness'};
IllSet.daylight.spd  = dCols;  IllSet.daylight.meta  = dMeta;
IllSet.smooth.spd    = sCols;  IllSet.smooth.meta    = sMeta;  IllSet.smooth.srcIdx = selSmooth(:);
IllSet.peaky.spd     = pCols;  IllSet.peaky.meta     = pMeta;  IllSet.peaky.srcIdx  = selPeaky(:);

% Write to data/illuminant_sets.mat (v7 keeps the file small and portable).
outFile = fullfile(dataDir,'illuminant_sets.mat');
save(outFile,'IllSet','-v7');
fprintf('Saved %s (%d daylight + %d smooth + %d peaky, CCT-matched).\n', ...
    outFile, cfg.M, cfg.M, cfg.M);
end

% ======================================================================
function pk = spd_peakiness(v)
% RMS spectral gradient of the max-normalised SPD (0 = smooth, large = peaky)
% Normalising by the peak makes the measure amplitude-independent (it depends
% on spectral SHAPE only); diff() approximates the wavelength derivative, so
% the RMS of diff() is small for slowly varying spectra and large for spiky
% multi-LED ones. This single scalar is what the smooth/peaky split is based on.
v = max(v(:),0); v = v ./ max(v);
pk = sqrt(mean(diff(v).^2));
end

function sel = select_stratified(pool, mCCT, mPk, miredTgt, extremum)
% Choose one source per CCT target from 'pool', CCT-stratified. Within a mired
% window around each target take the most extreme peakiness (min=smoothest,
% max=peakiest); fall back to the nearest-CCT candidate if the window is empty.
% Targets are processed sparsest-neighbourhood first so scarce CCTs get priority.
%   pool     : indices (into the 401 set) of the candidate half to draw from
%   mCCT,mPk : per-source CCT and peakiness for the whole 401 set
%   miredTgt : 1 x M mired targets (matches the daylight class)
%   extremum : 'min' -> smoothest, 'max' -> peakiest within the window
% Returns M x 1 selected 401-set indices, one per target, WITHOUT repeats.
M = numel(miredTgt);
% Work in mireds so "near in CCT" is symmetric; window is 3/4 of the target
% spacing, i.e. a bit narrower than half-way to the neighbouring target.
poolMired = 1e6 ./ mCCT(pool);
win = 0.75 * mean(abs(diff(miredTgt)));
% Count how many candidates fall within each target's window, then handle the
% targets with the FEWEST neighbours first: scarce CCTs claim their best source
% before commoner targets exhaust the pool (selection is greedy/without-repeat).
nearCnt = zeros(M,1);
for k=1:M, nearCnt(k) = sum(abs(poolMired - miredTgt(k)) <= win); end
[~,order] = sort(nearCnt,'ascend');
sel = nan(M,1);
for oo = 1:M
    k = order(oo);
    % candidates still unclaimed (exclude anything already selected)
    avail = find(~ismember(pool, sel(~isnan(sel))));   % indices into pool
    dM = abs(poolMired(avail) - miredTgt(k));
    inWin = avail(dM <= win);
    if isempty(inWin)
        % nothing in-window: take the single nearest-CCT candidate instead
        [~,mi] = min(dM); pick = avail(mi);
    else
        % in-window: pick the peakiness extremum asked for (smoothest/peakiest)
        pkw = mPk(pool(inWin));
        if strcmp(extremum,'min'), [~,mi]=min(pkw); else, [~,mi]=max(pkw); end
        pick = inWin(mi);
    end
    % map the pool-relative pick back to a 401-set index and record it
    sel(k) = pool(pick);
end
sel = sel(:);
end
