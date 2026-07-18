% ========================================================================
% compare_nonlinear_adapt_ucs.m
%
% Fair, common-ruler comparison of THREE adaptation models, all scored by the
% SAME CAM16-UCS distance:
%   M0  No-adapt      : observer's display XYZ, no correction
%   M2  CAT16 (linear): white-ratio diagonal gain in the M16 (sharpened) basis,
%                       mapping the observer's reproduced white -> reference white
%                       (this is von Kries in the sharpened cone basis)
%   M3  Non-linear    : the same white normalisation applied AFTER the
%                       Naka-Rushton compression (von Kries in the compressed
%                       response domain), then mapped back to XYZ
%   M5  Second-site   : white-anchored adaptation at the OPPONENT stage instead of
%                       the cone stage. The cone stage is anchored to the REFERENCE
%                       white (no per-observer cone gain, exactly as M0); the
%                       observer's own reproduced white is then SUBTRACTED in the
%                       post-compression opponent code (a,b), recentring chroma on
%                       that white. This is the classic second-site mechanism
%                       (Webster & Mollon).
%
% Why M3 is null but M5 need not be: in its working range the Naka-Rushton is
% essentially a power law, and a gain commutes with a power law -- (k*x)^p =
% k^p * x^p -- so putting the white gain before (M2) or after (M3) the
% compression is the SAME operation up to a constant. A SUBTRACTION after the
% power law does not commute with a gain before it, so M5 is a genuinely
% different mechanism rather than a reparameterisation of M2.
% Note a multiplicative gain is degenerate at the opponent stage: after the cone
% stage the adopting white sits at a = b = 0, so a white-anchored *gain* there
% would scale every chroma to zero. Subtraction is the meaningful operation.
%
% Common ruler U(x) = CAM16-UCS forward with adopting white = reference white,
% D=1. Reference (ground truth) = U(standard-observer world). Between-observer
% metamerism = median-over-surfaces UCS distance to the reference, per observer.
%
% Design identities (used for self-validation):
%   M0 == stored cam_d3c(D=1)   (common-white control)
%   M2 == stored cam_d3 (D=1)
% so the recomputation is checked against the validated 6.8 GB run without
% re-running it. Only M3 is genuinely new.
%
% The heavy N=10,000 analysis is NOT re-run: per-observer XYZ is reconstructed
% from the stored cone-fundamental family + display + world spectra.
%
% ------------------------------------------------------------------------
% CAM16-UCS FORWARD TRANSFORM (the "ruler", implemented in the helpers below)
%   Given a stimulus XYZ, an adopting white W, and a degree of adaptation D:
%   1. Sharpened cones: RGB = M16 * XYZ (M16 is the CAT16 sharpened basis).
%   2. Chromatic adaptation (CAT16): a per-channel von Kries gain
%        Dr = D*Yw/RGBw + (1-D)   ->   RGBc = Dr .* RGB
%      pulls the adopting white toward equal-energy; D=1 is full discounting,
%      D=0 is none.
%   3. Post-adaptation compression (Naka-Rushton / Hunt-Pointer-Estevez):
%        RGBa = sign*400*t/(27.13+t)+0.1,  t=(FL*|RGBc|/100)^0.42
%      a saturating nonlinearity giving the compressed cone responses.
%   4. Opponent + achromatic signals from RGBa: a = R-12G/11+B/11,
%        b = (R+G-2B)/9,  A = (2R+G+B/20-0.305)*Nbb.
%   5. Appearance correlates: lightness J from A/Aw, colourfulness M from a,b.
%   6. UCS compression: J' = 1.7 J /(1+0.007 J),  M' = ln(1+0.0228 M)/0.0228,
%      then the near-uniform Cartesian code  s3 = [J', M'cos h, M'sin h].
%   Between-observer distance is Euclidean in this [J';a';b'] space (=~ dE').
%
% WHY THE REFERENCE IS AT THE SAME D (the floor cancels)
%   At D<1 the shared adopting white is only partly discounted, so EVERY
%   observer's code carries the same global colour cast (the "adaptation-state
%   cast"). Encoding the standard-observer reference at the SAME D gives it the
%   identical cast; subtracting (taking distance) cancels this common mode. The
%   residual is purely the between-observer disagreement. The code checks this
%   explicitly: the reference-vs-itself distance ("floor") must be ~0 at each D.
%
% POOLING
%   Per (illuminant, SD, observer) we get one median-over-surfaces distance.
%   Illuminants are grouped into the three spectral classes; within a class the
%   per-observer distances are concatenated across illuminants (cat dim 2) and
%   the median over that pooled observer population is reported per SD and D.
%
% INPUTS (read from disk)
%   results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat
%       wls, params, display, cone_fund, illuminants, byIll (the main run)
%   data/CIE1931XYZ.txt          CIE 1931 2-deg CMFs (reference observer)
%   data/CIE_srf_cfi_1nm.csv     the 99 CIE test-colour reflectances
% OUTPUT
%   results/compare_nl_adapt_ucs.mat  (struct Res, CFG, classes) + printed
%       per-class tables of median metamerism vs D and the D0->D1 reduction.
%
% Run:  >> compare_nonlinear_adapt_ucs
% Output: results/compare_nl_adapt_ucs.mat  and printed reductions.
% ========================================================================
clearvars; close all; clc;

% ---- locate project root (three levels up from src/analysis/); fall back to a
% hard-coded path if the layout is unexpected, then add sources to the path. ---
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist(fullfile(projectRoot,'src'),'dir')
    projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
end
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
% the large main-run results file this script re-scores (must already exist)
inMat   = fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
assert(exist(inMat,'file')==2,'Missing results: %s',inMat);

% ---- config ------------------------------------------------------------
% SDsel selects which display-primary bandwidths (Gaussian SD in nm; smaller =
% narrower/peakier primaries) to evaluate; Dlist is the adaptation sweep.
CFG.SDsel    = [5 10 20];      % primary SDs to evaluate (intersected with sdList below)
CFG.Dlist    = [0 0.5 1];      % degree-of-adaptation sweep (own-white CAT16, full CIECAM16-UCS read-out)
CFG.testIll  = [];             % [] = all ok illuminants; else vector of indices (quick test)
CFG.maxIll   = [];             % cap #illuminants for a quick validation run ([]=all)
CFG.nayatani_beta = 0.75;      % M4: Nayatani-style power-law exponent (log+beta); beta=1 -> M2
CFG.validate = true;           % check M0/M2 against stored cam_d3c/cam_d3 at D=1
% CIECAM16 viewing condition (average surround): LA adapting luminance, Yb
% background luminance factor, c/Nc/F surround constants. Kept identical to the
% main run so the recomputation reproduces the stored values exactly.
cam = struct('LA',20,'Yb',20,'F',1.0,'c',0.69,'Nc',1.0);   % matches run_metamerism_analysis
% CAT16 sharpened-cone transform: XYZ -> sharpened RGB. Adaptation and the
% Naka-Rushton compression both act in this basis.
M16 = [ 0.401288  0.650173 -0.051461; ...
       -0.250268  1.204414  0.045854; ...
       -0.002079  0.048952  0.953127 ];
M16inv = inv(M16);

% ---- load (small top-level fields + byIll; keep only lean byIll fields) --
fprintf('Loading results (this loads the large .mat once)...\n');
S = load(inMat,'wls','params','display','cone_fund','illuminants','byIll');
wls    = S.wls(:);  nW = numel(wls);
sdList = S.params.sdList(:)';
% Pset{si}: nW x 3 display-primary SPDs at bandwidth si; Mdisp_set{si}: the
% display's RGB->XYZ matrix at that bandwidth.
Pset   = S.display.Pset;  Mdisp_set = S.display.Mdisp_set;
% cone_fund is stored as nObs x nW x 3; permute to 3 x nW x nObs so each page
% is one observer's LMS cone fundamentals (matrix-multiply friendly).
T_indiv = permute(S.cone_fund,[3 2 1]);            % 3 x nW x nObs
nObs   = size(T_indiv,3);
Ill    = S.illuminants;  nIll = numel(Ill);
% per-illuminant spectral class label ("daylight"/"smooth"/"peaky")
catVec = string({Ill.category});
% Copy only the small per-illuminant fields we need into a lean struct B, so the
% ~6.8 GB byIll can be dropped immediately afterwards. keepIdx = surviving
% surface indices, idxWhite = which surface is the white, XYZ_white = reference
% white XYZ. cam_d3/cam_d3c/cam_Dlist are kept only for the D=1 self-validation.
% lean copy of the per-illuminant fields we need, then drop the giant struct
B = struct('ok',{},'keepIdx',{},'idxWhite',{},'XYZ_white',{}, ...
           'cam_d3',{},'cam_d3c',{},'cam_Dlist',{});
for k = 1:nIll
    r = S.byIll(k);
    % ok flag: this illuminant produced usable results in the main run
    B(k).ok        = isfield(r,'ok') && ~isempty(r.ok) && r.ok;
    if ~B(k).ok, continue; end
    B(k).keepIdx   = r.keepIdx(:);
    B(k).idxWhite  = r.idxWhite;
    B(k).XYZ_white = r.XYZ_white(:);
    if CFG.validate && isfield(r,'cam_d3')
        B(k).cam_d3   = r.cam_d3;      % nSD x nObs x nD
        B(k).cam_d3c  = r.cam_d3c;
        B(k).cam_Dlist= r.cam_Dlist(:)';
    end
end
S = rmfield(S,'byIll'); clear r;   % free ~6.8 GB

% ---- T_xyz and reflectances on the wls grid ----------------------------
% Reference observer = CIE 1931 2-deg CMFs, resampled to the working grid (col 1
% is wavelength). This T_xyz encodes the ground-truth/standard observer.
Axyz  = readmatrix(fullfile(dataDir,'CIE1931XYZ.txt'));
T_xyz = interp1(Axyz(:,1),Axyz(:,2:4),wls,'pchip',0)';         % 3 x nW
% The 99 CIE test-colour reflectances, resampled and clamped to [0,1].
Rmat  = readmatrix(fullfile(dataDir,'CIE_srf_cfi_1nm.csv'));
R_all = min(max(interp1(Rmat(:,1),Rmat(:,2:end),wls,'pchip','extrap'),0),1);  % nW x 99

% ---- viewing-condition FL (scalar) -------------------------------------
% Luminance-level adaptation factor FL (a CIECAM16 constant). Computed once
% here for reuse; the helpers recompute the identical FL from cam.LA internally.
LA = cam.LA;  kF = 1/(5*LA+1);
FL = 0.2*kF^4*(5*LA) + 0.1*(1-kF^4)^2*(5*LA)^(1/3);

% ---- precompute invC per selected SD (illuminant-independent) -----------
% C = T_indiv * P is each observer's 3x3 cone-response-to-display-primary matrix
% (LMS produced by unit R,G,B). Its inverse invC turns a desired LMS into the
% display RGB that reproduces it for that observer -- the metameric match. It
% depends only on cones and primaries (not the world SPD), so cache it per SD.
CFG.SDsel = intersect(CFG.SDsel, sdList, 'stable');   % keep only SDs present in this run
nSDsel = numel(CFG.SDsel);
invC_set = cell(nSDsel,1);  siStored = zeros(nSDsel,1);
for j = 1:nSDsel
    % map each requested SD to its stored bandwidth index si
    si = find(sdList==CFG.SDsel(j),1);  siStored(j) = si;
    P  = Pset{si};                                   % nW x 3
    C  = pagemtimes(T_indiv, P);                     % 3 x 3 x nObs
    iC = zeros(3,3,nObs);
    % invert per observer (page); no batched inv over the 3rd dim in base MATLAB
    for o = 1:nObs, iC(:,:,o) = inv(C(:,:,o)); end
    invC_set{j} = iC;
end

% ---- illuminant list ----------------------------------------------------
% Evaluate only the illuminants flagged ok; CFG.testIll/CFG.maxIll optionally
% restrict this for quick validation runs.
okAll = find(arrayfun(@(k) B(k).ok, 1:nIll));
if ~isempty(CFG.testIll), okAll = intersect(okAll, CFG.testIll,'stable'); end
if ~isempty(CFG.maxIll), okAll = okAll(1:min(CFG.maxIll,numel(okAll))); end
fprintf('Evaluating %d illuminants x %d SDs x %d observers.\n', numel(okAll), nSDsel, nObs);

% storage: per illuminant, nSDsel x nObs distance-to-reference for each model
D0 = cell(nIll,1); D2 = cell(nIll,1); D3 = cell(nIll,1); D4 = cell(nIll,1); D4b = cell(nIll,1);
D5 = cell(nIll,1);
% D-sweep: own-white CAT16 at each degree of adaptation, full CIECAM16-UCS read-out
nD  = numel(CFG.Dlist);
Dsw = cell(nD,1); for iD=1:nD, Dsw{iD} = cell(nIll,1); end
% metamerism-free floor (pure adaptation-state cast) per illuminant x D
Bfloor    = nan(nIll,nD);
Bfloorcls = strings(nIll,1);
% ensures the (expensive) self-validation print runs on the first illuminant only
validatedOnce = false;

for k = okAll
    % per-illuminant bookkeeping: surviving surfaces, the white surface index,
    % and the reference (standard-observer) white XYZ.
    keepIdx  = B(k).keepIdx;
    idxWhite = B(k).idxWhite;
    Wref     = B(k).XYZ_white;                       % 3 x 1
    % world radiance = illuminant SPD times each surface reflectance, then the
    % standard-observer XYZ of every surface (the ground-truth colour signals).
    spd      = Ill(k).spd(:);
    SPD_world= spd .* R_all(:, keepIdx);             % nW x nSamp
    XYZ_world= T_xyz * SPD_world;                    % 3 x nSamp
    % per-model reference = standard observer's world encoded by THAT model
    % (adopting the reference white). Matches the paper's same-model reference
    % convention and removes any model self-bias, so each distance measures pure
    % between-observer spread. M0 and M2 share the linear reference.
    ref_s3   = cam16_s3(XYZ_world, Wref, 1, M16, cam);              % M0, M2
    ref_nl   = cam16_nl(XYZ_world, Wref, M16, cam);                 % M3
    ref_ny   = cam16_naya(XYZ_world, Wref, CFG.nayatani_beta, M16, cam);  % M4

    % D-matched references: the standard observer's world at each degree of
    % adaptation. Measuring the observer AND the reference at the SAME D cancels the
    % common-mode adaptation-state cast (the global colour shift from incomplete white
    % discounting), leaving only the between-observer metamerism at that D.
    ref_sw = cell(nD,1);
    for iD = 1:nD
        ref_sw{iD} = cam16_s3(XYZ_world, Wref, CFG.Dlist(iD), M16, cam);
        % floor = the reference stimulus vs itself at the same D -> must be ~0,
        % confirming the cast is removed (kept as a running sanity check).
        Bfloor(k,iD)  = median(dist_to_ref(ref_sw{iD}, ref_sw{iD}),'omitnan');
        Bfloorcls(k)  = catVec(k);
    end

    % each observer's actual cone excitations for every world surface
    LMS_world = pagemtimes(T_indiv, SPD_world);      % 3 x nSamp x nObs
    d0k = nan(nSDsel,nObs); d2k = nan(nSDsel,nObs); d3k = nan(nSDsel,nObs);
    d4k = nan(nSDsel,nObs); d4bk = nan(nSDsel,nObs);   % d4bk: beta-matched, white-correction OFF
    d5k = nan(nSDsel,nObs);                            % second-site (opponent) adaptation
    dDk = cell(nD,1); for iD=1:nD, dDk{iD}=nan(nSDsel,nObs); end   % D-sweep distances

    for j = 1:nSDsel
        Mdisp = Mdisp_set{siStored(j)};              % 3 x 3
        invC  = invC_set{j};                         % 3 x 3 x nObs
        % metameric match: display RGB that gives each observer the SAME LMS as
        % the world surface, then convert that RGB to XYZ (standard-observer XYZ
        % of the display match). This is what a given observer sees on screen.
        pRGB  = pagemtimes(invC, LMS_world);         % 3 x nSamp x nObs
        XYZ   = pagemtimes(Mdisp, pRGB);             % 3 x nSamp x nObs (no-adapt)
        % the observer's reproduced white (the display match to the white surface)
        Wown  = XYZ(:, idxWhite, :);                 % 3 x 1 x nObs

        s3_0 = cam16_s3(XYZ, Wref, 1, M16, cam);     % M0: adopt reference white (legacy no-adapt)
        s3_2 = cam16_s3(XYZ, Wown, 1, M16, cam);     % M2: linear CAT16 to own white (D=1)
        s3_3 = cam16_nl(XYZ, Wown, M16, cam);        % M3: same, white-normalised AFTER compression
        s3_4 = cam16_naya(XYZ, Wown, CFG.nayatani_beta, M16, cam);  % M4: Nayatani, white-corrected
        s3_4b= cam16_naya(XYZ, Wref, CFG.nayatani_beta, M16, cam);  % M4_0: same beta tone, NO white correction
        s3_5 = cam16_opp(XYZ, Wown, Wref, M16, cam);   % M5: second-site (opponent) adaptation

        d0k(j,:)  = dist_to_ref(s3_0,  ref_s3);
        d2k(j,:)  = dist_to_ref(s3_2,  ref_s3);
        d3k(j,:)  = dist_to_ref(s3_3,  ref_nl);
        d4k(j,:)  = dist_to_ref(s3_4,  ref_ny);
        d4bk(j,:) = dist_to_ref(s3_4b, ref_ny);
        % M5 shares M0/M2's reference: with Wown = Wref the opponent subtraction is
        % exactly zero, so cam16_opp reduces to cam16_s3(.,Wref,1) = the M0 encoding.
        d5k(j,:)  = dist_to_ref(s3_5,  ref_s3);

        % ---- degree-of-adaptation sweep: own-white CAT16 at D in {0,0.5,1},
        % full CIECAM16-UCS read-out, distance to the D-MATCHED reference so that
        % only between-observer metamerism (not the global cast) is measured.
        for iD = 1:nD
            s3_D = cam16_s3(XYZ, Wown, CFG.Dlist(iD), M16, cam);
            dDk{iD}(j,:) = dist_to_ref(s3_D, ref_sw{iD});
        end
    end
    % stash this illuminant's nSDsel x nObs distance matrices for every model
    D0{k}=d0k; D2{k}=d2k; D3{k}=d3k; D4{k}=d4k; D4b{k}=d4bk; D5{k}=d5k;
    for iD=1:nD, Dsw{iD}{k}=dDk{iD}; end

    % ---- self-validation on the first illuminant --------------------
    if CFG.validate && ~validatedOnce && ~isempty(B(k).cam_d3)
        iD1 = find(abs(B(k).cam_Dlist-1)<1e-9,1);
        e0=0; e2=0;
        for j=1:nSDsel
            st0 = B(k).cam_d3c(siStored(j),:,iD1);   % M0 == common-white D=1
            st2 = B(k).cam_d3 (siStored(j),:,iD1);   % M2 == own-white   D=1
            e0 = max(e0, max(abs(d0k(j,:)-st0)) / max(median(st0),eps));
            e2 = max(e2, max(abs(d2k(j,:)-st2)) / max(median(st2),eps));
        end
        fprintf('[validate ill %d] max rel err  M0 vs cam_d3c(D1)=%.2e   M2 vs cam_d3(D1)=%.2e\n',k,e0,e2);
        % M5 identity: feeding the REFERENCE white as the observer's own white makes
        % the opponent subtraction vanish, so M5 must reproduce the M0 encoding exactly.
        chk = cam16_opp(XYZ_world, Wref, Wref, M16, cam);
        fprintf('[validate ill %d] max abs err  M5(Wown=Wref) vs M0 encoding = %.2e\n', ...
                k, max(abs(chk(:)-ref_s3(:))));
        validatedOnce = true;
    end
    fprintf('  ill %d done.\n', k);
end

% ---- pool by illuminant class and report -------------------------------
% Group the evaluated illuminants into the three spectral classes and, within
% each class, POOL the per-observer distances across all illuminants of that
% class (cat along dim 2 = observers*illuminants) before taking medians. The
% pooled median per SD/model is the class's reported between-observer metamerism.
classes = intersect(["daylight","smooth","peaky"], unique(catVec),'stable');
Res = struct();
fprintf('\n=== Between-observer metamerism in CAM16-UCS (median over observers) ===\n');
for ci = 1:numel(classes)
    cls = classes(ci);
    % illuminants belonging to this class that were actually evaluated
    kk  = intersect(okAll, find(catVec==cls));
    if isempty(kk), continue; end
    % pool distances over this class: each p* is nSDsel x (nObs*nIll_in_class)
    p0=cat(2,D0{kk}); p2=cat(2,D2{kk}); p3=cat(2,D3{kk}); p4=cat(2,D4{kk}); p4b=cat(2,D4b{kk});
    p5=cat(2,D5{kk});
    % median over the pooled observer population, per SD (one value per row)
    m0=median(p0,2,'omitnan'); m2=median(p2,2,'omitnan'); m3=median(p3,2,'omitnan');
    m4=median(p4,2,'omitnan'); m4b=median(p4b,2,'omitnan'); m5=median(p5,2,'omitnan');
    Res.(cls).sd=CFG.SDsel(:); Res.(cls).m0=m0; Res.(cls).m2=m2; Res.(cls).m3=m3; Res.(cls).m4=m4; Res.(cls).m4b=m4b; Res.(cls).m5=m5;
    Res.(cls).p0=p0; Res.(cls).p2=p2; Res.(cls).p3=p3; Res.(cls).p4=p4; Res.(cls).p4b=p4b; Res.(cls).p5=p5;
    % ---- D-sweep pooled distributions + floor ----
    % Same pooling for the adaptation sweep: mD is nSDsel x nD, the median
    % pooled metamerism at each degree of adaptation D. This feeds the Fig-4 box
    % plots and the D0->D1 reduction reported below.
    Res.(cls).Dlist = CFG.Dlist(:)';
    pD = cell(1,nD); mD = zeros(nSDsel,nD);
    for iD=1:nD
        pD{iD} = cat(2, Dsw{iD}{kk});
        mD(:,iD)= median(pD{iD},2,'omitnan');
    end
    Res.(cls).pD = pD;  Res.(cls).mD = mD;
    % pooled floor per D: the reference-vs-itself distance, must be ~0 (confirms
    % the common adaptation-state cast cancels and only metamerism remains).
    Res.(cls).floor = median(Bfloor(kk,:),1,'omitnan');   % 1 x nD, must be ~0 (sanity)
    fprintf('\n[%s]  CIECAM16-UCS metamerism (median over observers), degree-of-adaptation sweep\n', cls);
    fprintf('  SD |'); for iD=1:nD, fprintf('  D=%-4g |', CFG.Dlist(iD)); end
    fprintf('  red(D1 vs D0)%%\n');
    for j=1:nSDsel
        fprintf('  %2g |', CFG.SDsel(j));
        for iD=1:nD, fprintf('  %5.3f  |', mD(j,iD)); end
        % reduction from no-adaptation (D0=first column) to full (D1=last column)
        rD = 100*(1 - mD(j,end)/mD(j,1));
        fprintf('     %5.1f\n', rD);
    end
    fprintf('  (floor / metamerism-free cast at each D: %s -- must be ~0)\n', ...
        mat2str(round(Res.(cls).floor,4)));
end
% Persist the pooled results (consumed by the Fig-4 box-plot script). -v7.3
% because Res can exceed the 2 GB limit of the older .mat format.
outMat = fullfile(projectRoot,'results','compare_nl_adapt_ucs.mat');
save(outMat,'Res','CFG','classes','-v7.3');
fprintf('\nSaved %s\n', outMat);

% ======================= local helpers =======================
function d = dist_to_ref(s3, ref_s3)
% median-over-surfaces Euclidean UCS distance to the reference, per observer.
% s3 is the observer's [J';a';b'] code (3 x N surfaces x O observers); ref_s3
% is the standard-observer reference (broadcast over the O dimension). The
% Euclidean norm over dim 1 is the CAM16-UCS colour difference dE' per surface;
% the median over surfaces summarises each observer as one number.
ds = s3 - ref_s3;                                    % 3 x N x O
d  = squeeze(median(sqrt(sum(ds.^2,1)), 2, 'omitnan'))';   % 1 x O
end

function s3 = cam16_s3(XYZ, W, D, M16, cam)
% Standard CIECAM16 forward -> CAM16-UCS [J';a';b'] (stage 3). 3 x N x O.
% XYZ: 3 x N x O ; W: 3 x 1 or 3 x 1 x O (per-page adopting white). D in [0,1].
% Faithful copy of ciecam16_stages_paged (returns s3 only). Reproduces the
% stored cam_d3 (own white) and cam_d3c (reference white) at any D.
LA=cam.LA; Yb=cam.Yb; c=cam.c; Nc=cam.Nc;
% Yw: luminance of the adopting white (per-page if W varies by observer).
Yw = W(2,1,:);
% FL: luminance-adaptation factor; n,z,Nbb,Ncb: background induction constants.
k  = 1/(5*LA+1);  FL = 0.2*k^4*(5*LA) + 0.1*(1-k^4)^2*(5*LA)^(1/3);
n  = max(Yb./Yw,1e-6); z=1.48+sqrt(n); Nbb=0.725.*(1./n).^0.2; Ncb=Nbb;
% CAT16 chromatic adaptation of the white: Dr is the von Kries gain that pulls
% the adopting white toward equal-energy (D=1 full, D=0 none); RGBcw is the
% adapted white in the sharpened basis.
RGBw = pagemtimes(M16,W);  Dr = D.*Yw./RGBw + (1-D);  RGBcw = Dr.*RGBw;
% Naka-Rushton compression of the adapted white -> its achromatic signal Aw
% (the reference against which stimulus lightness J is scaled).
tw   = (FL.*abs(RGBcw)/100).^0.42;  RGBaw = sign(RGBcw).*400.*tw./(27.13+tw)+0.1;
Aw   = (2*RGBaw(1,1,:)+RGBaw(2,1,:)+RGBaw(3,1,:)/20-0.305).*Nbb;
% Same adaptation (Dr) and compression applied to the stimulus.
RGB  = pagemtimes(M16,XYZ);  RGBc = Dr.*RGB;
t    = (FL.*abs(RGBc)/100).^0.42;  RGBa = sign(RGBc).*400.*t./(27.13+t)+0.1;
% Shared back end -> CAM16-UCS [J';a';b'].
s3   = ucs_from_compressed(RGBa, Aw, Nbb, Ncb, Nc, c, z, n, FL);
end

function s3 = cam16_nl(XYZ, W, M16, cam)
% Non-linear adaptation variant of cam16_s3: identical adopting white and
% viewing constants, but the white normalisation is applied in the Naka-Rushton
% COMPRESSED response domain (von Kries after compression) instead of the linear
% CAT16 gain before it. The white target is the same as linear CAT16 at D=1
% ([Yw;Yw;Yw] in M16), so cam16_nl and cam16_s3(.,W,1,..) differ ONLY in whether
% the observer's white is matched before or after the compressive nonlinearity.
LA=cam.LA; Yb=cam.Yb; c=cam.c; Nc=cam.Nc;
Yw = W(2,1,:);
k  = 1/(5*LA+1);  FL = 0.2*k^4*(5*LA) + 0.1*(1-k^4)^2*(5*LA)^(1/3);
n  = max(Yb./Yw,1e-6); z=1.48+sqrt(n); Nbb=0.725.*(1./n).^0.2; Ncb=Nbb;
RGBw = pagemtimes(M16,W);                 % observer's own reproduced white (M16)
aw   = nr_compress(RGBw, FL);             % its compressed response
Wtar = repmat(Yw,3,1,1);                  % CAT16 D=1 white target [Yw;Yw;Yw]
tgt  = nr_compress(Wtar, FL);
gain = tgt ./ aw;                         % compressed-domain von Kries gain (3 x 1 x O)
Aw   = (2*tgt(1,1,:)+tgt(2,1,:)+tgt(3,1,:)/20-0.305).*Nbb;   % same Aw as linear D=1
RGB  = pagemtimes(M16,XYZ);
RGBa = nr_compress(RGB, FL) .* gain;      % compress first, adapt after
s3   = ucs_from_compressed(RGBa, Aw, Nbb, Ncb, Nc, c, z, n, FL);
end

function s3 = cam16_naya(XYZ, W, beta, M16, cam)
% Nayatani-style power-law (log+beta) adaptation, read out through the same
% CAM16-UCS ruler. In the M16 basis the observer's white is mapped to the CAT16
% D=1 target [Yw;Yw;Yw] but with a nonlinear exponent on the deviation from white:
%   rho_hat_i = Yw * sign(rho_i) * (|rho_i| / rho_w_i)^beta
% beta=1 recovers linear CAT16 (M2); beta<1 is compressive (Nayatani-type
% incomplete adaptation). Equivalent to exp(beta*log(rho/rho_w)) i.e. the log+beta
% mechanism. The white maps to [Yw;Yw;Yw] for any beta, so Aw matches M2/M3.
LA=cam.LA; Yb=cam.Yb; c=cam.c; Nc=cam.Nc;
Yw = W(2,1,:);
k  = 1/(5*LA+1);  FL = 0.2*k^4*(5*LA) + 0.1*(1-k^4)^2*(5*LA)^(1/3);
n  = max(Yb./Yw,1e-6); z=1.48+sqrt(n); Nbb=0.725.*(1./n).^0.2; Ncb=Nbb;
RGBw = pagemtimes(M16,W);                                     % 3 x 1 x O
RGB  = pagemtimes(M16,XYZ);                                   % 3 x N x O
rho_hat = Yw .* sign(RGB) .* (abs(RGB)./RGBw).^beta;         % power-law adaptation
tgt  = repmat(Yw,3,1,1);                                      % adapted white = [Yw;Yw;Yw]
tw   = (FL.*abs(tgt)/100).^0.42;  RGBaw = sign(tgt).*400.*tw./(27.13+tw)+0.1;
Aw   = (2*RGBaw(1,1,:)+RGBaw(2,1,:)+RGBaw(3,1,:)/20-0.305).*Nbb;
t    = (FL.*abs(rho_hat)/100).^0.42;  RGBa = sign(rho_hat).*400.*t./(27.13+t)+0.1;
s3   = ucs_from_compressed(RGBa, Aw, Nbb, Ncb, Nc, c, z, n, FL);
end

function s3 = cam16_opp(XYZ, Wown, Wref, M16, cam)
% M5: second-site (opponent-stage) white-anchored adaptation.
% Cone stage is anchored to the REFERENCE white (identical to M0 -- no per-observer
% cone gain). The per-observer correction happens at the opponent stage instead:
% the observer's own reproduced white is pushed through the same front end and its
% opponent response (aw,bw) is SUBTRACTED from the stimulus, recentring the
% chromatic channels on that white. The achromatic normalisation (Aw, Yw) is left
% at the reference, so this isolates second-site CHROMATIC adaptation.
% With Wown = Wref the subtraction is exactly zero and this reduces to M0.
LA=cam.LA; Yb=cam.Yb; c=cam.c; Nc=cam.Nc;
Yw = Wref(2,1,:);
k  = 1/(5*LA+1);  FL = 0.2*k^4*(5*LA) + 0.1*(1-k^4)^2*(5*LA)^(1/3);
n  = max(Yb./Yw,1e-6); z=1.48+sqrt(n); Nbb=0.725.*(1./n).^0.2; Ncb=Nbb;
RGBwr = pagemtimes(M16,Wref);  Dr = Yw./RGBwr;          % CAT16 to the reference white, D=1
RGBaw = nr_compress(Dr.*RGBwr, FL);
Aw    = (2*RGBaw(1,1,:)+RGBaw(2,1,:)+RGBaw(3,1,:)/20-0.305).*Nbb;
RGBa  = nr_compress(Dr.*pagemtimes(M16,XYZ),  FL);      % stimulus
RGBaW = nr_compress(Dr.*pagemtimes(M16,Wown), FL);      % observer's own white, same front end
s3    = ucs_from_compressed(RGBa, Aw, Nbb, Ncb, Nc, c, z, n, FL, RGBaW);
end

function s3 = ucs_from_compressed(RGBa, Aw, Nbb, Ncb, Nc, c, z, n, FL, RGBaW)
% Shared CAM16 back end: compressed cone responses -> CAM16-UCS [J';a';b'].
% Optional RGBaW: compressed responses of the observer's own white. When given,
% its opponent response is subtracted from (a,b) -- the M5 second-site stage.
% Opponent signals a (red-green) and b (yellow-blue) from the compressed cones.
R1=RGBa(1,:,:); R2=RGBa(2,:,:); R3=RGBa(3,:,:);
a = R1 - 12*R2/11 + R3/11;
b = (R1 + R2 - 2*R3)/9;
% M5 second-site stage: subtract the observer's own white opponent code so the
% chromatic channels are recentred on that white (only when RGBaW is supplied).
if nargin >= 10 && ~isempty(RGBaW)
    W1=RGBaW(1,:,:); W2=RGBaW(2,:,:); W3=RGBaW(3,:,:);
    a = a - (W1 - 12*W2/11 + W3/11);
    b = b - (W1 + W2 - 2*W3)/9;
end
% Achromatic signal A, then lightness J as A relative to the white's Aw.
A = (2*R1 + R2 + R3/20 - 0.305).*Nbb;
J = 100.*sign(A).*(abs(A)./Aw).^(c.*z);
% Hue angle h and its eccentricity factor et; temporary quantity t drives
% chroma C, converted to colourfulness M.
h = atan2d(b,a);  et = 0.25*(cos(deg2rad(h)+2)+3.8);
denC = R1 + R2 + (21/20)*R3;
tC = (50000/13*Nc.*Ncb).*et.*hypot(a,b)./max(denC,1e-9);
Cc = (max(tC,0)).^0.9 .* sqrt(max(J,0)/100) .* (1.64-0.29.^n).^0.73;
Mc = Cc.*FL^0.25;
% CAM16-UCS uniformity compression of lightness (J') and colourfulness (M'),
% then the Cartesian code [J', a'=M'cos h, b'=M'sin h] used for the dE' metric.
Jp = (1+100*0.007).*J./(1+0.007.*J);
Mp = (1/0.0228)*log(1+0.0228.*Mc);
s3 = cat(1, Jp, Mp.*cosd(h), Mp.*sind(h));
end

function a = nr_compress(v, FL)
% Naka-Rushton / Hunt-Pointer-Estevez post-adaptation compression (the same
% saturating nonlinearity used inline in cam16_s3), factored out so the
% nonlinear/opponent variants can apply it explicitly.
t = (FL.*abs(v)/100).^0.42;
a = sign(v).*400.*t./(27.13+t) + 0.1;
end
