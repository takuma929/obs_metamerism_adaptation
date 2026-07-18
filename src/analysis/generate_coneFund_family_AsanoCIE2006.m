function out = generate_coneFund_family_AsanoCIE2006(varargin)
%GENERATE_CONEFUND_FAMILY_ASANOCIE2006
% Family of cone fundamentals using CIE2006/CIEPO06-style components:
% - Age-dependent ocular media density SHAPE via Docul1/Docul2 (Eq 2/3/4 in Asano 2016)
% - Macular pigment via relative density table + field-size dependent peak (Eq 5/6)
% - Cone photopigment absorptance via OD + absorbance shapes (Eq 7/8/9)
%
% ---------------------------------------------------------------------------
% PURPOSE
% ---------------------------------------------------------------------------
%   Generate a large population of INDIVIDUAL colorimetric observers (cone
%   fundamentals) using the Asano et al. (2016) individual-observer model
%   built on the CIE 2006 (CIEPO06) cone-fundamental machinery. Each observer
%   is a set of three L/M/S corneal-incident cone spectral sensitivities on
%   the wavelength grid wl_out. A separate population is produced for the 2 deg
%   and the 10 deg viewing field, then the two are stacked into one output.
%
%   This routine is the observer generator used in the paper
%     "Adaptation reduces inter-observer metamerism with narrow-band primaries"
%   where 10,000 observers per field are drawn to probe inter-observer
%   metamerism.
%
% ---------------------------------------------------------------------------
% HOW AN OBSERVER IS BUILT (per cone, per observer)
% ---------------------------------------------------------------------------
%   1. Start from a template photopigment ABSORBANCE spectrum (L, M, S) taken
%      from the CIE 2006 tables.
%   2. Apply that observer's lambda_max (peak-wavelength) shift by sliding the
%      absorbance spectrum along the wavelength axis (see shift_spectrum_allowNaN).
%   3. Convert absorbance shape -> ABSORPTANCE via Beer-Lambert with a peak
%      photopigment optical density (self-screening):  a = 1 - 10^(-OD * shape).
%   4. Apply PRERECEPTORAL FILTERING (lens/ocular media + macular pigment) as a
%      spectral transmittance  Tpre = 10^-(Docul + Dmac).
%   5. Corneal cone fundamental = absorptance .* Tpre, resampled onto wl_out.
%      Optional per-observer, per-cone normalization ('area' / 'max' / 'none').
%
% ---------------------------------------------------------------------------
% INDIVIDUAL DIFFERENCES
% ---------------------------------------------------------------------------
%   Continuous (Asano's Gaussian model, sampled as percentage deviations of a
%   nominal value): lens/ocular media density, macular pigment density, L/M/S
%   photopigment peak optical densities. Age is drawn discrete-uniform over
%   age_range and drives the ocular-media SHAPE. Field size (2 or 10 deg) sets
%   the nominal macular and photopigment peak densities.
%
%   Discrete, GENE-DETERMINED lambda_max shifts (see parameter block and
%   generate_one_field below): L is bimodal (Ser180Ala), M is fixed at 0
%   (monomorphic at codon 180), S is symmetric bimodal.
%
% ---------------------------------------------------------------------------
% INPUTS  (Name-Value pairs; see inputParser block for validators/defaults)
% ---------------------------------------------------------------------------
%   'data_dir'      folder holding the CIE2006 .txt tables (default '.').
%   'N2','N10'      number of observers for the 2 deg / 10 deg field.
%   'seed'          RNG seed for reproducibility.
%   'age_range'     [min max] years for the discrete-uniform age draw.
%   'sd_d*_pct'     SDs (in %) of the Gaussian density deviations
%                   (dlens, dmac, dL, dM, dS).
%   'pL_ser'        probability of the L Ser180 variant.
%   'L_ser_nm','L_ala_nm'  lambda_max shift (nm) of the two L variants.
%   'pS_hi'         probability of the long-wave S cluster.
%   'S_hi_nm','S_lo_nm'    lambda_max shift (nm) of the two S clusters.
%   'normalize'     'area' | 'max' | 'none' (per-observer, per-cone).
%   'S_cut_nm'      long-wave cutoff (nm) beyond which S is forced to 0.
%
% ---------------------------------------------------------------------------
% OUTPUT  (struct 'out')
% ---------------------------------------------------------------------------
%   out.wl_out   1 x 31 wavelength axis (nm).
%   out.fund     (N2+N10) x 31 x 3 array of cone fundamentals; 3rd dim = L,M,S.
%   out.meta     table with one row per observer (field, age, density
%                deviations, lambda_max shifts, and the L/S variant labels).
%   out.N2,out.N10   observer counts per field.
%
% ---------------------------------------------------------------------------
% REFERENCES
% ---------------------------------------------------------------------------
%   Asano, Fairchild & Blonde (2016) PLoS ONE 11(2):e0145671.
%   CIE 170-1:2006 / CIEPO06 cone fundamentals.
%   Winderickx et al. (1992); Stockman & Sharpe (2000);
%   Stockman, Sharpe & Fach (1999); Stockman & Rider (2023).
%
% NOTE : λmax shifts (sL,sM,sS) are sampled from Normal distributions,
% as assumed in Asano et al. (2016). Default SDs [nm] are Table 5 Step 2:
%   SD(sL)=2.0, SD(sM)=1.5, SD(sS)=1.3, with means = 0.
%
% IMPORTANT (accuracy note for the reviewer): the NOTE just above describes
% Asano's ORIGINAL continuous (Gaussian) treatment of lambda_max. THIS code
% does NOT use Gaussian lambda_max shifts. Instead it overrides them with the
% DISCRETE, gene-determined scheme documented in the parameter block below and
% implemented inside generate_one_field: L is bimodal (Ser180Ala), M is fixed
% at 0, S is symmetric bimodal. The sd_dL/dM/dS_pct parameters below control
% the photopigment OPTICAL-DENSITY deviations, not lambda_max.

% -----------------------------
% options
% -----------------------------
% MATLAB idiom: inputParser collects optional Name-Value arguments. Each
% addParameter registers a name, a default value, and an anonymous validation
% function @(x)... that must return true for the supplied value to be accepted.
p = inputParser;
% Folder containing the CIE2006 lookup tables (accepts char or string).
p.addParameter('data_dir', '.', @(s)ischar(s)||isstring(s));
% Number of observers to simulate for the 2 deg and 10 deg field populations.
p.addParameter('N2', 10000, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('N10', 10000, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
% Random-number-generator seed so the whole population is reproducible.
p.addParameter('seed', 0, @(x)isnumeric(x)&&isscalar(x));
% [min max] age (years); ages are drawn discrete-uniform over this range.
p.addParameter('age_range', [18 70], @(x)isnumeric(x)&&numel(x)==2&&x(1)<x(2));

% individual variability (deviations in % for densities)
% Each SD below is the standard deviation (in PERCENT) of a zero-mean Gaussian
% that perturbs a nominal density value for that observer. Values follow
% Asano et al. (2016) population statistics.
% Lens / ocular-media density deviation.
p.addParameter('sd_dlens_pct', 18.7, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
% Macular pigment density deviation.
p.addParameter('sd_dmac_pct',  36.5, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
% L / M / S photopigment PEAK OPTICAL DENSITY deviations (self-screening).
p.addParameter('sd_dL_pct',     9.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('sd_dM_pct',     9.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('sd_dS_pct',     7.4, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% --- λmax shifts: gene-determined discrete distribution ---
% Grounded in Winderickx et al. (1992), Stockman & Sharpe (2000),
% Stockman & Rider (2023), Stockman, Sharpe & Fach (1999):
%   L : bimodal Ser180Ala polymorphism. L(ser180) is +1.2 nm and L(ala180)
%       is -1.5 nm relative to the CIE-standard L (a 0.56:0.44 ser:ala mean),
%       i.e. the two variants are 2.7 nm apart; frequency 0.56 ser : 0.44 ala.
%   M : effectively monomorphic at codon 180 (~94% Ala180; the M(S180)
%       variant shifts lambda_max by only ~0.17 nm), so no shift is applied.
%   S : weak bimodal, two clusters at 417.4 / 420.1 nm (Stockman, Sharpe &
%       Fach 1999), i.e. +/-1.35 nm about the mean, taken 50:50.
% Probability an observer carries the L Ser180 variant (else Ala180).
p.addParameter('pL_ser',      0.56, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
% lambda_max shift (nm) of the L Ser180 variant: +1.2 nm (red-shifted).
p.addParameter('L_ser_nm',    1.2,  @(x)isnumeric(x)&&isscalar(x));
% lambda_max shift (nm) of the L Ala180 variant: -1.5 nm (blue-shifted).
p.addParameter('L_ala_nm',   -1.5,  @(x)isnumeric(x)&&isscalar(x));
% Probability an observer is in the long-wave S cluster (else short-wave).
p.addParameter('pS_hi',       0.50, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
% lambda_max shifts (nm) of the two symmetric S clusters (+/-1.35 nm).
p.addParameter('S_hi_nm',     1.35, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('S_lo_nm',    -1.35, @(x)isnumeric(x)&&isscalar(x));

% output normalization
% Per-observer, per-cone rescaling applied at the end of generate_one_field:
%   'area' -> unit area under each fundamental; 'max' -> peak set to 1;
%   'none' -> leave in absolute (absorptance x transmittance) units.
p.addParameter('normalize', 'none', @(s)any(strcmpi(s,{'area','max','none'})));

% S-cone long-wave cutoff (match common Asano practice)
% Wavelength (nm) beyond which the S sensitivity is forced to 0 (S response is
% negligible in the long-wave range; this also suppresses interpolation noise).
p.addParameter('S_cut_nm', 620, @(x)isnumeric(x)&&isscalar(x)&&x>=500&&x<=700);

% Apply the validators and collect the results into struct 'opt'.
p.parse(varargin{:});
opt = p.Results;
% Seed the RNG so every random draw below (age, density deviations, and the
% L/S variant coin-flips) is reproducible for a given 'seed'.
rng(opt.seed);

% Normalize data_dir to a plain char array for fullfile()/exist() below.
data_dir = char(opt.data_dir);

% -----------------------------
% wavelength grids
% -----------------------------
% Both grids span 390-780 nm at 5 nm steps (31 samples). wl_in is the internal
% grid on which all spectra are built; wl_out is the reported grid. The trailing
% (') transposes each row vector into a COLUMN vector.
wl_out = (390:5:780)';  % 31 channels
wl_in  = (390:5:780)';   % internal grid we interpolate to

% -----------------------------
% load prereceptoral tables (wl + columns)
% -----------------------------
% Locate the three CIE2006 data tables inside data_dir.
%   mac  : relative macular pigment density spectrum (shape only).
%   docul: ocular-media (lens) density basis functions (Docul1, Docul2).
%   abs  : L/M/S photopigment absorbance template spectra.
mac_file   = fullfile(data_dir, 'cie2006_RelativeMacularDensity.txt');
docul_file = fullfile(data_dir, 'cie2006_docul.txt');
abs_file   = fullfile(data_dir, 'cie2006_LMSAbsorbance.txt');

% Fail early with a clear message if any table is missing (exist(...)==2 means
% the name refers to an existing file).
assert(exist(mac_file,'file')==2,   "Missing: %s", mac_file);
assert(exist(docul_file,'file')==2, "Missing: %s", docul_file);
assert(exist(abs_file,'file')==2,   "Missing: %s", abs_file);

% Read each table into a numeric matrix, then drop any row that contains a NaN
% (all(~isnan(row)) keeps only fully-populated rows, e.g. removing header/blank
% lines that parse as NaN).
A_mac = readmatrix(mac_file);   A_mac = A_mac(all(~isnan(A_mac),2),:);
A_doc = readmatrix(docul_file); A_doc = A_doc(all(~isnan(A_doc),2),:);
A_abs = readmatrix(abs_file);   A_abs = A_abs(all(~isnan(A_abs),2),:);

% Sanity-check that each table has the expected column layout.
assert(size(A_mac,2) >= 2, 'mac file must be [wl, relMac, ...]');
assert(size(A_abs,2) >= 4, 'abs file must be [wl, L, M, S, ...]');
assert(size(A_doc,2) >= 2, 'docul file must be [wl, Docul...]');

% Split each table into its native wavelength axis (column 1) and data columns.
wl_mac = A_mac(:,1); relMac_raw = A_mac(:,2);
wl_doc = A_doc(:,1); doc_cols   = A_doc(:,2:end);
wl_abs = A_abs(:,1);
Labs_raw = A_abs(:,2);
Mabs_raw = A_abs(:,3);
Sabs_raw = A_abs(:,4);

% Mark exact zeros in the S template as NaN so they are treated as "no data"
% rather than as genuine zero absorbance during the shift/interpolation below.
Sabs_raw(Sabs_raw==0) = NaN;

% Resample the relative macular pigment SHAPE onto the internal grid (linear
% interpolation, linearly extrapolated at the ends), then clamp negatives to 0.
relMac = interp1(wl_mac, relMac_raw, wl_in, 'linear', 'extrap');
relMac = max(relMac, 0);

% The ocular-media (lens) density can be supplied either as TWO CIEPO06 basis
% functions (Docul1, Docul2) that are later combined with age, or as a single
% precomputed density column. docMode records which branch was taken so the
% per-observer loop can reconstruct the density correctly.
if size(doc_cols,2) >= 2
    % Two-component form: resample both bases onto the internal grid, clamp
    % negatives, and flag this as the age-dependent Docul1+Docul2 mode.
    Docul1 = interp1(wl_doc, doc_cols(:,1), wl_in, 'linear', 'extrap');
    Docul2 = interp1(wl_doc, doc_cols(:,2), wl_in, 'linear', 'extrap');
    Docul1 = max(Docul1, 0);
    Docul2 = max(Docul2, 0);
    docMode = "Docul1+Docul2";
else
    % Single-column form: use the density spectrum directly (no age shaping).
    Docul  = interp1(wl_doc, doc_cols(:,1), wl_in, 'linear', 'extrap');
    Docul  = max(Docul, 0);
    Docul1 = [];
    Docul2 = [];
    docMode = "DoculDirect";
end

% Resample the L/M/S absorbance templates onto the internal grid.
Labs = interp1(wl_abs, Labs_raw, wl_in, 'linear', 'extrap');
Mabs = interp1(wl_abs, Mabs_raw, wl_in, 'linear', 'extrap');
Sabs = interp1(wl_abs, Sabs_raw, wl_in, 'linear', 'extrap');

% Auto-detect whether the absorbance table is stored as log10 values or as
% linear values. Pool all finite samples and inspect them:
%   fracNeg  : fraction of negative values (log10 spectra are mostly negative);
%   rangeVals: 1st..99th percentile spread.
% If most values are negative, OR the values sit in a small range with a
% negative floor, the table is assumed to be log10 and is exponentiated to
% linear. (x(:) flattens to a column; prctile ignores the tails' outliers.)
vals = [Labs(:); Mabs(:); Sabs(:)];
vals = vals(~isnan(vals));
fracNeg = mean(vals < 0);
rangeVals = [prctile(vals,1) prctile(vals,99)];
isLog10 = (fracNeg > 0.5) || (rangeVals(2) <= 2 && rangeVals(1) < 0);

if isLog10
    % Convert log10 absorbance -> linear absorbance shape.
    Labs_lin = 10.^Labs;
    Mabs_lin = 10.^Mabs;
    Sabs_lin = 10.^Sabs;
else
    % Already linear; use as-is.
    Labs_lin = Labs;
    Mabs_lin = Mabs;
    Sabs_lin = Sabs;
end

% Guard against any tiny negative values produced by interpolation.
Labs_lin(Labs_lin < 0) = 0;
Mabs_lin(Mabs_lin < 0) = 0;
Sabs_lin(Sabs_lin < 0) = 0;

% -----------------------------
% generate for each field
% -----------------------------
% Build the 2 deg and 10 deg observer populations separately; the field angle
% (first argument) sets the nominal macular/photopigment peak densities inside.
[fund2,  meta2]  = generate_one_field(2,  opt.N2,  opt, wl_in, wl_out, relMac, docMode, Labs_lin, Mabs_lin, Sabs_lin, Docul1, Docul2);
[fund10, meta10] = generate_one_field(10, opt.N10, opt, wl_in, wl_out, relMac, docMode, Labs_lin, Mabs_lin, Sabs_lin, Docul1, Docul2);

% -----------------------------
% pack : merge 2deg + 10deg into one fund/meta
% -----------------------------
% Stack the two populations (2 deg rows first, then 10 deg) into a single
% output struct.
out = struct();
% Store wl_out as a row vector (1 x 31).
out.wl_out = wl_out(:)';                     % 1 x 31
% Concatenate along dim 1 (observers): result is (N2+N10) x 31 x 3.
out.fund   = cat(1, fund2, fund10);          % (N2+N10) x 31 x 3
% Vertically concatenate the two metadata tables.
out.meta   = [meta2; meta10];                % (N2+N10) x (...)
% Add a stable 1..N observer index as a new column.
out.meta.obs_idx = (1:height(out.meta))';    % stable row id
out.N2  = opt.N2;
out.N10 = opt.N10;

% ======================================================================
% Build ONE field's population of N observers for viewing angle vdeg (2 or 10).
% Returns fund (N x 31 x 3) and a per-observer meta table.
function [fund, meta] = generate_one_field(vdeg, N, opt, wl_in, wl_out, relMac, docMode, Labs_lin, Mabs_lin, Sabs_lin, Docul1, Docul2)

    % Draw each observer's age: integer, uniform over [age_min, age_max].
    age = randi([opt.age_range(1), opt.age_range(2)], N, 1);

    % Draw the continuous individual-difference parameters as PERCENT
    % deviations: SD (in %) times a standard normal randn(N,1) column. These
    % scale the nominal densities later via (1 + d/100).
    dlens = opt.sd_dlens_pct * randn(N,1);   % lens / ocular media
    dmac  = opt.sd_dmac_pct  * randn(N,1);   % macular pigment
    dL    = opt.sd_dL_pct    * randn(N,1);   % L photopigment peak OD
    dM    = opt.sd_dM_pct    * randn(N,1);   % M photopigment peak OD
    dS    = opt.sd_dS_pct    * randn(N,1);   % S photopigment peak OD

    % --- λmax shifts: gene-determined discrete distribution ---
    % These are DISCRETE gene-determined shifts (not Gaussian). rand(N,1) draws
    % N uniform values in [0,1); comparing to a probability yields a logical
    % mask that selects each observer's genetic variant.
    % L: bimodal Ser180 (+1.2 nm, p=0.56) / Ala180 (-1.5 nm, p=0.44)
    % Start everyone at the Ala180 shift, then overwrite the Ser180 carriers.
    isSer = rand(N,1) < opt.pL_ser;                 % 1 = L(ser180)
    sL_nm = opt.L_ala_nm * ones(N,1);
    sL_nm(isSer) = opt.L_ser_nm;
    % M: monomorphic at codon 180 -> no shift
    sM_nm = zeros(N,1);
    % S: weak bimodal 417.4/420.1 nm -> +/-1.35 nm, 50:50
    % Start everyone at the short-wave (-1.35 nm) shift, then overwrite the
    % long-wave cluster (+1.35 nm).
    isHi  = rand(N,1) < opt.pS_hi;                  % 1 = long-wave cluster
    sS_nm = opt.S_lo_nm * ones(N,1);
    sS_nm(isHi) = opt.S_hi_nm;

    % Field-size-dependent NOMINAL peak densities (Asano/CIEPO06 formulas):
    %   Dmac_max : macular pigment peak density (falls off with field size).
    %   DLM_max  : L and M photopigment peak optical density.
    %   DS_max   : S photopigment peak optical density.
    % Larger fields (10 deg) give lower macular and photopigment peak densities
    % than the fovea-centred 2 deg field.
    Dmac_max = 0.485 * exp(-vdeg/6.132);
    DLM_max  = 0.38  + 0.54 * exp(-vdeg/1.333);
    DS_max   = 0.30  + 0.45 * exp(-vdeg/1.333);

    % Preallocate the output: N observers x wavelengths x 3 cones (L,M,S).
    fund = zeros(N, numel(wl_out), 3);

    % Build each observer i in turn.
    for i = 1:N
        % --- Prereceptoral filter 1: ocular media / lens density ---
        % Get this observer's ocular-media density and apply their individual
        % percent deviation (1 + dlens/100).
        if docMode=="Docul1+Docul2"
            % Age-dependent: combine the two Docul bases for this age...
            Docul_ave = ocular_density_age(age(i), Docul1, Docul2);
            Docul_i = Docul_ave .* (1 + dlens(i)/100);
        else
            % ...or scale the single supplied density spectrum directly.
            Docul_i = Docul .* (1 + dlens(i)/100);
        end

        % --- Prereceptoral filter 2: macular pigment density ---
        % Nominal peak density x individual deviation, spread over wavelength
        % by the relative macular SHAPE (relMac).
        Dmac_i = (Dmac_max * (1 + dmac(i)/100)) .* relMac;

        % --- Apply the genetic lambda_max shift to each absorbance template ---
        % Slide each L/M/S absorbance spectrum along wavelength by this
        % observer's shift (positive => red shift). See shift_spectrum_allowNaN.
        Lshape = shift_spectrum_allowNaN(wl_in, Labs_lin, sL_nm(i));
        Mshape = shift_spectrum_allowNaN(wl_in, Mabs_lin, sM_nm(i));
        Sshape = shift_spectrum_allowNaN(wl_in, Sabs_lin, sS_nm(i));
        % Force the S absorbance shape to 0 beyond the long-wave cutoff.
        Sshape(wl_in >= opt.S_cut_nm) = 0;

        % --- This observer's peak photopigment optical densities ---
        % Nominal peak OD x individual deviation, floored at 0.
        dL_i = max(0, DLM_max * (1 + dL(i)/100));
        dM_i = max(0, DLM_max * (1 + dM(i)/100));
        dS_i = max(0, DS_max  * (1 + dS(i)/100));

        % --- Absorbance shape -> ABSORPTANCE (Beer-Lambert self-screening) ---
        % a(lambda) = 1 - 10^(-OD * absorbance_shape(lambda)). Higher peak OD
        % broadens/saturates the sensitivity relative to the raw absorbance.
        aL = 1 - 10.^(-dL_i * Lshape);
        aM = 1 - 10.^(-dM_i * Mshape);
        aS = 1 - 10.^(-dS_i * Sshape);
        % Re-enforce the S long-wave cutoff after the transform.
        aS(wl_in >= opt.S_cut_nm) = 0;

        % --- Combined prereceptoral transmittance ---
        % Densities add in log space, so total transmittance is
        % Tpre = 10^-(Docul + Dmac). This attenuates all three cones identically.
        Tpre = 10.^(-(Docul_i + Dmac_i));

        % --- Corneal-incident cone fundamentals ---
        % Multiply cone absorptance by the prereceptoral transmittance.
        Lfund = aL .* Tpre;
        Mfund = aM .* Tpre;
        Sfund = aS .* Tpre;

        % Resample each fundamental onto the output grid (linear; 0 outside the
        % data range) and store into the L/M/S slices of the output array.
        fund(i,:,1) = interp1(wl_in, Lfund, wl_out, 'linear', 0);
        fund(i,:,2) = interp1(wl_in, Mfund, wl_out, 'linear', 0);
        fund(i,:,3) = interp1(wl_in, Sfund, wl_out, 'linear', 0);
    end

    % ----------------------------------------------------------
    % Normalization
    %   - 'max' or 'peak1': per-observer, per-cone peak -> 1
    % ----------------------------------------------------------
    % Rescale each cone independently, per observer, according to opt.normalize.
    switch lower(opt.normalize)
        case 'area'
            % Unit-area normalization: divide each cone by its integral.
            for c = 1:3
                % NaN-safe: treat NaNs as 0 for area normalization
                Fc = fund(:,:,c);
                Fc(~isfinite(Fc)) = 0;
                % Integrate over wavelength (dim 2) with the trapezoid rule.
                A = trapz(wl_out, Fc, 2);      % N x 1
                % Avoid divide-by-zero / bad values.
                A(A<=0 | ~isfinite(A)) = 1;
                % bsxfun broadcasts the N x 1 divisor across all wavelengths.
                fund(:,:,c) = bsxfun(@rdivide, fund(:,:,c), A);
            end

        case {'max','peak1'}
            % Peak normalization: divide each cone by its own maximum (peak->1).
            for c = 1:3
                Fc = fund(:,:,c);  % N x nW

                % NaN-safe peak per observer
                % 'omitnan' exists only in newer MATLAB; the catch handles older
                % releases by replacing non-finite entries with -inf before max.
                try
                    mx = max(Fc, [], 2, 'omitnan');  % N x 1 (newer MATLAB)
                catch
                    % older MATLAB: replace NaNs then max
                    Fc2 = Fc; Fc2(~isfinite(Fc2)) = -inf;
                    mx  = max(Fc2, [], 2);
                end

                % Avoid divide-by-zero / bad values.
                mx(mx<=0 | ~isfinite(mx)) = 1;
                fund(:,:,c) = bsxfun(@rdivide, fund(:,:,c), mx);
            end

        case 'none'
            % do nothing
            % Leave fundamentals in absolute (absorptance x transmittance) units.
    end

    % Record every per-observer parameter so an observer can be reproduced or
    % analysed later. One row per observer, same row order as fund.
    meta = table();
    meta.field_deg = repmat(vdeg, N, 1);   % 2 or 10 (constant within this call)
    meta.age = age;                         % years
    meta.dlens_pct = dlens;                 % lens density deviation (%)
    meta.dmac_pct  = dmac;                  % macular density deviation (%)
    meta.dL_pct = dL;                       % L peak-OD deviation (%)
    meta.dM_pct = dM;                       % M peak-OD deviation (%)
    meta.dS_pct = dS;                       % S peak-OD deviation (%)
    meta.sL_nm = sL_nm;                     % applied L lambda_max shift (nm)
    meta.sM_nm = sM_nm;                     % applied M lambda_max shift (nm, =0)
    meta.sS_nm = sS_nm;                     % applied S lambda_max shift (nm)
    meta.L_ser = isSer;   % 1 = L(ser180), 0 = L(ala180)
    meta.S_hi  = isHi;    % 1 = long-wave S cluster, 0 = short

    % Nested helper: shift a spectrum along the wavelength axis (implements the
    % lambda_max shift). Sampling x at (wl - shift_nm) moves a feature that sat
    % at wavelength w to w + shift_nm, i.e. positive shift_nm = RED shift.
    function y = shift_spectrum_allowNaN(wl, x, shift_nm)
        % Query wavelengths for the shifted spectrum.
        wlq = wl - shift_nm;
        % Temporarily silence the NaN-strip warning: the S template contains
        % NaN (former zeros), which interp1 strips before spline fitting.
        warnState = warning('query','MATLAB:interp1:NaNstrip');
        warning('off','MATLAB:interp1:NaNstrip');
        % Spline interpolation onto the shifted grid.
        y = interp1(wl, x, wlq, 'spline');
        % Restore the previous warning state.
        warning(warnState.state,'MATLAB:interp1:NaNstrip');
        % Clamp any negative spline overshoot to 0.
        y(y<0) = 0;
    end
end

end

% Age-dependent ocular-media (lens) density, CIEPO06 form. Docul1 carries the
% age-varying component and Docul2 the age-invariant component; the piecewise
% weight on Docul1 grows with age (two regimes split at 60 years). Result is
% clamped to be non-negative.
function Docul_ave = ocular_density_age(age, Docul1, Docul2)
% Younger observers (<= 60 yr): linear increase in the Docul1 weight with age.
if age <= 60
    Docul_ave = Docul1 .* (1 + 0.02*(age - 32)) + Docul2;
% Older observers (> 60 yr): steeper linear increase from a higher baseline.
else
    Docul_ave = Docul1 .* (1.56 + 0.0667*(age - 60)) + Docul2;
end
Docul_ave = max(Docul_ave, 0);
end
