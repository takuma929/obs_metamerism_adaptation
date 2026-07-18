% ========================================================================
% check_uv_saturation.m
%
% PURPOSE
%   Sanity check on stimulus saturation: report the CIE 1976 u'v' chromaticity
%   span of the retained CFI surfaces relative to (A) the D65 white and (B) each
%   illuminant's own adopted white. Confirms the manuscript claim that surfaces
%   sit within ~0.14 u'v' of the adopted white (i.e. not over-saturated).
%
% INPUT (.mat / data)
%   results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat
%     (wls, illuminants.spd, byIll.ok/.keepIdx)
%   data/CIE1931XYZ.txt        reference 2 deg CMFs
%   data/CIE_srf_cfi_1nm.csv   CIE 99 CFI surface reflectances (1 nm)
%   data/D65.mat               D65 relative SPD
%
% PRINTS (nothing saved)
%   - number of surfaces retained under EVERY illuminant (intersection set)
%   - max and median surface u'v' distance from the D65 white
%   - max and median surface u'v' distance from each illuminant's adopted white
%
% METHOD
%   For each surface, XYZ = CMFs * (SPD .* reflectance); convert to u'v' and
%   take the Euclidean chromaticity distance to the corresponding white's u'v'.
% ========================================================================

% Check the u'v' saturation span of the retained (76) surfaces, and the reference
% white it is measured from. Manuscript states "u'v' ~ 0.14 from the adopted white".
clearvars; clc;
projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
inMat = fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
% Only the wavelength grid, illuminant SPDs and per-illuminant kept surfaces.
S = load(inMat,'wls','illuminants','byIll');
wls = S.wls(:);
% Reference CIE 1931 2 deg CMFs on the working wavelength grid.
Axyz  = readmatrix(fullfile(dataDir,'CIE1931XYZ.txt'));
T_xyz = interp1(Axyz(:,1),Axyz(:,2:4),wls,'pchip',0)';       % 3 x nW
% CIE 99 CFI surface reflectances, interpolated and clamped to [0,1].
Rmat  = readmatrix(fullfile(dataDir,'CIE_srf_cfi_1nm.csv'));
R_all = min(max(interp1(Rmat(:,1),Rmat(:,2:end),wls,'pchip','extrap'),0),1);   % nW x 99

% common retained set (intersection across illuminants)
% Build the set of surfaces kept under EVERY converged illuminant.
nIll=numel(S.illuminants); keep=[];
for k=1:nIll
    r=S.byIll(k);
    if isfield(r,'ok') && ~isempty(r.ok) && r.ok
        % First ok illuminant seeds the set; later ones intersect it.
        if isempty(keep), keep=r.keepIdx(:); else, keep=intersect(keep,r.keepIdx(:)); end
    end
end
fprintf('common retained surfaces: %d\n', numel(keep));

% Anonymous CIE 1976 u'v' chromaticity from XYZ (columns = colours).
upvp = @(XYZ) deal(4*XYZ(1,:)./(XYZ(1,:)+15*XYZ(2,:)+3*XYZ(3,:)), ...
                   9*XYZ(2,:)./(XYZ(1,:)+15*XYZ(2,:)+3*XYZ(3,:)));

% ---- reference A: D65 ----
% Load and resample the D65 SPD (clamped non-negative).
Dd = load(fullfile(dataDir,'D65.mat'));
d65 = interp1(Dd.wls(:), Dd.d65(:), wls, 'pchip','extrap'); d65=max(d65,0);
% u'v' of the D65 white itself.
XYZ_w65 = T_xyz*d65;
[uw,vw] = upvp(XYZ_w65);
% u'v' of each retained surface seen under D65.
XYZ_s65 = T_xyz*(d65.*R_all(:,keep));            % 3 x 76
[us,vs] = upvp(XYZ_s65);
% Chromaticity distance of each surface from the D65 white.
d_from_D65 = hypot(us-uw, vs-vw);
fprintf('[D65 ref] surface u''v'' distance from D65 white: max=%.3f  median=%.3f  (n=%d)\n', ...
    max(d_from_D65), median(d_from_D65), numel(d_from_D65));

% ---- reference B: each illuminant's own adopted white (max over all illum) ----
% Accumulate surface-to-white distances across every illuminant; track the max.
maxd=0; alld=[];
for k=1:nIll
    r=S.byIll(k);
    if ~(isfield(r,'ok')&&~isempty(r.ok)&&r.ok), continue; end
    % This illuminant's white u'v'.
    spd=S.illuminants(k).spd(:);
    XYZw=T_xyz*spd; [uw2,vw2]=upvp(XYZw);
    % Retained surfaces under this illuminant, in u'v'.
    XYZs=T_xyz*(spd.*R_all(:,keep)); [us2,vs2]=upvp(XYZs);
    % Distances from this white; append to the pooled list and update the max.
    dd=hypot(us2-uw2, vs2-vw2); alld=[alld; dd(:)];
    maxd=max(maxd,max(dd));
end
fprintf('[adopted-white ref] surface u''v'' distance from each illuminant white: max=%.3f  median=%.3f\n', ...
    maxd, median(alld));
