% First developed by ACS, modified by TM.
function out = compute_tm30(spd, wl, D)
% compute_tm30  ANSI/IES TM-30-18 colour rendition metrics (pure MATLAB).
%   out = compute_tm30(spd) computes TM-30-18 quantities for the test light
%   source whose spectral power distribution is spd, sampled at 380:5:780 nm.
%
%   out = compute_tm30(spd, wl) resamples spd from wl onto 380:5:780 nm first.
%   out = compute_tm30(spd, wl, D) reuses pre-loaded static data D (see
%   load_tm30_data) to avoid repeated disk reads in a loop.
%
%   Returned struct fields:
%     CCT     correlated colour temperature (K), Ohno (2013)
%     Duv     distance from the Planckian locus
%     Rf      fidelity index R_f
%     Rg      gamut index R_g
%     Rcs     16x1 local chroma shifts (%), Rcs(1) is the red-bin Rcs,h1
%     Rfs     16x1 local fidelity indices
%     Rhs     16x1 local hue shifts
%     bins    99x1 hue-bin assignment (0..15) of the reference samples
%
%   This is a faithful reimplementation of colour-science's
%   colour_fidelity_index_ANSIIESTM3018 (CIE 2017 CFI + TM-30 gamut/chroma),
%   validated to match it on FL2 (Rf=70.12), D65/A (Rf=Rg=100) and the 401
%   Houser illuminants. See functions/tm30/validate_tm30.m.

    % Load the static reference data once (or reuse a caller-supplied D).
    if nargin < 3 || isempty(D), D = load_tm30_data(); end
    wlStd = D.wl;
    % Resample the test SPD onto the standard 380:5:780 grid if wl differs.
    if nargin >= 2 && ~isempty(wl) && (numel(wl)~=numel(wlStd) || any(abs(wl(:)-wlStd)>1e-9))
        spd = interp1(wl(:), spd(:), wlStd, 'linear', 0);
    end
    % Clamp negatives and force a column vector.
    spd = max(spd(:), 0);

    % --- 1. CCT / Duv from the 2-deg test chromaticity (Ohno 2013) ---
    XYZt = D.cmf2(:,2:4)' * spd;                       % 3x1
    den  = XYZt(1) + 15*XYZt(2) + 3*XYZt(3);
    uv   = [4*XYZt(1)/den, 6*XYZt(2)/den];
    [CCT, Duv] = uv_to_cct_ohno2013(uv, D.cmf2, 1000, 25000, 1.001);

    % --- 2. reference illuminant SPD ---
    sdRef = tm30RefSPD(CCT, wlStd, D);

    % --- 3. TCS colorimetry (10-deg) under test and reference ---
    [Jt, at, bt, ~]      = tcsCAM(spd,  D);
    [Jr, ar, br, hRef]   = tcsCAM(sdRef, D);
    JpapbpT = [Jt, at, bt];
    JpapbpR = [Jr, ar, br];

    % --- 4. fidelity ---
    dE = sqrt(sum((JpapbpT - JpapbpR).^2, 2));         % 99x1
    Rf = deltaEToRf(mean(dE));

    % --- 5. hue bins from the reference sample hue angles ---
    bins = floor(hRef / 22.5);                          % 0..15
    avgT = zeros(16,2); avgR = zeros(16,2); binDE = zeros(16,1);
    for kbin = 0:15
        m = bins == kbin;
        avgT(kbin+1,:) = mean([at(m) bt(m)], 1);
        avgR(kbin+1,:) = mean([ar(m) br(m)], 1);
        binDE(kbin+1)  = mean(dE(m));
    end

    % --- 6. gamut index ---
    Rg = 100 * shoelace(avgT) / shoelace(avgR);

    % --- 7. local chroma / hue shifts ---
    ang    = (22.5*(0:15) + 11.25)' * pi/180;
    normsR = sqrt(sum(avgR.^2, 2));
    da = avgT(:,1) - avgR(:,1);
    db = avgT(:,2) - avgR(:,2);
    Rcs = 100 * (da.*cos(ang) + db.*sin(ang)) ./ normsR;
    Rhs =       (-da.*sin(ang) + db.*cos(ang)) ./ normsR;
    Rfs = deltaEToRf(binDE);

    out = struct('CCT',CCT,'Duv',Duv,'Rf',Rf,'Rg',Rg, ...
                 'Rcs',Rcs,'Rfs',Rfs,'Rhs',Rhs,'bins',bins);
end

% ------------------------------------------------------------------ helpers
function sd = tm30RefSPD(CCT, wl, D)
% Reference illuminant per CIE 2017 / TM-30: Planckian below 4000 K, CIE
% daylight above 5000 K, and a normalised Planckian/daylight blend between.
    if CCT < 4000
        sd = planckian_spd(CCT, wl);
    elseif CCT > 5000
        sd = cie_daylight_spd(CCT, wl, D.dayBasis);
    else
        sp = planckian_spd(CCT, wl);
        sd_ = cie_daylight_spd(CCT, wl, D.dayBasis);
        yb = D.cmf2(:,3);                       % 2-deg ybar for the blend norm
        sp  = sp  / (yb' * sp);
        sd_ = sd_ / (yb' * sd_);
        m = (CCT - 4000) / 1000;
        sd = sp + m*(sd_ - sp);                 % linstep_function(m, sp, sd_)
    end
end

function [Jp, ap, bp, h] = tcsCAM(illum, D)
% CAM02-UCS coordinates of the 99 TCS under illuminant illum (10-deg).
    cmf = D.cmf10(:,2:4);                        % 81x3
    Xw  = cmf' * illum;                          % 3x1 white tristimulus
    k   = 100 / Xw(2);
    XYZw = (k * Xw)';                            % 1x3, Yw=100
    ills = k * illum;                            % scaled illuminant
    XYZ = (D.tcs .* ills)' * cmf;                % 99x3
    [J, M, h] = ciecam02JMh(XYZ, XYZw);
    % CAM02-UCS
    Jp = (1 + 100*0.007) * J ./ (1 + 0.007*J);
    Mp = (1/0.0228) * log(1 + 0.0228*M);
    ap = Mp .* cos(h*pi/180);
    bp = Mp .* sin(h*pi/180);
end

function [J, M, h] = ciecam02JMh(XYZ, XYZw)
% CIECAM02 forward model (Average surround, L_A=100, Y_b=20,
% discount_illuminant = true) returning lightness J, colourfulness M and
% hue angle h (deg). XYZ is Nx3 in a 0..100 scale, XYZw is 1x3 (Yw=100).
    LA = 100; Yb = 20; F = 1.0; c = 0.69; Nc = 1.0;
    MCAT02 = [ 0.7328  0.4296 -0.1624; ...
              -0.7036  1.6975  0.0061; ...
               0.0030  0.0136  0.9834];
    MHPE   = [ 0.38971 0.68898 -0.07868; ...
              -0.22981 1.18340  0.04641; ...
               0.00000 0.00000  1.00000];
    Minv = MHPE / MCAT02;                        % MHPE * inv(MCAT02)

    RGB  = XYZ  * MCAT02';                        % Nx3
    RGBw = XYZw * MCAT02';                        % 1x3
    Yw = XYZw(2);
    Dd = 1.0;                                     % discount_illuminant = true
    Dfac = Yw*Dd ./ RGBw + 1 - Dd;               % 1x3
    RGBc  = RGB  .* Dfac;
    RGBwc = RGBw .* Dfac;

    n = Yb/Yw; z = 1.48 + sqrt(n);
    Nbb = 0.725*(1/n)^0.2; Ncb = Nbb;
    kk = 1/(5*LA+1);
    FL = 0.2*kk^4*(5*LA) + 0.1*(1-kk^4)^2*(5*LA)^(1/3);

    RGBp  = RGBc  * Minv';                        % Nx3 in HPE space
    RGBwp = RGBwc * Minv';                        % 1x3

    pa  = postAdapt(RGBp,  FL);
    pwa = postAdapt(RGBwp, FL);

    a = pa(:,1) - 12/11*pa(:,2) + 1/11*pa(:,3);
    b = (1/9)*(pa(:,1) + pa(:,2) - 2*pa(:,3));
    h = atan2(b, a) * 180/pi; h(h<0) = h(h<0) + 360;
    et = 0.25*(cos(h*pi/180 + 2) + 3.8);

    A  = (2*pa(:,1)  + pa(:,2)  + (1/20)*pa(:,3)  - 0.305) * Nbb;
    Aw = (2*pwa(1)   + pwa(2)   + (1/20)*pwa(3)   - 0.305) * Nbb;
    J  = 100*(A/Aw).^(c*z);
    t  = (50000/13*Nc*Ncb .* et .* sqrt(a.^2+b.^2)) ./ ...
         (pa(:,1) + pa(:,2) + (21/20)*pa(:,3));
    C  = t.^0.9 .* sqrt(J/100) .* (1.64 - 0.29^n)^0.73;
    M  = C * FL^0.25;
end

function pa = postAdapt(RGB, FL)
% CIECAM02 nonlinear post-adaptation (cone) response compression.
    x = (FL * abs(RGB) / 100).^0.42;
    pa = 400 * sign(RGB) .* x ./ (27.13 + x) + 0.1;
end

function R = deltaEToRf(dE)
% Convert a colour difference dE into a TM-30 fidelity score (soft 0..100).
    R = 10 * log(1 + exp((100 - 6.73*dE) / 10));
end

function A = shoelace(P)
% Signed polygon area of the 16 hue-bin (a',b') averages.
    N = size(P,1);
    A = 0;
    for i = 1:N
        j = mod(i, N) + 1;
        A = A + (P(i,1)*P(j,2) - P(i,2)*P(j,1)) / 2;
    end
end
