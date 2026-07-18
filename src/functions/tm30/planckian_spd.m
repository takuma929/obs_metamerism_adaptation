% First developed by ACS, modified by TM.
function spd = planckian_spd(T, wl)
% planckian_spd  Planckian (blackbody) spectral radiance at temperature T.
%   spd = planckian_spd(T, wl) returns the blackbody spectral power at the
%   wavelengths wl (nm) for absolute temperature T (K), as a column vector.
%
%   Uses Planck's law with the colorimetric constants adopted by CIE / the
%   colour-science library (c1 = 3.741771e-16 W*m^2, c2 = 1.4388e-2 m*K,
%   n = 1). The absolute scale is irrelevant for TM-30 because every
%   illuminant is renormalised to Y = 100 before use; only the spectral
%   shape matters.
%
%   Matches colour.colorimetry.sd_blackbody.

    % First and second radiation constants; refractive index n = 1 (vacuum).
    c1 = 3.741771e-16;      % 2*pi*h*c^2
    c2 = 1.4388e-2;         % h*c/k_B
    n  = 1;

    % Planck's law evaluated per wavelength (wl converted nm -> m).
    l = wl(:) * 1e-9;       % nm -> m
    spd = (c1 * n^-2 .* l.^-5 / pi) ./ (exp(c2 ./ (n .* l .* T)) - 1);
    spd = spd * 1e-9;       % per-nm scaling as in colour (cosmetic; cancels)
end
