% First developed by ACS, modified by TM.
function spd = cie_daylight_spd(cct, wl, basis)
% cie_daylight_spd  CIE daylight illuminant spectral power distribution.
%   spd = cie_daylight_spd(cct, wl, basis) returns the CIE D-series daylight
%   SPD for correlated colour temperature cct (K), sampled at wavelengths wl
%   (nm), as a column vector.
%
%   basis is the [wl S0 S1 S2] matrix of daylight component functions read
%   from data/cie_daylight_basis_5nm.csv. If its wavelength column differs
%   from wl the components are linearly interpolated onto wl.
%
%   Reproduces colour.colorimetry.sd_CIE_illuminant_D_series /
%   colour.temperature.CCT_to_xy_CIE_D.

    % --- daylight locus chromaticity (CCT_to_xy_CIE_D) ---
    if cct <= 7000
        xd = -4.6070e9/cct^3 + 2.9678e6/cct^2 + 0.09911e3/cct + 0.244063;
    else
        xd = -2.0064e9/cct^3 + 1.9018e6/cct^2 + 0.24748e3/cct + 0.237040;
    end
    yd = -3.000*xd^2 + 2.870*xd - 0.275;

    % --- component weights ---
    M  = 0.0241 + 0.2562*xd - 0.7341*yd;
    M1 = (-1.3515 -  1.7703*xd +  5.9114*yd) / M;
    M2 = ( 0.0300 - 31.4424*xd + 30.0717*yd) / M;

    % Split the basis into its wavelength column and S0,S1,S2 components.
    wlB = basis(:,1);
    S0 = basis(:,2); S1 = basis(:,3); S2 = basis(:,4);
    % Interpolate the components onto wl only if the grids differ.
    if numel(wlB) ~= numel(wl) || any(abs(wlB(:)-wl(:)) > 1e-9)
        S0 = interp1(wlB, S0, wl(:), 'linear', 'extrap');
        S1 = interp1(wlB, S1, wl(:), 'linear', 'extrap');
        S2 = interp1(wlB, S2, wl(:), 'linear', 'extrap');
    end

    % Reconstruct the daylight SPD and return it as a column vector.
    spd = S0 + M1*S1 + M2*S2;
    spd = spd(:);
end
