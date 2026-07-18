function daylightspd = get_daylightspd(CCT,wls)
% get_daylightspd  Reconstruct a CIE daylight SPD at a given CCT.
%   daylightspd = get_daylightspd(CCT, wls) returns the relative spectral
%   power distribution of the CIE D-series daylight illuminant for correlated
%   colour temperature CCT (K), sampled at the wavelengths wls (nm).
%
%   INPUTS
%     CCT : correlated colour temperature in kelvin (scalar).
%     wls : column vector of wavelengths (nm) at which to evaluate the SPD.
%
%   OUTPUT
%     daylightspd : column vector (numel(wls)x1) of relative spectral power.
%
%   METHOD
%     Classic CIE reconstruction daylightspd = S0 + M1*S1 + M2*S2, where the
%     three daylight component functions S0,S1,S2 are read from data/DSPD.mat
%     (columns 2,3,4; column 1 is wavelength) and interpolated onto wls. The
%     daylight-locus chromaticity (xd,yd) and the weights M1,M2 follow the
%     standard CIE polynomials in CCT.

    % Load the daylight basis: DSPD = [wl S0 S1 S2].
    load(fullfile('data','DSPD'),'DSPD')

    % Linearly interpolate the S0,S1,S2 components onto the requested wls;
    % keep wls as the first column so DSPD_interpolated stays [wl S0 S1 S2].
    %linearly interpolate DSPD
    DSPD_interpolated = [wls,interp1(DSPD(:,1),DSPD(:,[2 3 4]),wls','linear')];

    % Daylight-locus x chromaticity xd: two CCT ranges use different
    % coefficients (<=7000 K vs >7000 K), per the CIE daylight definition.
    % calculate x_d,y_d based on CCT color temperature
    if CCT <= 7000
        xd = .244063 + .09911*(1e3/CCT) + 2.9678*(1e6/(CCT^2)) - 4.6070*(1e9/(CCT^3));
    else 
        xd = .237040 + .24748*(1e3/CCT) + 1.9018*(1e6/CCT^2) - 2.0064*(1e9/CCT^3);
    end
    
    % Daylight-locus y chromaticity yd as a quadratic in xd.
    yd = -3.000*xd^2 + 2.870*xd - 0.275;

    % Component weights M1,M2 from (xd,yd) via the CIE normalisation term M.
    % calculate relatative SPD
    M = 0.0241 + 0.2562*xd - 0.7341*yd;
    M1 = (-1.3515 - 1.7703*xd + 5.9114*yd)/M;
    M2 = (0.03000 - 31.4424*xd + 30.0717*yd)/M;

    % Reconstruct the SPD: S0 + M1*S1 + M2*S2 (columns 2,3,4 of the basis).
    daylightspd = DSPD_interpolated(:,2) + M1*DSPD_interpolated(:,3) + M2*DSPD_interpolated(:,4);
end