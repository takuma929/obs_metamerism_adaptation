% First developed by ACS, modified by TM.
function [CCT, Duv] = uv_to_cct_ohno2013(uv, cmf2, start, stop, spacing)
% uv_to_cct_ohno2013  Correlated colour temperature via the Ohno (2013) method.
%   [CCT, Duv] = uv_to_cct_ohno2013(uv, cmf2) returns the correlated colour
%   temperature (K) and distance from the Planckian locus (Duv) for the
%   CIE 1960 UCS coordinates uv = [u v].
%
%   cmf2 is the [wl xbar ybar zbar] CIE 1931 2-deg observer used to build the
%   Planckian locus (data/tm30_cmf_2deg_5nm.csv). start/stop/spacing are the
%   Ohno table bounds/multiplier; TM-30 uses 1000, 25000, 1.001.
%
%   Faithful port of colour.temperature.uv_to_CCT_Ohno2013 (cascade-expanded
%   Planckian table + triangular/parabolic solution).

    if nargin < 3 || isempty(start),   start   = 1000;  end
    if nargin < 4 || isempty(stop),    stop    = 25000; end
    if nargin < 5 || isempty(spacing), spacing = 1.001; end

    CCT_MIN = 1000; CCT_MAX = 100000;

    % --- cascade-expanded temperature list (matches colour) ---
    Ti = [start, start+1];
    next_ti = start + 1;
    next_spacing = spacing;
    while true
        next_ti = next_ti * next_spacing;
        if next_ti >= stop, break; end
        Ti(end+1) = next_ti; %#ok<AGROW>
        D = (next_ti - CCT_MIN) / (CCT_MAX - CCT_MIN);
        D = min(max(D,0),1);
        next_spacing = spacing*(1-D) + (1 + (spacing-1)/10)*D;
    end
    Ti = [Ti, stop-1, stop];

    % --- Planckian locus uv at each temperature (2-deg observer) ---
    wl = cmf2(:,1); XYZbar = cmf2(:,2:4);
    [uL, vL] = deal(zeros(numel(Ti),1));
    for i = 1:numel(Ti)
        Sb = planckian_spd(Ti(i), wl);
        XYZ = XYZbar' * Sb;                 % 3x1 (dlambda constant cancels)
        den = XYZ(1) + 15*XYZ(2) + 3*XYZ(3);
        uL(i) = 4*XYZ(1)/den;
        vL(i) = 6*XYZ(2)/den;
    end

    % --- nearest locus point ---
    % Index of the locus sample closest to uv; clamp off the ends so idx-1 and
    % idx+1 neighbours always exist.
    d = hypot(uL - uv(1), vL - uv(2));
    [~, idx] = min(d);
    if idx == 1,           idx = 2;
    elseif idx == numel(Ti), idx = numel(Ti) - 1; end

    % Previous / current / next locus points (temperature, uv, distance).
    Tip = Ti(idx-1); uip = uL(idx-1); vip = vL(idx-1); dip = d(idx-1);
    Tii = Ti(idx);                                     di  = d(idx);
    Tin = Ti(idx+1); uin = uL(idx+1); vin = vL(idx+1); din = d(idx+1);

    % --- triangular solution ---
    l = hypot(uin - uip, vin - vip);
    x = (dip^2 - din^2 + l^2) / (2*l);
    T_t = Tip + (Tin - Tip) * (x/l);
    vtx = vip + (vin - vip) * (x/l);
    sgn = sign(uv(2) - vtx);
    D_uv_t = sqrt(max(dip^2 - x^2, 0)) * sgn;

    % --- parabolic solution ---
    X = (Tin - Tii) * (Tip - Tin) * (Tii - Tip);
    a = (Tip*(din - di)  + Tii*(dip - din) + Tin*(di  - dip)) / X;
    b = -(Tip^2*(din - di) + Tii^2*(dip - din) + Tin^2*(di - dip)) / X;
    c = -(dip*(Tin - Tii)*Tii*Tin + di*(Tip - Tin)*Tip*Tin + din*(Tii - Tip)*Tip*Tii) / X;
    T_p = -b / (2*a);
    D_uv_p = (a*T_p^2 + b*T_p + c) * sgn;

    % Ohno's rule: use the parabolic estimate when the point lies well off the
    % locus (|Duv| >= 0.002), otherwise the triangular estimate is accurate.
    if abs(D_uv_t) >= 0.002
        CCT = T_p;  Duv = D_uv_p;
    else
        CCT = T_t;  Duv = D_uv_t;
    end
end
