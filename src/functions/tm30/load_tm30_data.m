function D = load_tm30_data(~)
% load_tm30_data  Load and cache the static reference data for TM-30-18.
%   D = load_tm30_data() returns a struct with fields:
%     wl       81x1  wavelengths 380:5:780 (nm)
%     tcs      81x99 CIE 2017 / TM-30 test colour sample reflectances
%     cmf2     81x4  [wl xbar ybar zbar] CIE 1931 2-deg observer (CCT step)
%     cmf10    81x4  [wl xbar ybar zbar] CIE 1964 10-deg observer (everything else)
%     dayBasis 81x4  [wl S0 S1 S2] CIE daylight component functions
%
%   The arrays are read from data/tm30_refdata.mat (exported from the
%   validated photoSim_revision TM-30 implementation). The result is cached
%   across calls. The optional argument is accepted for backwards
%   compatibility and ignored.

    % Return the cached struct on any call after the first.
    persistent CACHE
    if ~isempty(CACHE), D = CACHE; return; end

    % Locate the .mat on the MATLAB path (expected under ./data).
    f = which('tm30_refdata.mat');
    if isempty(f)
        error('load_tm30_data:missing', ...
            'tm30_refdata.mat not found on the MATLAB path (expected under ./data).');
    end
    % Load the tm30ref struct (wl, tcs, cmf2, cmf10, dayBasis).
    S = load(f, 'tm30ref');
    D = S.tm30ref;

    % Guard: all TM-30 quantities must share the 380:5:780 nm (81-point) grid.
    assert(numel(D.wl)==81 && abs(D.wl(1)-380)<1e-9 && abs(D.wl(end)-780)<1e-9, ...
        'TM-30 data must be on the 380:5:780 grid.');
    % Cache for subsequent calls.
    CACHE = D;
end
