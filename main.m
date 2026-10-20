function main(REBUILD)
% MAIN  One-click regeneration of all results and figures for
%       "Adaptation reduces inter-observer metamerism with narrow-band primaries".
%
%   >> main            % rebuild results/*.mat, then redraw every figure
%   >> main(false)     % reuse existing results/*.mat, only redraw figures
%
%   What the study computes (one-paragraph overview for the reader)
%     We generate 10,000 individual colorimetric observers (Asano individual-
%     observer model on the CIE 2006 basis, with gene-determined DISCRETE
%     lambda_max shifts). For each surface reflectance under each test
%     illuminant we compute every observer's LMS cone excitations to the
%     "world" spectrum, then solve the RGB drives that make a virtual Gaussian-
%     primary display reproduce THAT observer's LMS (a per-observer metameric
%     match). Because a display driven to match one observer is not a match for
%     another, the reproduced stimuli disagree once read out in a common space;
%     that disagreement, scored as CAM16-UCS Delta-E', is "observer metamerism".
%     We sweep the display primary bandwidth (sigma = 5/10/20 nm) and the degree
%     of chromatic adaptation D (0, 0.5, 1; CIECAM16/CAT16) across three
%     CCT-matched illuminant classes (daylight, smooth, peaky) to show how much
%     of the metamerism chromatic adaptation removes.
%
%   Stages
%     (1) run_metamerism_analysis  - simulate 10,000 Asano observers (with
%         gene-determined discrete lambda_max shifts), sweep primary bandwidth
%         and illuminant, and quantify observer metamerism before/after
%         adaptation via the CIECAM16 / CAT16 framework. Writes
%         results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat.
%     (2) compare_nonlinear_adapt_ucs - pool the CAM16-UCS metamerism over the
%         degree-of-adaptation sweep. Writes results/compare_nl_adapt_ucs.mat
%         (used by the Fig. 4 box plots).
%     (2b) adapting_white_variants - rerun the D = 1 scoring with the adapting
%         white defined as the mean LMS or the mean chromaticity instead of the
%         brightest surface. Writes results/adapting_white_variants.mat (Fig. 7).
%     (3) redraw every manuscript figure into figs/.
%
%   Requirements
%     - MATLAB (developed with R2025b), Statistics Toolbox (signrank/friedman).
%     - data/  (bundled inputs: CIE tables, reflectances, illuminants).
%
%   Notes
%     Each analysis/plotting script clears the base workspace and derives its
%     own paths from pwd, so they are launched with evalin('base', ...) from
%     the project root to keep this wrapper's variables isolated.

if nargin < 1 || isempty(REBUILD), REBUILD = true; end

% ------------------------------------------------------------------- paths
projectRoot = fileparts(mfilename('fullpath'));
cd(projectRoot);
addpath(genpath(fullfile(projectRoot, 'src')));
for d = {'results','figs'}
    if ~exist(fullfile(projectRoot, d{1}), 'dir'); mkdir(fullfile(projectRoot, d{1})); end
end

% ------------------------------------------------ figure style parameters
if ~exist(fullfile(projectRoot,'data','fig_parameters.mat'), 'file')
    fprintf('Generating data/fig_parameters.mat ...\n');
    evalin('base', 'save_fig_parameters');
end

% ------------------------------------------- 0. build illuminant classes
% daylight + artificial (CCT-matched) reference sets used by the analysis
if REBUILD || ~exist(fullfile(projectRoot,'data','illuminant_sets.mat'), 'file')
    fprintf('\n=== [0/3] Building illuminant classes (daylight + artificial) ===\n');
    evalin('base', 'build_illuminant_sets');
end

% ------------------------------------------------------------ 1. analysis
if REBUILD
    fprintf('\n=== [1/3] Running metamerism analysis (~30-35 min, 10000 observers) ===\n');
    evalin('base', 'run_metamerism_analysis');
    fprintf('\n=== [2/3] Pooling CAM16-UCS metamerism over the adaptation sweep ===\n');
    evalin('base', 'compare_nonlinear_adapt_ucs');
    fprintf('\n=== [2b/3] Adapting-white definitions (mean LMS / mean chromaticity) ===\n');
    evalin('base', 'adapting_white_variants');
else
    fprintf('\n=== [1-2/3] Skipped analysis (reusing results/*.mat) ===\n');
end

% Figure 7 also needs the optimised-white upper bound; that search is slow
% (per-observer Nelder-Mead), so it is only run when its result is missing.
if ~exist(fullfile(projectRoot,'results','optimal_white_upperlimit.mat'), 'file')
    fprintf('\n=== Optimised adapting white (upper bound; slow) ===\n');
    evalin('base', 'optimal_white_upperlimit');
end

% ------------------------------------------------------------- 2. figures
fprintf('\n=== [3/3] Drawing figures into figs/ ===\n');
% start from a clean figs/ so stale figures from earlier runs never linger.
% Each script writes into a per-figure subfolder (figs/fig1 ... figs/fig7);
% clear those subfolders as well as any stray top-level PDFs/PNGs.
figRoot = fullfile(projectRoot,'figs');
delete(fullfile(figRoot,'*.pdf'));
delete(fullfile(figRoot,'*.png'));
subs = dir(fullfile(figRoot,'fig*'));
for s = 1:numel(subs)
    if subs(s).isdir
        delete(fullfile(figRoot,subs(s).name,'*.pdf'));
        delete(fullfile(figRoot,subs(s).name,'*.png'));
    end
end
% One script per manuscript figure; each is self-contained and can also be
% run on its own. NOTE: the script file names keep their original numbering,
% which no longer matches the manuscript figure order (the figures were
% renumbered). The mapping script -> manuscript figure is given below.
figures = { ...
    'plot_fig1_observers', ...        % Fig 1: individual-observer cone fundamentals
    'plot_fig2_stimuli', ...          % Fig 2: reflectances, illuminant classes, primaries, gamut, sRGB
    'plot_fig3_schematic', ...        % Fig 3: metamerism-quantification schematic (reproduced spectra)
    'plot_fig4_degree_of_metamerism', ... % Fig 4: metamerism vs degree of adaptation (box plots + reduction)
    'plot_fig5_ellipses', ...         % Fig 5: CAM16-UCS (a',b') ellipses before/after adaptation
    'plot_fig6_diagonality', ...      % Fig 6: cross-observer map Γ distribution
    'plot_fig7_adapting_white'};      % Fig 7: adapting-white definitions (bar graph)

for i = 1:numel(figures)
    fprintf('  -> %s\n', figures{i});
    try
        evalin('base', figures{i});
    catch ME
        warning('main:figureFailed', '  %s failed: %s', figures{i}, ME.message);
    end
end

fprintf('\nDone. Figures written to %s\n', fullfile(projectRoot,'figs'));
end
