%--------------------------------------------------------------------------
% save_fig_parameters.m
% Figure Parameter Initialization Script
%
% PURPOSE
%   Define a single, shared set of figure-style parameters (column widths,
%   font sizes, font name) and save them to disk so that every figure-making
%   script in the project can load the SAME values. Centralising the style
%   here guarantees a consistent look across all figures in the manuscript
%   (e.g. all panels use the identical two-column width and font).
%
% METHOD
%   The parameters are collected into a single struct 'figp' and written to a
%   .mat file. Widths are given in centimetres to match the journal's column
%   specification (Optics Express two-column layout); a one-column width is
%   derived as half the two-column width. Downstream scripts do
%   load('data/fig_parameters.mat') and then size axes/fonts from figp.
%
% INPUTS
%   (none) -- all values are hard-coded below.
%
% OUTPUTS
%   - Ensures the output directories exist (creates 'figs/' if absent).
%   - data/fig_parameters.mat  containing the struct 'figp' with fields:
%       figp.twocolumn      two-column figure width  [cm]
%       figp.onecolumn      one-column figure width  [cm]
%       figp.fontsize       general text font size   [pt]
%       figp.fontsize_axis  axis-label font size     [pt]
%       figp.fontname       font family name (string)
%
% HOW TO RUN
%   Run from the project root (so that pwd/data and pwd/figs resolve):
%     >> save_fig_parameters
%
% Author: TM, 2025
%--------------------------------------------------------------------------

% Start from a clean workspace and close any open figures so nothing stale
% leaks into the saved parameters.
clearvars; close all;

disp('Saving figure parameters...')

% Create the output figure directory (relative to the current working
% directory) if it is not already present, so later save() calls do not fail.
% Create directories for saving figures and data if they do not already exist
if ~exist(fullfile(pwd,'figs'),'dir')
    mkdir(fullfile(pwd,'figs'))
end

% Figure widths in centimetres, matching the journal column layout. The
% one-column width is exactly half the two-column width so single-panel
% figures line up with a column of the double-width ones.
% Set figure width (in centimeters) for two-column and one-column formats
figp.twocolumn = 13.04;              % Standard two-column width
figp.onecolumn = figp.twocolumn/2;  % Standard one-column width

% Font sizes in points: a smaller size for general text/annotations and a
% slightly larger size for axis labels (kept distinct for readability).
% Set font sizes for general text and axes
figp.fontsize      = 7;   % General font size
figp.fontsize_axis = 8;   % Axis label font size

% Single font family applied to every figure for a uniform typeface.
% Set font name for all figures
figp.fontname = 'Arial';

% Persist the whole 'figp' struct to data/fig_parameters.mat; other scripts
% load this file to reuse the identical style settings.
% Save the figure parameter structure for use in other scripts
save(fullfile(pwd,'data','fig_parameters'), 'figp')

disp('Done.')
