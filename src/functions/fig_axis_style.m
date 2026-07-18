function fig_axis_style(ax, figp)
% fig_axis_style  Shared axis look for all paper figures (light-grey panel,
% faint horizontal grid, Arial at the project font size). Matches the style
% used across plot_fig*.m so every panel reads as one system.
%
%   fig_axis_style(ax, figp) restyles the axes in place; nothing is returned.
%
%   INPUTS
%     ax   : handle to the axes to restyle.
%     figp : project figure-parameter struct; only figp.fontsize is used here.
%
%   OUTPUT
%     none (ax is modified in place).

% Arial text at the project font size on a near-white (0.97 grey) panel.
ax.FontName = 'Arial'; ax.Color = ones(3,1)*0.97; ax.FontSize = figp.fontsize;
% Black axis lines, thin (0.5 pt) rules.
ax.XColor = 'k'; ax.YColor = 'k'; ax.LineWidth = 0.5;
% Faint horizontal-only grid (y-grid on, x-grid off, low alpha) and no box.
grid(ax,'on'); ax.YGrid = 'on'; ax.XGrid = 'off'; ax.GridAlpha = 0.1; box(ax,'off');
end
