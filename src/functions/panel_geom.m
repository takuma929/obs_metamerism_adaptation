function [figW, axPos, figH] = panel_geom(nPerRow, figp, aspect)
% panel_geom  Consistent panel geometry so N panels tile to the two-column width.
%   [figW, axPos, figH] = panel_geom(nPerRow, figp, aspect)
%     nPerRow : how many of these panels sit side-by-side in the final figure
%               (2 -> each is twocolumn/2, 3 -> twocolumn/3, 1 -> full width)
%     figp    : project figure-parameter struct (uses figp.twocolumn)
%     aspect  : axes height / axes width (default 0.72)
%   Returns the figure width figW (cm), the axes position axPos = [l b w h] (cm),
%   and the figure height figH (cm) WITHOUT a legend (legend_above adds its own).
%   Margins are fixed in cm so tick labels and axis labels are never clipped and
%   all panels share the same plotting-box aspect regardless of column fraction.

% Default plotting-box aspect (height/width) when not supplied.
if nargin < 3 || isempty(aspect), aspect = 0.72; end
% Figure width = two-column width divided into nPerRow equal panels.
figW = figp.twocolumn / nPerRow;

% Fixed margins in cm (reserve room for labels so nothing is clipped).
axL = 1.05;      % left margin (y-label)
axB = 0.82;      % bottom margin (x-label)
axR = 0.20;      % right margin
axT = 0.18;      % top margin

% Axes width fills what the L/R margins leave; height follows the aspect.
axW = figW - axL - axR;
axH = axW * aspect;
% Figure height = bottom margin + axes height + top margin (no legend).
figH = axB + axH + axT;
% Axes position [left bottom width height] in cm.
axPos = [axL, axB, axW, axH];
end
