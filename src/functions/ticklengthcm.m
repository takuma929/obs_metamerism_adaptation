function ticklengthcm(axh,cm)
% this code was downloaded from code exchange (01/09/2024 TM)
% https://uk.mathworks.com/matlabcentral/fileexchange/55549-ticklengthcm-axh-cm

% ticklengthcm allows you to specify the length of Axes ticks in
% centimeters. In MATLAB, TickLength of Axes object is specified in
% relaiton to the longest axes. When you have mutliple axes of different
% sizes in one figure, then TickLenght for each Axes objects will be
% different. ticklengthcm allows you to have tick marks with an uniform
% length throughout the figure.
%
% ticklengthcm(axh,cm)
%
% INPUT ARGUMENTS
% axh       an Axes, ColorBar or a graphics array containg Axes handles
%
% cm        Legnth of ticks in centimeters (scalar)
%
%
% LIMITATION
% only supports 2D plots
% How to get the length of the longest axes in 3D plot?
%
% DOCUMENTATION on TickLength propertis of Axes
%
%   [0.01 0.025] (default) | two element vector
%
%   Tick mark length, specified as a two-element vector of the form
%   [2Dlength 3Dlength]. The first element is the tick mark length in 2-D
%   views and the second element is the tick mark length in 3-D views.
%   Specify the values in units normalized relative to the longest of the
%   visible x-axis, y-axis, or z-axis lines.
%
% See also
% doc Axes properties
% setupFigure, putpanelID, savefigplus
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 11-Jan-2017 15:45:14


% Validate inputs: axh is a scalar axes/colorbar (or a graphics array), and
% cm is a non-negative scalar length.
p = inputParser;
p.addRequired('axh',@(x) (isscalar(x) && ...
    all(isgraphics(x,'axes') | isgraphics(x,'colorbar'))) ...
    || isa(x,'matlab.graphics.Graphics'));
p.addRequired('cm', @(x) isscalar(x) && x >= 0);
p.parse(axh,cm);

% Setting TickLength via a scalar requires R2014b (8.4.0) or later.
assert(~verLessThan('matlab','8.4.0'),...
    'ticklength only works with MATLAB R2014b or later.');

% Apply to each handle in the (possibly array) input.
for i = 1:numel(axh)
    if isgraphics(axh(i),'axes') || isgraphics(axh(i),'colorbar')

        % Measure the axes in cm so the requested length is absolute.
        oriignalunits = axh(i).Units;
        axh(i).Units = 'centimeters';

        pos = axh(i).Position;

        % MATLAB normalises TickLength to the longest visible axis; find it.
        if pos(3) > pos(4)
            longest = pos(3);
        else
            longest = pos(4);
        end

        % Convert the desired cm length into the normalised fraction.
        newlength = cm/longest; % longest may be wrong for some cases

        % Axes take a [2D 3D] pair (keep the existing 3D value); colorbars
        % take a scalar.
        if isgraphics(axh(i),'axes')
            axh(i).TickLength = [newlength,axh(i).TickLength(2)];
        elseif isgraphics(axh(i),'colorbar')
            % must be scalar
            axh(i).TickLength = newlength;

        end

        % Restore the handle's original units.
        axh(i).Units = oriignalunits;
    end
end

end