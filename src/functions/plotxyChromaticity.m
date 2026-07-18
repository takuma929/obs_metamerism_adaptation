function plotxyChromaticity(varargin)
%%%%%%% Note: This is the code included in PsychToolbox3 (http://psychtoolbox.org/) %%%%%%%

%plotChromaticity Plot color reproduction on chromaticity diagram
%
%   plotChromaticity(colorTable) plots the measured and reference colors of
%   an esfrChart captured in a colorTable on a chromaticity diagram.
%
%   plotChromaticity(___, Name, Value, ___) plots the measured and reference colors of
%   an esfrChart with additional parameters controlling aspects of the
%   display.
%
%   Parameters are:
%   'alpha'   :   alpha value for the chromaticity diagram.  Default is 0.5.

%   'displayROIIndex'   :   Logical controlling whether color patch indices are 
%                           overlaid or not. Default is true.
%
%   'Parent'            :   Handle of an axes that specifies the parent of the 
%                           image object created by displayColorPatch.
%
%   plotChromaticity() plots an empty chromaticity diagram.
%
%   plotChromaticity(Name, Value) plots an empty chromaticity diagram with
%   additional parameter 'Parent' which specifies a handle of a parent axes
%   of the plot object.
%
%   Class Support
%   -------------
%   colorTable is a table which can be computed by measureColor function.
%
%   Notes
%   -----
%   1.  Numbering convention for color patches match the displayed numbers on 
%       an esfrChart using displayChart function.
%
%   Example
%   -------
%   % This example shows the procedure for measuring color
%   % and plotting the results
%
%   I = imread('eSFRTestImage.jpg');
%   I = rgb2lin(I);
%   chart = esfrChart(I);
%   figure
%   displayChart(chart);
%   colorTable = measureColor(chart);
%
%   % Plot results
%   figure
%   plotChromaticity(colorTable)
%
%   See also esfrChart, measureColor, displayColorPatch

%   Copyright 2017-2019 The MathWorks, Inc.


% METHOD (as used in this project): draws a CIE 1931 xy chromaticity diagram
% as a backdrop -- the spectral locus outline, the purple (non-spectral) line,
% and an interpolated sRGB fill inside the locus. If a colorTable is supplied,
% it also overlays reference/measured chromaticities and arrows between them.
% Sub-functions below (parseInputs, validate*) handle name/value parsing.

% Accept 0..5 positional/name-value inputs; reject gpuArray inputs.
narginchk(0,5);
matlab.images.internal.errorIfgpuArray(varargin{:});

% Parse inputs: 2-arg form is the ('Parent', ax) empty-diagram case; any
% other non-zero nargin carries a colorTable plus optional name/value pairs.
if nargin > 0
    if nargin ~= 2
        options = parseInputs(varargin{:});
        colorTable = options.colorTable;
        displayROIIndex = options.displayROIIndex;
        parentAxis = options.Parent;
    elseif nargin == 2
        options = parseInputParent(varargin{:});
        parentAxis = options.Parent;
    end
else
    parentAxis = [];
end

% Load spectral locus xy values at 1-nm intervals
dat = load(fullfile(toolboxdir('images'),'images','+images','+internal', ...
    '+testchart','locus.mat'));
locus = dat.locus;
plotLineWidth = 1;

if isempty(parentAxis)
    h = gcf;
    set(h,'Name','Chromaticity Diagram');
    plot(locus(:,1),locus(:,2),'k','LineWidth',plotLineWidth,'HandleVisibility','off');
else
    plot(parentAxis, locus(:,1),locus(:,2),'k','LineWidth',plotLineWidth,'HandleVisibility','off');
end

% Fix the aspect ratio and axis limits so the diagram is undistorted.
plotAxis = gca;
set(plotAxis,'DataAspectRatioMode','Manual');
grid on
hold on
axis([0.0 0.85 0.0 0.85])
xlabel('x')
ylabel('y')

% Close the horseshoe: straight line from the last to the first locus point.
% plot the non-spectral locus
plot(plotAxis, [locus(1,1) locus(end,1)], [locus(1,2) locus(end,2)],'k','LineWidth',plotLineWidth,'HandleVisibility','off')
% chromaticity coordinates of spectrum locus
x = [ 0.175596 0.172787 0.170806 0.170085 0.160343 0.146958 0.139149 ...
    0.133536 0.126688 0.115830 0.109616 0.099146 0.091310 0.078130 ...
    0.068717 0.054675 0.040763 0.027497 0.016270 0.008169 0.004876 ...
    0.003983 0.003859 0.004646 0.007988 0.013870 0.022244 0.027273 ...
    0.032820 0.038851 0.045327 0.052175 0.059323 0.066713 0.074299 ...
    0.089937 0.114155 0.138695 0.154714 0.192865 0.229607 0.265760 ...
    0.301588 0.337346 0.373083 0.408717 0.444043 0.478755 0.512467 ...
    0.544767 0.575132 0.602914 0.627018 0.648215 0.665746 0.680061 ...
    0.691487 0.700589 0.707901 0.714015 0.719017 0.723016 0.734674 ]';
y = [ 0.005295 0.004800 0.005472 0.005976 0.014496 0.026643 0.035211 ...
    0.042704 0.053441 0.073601 0.086866 0.112037 0.132737 0.170464 ...
    0.200773 0.254155 0.317049 0.387997 0.463035 0.538504 0.587196 ...
    0.610526 0.654897 0.675970 0.715407 0.750246 0.779682 0.792153 ...
    0.802971 0.812059 0.819430 0.825200 0.829460 0.832306 0.833833 ...
    0.833316 0.826231 0.814796 0.805884 0.781648 0.754347 0.724342 ...
    0.692326 0.658867 0.624470 0.589626 0.554734 0.520222 0.486611 ...
    0.454454 0.424252 0.396516 0.372510 0.351413 0.334028 0.319765 ...
    0.308359 0.299317 0.292044 0.285945 0.280951 0.276964 0.265326 ]';
% Build the coloured fill by sweeping the locus in wedges from the white
% point e=(1/3,1/3) outward. Each wedge is subdivided into `steps` quad
% patches; xy4rgb accumulates [x y R G B] for the four vertices of each quad.
N = length(x);
i = 1;
e = 1/3;
steps = 25;
xy4rgb = zeros(N*steps*4,5,'double');
for w = 1:N                                     % wavelength
    w2 = mod(w,N)+1;
    a1 = atan2(y(w) - e,x(w) - e);              % start angle
    a2 = atan2(y(w2) - e,x(w2) - e);            % end angle
    r1 = ((x(w) - e)^2 + (y(w) - e)^2)^0.5;     % start radius
    r2 = ((x(w2) - e)^2 + (y(w2) - e)^2)^0.5;   % end radius
    for c = 1:steps                             % colorfulness
        % patch polygon
        xyz(1,1) = e+r1*cos(a1)*c/steps;
        xyz(1,2) = e+r1*sin(a1)*c/steps;
        xyz(1,3) = 1 - xyz(1,1) - xyz(1,2);
        xyz(2,1) = e+r1*cos(a1)*(c-1)/steps;
        xyz(2,2) = e+r1*sin(a1)*(c-1)/steps;
        xyz(2,3) = 1 - xyz(2,1) - xyz(2,2);
        xyz(3,1) = e+r2*cos(a2)*(c-1)/steps;
        xyz(3,2) = e+r2*sin(a2)*(c-1)/steps;
        xyz(3,3) = 1 - xyz(3,1) - xyz(3,2);
        xyz(4,1) = e+r2*cos(a2)*c/steps;
        xyz(4,2) = e+r2*sin(a2)*c/steps;
        xyz(4,3) = 1 - xyz(4,1) - xyz(4,2);
        % compute sRGB for vertices
        rgb = images.internal.testchart.xyz2srgb(xyz');
        % store the results
        xy4rgb(i:i+3,1:2) = xyz(:,1:2);
        xy4rgb(i:i+3,3:5) = rgb';
        i = i + 4;
    end
end
% Render each stored quad as a patch with per-vertex interpolated sRGB.
[rows, ~] = size(xy4rgb);
f = [1 2 3 4];
v = zeros(4,3,'double');
for i = 1:4:rows
    v(:,1:2) = xy4rgb(i:i+3,1:2);
    patch('Vertices',v, 'Faces',f, 'EdgeColor','none', ...
        'FaceVertexCData',xy4rgb(i:i+3,3:5),'FaceColor','interp','FaceAlpha',0.3,'HandleVisibility','off')
end

% When a colorTable was supplied, overlay reference vs measured patches.
if nargin > 0 && nargin ~= 2
    hold on;

    numColorPatches = 16;

    % Reference colours: Lab -> XYZ -> xy chromaticity (normalise by sum).
    XYZ_ref = lab2xyz([colorTable.Reference_L colorTable.Reference_a colorTable.Reference_b]);
    xyz_ref = XYZ_ref./sum(XYZ_ref,2);

    % Measured colours: RGB -> XYZ -> xy chromaticity (normalise by sum).
    XYZ_measured = rgb2xyz([colorTable.Measured_R colorTable.Measured_G colorTable.Measured_B]);
    xyz_measured = XYZ_measured./sum(XYZ_measured,2);

    % Marker colours: reference points red (1:16), measured points green (17:32).
    scatter_c = zeros(32,3);
    scatter_c(1:16,:) = repmat([1 0 0],16,1);
    scatter_c(17:end,:) = repmat([0 1 0],16,1);

    % Scatter all reference then measured chromaticities.
    scatter(plotAxis, [xyz_ref(:,1);xyz_measured(:,1)],[xyz_ref(:,2);xyz_measured(:,2)],50, scatter_c, 'filled','HandleVisibility','off');hold on;
    % Draw an arrow from each reference point to its measured counterpart.
    pt_txt = cell(numColorPatches,1);
    for j=1:numColorPatches
        pt_txt{j} = num2str(j);
        p1 = [xyz_ref(j,1) xyz_ref(j,2)];
        p2= [xyz_measured(j,1) xyz_measured(j,2)];
        dp = p2-p1;
        quiver(plotAxis, p1(1),p1(2),dp(1),dp(2),0,'k','LineWidth',3, 'MaxHeadSize',0.9,'HandleVisibility','off');hold on;
    end
    % Optionally label each reference point with its patch index.
    if displayROIIndex
        text(xyz_ref(:,1),xyz_ref(:,2),pt_txt,'FontSize',12,'FontWeight','bold','Color',[1 1 1],'HandleVisibility','off');    
    end
end
hold off;
end

function options = parseInputs(varargin)
% Parse the colorTable form: required colorTable + optional displayROIIndex/Parent.

parser = inputParser();
parser.addRequired('colorTable',@validateTable);
parser.addParameter('displayROIIndex',true,@validateDisplayFlag);
parser.addParameter('Parent',[],@validateParentAxis);

parser.parse(varargin{:});
options = parser.Results;
end

function options = parseInputParent(varargin)
% Parse the 2-arg empty-diagram form: only the 'Parent' axes name/value pair.
parser = inputParser();
parser.addParameter('Parent',[],@validateParentAxis);
parser.parse(varargin{:});
options = parser.Results;
end

function validateTable(colorTable)
% Require colorTable to be a non-empty table.
    validateattributes(colorTable,{'table'},{'nonempty'},mfilename,'colorTable',1);
end

function validateDisplayFlag(flag)
% Require displayROIIndex to be a scalar logical.
    supportedClasses = {'logical'};
    attributes = {'nonempty','finite','nonsparse','scalar','nonnan'};
    validateattributes(flag,supportedClasses,attributes,...
        mfilename);
end

function validateParentAxis(Parent)
% Require Parent to be a valid axes handle.
validateattributes(Parent, {'matlab.graphics.axis.Axes'},{'nonempty','nonsparse'});
end
