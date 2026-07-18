% ========================================================================
% plot_fig6_ellipses.m   —   Figure 6 (CAM16-UCS a'b' covariance ellipses)
%
% NOTE ON FIGURE NUMBERING: this SCRIPT is named "fig6", but it produces the
% MANUSCRIPT's Figure 5 (script and manuscript numbers differ because the
% figures were renumbered during revision). Output PDFs go to figs/fig6/.
%
% WHAT THIS FIGURE SHOWS (manuscript Figure 5): the CAM16-UCS (a',b') chromatic
% plane with, per representative surface, the 95% between-observer covariance
% ellipse, shown at D = 0 / 0.5 / 1 for narrow (SD=5) and broad (SD=20)
% primaries. Ellipses shrink toward the adopted white (origin) as adaptation
% increases; blue surfaces carry the largest residual spread.
%
% Observer metamerism in the CAM16-UCS chromatic plane (a',b'): 95% covariance
% ellipses of the observer distribution for representative surfaces, as a function
% of the degree of chromatic adaptation D, scored in the SAME CAM16-UCS space:
%   D = 0 (no adaptation) / D = 0.5 (partial) / D = 1 (complete)
% One figure per illuminant CLASS (daylight/smooth/peaky) and an 'all' group,
% at narrow (SD=5) and broad (SD=20) primaries; each figure has three panels.
% Observer spread is pooled over all illuminants of the class as DEVIATIONS from
% the standard-observer reference at the SAME D (so the global adaptation-state
% cast cancels), then placed at the surface's fully-adapted reference a'b' so the
% ellipse centres align across panels. Increasing D shrinks the ellipses.
%
% Variants of the adaptation's functional form are not drawn: a gain commutes
% with the near-power-law Naka-Rushton compression, so applying the white
% normalisation after the compression is the same operation up to a constant and
% the ellipses are super-imposabe. Reported in the text, not as extra panels.
%
% Per-observer XYZ is reconstructed from the stored cone-fundamental family (no
% re-run of the heavy analysis). Observers are subsampled for the covariance
% estimate (CFG.nObsSub).
%
% Run section by section: SETUP first, then any group section.
% Output: figs/fig6/fig6_<class>_uvEllipse_SD{05,20}.pdf
% ========================================================================

%% ============================ SETUP ============================
clearvars; close all; clc;
% Project root inferred from this file (.../src/plotting/plot_fig6_...);
% three fileparts() strip filename -> plotting -> src -> repo root.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist(fullfile(projectRoot,'src'),'dir')
    projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
end
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
figDir  = fullfile(projectRoot,'figs','fig6'); if ~exist(figDir,'dir'); mkdir(figDir); end
% figp: shared figure-style parameters (fonts, column widths).
T = load(fullfile(dataDir,'fig_parameters.mat')); figp = T.figp;
% Main results file (the heavy metamerism analysis); loaded selectively below.
inMat = fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
assert(exist(inMat,'file')==2,'Missing results: %s',inMat);

% ----- Configuration -----
CFG.SDsel   = [5 20];        % narrow / broad primaries
CFG.Dlist   = [0 0.5 1];     % degree-of-adaptation panels (own white, D-matched ref)
CFG.nObsSub = 10000;         % observers used for the covariance ellipse (full sample)
CFG.nSurf   = 30;            % well-spread surfaces drawn
% CIECAM16 viewing-condition parameters:
%   LA = adapting luminance, Yb = background relative luminance,
%   F = surround factor, c = impact of surround, Nc = chromatic induction.
cam = struct('LA',20,'Yb',20,'F',1.0,'c',0.69,'Nc',1.0);
% CAT16 chromatic-adaptation matrix (XYZ -> sharpened R,G,B cone-like space).
M16 = [ 0.401288 0.650173 -0.051461; -0.250268 1.204414 0.045854; -0.002079 0.048952 0.953127 ];

fprintf('Loading results...\n');
% Load only the fields needed here (the full file is large).
Sres = load(inMat,'wls','params','display','cone_fund','illuminants','byIll');
wls    = Sres.wls(:);  nW = numel(wls);         % wavelength sampling grid
sdList = Sres.params.sdList(:)';                % all primary bandwidths analysed
% Display model: Pset{si} = 3x3 primaries matrix, Mdisp_set{si} = RGB->XYZ, per SD.
Pset   = Sres.display.Pset;  Mdisp_set = Sres.display.Mdisp_set;
Ill    = Sres.illuminants;  nIll = numel(Ill);
catVec = string({Ill.category});               % class label per illuminant
% Subsample observers evenly across the stored family (here: full 10000).
obsIdx = round(linspace(1, size(Sres.cone_fund,1), CFG.nObsSub));
% cone_fund is stored obs x wl x LMS; permute to 3(LMS) x nW x nObs for pagemtimes.
T_indiv = permute(Sres.cone_fund(obsIdx,:,:),[3 2 1]);       % 3 x nW x nObsSub
nObs = size(T_indiv,3);
% Per-illuminant bookkeeping copied out of byIll:
%   ok        = illuminant produced usable matches,
%   keepIdx   = indices of surfaces retained (in-gamut etc.),
%   idxWhite  = which retained sample is the white/reference,
%   XYZ_white = reference (standard-observer) white XYZ.
B = struct('ok',{},'keepIdx',{},'idxWhite',{},'XYZ_white',{});
for k=1:nIll
    r = Sres.byIll(k);
    B(k).ok = isfield(r,'ok') && ~isempty(r.ok) && r.ok;
    if ~B(k).ok, continue; end
    B(k).keepIdx=r.keepIdx(:); B(k).idxWhite=r.idxWhite; B(k).XYZ_white=r.XYZ_white(:);
end
% Drop the bulky byIll now that the needed fields are cached in B.
Sres = rmfield(Sres,'byIll'); clear r;

% CIE 1931 2 deg colour-matching functions (the standard-observer reference),
% resampled to wls; ' -> 3 x nW. interp1 ...,0 sets out-of-range to 0.
Axyz  = readmatrix(fullfile(dataDir,'CIE1931XYZ.txt'));
T_xyz = interp1(Axyz(:,1),Axyz(:,2:4),wls,'pchip',0)';       % 3 x nW
% Surface reflectance set (columns = surfaces); clamped to [0,1] after resample.
Rmat  = readmatrix(fullfile(dataDir,'CIE_srf_cfi_1nm.csv'));
R_all = min(max(interp1(Rmat(:,1),Rmat(:,2:end),wls,'pchip','extrap'),0),1);   % nW x 99
nG    = size(R_all,2);
% surface sRGB under D65 for ellipse colouring (illuminant-independent)
Dd    = load(fullfile(dataDir,'D65.mat'));
d65   = interp1(Dd.wls(:), Dd.d65(:), wls, 'pchip', 'extrap');
% Normalise D65 so a perfect diffuser gives Y = 100 (white-referenced XYZ).
d65   = max(d65,0) * (100/max(T_xyz(2,:)*max(d65,0),eps));
% Each surface's XYZ under D65, then converted to an sRGB swatch for colouring.
XYZ65 = T_xyz * (d65 .* R_all);
srgbG = xyz_to_srgb_d65(XYZ65);                              % nG x 3
% Precompute, per selected SD, the inverse of each observer's display response
% matrix C = T_indiv * Pset (LMS produced by the display primaries). invC maps a
% target LMS back to display RGB, i.e. it solves each observer's metameric match.
nSDsel=numel(CFG.SDsel); invC_set=cell(nSDsel,1); siStored=zeros(nSDsel,1);
for j=1:nSDsel
    % Locate the requested SD within the stored sdList; remember its index.
    si=find(sdList==CFG.SDsel(j),1); siStored(j)=si;
    C=pagemtimes(T_indiv,Pset{si}); iC=zeros(3,3,nObs);
    % Per-observer 3x3 inverse (loop because inv is not paged).
    for o=1:nObs, iC(:,:,o)=inv(C(:,:,o)); end
    invC_set{j}=iC;
end
% Indices of all usable illuminants (used to build class subsets below).
okAll = find(arrayfun(@(k)B(k).ok,1:nIll));
fprintf('Setup done. Run a group section (daylight / smooth / peaky / all).\n');

% Each section renders one illuminant class by passing the subset of usable
% illuminant indices whose category matches. Run SETUP once, then any section.

%% ============================ daylight ============================
render_ellipse_group('daylight', intersect(okAll,find(catVec=="daylight")), ...
    CFG,cam,M16,T_indiv,T_xyz,R_all,Ill,B,Pset,Mdisp_set,invC_set,siStored,srgbG,nG,figDir,figp);

%% ============================= smooth =============================
render_ellipse_group('smooth', intersect(okAll,find(catVec=="smooth")), ...
    CFG,cam,M16,T_indiv,T_xyz,R_all,Ill,B,Pset,Mdisp_set,invC_set,siStored,srgbG,nG,figDir,figp);

%% ============================= peaky ==============================
render_ellipse_group('peaky', intersect(okAll,find(catVec=="peaky")), ...
    CFG,cam,M16,T_indiv,T_xyz,R_all,Ill,B,Pset,Mdisp_set,invC_set,siStored,srgbG,nG,figDir,figp);

%% ======================== all illuminants =========================
% 'all' pools every usable illuminant regardless of class.
render_ellipse_group('all', okAll, ...
    CFG,cam,M16,T_indiv,T_xyz,R_all,Ill,B,Pset,Mdisp_set,invC_set,siStored,srgbG,nG,figDir,figp);

% ===================== local helpers =====================
function render_ellipse_group(tag, idxC, CFG,cam,M16,T_indiv,T_xyz,R_all,Ill,B,Pset,Mdisp_set,invC_set,siStored,srgbG,nG,figDir,figp)
if isempty(idxC), warning('plot_fig6:noIll','no illuminants for "%s"',tag); return; end
% Panel order: 1 = D=0, 2 = D=0.5, 3 = D=1. Not titled on the figure -- the panels
% are identified in the caption.
Dlist = CFG.Dlist(:)'; nDp = numel(Dlist);
% One figure per selected SD (narrow / broad); each has nDp panels (one per D).
for j = 1:numel(CFG.SDsel)
    % Display RGB->XYZ and per-observer inverse response for this SD.
    SD = CFG.SDsel(j); Mdisp = Mdisp_set{siStored(j)}; invC = invC_set{j};
    % dev{iD}{g}: 2 x (obs*illum) accumulator of a'b' deviations for surface g at D level iD.
    dev = repmat({cell(nG,1)},1,nDp);                    % pooled deviations per D
    % Running sum/count of each surface's reference a'b' position (for the centre).
    refSum = zeros(nG,2); refN = zeros(nG,1);
    for k = idxC
        keepIdx=B(k).keepIdx; idxWhite=B(k).idxWhite; Wref=B(k).XYZ_white;
        % Colour signal = illuminant SPD times each retained surface reflectance.
        SPD = Ill(k).spd(:) .* R_all(:,keepIdx);             % nW x nSamp
        % Reference (standard-observer) tristimulus of the surfaces and the white.
        XYZw = T_xyz*SPD;                                    % 3 x nSamp (world/reference XYZ)
        % Each observer's cone response to the same signals.
        LMS  = pagemtimes(T_indiv, SPD);                     % 3 x nSamp x nObs
        % Observer's metameric DISPLAY match: invC turns the desired LMS into
        % display RGB, then Mdisp turns that RGB back into XYZ (what each observer
        % would actually set on the display to match the reference stimulus).
        XYZ  = pagemtimes(Mdisp, pagemtimes(invC, LMS));     % 3 x nSamp x nObs
        % Each observer's reproduced white = their match at the white sample.
        Wown = XYZ(:,idxWhite,:);
        % Placement centre: reference a'b' at full adaptation (D=1), so all panels
        % share one centre per surface and only the SPREAD changes with D.
        rab1 = getab(cam16_s3(XYZw, Wref, 1, M16, cam));     % D=1 ref -> common placement
        refD = cell(nDp,1); obsD = cell(nDp,1);
        for iD = 1:nDp
            % Reference a'b' at this D (adopted white = reference white Wref).
            refD{iD} = getab(cam16_s3(XYZw, Wref, Dlist(iD), M16, cam));   % 2 x nSamp
            % Observer a'b' at this D (each adapts to their OWN reproduced white).
            obsD{iD} = getab(cam16_s3(XYZ,  Wown, Dlist(iD), M16, cam));   % 2 x nSamp x nObs
        end
        gid = keepIdx;
        for l = 1:numel(gid)
            g = gid(l);
            for iD = 1:nDp
                % Deviation of every observer from the D-matched reference for this
                % surface; appended across illuminants -> pooled between-observer spread.
                dev{iD}{g} = [dev{iD}{g}, squeeze(obsD{iD}(:,l,:)) - refD{iD}(:,l)];
            end
            % Accumulate the shared (D=1) reference position for this surface.
            refSum(g,:) = refSum(g,:) + rab1(:,l)'; refN(g)=refN(g)+1;
        end
    end
    % Mean reference a'b' per present surface = ellipse centre for the panels.
    refPos = nan(nG,2); pres = find(refN>0); refPos(pres,:) = refSum(pres,:)./refN(pres);
    % Choose up to nSurf well-separated surfaces (>=0.5 apart) so ellipses don't overlap.
    loc = subsample_by_distance_farthest(refPos(pres,:), 0.5, CFG.nSurf); keptG = pres(loc);

    % Build the 95% covariance ellipse for each kept surface at every D, and
    % collect all vertices to compute one shared square axis range.
    E = cell(nDp,1); allx=[]; ally=[];
    for m=1:nDp
        [Xe,Ye,Col,Ar] = ellipses_for(keptG, refPos, dev{m}, srgbG);
        E{m}=struct('X',{Xe},'Y',{Ye},'C',{Col},'A',{Ar});
        allx=[allx; vertcat(Xe{:})]; ally=[ally; vertcat(Ye{:})];
    end
    % Square, centred limits covering all ellipses across all D panels (+1 margin).
    allx=allx(isfinite(allx)); ally=ally(isfinite(ally));
    cx=mean([min(allx) max(allx)]); cy=mean([min(ally) max(ally)]);
    half=max(max(allx)-min(allx), max(ally)-min(ally))/2 + 1;
    xlm=[cx-half cx+half]; ylm=[cy-half cy+half];

    fig=figure('Color','w','ToolBar','none','MenuBar','none');
    tl=tiledlayout(fig,1,nDp,'Padding','compact','TileSpacing','tight');
    for m=1:nDp
        ax=nexttile(tl); hold(ax,'on');
        % Draw largest ellipses first so small ones stay visible on top.
        [~,ord]=sort(E{m}.A,'descend','MissingPlacement','last');
        fillCols=desaturate_rgb(E{m}.C,0.85);   % keep most of the surface's saturation
        edgeCols=max(E{m}.C*0.72,0);   % darker outline: pale surfaces washed out at full sRGB
        for kk=1:numel(ord)
            ii=ord(kk); if isempty(E{m}.X{ii})||~all(isfinite(E{m}.X{ii})), continue; end
            % Translucent fill in the surface colour.
            patch('Parent',ax,'XData',E{m}.X{ii},'YData',E{m}.Y{ii}, ...
                'FaceColor',fillCols(ii,:),'FaceAlpha',0.20,'EdgeColor','none');
            % Darker outline of the same ellipse.
            plot(ax,E{m}.X{ii},E{m}.Y{ii},'-','LineWidth',0.6,'Color',edgeCols(ii,:));
        end
        % adopted white point: maps to the CAM16-UCS achromatic origin (a'=b'=0)
        plot(ax,0,0,'+','Color','k','MarkerSize',8,'LineWidth',1.0);
        xlim(ax,xlm); ylim(ax,ylm); axis(ax,'square');
        xticks(ax,-100:20:100); yticks(ax,-100:20:100);   % 0-centred, 20-unit spacing
        % Only the leftmost panel carries the b' label; every panel gets a'.
        if m==1, ylabel(ax,"b'",'FontWeight','Bold'); end
        xlabel(ax,"a'",'FontWeight','Bold');
        ax.FontName='Arial'; ax.Color=ones(3,1)*0.97; ax.FontSize=figp.fontsize;
        ax.XColor='k'; ax.YColor='k'; ax.LineWidth=0.5;
        grid(ax,'on'); ax.GridAlpha=0.08; box(ax,'off'); ticklengthcm(ax,0.08);
    end
    % --- Absolute sizing: the nDp panels fill the journal two-column width -----
    fig.Units='centimeters'; fig.InvertHardcopy='off';
    panelScale = 1.0;                                 % 3 panels exactly fill the two-column width
    panelW = panelScale*figp.twocolumn/3;             % width per panel
    fig.Position=[10,10, nDp*panelW, panelScale*figp.twocolumn/3.1];   % height scaled to match
    % Let the renderer settle before exporting the vector PDF.
    pause(0.5)
    exportgraphics(fig, fullfile(figDir,sprintf('fig6_%s_uvEllipse_SD%02d.pdf',tag,SD)),'ContentType','vector');
    close(fig);
end
fprintf('Fig 6 group "%s" done (%d SDs x %d D levels).\n', tag, numel(CFG.SDsel), nDp);
end

% Extract the two chromatic coordinates (a',b') from a stacked [J' a' b'] array;
% works for 2-D (samples) or 3-D (samples x observers) inputs.
function ab = getab(s3)   % rows 2:3 (a',b') of a CAM16-UCS array
ab = s3(2:3,:,:);
end

% CIECAM16 forward transform -> CAM16-UCS [J' a' b'] using the STANDARD (CIE)
% chromatic-adaptation form with degree D (CAT16 white discounting).
% XYZ: 3 x ... stimulus tristimulus; W: 3 x 1 x ... adopted white; D in [0,1].
function s3 = cam16_s3(XYZ, W, D, M16, cam)
% Unpack viewing conditions; Yw = luminance of the adopted white.
LA=cam.LA; Yb=cam.Yb; c=cam.c; Nc=cam.Nc; Yw=W(2,1,:);
% Luminance-level adaptation factor FL (controls the nonlinearity strength).
k=1/(5*LA+1); FL=0.2*k^4*(5*LA)+0.1*(1-k^4)^2*(5*LA)^(1/3);
% Background induction terms: n = relative background, z = base-exponent,
% Nbb/Ncb = brightness/chromatic background induction factors.
n=max(Yb./Yw,1e-6); z=1.48+sqrt(n); Nbb=0.725.*(1./n).^0.2; Ncb=Nbb;
% CAT16: transform white to sharpened RGB, form per-channel adaptation gain Dr
% (D=0 -> no discounting Dr=1; D=1 -> full discounting to equal-energy), apply it.
RGBw=pagemtimes(M16,W); Dr=D.*Yw./RGBw+(1-D); RGBcw=Dr.*RGBw;
% Post-adaptation nonlinear (Naka-Rushton-type) compression of the white.
tw=(FL.*abs(RGBcw)/100).^0.42; RGBaw=sign(RGBcw).*400.*tw./(27.13+tw)+0.1;
% Achromatic response of the white (reference for lightness J).
Aw=(2*RGBaw(1,1,:)+RGBaw(2,1,:)+RGBaw(3,1,:)/20-0.305).*Nbb;
% Same adaptation gain Dr and compression applied to the stimulus.
RGB=pagemtimes(M16,XYZ); RGBc=Dr.*RGB;
t=(FL.*abs(RGBc)/100).^0.42; RGBa=sign(RGBc).*400.*t./(27.13+t)+0.1;
% Convert compressed cone signals to CAM16-UCS coordinates.
s3=ucs_from_compressed(RGBa,Aw,Nbb,Ncb,Nc,c,z,n,FL);
end

% ALTERNATIVE forward form (post-compression gain instead of pre-compression
% CAT16). NOT used to draw the figure -- kept for the algebraic-identity check
% discussed in the header: normalising after the near-power-law compression is
% equivalent to the standard form up to a constant.
function s3 = cam16_nl(XYZ, W, M16, cam)
LA=cam.LA; Yb=cam.Yb; c=cam.c; Nc=cam.Nc; Yw=W(2,1,:);
k=1/(5*LA+1); FL=0.2*k^4*(5*LA)+0.1*(1-k^4)^2*(5*LA)^(1/3);
n=max(Yb./Yw,1e-6); z=1.48+sqrt(n); Nbb=0.725.*(1./n).^0.2; Ncb=Nbb;
% Compress first, then normalise by the compressed white (gain applied AFTER).
RGBw=pagemtimes(M16,W); aw=nr_compress(RGBw,FL);
tgt=nr_compress(repmat(Yw,3,1,1),FL); gain=tgt./aw;
Aw=(2*tgt(1,1,:)+tgt(2,1,:)+tgt(3,1,:)/20-0.305).*Nbb;
RGBa=nr_compress(pagemtimes(M16,XYZ),FL).*gain;
s3=ucs_from_compressed(RGBa,Aw,Nbb,Ncb,Nc,c,z,n,FL);
end

% Map compressed cone signals RGBa to CAM16-UCS [J' a' b'].
% Aw = white achromatic response; the remaining args are the induction/surround
% factors computed in the forward transform above.
function s3 = ucs_from_compressed(RGBa,Aw,Nbb,Ncb,Nc,c,z,n,FL)
R1=RGBa(1,:,:); R2=RGBa(2,:,:); R3=RGBa(3,:,:);
% Opponent chromatic signals a,b and achromatic response A.
a=R1-12*R2/11+R3/11; b=(R1+R2-2*R3)/9; A=(2*R1+R2+R3/20-0.305).*Nbb;
% Lightness J (relative to the white's achromatic response Aw).
J=100.*sign(A).*(abs(A)./Aw).^(c.*z);
% Hue angle h and its eccentricity weighting et.
h=atan2d(b,a); et=0.25*(cos(deg2rad(h)+2)+3.8);
% Temporary magnitude t -> chroma C -> colourfulness M (CIECAM16 chain).
denC=R1+R2+(21/20)*R3; tC=(50000/13*Nc.*Ncb).*et.*hypot(a,b)./max(denC,1e-9);
Cc=(max(tC,0)).^0.9.*sqrt(max(J,0)/100).*(1.64-0.29.^n).^0.73; Mc=Cc.*FL^0.25;
% CAM16-UCS compression of lightness (J') and colourfulness (M'); a'=M'cos h,
% b'=M'sin h. These are the perceptually-uniform coordinates dE' is measured in.
Jp=(1+100*0.007).*J./(1+0.007.*J); Mp=(1/0.0228)*log(1+0.0228.*Mc);
s3=cat(1,Jp,Mp.*cosd(h),Mp.*sind(h));
end

% Naka-Rushton hyperbolic compression used by the alternative forward form
% (cam16_nl); FL is the luminance-adaptation factor. Same nonlinearity applied
% inline inside cam16_s3.
function a = nr_compress(v,FL)
t=(FL.*abs(v)/100).^0.42; a=sign(v).*400.*t./(27.13+t)+0.1;
end

% Build a 95% covariance ellipse polygon for each kept surface. Each surface's
% cloud = its reference centre refPos(g,:) plus the pooled observer deviations
% devCell{g}. Also returns the surface sRGB colour and the ellipse area.
function [Xe,Ye,col,area] = ellipses_for(keptG, refPos, devCell, srgbG)
n=numel(keptG); Xe=cell(n,1); Ye=cell(n,1); col=nan(n,3); area=nan(n,1);
for ii=1:n
    g=keptG(ii); col(ii,:)=srgbG(g,:);
    if isempty(devCell{g})||any(~isfinite(refPos(g,:))), continue; end
    % Reconstruct the observer point cloud (centre + deviations, one row each).
    P = refPos(g,:) + devCell{g}.';
    [xe,ye,ok]=ellipse_poly(P,0.95,120);
    % polyarea gives the ellipse size, used to order drawing (large first).
    if ok, Xe{ii}=xe(:); Ye{ii}=ye(:); area(ii)=polyarea(xe,ye); end
end
end

% Farthest-point sampling: greedily pick up to maxKeep points that are mutually
% well separated (stop once the next-farthest is closer than minDist). Keeps the
% displayed surfaces spread out so their ellipses do not pile up.
function idxKeep = subsample_by_distance_farthest(P,minDist,maxKeep)
ok=all(isfinite(P),2); idxAll=find(ok); P2=P(ok,:); n=size(P2,1);
if n==0, idxKeep=[]; return; end
% Seed with the point farthest from the centroid.
mu=mean(P2,1); d2=sum((P2-mu).^2,2); [~,i0]=max(d2); sel=i0; md2=minDist^2;
while numel(sel)<min(maxKeep,n)
    bi=-1; bd=-inf;
    % Find the unselected point whose nearest selected neighbour is farthest away.
    for i=1:n
        if any(sel==i), continue; end
        di=inf; for s=sel, v=P2(i,:)-P2(s,:); di=min(di,v(1)^2+v(2)^2); end
        if di>bd, bd=di; bi=i; end
    end
    % Stop if the best candidate is still within the minimum spacing.
    if bi<0||bd<md2, break; end
    sel(end+1)=bi; %#ok<AGROW>
end
% Map back to original (pre-mask) indices.
idxKeep=idxAll(sel);
end

% Covariance-error ellipse of a 2-D point set P at confidence `conf`, returned as
% an nPts polygon. Standard construction: eigen-decompose the covariance, scale
% axes by sqrt(chi-square quantile), rotate by eigenvectors, translate to mean.
function [xe,ye,ok]=ellipse_poly(P,conf,nPts)
ok=true; P=P(all(isfinite(P),2),:);
% Need at least a few points for a meaningful covariance.
if size(P,1)<5, xe=nan; ye=nan; ok=false; return; end
% Mean (centre), covariance C, and the chi-square(2 dof) scale for `conf`.
mu=mean(P,1); C=cov(P); q=-2*log(max(1-conf,eps));
% Eigenvectors V give the axis directions; radii r = sqrt(q)*sqrt(eigenvalues).
[V,Dg]=eig(C); r=sqrt(q)*sqrt(max(diag(Dg),0));
% Parametric unit circle -> scaled -> rotated -> shifted to the mean.
t=linspace(0,2*pi,nPts); xy=(V*diag(r)*[cos(t);sin(t)])+mu(:);
xe=xy(1,:); ye=xy(2,:);
end

% Blend an sRGB colour toward its luminance-matched grey, keeping fraction
% keepSat of the original saturation (used for the pale ellipse fills).
function rgb2=desaturate_rgb(rgb,keepSat)
keepSat=max(0,min(1,keepSat));
% Rec.709 luminance as the grey anchor.
Y=0.2126*rgb(:,1)+0.7152*rgb(:,2)+0.0722*rgb(:,3); g=[Y Y Y];
rgb2=min(max(g+keepSat*(rgb-g),0),1);
end

% Convert XYZ to display-ready sRGB assuming a D65 white (used only to colour
% the surface swatches, not for any colorimetry).
function rgb = xyz_to_srgb_d65(XYZ)
% Normalise by the brightest Y so the whole set fits the display range.
Yw=max(XYZ(2,:)); if isempty(Yw)||Yw<=0, Yw=100; end
XYZn=XYZ/Yw;
% XYZ (D65) -> linear sRGB, clamped to [0,1].
M=[3.2406 -1.5372 -0.4986; -0.9689 1.8758 0.0415; 0.0557 -0.2040 1.0570];
lin=max(min((M*XYZn)',1),0); a=0.055; thr=0.0031308;
% sRGB opto-electronic transfer function (gamma): linear segment below thr,
% power-law above.
rgb=zeros(size(lin)); lo=lin<=thr; rgb(lo)=12.92*lin(lo); rgb(~lo)=(1+a)*lin(~lo).^(1/2.4)-a;
rgb=max(min(rgb,1),0);
end
