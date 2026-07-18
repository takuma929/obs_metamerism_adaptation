% ========================================================================
% stats_ellipse_and_daylight.m
%
% PURPOSE
%   (1) Between-observer ellipse area vs degree of adaptation (Fig. 5): for
%       every surface, form the 95% covariance ellipse of the CAM16-UCS (a',b')
%       observer deviations from the CIE 1931 reference at each D, and test that
%       its area shrinks with adaptation (paired signrank D=0 vs D=1 / D=0.5).
%   (2) Per-observer daylight D=0 vs D=1 test on the full 10,000 observers,
%       averaged over the 12 daylight phases (CCTs), using the full CAM16-UCS
%       distance (J'a'b').
%   (3) Reports the data counts feeding each ellipse.
%
% INPUT (.mat)
%   results/results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat, fields:
%     wls, params.sdList, display.Pset/.Mdisp_set (per-bandwidth display
%     primaries + display->XYZ matrix), cone_fund (obs x wl x 3 individual cone
%     fundamentals), illuminants (spd + category), byIll (per-illuminant kept
%     surfaces keepIdx, white index idxWhite, reference white XYZ_white).
%   data/CIE1931XYZ.txt   reference 2 deg colour-matching functions.
%   data/CIE_srf_cfi_1nm.csv  CIE 99 CFI surface reflectances (1 nm).
%
% PRINTS
%   observer/illuminant/point counts; median ellipse area at each D; paired
%   signrank D0-vs-D1 and D0-vs-D0.5 with z and area ratio; daylight per-observer
%   D0-vs-D1 median dE', signrank p/z and Cohen's dz.
%
% METHOD
%   Each observer's metameric display match is the display XYZ that reproduces,
%   under this observer's cone fundamentals, the same LMS the surface produces:
%   XYZ = Mdisp * inv(C) * LMS with C = T_indiv * Pset (per-observer display
%   response). Reference and observer are both put through CIECAM16/CAT16 at the
%   same D (cam16_s3 helper) and compared in CAM16-UCS. Observer metamerism =
%   spread of these matches across observers; ellipse area quantifies it in the
%   (a',b') plane. NB the paper reports SD = 5 nm (narrow primaries) here.
% ========================================================================

% (1) Ellipse area vs degree of adaptation (Fig 5): paired test D=0 vs D=1.
% (2) Per-observer daylight D=0 vs D=1 test (n=10,000, averaged over 12 phases).
% (3) Data counts per ellipse.
clearvars; clc;
projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
inMat = fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
% CIECAM16 viewing-condition parameters (adapting luminance LA, background Yb,
% surround factors F/c/Nc for an average surround).
cam = struct('LA',20,'Yb',20,'F',1.0,'c',0.69,'Nc',1.0);
% CAT16 chromatic-adaptation transform matrix (XYZ -> sharpened RGB).
M16 = [0.401288 0.650173 -0.051461; -0.250268 1.204414 0.045854; -0.002079 0.048952 0.953127];
% Adaptation degrees to evaluate; SDsel = analysed primary bandwidth (nm).
Dlist=[0 0.5 1]; SDsel=5;

% Load simulation inputs/results.
Sres = load(inMat,'wls','params','display','cone_fund','illuminants','byIll');
wls=Sres.wls(:); sdList=Sres.params.sdList(:)';
% Per-bandwidth display primaries (Pset) and display->XYZ matrices (Mdisp_set).
Pset=Sres.display.Pset; Mdisp_set=Sres.display.Mdisp_set;
Ill=Sres.illuminants; nIll=numel(Ill);
% Individual cone fundamentals reshaped to 3 x nW x nObs (all 10,000 observers).
T_indiv=permute(Sres.cone_fund,[3 2 1]); nObs=size(T_indiv,3);   % FULL 10,000
% Per-illuminant bookkeeping: whether it converged (ok), which surfaces were
% retained (keepIdx), the display-white surface index and its reference XYZ.
B=struct('ok',{},'keepIdx',{},'idxWhite',{},'XYZ_white',{});
for k=1:nIll
    r=Sres.byIll(k); B(k).ok=isfield(r,'ok')&&~isempty(r.ok)&&r.ok;
    if B(k).ok, B(k).keepIdx=r.keepIdx(:); B(k).idxWhite=r.idxWhite; B(k).XYZ_white=r.XYZ_white(:); end
end
% Reference CIE 1931 2 deg CMFs, interpolated to the working wavelength grid.
Axyz=readmatrix(fullfile(dataDir,'CIE1931XYZ.txt'));
T_xyz=interp1(Axyz(:,1),Axyz(:,2:4),wls,'pchip',0)';
% CIE 99 CFI surface reflectances, interpolated and clamped to [0,1].
Rmat=readmatrix(fullfile(dataDir,'CIE_srf_cfi_1nm.csv'));
R_all=min(max(interp1(Rmat(:,1),Rmat(:,2:end),wls,'pchip','extrap'),0),1);
nG=size(R_all,2);
% Select the display model for the analysed bandwidth (SDsel nm).
si=find(sdList==SDsel,1); Mdisp=Mdisp_set{si};
% Per-observer display response C = T_indiv * Pset, and its inverse (used to
% solve for the display XYZ that metamerically matches each surface).
C=pagemtimes(T_indiv,Pset{si}); invC=zeros(3,3,nObs);
for o=1:nObs, invC(:,:,o)=inv(C(:,:,o)); end
% Indices of the illuminants that converged.
okAll=find(arrayfun(@(k)B(k).ok,1:nIll));
fprintf('n observers = %d, n illuminants (all, ok) = %d\n', nObs, numel(okAll));
% Each pooled ellipse draws one point per (observer x illuminant).
fprintf('points per ellipse (all group) = %d x %d = %d\n', nObs, numel(okAll), nObs*numel(okAll));

% gather deviations per D, pooled over all illuminants
% dev{iD}{g} accumulates the (a',b') observer-minus-reference deviations for
% surface g at adaptation degree D(iD), across all illuminants (2 x nPoints).
nD=numel(Dlist);
dev=repmat({cell(nG,1)},1,nD); refN=zeros(nG,1);
for k=okAll
    % Kept surfaces, display-white index, and reference white for this illuminant.
    keepIdx=B(k).keepIdx; idxWhite=B(k).idxWhite; Wref=B(k).XYZ_white;
    % Colour signals of each retained surface: SPD = illuminant * reflectance.
    % XYZw = reference tristimulus of those surfaces under the CIE 1931 CMFs.
    SPD=Ill(k).spd(:).*R_all(:,keepIdx); XYZw=T_xyz*SPD;
    % Per-observer metameric display match: LMS -> display XYZ via Mdisp*inv(C).
    LMS=pagemtimes(T_indiv,SPD); XYZ=pagemtimes(Mdisp,pagemtimes(invC,LMS));
    % Each observer's adopted (display) white = reproduced display-white surface.
    Wown=XYZ(:,idxWhite,:);
    for iD=1:nD
        % CAM16-UCS (a',b') of reference and observers at this adaptation degree.
        refD=getab(cam16_s3(XYZw,Wref,Dlist(iD),M16,cam));
        obsD=getab(cam16_s3(XYZ, Wown,Dlist(iD),M16,cam));
        for l=1:numel(keepIdx)
            % Append this illuminant's per-observer deviations for surface g.
            g=keepIdx(l);
            dev{iD}{g}=[dev{iD}{g}, squeeze(obsD(:,l,:))-refD(:,l)];
        end
    end
    % Count how many illuminants contributed to each surface.
    for l=1:numel(keepIdx), refN(keepIdx(l))=refN(keepIdx(l))+1; end
end
% Surfaces that appeared under at least one illuminant.
pres=find(refN>0);
% 95% covariance-ellipse area per surface (rows) and adaptation degree (cols).
area=nan(nG,nD);
for iD=1:nD
    for g=pres'
        if isempty(dev{iD}{g}), continue; end
        % Fit the 95% ellipse to the (a',b') deviation cloud (120-vertex polygon).
        [xe,ye,ok]=ellipse_poly(dev{iD}{g}.',0.95,120);
        % Ellipse area = area of the fitted polygon.
        if ok, area(g,iD)=polyarea(xe,ye); end
    end
end
% Keep only surfaces with a valid ellipse at every D (paired comparison).
val=find(all(isfinite(area),2));
fprintf('\n--- Ellipse area vs D (n=%d surfaces) ---\n', numel(val));
% Median ellipse area at each adaptation degree.
for iD=1:nD, fprintf('  D=%g: median area = %.1f\n', Dlist(iD), median(area(val,iD))); end
% Paired signrank of ellipse area, unadapted (col 1, D=0) vs fully adapted (col 3, D=1).
[p01,~,s01]=signrank(area(val,1),area(val,3));
fprintf('  paired signrank D=0 vs D=1: p=%.2e, z=%.2f, median ratio D1/D0 = %.2f\n', ...
    p01, s01.zval, median(area(val,3))./median(area(val,1)));
% Paired signrank D=0 vs partial adaptation D=0.5 (col 2).
[p005,~,s005]=signrank(area(val,1),area(val,2));
fprintf('  paired signrank D=0 vs D=0.5: p=%.2e, z=%.2f\n', p005, s005.zval);

% (2) per-observer daylight test (average over 12 daylight phases), n=10,000
% Restrict to the daylight illuminants (the 12 daylight-class CCTs/phases).
catVec=string({Ill.category}); dl=intersect(okAll,find(catVec=="daylight"));
fprintf('\n--- daylight per-observer D=0 vs D=1 (n=%d obs, avg over %d phases), SD=5 ---\n', nObs, numel(dl));
% dE(o,iD): observer o's metamerism at degree D(iD), summed over phases (then averaged).
dE=zeros(nObs,nD);
cnt=0;
for k=dl
    keepIdx=B(k).keepIdx; idxWhite=B(k).idxWhite; Wref=B(k).XYZ_white;
    % Same colour-signal and metameric-match construction as part (1).
    SPD=Ill(k).spd(:).*R_all(:,keepIdx); XYZw=T_xyz*SPD;
    LMS=pagemtimes(T_indiv,SPD); XYZ=pagemtimes(Mdisp,pagemtimes(invC,LMS));
    Wown=XYZ(:,idxWhite,:);
    for iD=1:nD
        % Full CAM16-UCS coordinates (J',a',b') of reference and observers.
        ref=cam16_s3(XYZw,Wref,Dlist(iD),M16,cam);
        obs=cam16_s3(XYZ, Wown,Dlist(iD),M16,cam);
        % Per-observer median CAM16-UCS distance over surfaces (Euclidean dE').
        dd=squeeze(median(sqrt(sum((obs-ref).^2,1)),2,'omitnan'));  % per obs median over surf
        % Accumulate across daylight phases.
        dE(:,iD)=dE(:,iD)+dd(:);
    end
    cnt=cnt+1;
end
% Average the accumulated dE' over the daylight phases.
dE=dE/cnt;
% Paired per-observer signrank D=0 (col 1) vs D=1 (col 3), with Cohen's dz.
[pd,~,sd_]=signrank(dE(:,1),dE(:,3)); dz=mean(dE(:,1)-dE(:,3))/std(dE(:,1)-dE(:,3));
fprintf('  median D0=%.2f D1=%.2f, signrank p=%.2e, z=%.1f, dz=%.2f\n', ...
    median(dE(:,1)), median(dE(:,3)), pd, sd_.zval, dz);

% --- local helpers ---
% getab: keep only the (a',b') chroma coordinates from CAM16-UCS J'a'b' (drop J').
function ab=getab(s3), ab=s3(2:3,:,:); end
% cam16_s3: CIECAM16 forward transform to CAM16-UCS J'a'b'.
%   XYZ (3 x nSurf x nObs) test tristimuli, W adopted white, D adaptation degree,
%   M16 CAT16 matrix, cam viewing conditions. Returns [J'; a'; b'] (page-wise).
function s3=cam16_s3(XYZ,W,D,M16,cam)
% Unpack viewing conditions; Yw = luminance of the adopted white.
LA=cam.LA;Yb=cam.Yb;c=cam.c;Nc=cam.Nc;Yw=W(2,1,:);
% Luminance-level adaptation factor FL.
k=1/(5*LA+1);FL=0.2*k^4*(5*LA)+0.1*(1-k^4)^2*(5*LA)^(1/3);
% Background induction factors: n (relative bg luminance), z, Nbb=Ncb.
n=max(Yb./Yw,1e-6);z=1.48+sqrt(n);Nbb=0.725.*(1./n).^0.2;Ncb=Nbb;
% CAT16 adaptation of the white: Dr = per-channel von Kries gain at degree D.
RGBw=pagemtimes(M16,W);Dr=D.*Yw./RGBw+(1-D);RGBcw=Dr.*RGBw;
% Post-adaptation nonlinear (Hunt-Pointer-Estevez) response of the white.
tw=(FL.*abs(RGBcw)/100).^0.42;RGBaw=sign(RGBcw).*400.*tw./(27.13+tw)+0.1;
% Achromatic response of the white, Aw (denominator of lightness J).
Aw=(2*RGBaw(1,1,:)+RGBaw(2,1,:)+RGBaw(3,1,:)/20-0.305).*Nbb;
% Same adaptation + nonlinearity applied to the test stimuli.
RGB=pagemtimes(M16,XYZ);RGBc=Dr.*RGB;
t=(FL.*abs(RGBc)/100).^0.42;RGBa=sign(RGBc).*400.*t./(27.13+t)+0.1;
R1=RGBa(1,:,:);R2=RGBa(2,:,:);R3=RGBa(3,:,:);
% Opponent signals a,b and achromatic response Ac of the test stimuli.
a=R1-12*R2/11+R3/11;b=(R1+R2-2*R3)/9;Ac=(2*R1+R2+R3/20-0.305).*Nbb;
% Lightness J (relative to white's Aw).
J=100.*sign(Ac).*(abs(Ac)./Aw).^(c.*z);
% Hue angle hh and eccentricity factor et.
hh=atan2d(b,a);et=0.25*(cos(deg2rad(hh)+2)+3.8);
% Chroma C and colourfulness M.
denC=R1+R2+(21/20)*R3;tC=(50000/13*Nc.*Ncb).*et.*hypot(a,b)./max(denC,1e-9);
Cc=(max(tC,0)).^0.9.*sqrt(max(J,0)/100).*(1.64-0.29.^n).^0.73;Mc=Cc.*FL^0.25;
% UCS compression: J' from J, M' from M; a'=M'cos(h), b'=M'sin(h).
Jp=(1+100*0.007).*J./(1+0.007.*J);Mp=(1/0.0228)*log(1+0.0228.*Mc);
s3=cat(1,Jp,Mp.*cosd(hh),Mp.*sind(hh));
end
% ellipse_poly: sample the conf-level covariance ellipse of 2-D points P
%   (rows = observations). Returns polygon vertices (xe,ye) and ok flag.
function [xe,ye,ok]=ellipse_poly(P,conf,nPts)
% Drop non-finite rows; need >=5 points for a meaningful covariance.
ok=true;P=P(all(isfinite(P),2),:);
if size(P,1)<5,xe=nan;ye=nan;ok=false;return;end
% Mean, covariance, and chi-square radius^2 for the requested confidence.
mu=mean(P,1);C=cov(P);q=-2*log(max(1-conf,eps));
% Eigen-decompose C: eigenvectors give axes, semi-axis lengths = sqrt(q*eig).
[V,Dg]=eig(C);r=sqrt(q)*sqrt(max(diag(Dg),0));
% Trace the ellipse: rotate/scale a unit circle and shift to the mean.
t=linspace(0,2*pi,nPts);xy=(V*diag(r)*[cos(t);sin(t)])+mu(:);
xe=xy(1,:);ye=xy(2,:);
end
