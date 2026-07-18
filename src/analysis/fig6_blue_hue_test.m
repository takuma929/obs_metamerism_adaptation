% ========================================================================
% fig6_blue_hue_test.m
% Test whether between-observer metamerism ellipses (CAM16-UCS a'b') are
% larger for BLUE surfaces, at the unadapted state D = 0 and narrow primaries
% (SD = 5 nm), pooled over all illuminants -- the regime of Fig. 6.
%
% For every surface we form the 95% between-observer covariance ellipse of the
% observer deviations from the standard-observer reference (same as Fig. 6),
% take its area, and its hue h = atan2(b',a') and blue-axis coordinate b' at
% the fully-adapted (D=1) reference placement. We then test:
%   (i)   Spearman correlation between ellipse area and b'  (blue = b'<0)
%   (ii)  circular-linear correlation between hue and log-area
%   (iii) blue hue sector (h in [230,290] deg) vs the rest (rank-sum)
% Output: prints manuscript-ready statistics.
% ========================================================================
clearvars; close all; clc;
projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
inMat = fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');

% Analysis configuration and CIECAM16 constants.
CFG.SDsel   = 5;             % narrow primaries only
CFG.nObsSub = 1500;          % observers subsampled (ellipse shape is stable)
% CIECAM16 average-surround viewing conditions.
cam = struct('LA',20,'Yb',20,'F',1.0,'c',0.69,'Nc',1.0);
% CAT16 adaptation matrix (XYZ -> sharpened RGB).
M16 = [ 0.401288 0.650173 -0.051461; -0.250268 1.204414 0.045854; -0.002079 0.048952 0.953127 ];

% Load simulation inputs/results.
Sres = load(inMat,'wls','params','display','cone_fund','illuminants','byIll');
wls    = Sres.wls(:);  nW = numel(wls);
sdList = Sres.params.sdList(:)';
% Per-bandwidth display primaries and display->XYZ matrices.
Pset   = Sres.display.Pset;  Mdisp_set = Sres.display.Mdisp_set;
Ill    = Sres.illuminants;  nIll = numel(Ill);
catVec = string({Ill.category});
% Evenly subsample observers, reshaped to 3 x nW x nObs cone fundamentals.
obsIdx = round(linspace(1, size(Sres.cone_fund,1), CFG.nObsSub));
T_indiv = permute(Sres.cone_fund(obsIdx,:,:),[3 2 1]);   % 3 x nW x nObs
nObs = size(T_indiv,3);
% Per-illuminant bookkeeping (converged?, kept surfaces, white index/XYZ).
B = struct('ok',{},'keepIdx',{},'idxWhite',{},'XYZ_white',{});
for k=1:nIll
    r = Sres.byIll(k);
    B(k).ok = isfield(r,'ok') && ~isempty(r.ok) && r.ok;
    if ~B(k).ok, continue; end
    B(k).keepIdx=r.keepIdx(:); B(k).idxWhite=r.idxWhite; B(k).XYZ_white=r.XYZ_white(:);
end
% Reference CIE 1931 2 deg CMFs on the working grid.
Axyz  = readmatrix(fullfile(dataDir,'CIE1931XYZ.txt'));
T_xyz = interp1(Axyz(:,1),Axyz(:,2:4),wls,'pchip',0)';
% CIE 99 CFI reflectances, interpolated and clamped to [0,1].
Rmat  = readmatrix(fullfile(dataDir,'CIE_srf_cfi_1nm.csv'));
R_all = min(max(interp1(Rmat(:,1),Rmat(:,2:end),wls,'pchip','extrap'),0),1);
nG    = size(R_all,2);

% Display model for the narrow bandwidth; per-observer display response C and inv(C).
si = find(sdList==CFG.SDsel,1); Mdisp = Mdisp_set{si};
C=pagemtimes(T_indiv,Pset{si}); invC=zeros(3,3,nObs);
for o=1:nObs, invC(:,:,o)=inv(C(:,:,o)); end
% Indices of converged illuminants.
okAll = find(arrayfun(@(k)B(k).ok,1:nIll));

% -------- gather deviations at D=0 across ALL illuminants --------
% dev{g}: (a',b') observer-minus-reference deviations for surface g at D=0.
% refSum/refN accumulate the D=1 reference (a',b') position used to place hue.
dev = cell(nG,1); refSum=zeros(nG,2); refN=zeros(nG,1);
for k = okAll
    keepIdx=B(k).keepIdx; idxWhite=B(k).idxWhite; Wref=B(k).XYZ_white;
    % Colour signal and per-observer metameric display match (see other scripts).
    SPD = Ill(k).spd(:) .* R_all(:,keepIdx);
    XYZw = T_xyz*SPD;
    LMS  = pagemtimes(T_indiv, SPD);
    XYZ  = pagemtimes(Mdisp, pagemtimes(invC, LMS));
    Wown = XYZ(:,idxWhite,:);
    % Reference hue anchor at full adaptation, and D=0 reference/observer (a',b').
    rab1 = getab(cam16_s3(XYZw, Wref, 1, M16, cam));          % D=1 placement (hue)
    refD0 = getab(cam16_s3(XYZw, Wref, 0, M16, cam));         % D=0 reference
    obsD0 = getab(cam16_s3(XYZ,  Wown, 0, M16, cam));         % D=0 observers
    for l = 1:numel(keepIdx)
        % Accumulate deviations and the mean D=1 reference position per surface.
        g = keepIdx(l);
        dev{g} = [dev{g}, squeeze(obsD0(:,l,:)) - refD0(:,l)];
        refSum(g,:) = refSum(g,:) + rab1(:,l)'; refN(g)=refN(g)+1;
    end
end
% Mean D=1 reference (a',b') per surface (its nominal hue location).
refPos = nan(nG,2); pres=find(refN>0); refPos(pres,:)=refSum(pres,:)./refN(pres);

% -------- per-surface ellipse area + hue --------
% area = ellipse area; apr/bpr = reference a'/b'; hue = reference hue angle (deg).
area=nan(nG,1); apr=nan(nG,1); bpr=nan(nG,1); hue=nan(nG,1);
for g=pres'
    if isempty(dev{g}), continue; end
    % Place the deviation cloud at the surface's reference position, fit ellipse.
    P = refPos(g,:) + dev{g}.';
    [xe,ye,ok]=ellipse_poly(P,0.95,120);
    if ~ok, continue; end
    % Metamerism magnitude = ellipse area; record reference chroma and hue.
    area(g)=polyarea(xe,ye);
    apr(g)=refPos(g,1); bpr(g)=refPos(g,2); hue(g)=mod(atan2d(bpr(g),apr(g)),360);
end
% Keep surfaces with a valid ellipse and reference b'.
val = find(isfinite(area) & isfinite(bpr));
A=area(val); ap=apr(val); bp=bpr(val); h=hue(val); N=numel(val);
fprintf('\n==== Fig 6 blue-hue test: D=0, SD=5 nm, pooled all illuminants ====\n');
fprintf('surfaces analysed: N = %d\n', N);

% (i) Spearman area vs b'
% Rank correlation of ellipse area with the reference b' axis (b'<0 = bluer).
[rho_b,p_b] = corr(A, bp, 'type','Spearman');
fprintf('(i) Spearman(area, b'''')      rho = %+.3f, p = %.2e  (blue = b''<0)\n', rho_b, p_b);
% area vs -b' (blueness) for sign clarity
% Same test against blueness (-b'); positive rho means bluer -> larger ellipse.
[rho_bl,p_bl] = corr(A, -bp, 'type','Spearman');
fprintf('    Spearman(area, blueness=-b'''') rho = %+.3f, p = %.2e\n', rho_bl, p_bl);

% (ii) circular-linear correlation between hue and log-area (Mardia)
% Mardia's circular-linear correlation between hue (circular) and log-area
% (linear): combine correlations of log-area with cos/sin of hue.
la = log(A);
c1 = corr(la, cosd(h)); s1 = corr(la, sind(h)); rcs = corr(cosd(h), sind(h));
Rcl = sqrt( (c1^2 + s1^2 - 2*c1*s1*rcs) / (1 - rcs^2) );
% F approximation of the circular-linear correlation for significance.
Fcl = (N-3)/2 * Rcl^2/(1-Rcl^2); p_cl = 1 - fcdf(Fcl, 2, N-3);
fprintf('(ii) circular-linear corr(hue, log-area) R = %.3f, p = %.2e\n', Rcl, p_cl);
% hue of the largest ellipses
% Where the biggest ellipses sit in hue (expected in the blue region).
[~,ord]=sort(A,'descend');
fprintf('     hue of 5 largest ellipses: %s deg\n', num2str(round(h(ord(1:5))')));
% Area-weighted resultant mean hue.
fprintf('     mean hue weighted by area: %.0f deg\n', mod(atan2d(sum(A.*sind(h)),sum(A.*cosd(h))),360));

% (iii) blue sector vs rest
% Split surfaces into a blue hue sector (200-300 deg) and the rest.
isBlue = h>=200 & h<=300;
fprintf('(iii) blue sector h in [200,300] deg: n=%d (rest n=%d)\n', sum(isBlue), sum(~isBlue));
% Median ellipse area in each group and their ratio.
med_b = median(A(isBlue)); med_r = median(A(~isBlue));
% Wilcoxon rank-sum (Mann-Whitney) test that blue areas exceed the rest.
[p_rs,~,st] = ranksum(A(isBlue), A(~isBlue));
fprintf('     median area  blue = %.1f,  rest = %.1f,  ratio = %.2f\n', med_b, med_r, med_b/med_r);
fprintf('     Wilcoxon rank-sum p = %.2e, z = %.2f\n', p_rs, st.zval);

% report area range
% Overall spread of ellipse areas.
fprintf('area range [%.1f, %.1f], median %.1f\n', min(A),max(A),median(A));
fprintf('=================================================================\n');

% Save per-surface arrays and test statistics for the figure/manuscript.
save(fullfile(projectRoot,'results','fig6_blue_hue_test.mat'), ...
    'A','ap','bp','h','val','rho_b','p_b','Rcl','p_cl','p_rs','med_b','med_r');

% ===================== local helpers (copied from plot_fig6) =====================
% getab: drop lightness J', keep the (a',b') chroma coordinates.
function ab = getab(s3), ab = s3(2:3,:,:); end
% cam16_s3: CIECAM16 forward transform to CAM16-UCS J'a'b'.
%   XYZ test tristimuli, W adopted white, D adaptation degree, M16 CAT16 matrix,
%   cam viewing conditions. Returns [J'; a'; b'] page-wise. Same math as the
%   other analysis scripts (see stats_ellipse_and_daylight.m for step comments).
function s3 = cam16_s3(XYZ, W, D, M16, cam)
% Viewing conditions; Yw = white luminance.
LA=cam.LA; Yb=cam.Yb; c=cam.c; Nc=cam.Nc; Yw=W(2,1,:);
% Luminance adaptation factor FL.
k=1/(5*LA+1); FL=0.2*k^4*(5*LA)+0.1*(1-k^4)^2*(5*LA)^(1/3);
% Background induction factors n, z, Nbb=Ncb.
n=max(Yb./Yw,1e-6); z=1.48+sqrt(n); Nbb=0.725.*(1./n).^0.2; Ncb=Nbb;
% CAT16 von Kries gain Dr and adapted/nonlinear white response.
RGBw=pagemtimes(M16,W); Dr=D.*Yw./RGBw+(1-D); RGBcw=Dr.*RGBw;
tw=(FL.*abs(RGBcw)/100).^0.42; RGBaw=sign(RGBcw).*400.*tw./(27.13+tw)+0.1;
% White achromatic response Aw.
Aw=(2*RGBaw(1,1,:)+RGBaw(2,1,:)+RGBaw(3,1,:)/20-0.305).*Nbb;
% Adapt + nonlinearity for the test stimuli.
RGB=pagemtimes(M16,XYZ); RGBc=Dr.*RGB;
t=(FL.*abs(RGBc)/100).^0.42; RGBa=sign(RGBc).*400.*t./(27.13+t)+0.1;
R1=RGBa(1,:,:); R2=RGBa(2,:,:); R3=RGBa(3,:,:);
% Opponent signals a,b and achromatic response Ac.
a=R1-12*R2/11+R3/11; b=(R1+R2-2*R3)/9; Ac=(2*R1+R2+R3/20-0.305).*Nbb;
% Lightness J.
J=100.*sign(Ac).*(abs(Ac)./Aw).^(c.*z);
% Hue hh and eccentricity et.
hh=atan2d(b,a); et=0.25*(cos(deg2rad(hh)+2)+3.8);
% Chroma C and colourfulness M.
denC=R1+R2+(21/20)*R3; tC=(50000/13*Nc.*Ncb).*et.*hypot(a,b)./max(denC,1e-9);
Cc=(max(tC,0)).^0.9.*sqrt(max(J,0)/100).*(1.64-0.29.^n).^0.73; Mc=Cc.*FL^0.25;
% UCS compression to J', M'; a'=M'cos h, b'=M'sin h.
Jp=(1+100*0.007).*J./(1+0.007.*J); Mp=(1/0.0228)*log(1+0.0228.*Mc);
s3=cat(1,Jp,Mp.*cosd(hh),Mp.*sind(hh));
end
% ellipse_poly: conf-level covariance ellipse of 2-D points P as an nPts polygon.
function [xe,ye,ok]=ellipse_poly(P,conf,nPts)
% Drop non-finite rows; need >=5 points for a covariance estimate.
ok=true; P=P(all(isfinite(P),2),:);
if size(P,1)<5, xe=nan; ye=nan; ok=false; return; end
% Mean, covariance, chi-square radius^2 for the confidence level.
mu=mean(P,1); C=cov(P); q=-2*log(max(1-conf,eps));
% Eigen-axes and semi-axis lengths.
[V,Dg]=eig(C); r=sqrt(q)*sqrt(max(diag(Dg),0));
% Trace unit circle -> ellipse, shifted to the mean.
t=linspace(0,2*pi,nPts); xy=(V*diag(r)*[cos(t);sin(t)])+mu(:);
xe=xy(1,:); ye=xy(2,:);
end
