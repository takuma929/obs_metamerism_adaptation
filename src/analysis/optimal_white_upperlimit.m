% ========================================================================
% optimal_white_upperlimit.m
%
% Upper limit of metamerism reduction from the choice of adapting white.
% Throughout the paper the adapting white is the reproduced highest-luminance
% (display-white) surface. Here we ask how much better a NUMERICALLY OPTIMISED
% adapting white could do: for every observer we search the adapting white W
% that minimises that observer's median CAM16-UCS metamerism to the standard
% observer at full adaptation (D = 1), then compare the resulting reduction to
% the white-anchored D = 1 reduction. The reference (standard observer) is left
% at D = 1 anchored to its own white, exactly as in the main analysis.
%
% Reports, POOLED over all illuminant classes, the median metamerism reduction
% at D = 1 for each primary bandwidth (sigma = 5, 10, 20 nm):
%     - white-anchored  (highest-luminance surface)
%     - optimised white (numerical upper limit)
%
% Reduction per observer-illuminant = 100 * (1 - d / d0), d0 = unadapted (D=0).
% Subsamples observers for the search (nSub); the median is insensitive to it.
% ========================================================================
clearvars; clc;
projectRoot = '/Users/takuma/Documents/GitHub/observer_metamerism_adaptation';
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
inMat = fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
% CIECAM16 average-surround viewing conditions and CAT16 matrix.
cam = struct('LA',20,'Yb',20,'F',1.0,'c',0.69,'Nc',1.0);
M16 = [0.401288 0.650173 -0.051461; -0.250268 1.204414 0.045854; -0.002079 0.048952 0.953127];
% Primary bandwidths to evaluate (nm) and observer subsample size.
SDsel = [5 10 20];
nSub  = 1500;                          % observers subsampled for the search
% Nelder-Mead options for the per-observer white search (capped iterations).
optset = optimset('Display','off','MaxFunEvals',80,'MaxIter',80,'TolX',1e-3,'TolFun',1e-3);

% Load simulation inputs/results.
Sres = load(inMat,'wls','params','display','cone_fund','illuminants','byIll');
wls=Sres.wls(:); sdList=Sres.params.sdList(:)';
% Per-bandwidth display primaries and display->XYZ matrices.
Pset=Sres.display.Pset; Mdisp_set=Sres.display.Mdisp_set;
Ill=Sres.illuminants; nIll=numel(Ill);
% Evenly subsampled observers reshaped to 3 x nW x nObs cone fundamentals.
obsIdx = round(linspace(1,size(Sres.cone_fund,1),nSub));
T_indiv = permute(Sres.cone_fund(obsIdx,:,:),[3 2 1]); nObs=size(T_indiv,3);
% Per-illuminant bookkeeping (converged?, kept surfaces, white index/XYZ).
B=struct('ok',{},'keepIdx',{},'idxWhite',{},'XYZ_white',{});
for k=1:nIll
    r=Sres.byIll(k); B(k).ok=isfield(r,'ok')&&~isempty(r.ok)&&r.ok;
    if B(k).ok, B(k).keepIdx=r.keepIdx(:); B(k).idxWhite=r.idxWhite; B(k).XYZ_white=r.XYZ_white(:); end
end
% Reference CIE 1931 CMFs and CIE 99 CFI reflectances on the working grid.
Axyz=readmatrix(fullfile(dataDir,'CIE1931XYZ.txt'));
T_xyz=interp1(Axyz(:,1),Axyz(:,2:4),wls,'pchip',0)';
Rmat=readmatrix(fullfile(dataDir,'CIE_srf_cfi_1nm.csv'));
R_all=min(max(interp1(Rmat(:,1),Rmat(:,2:end),wls,'pchip','extrap'),0),1);
% Converged illuminants.
okAll=find(arrayfun(@(k)B(k).ok,1:nIll));
fprintf('observers (subsampled) = %d, illuminants = %d\n', nObs, numel(okAll));

nSD=numel(SDsel);
redD1  = cell(nSD,1);   % pooled per-(obs,illum) reduction, white-anchored
redOpt = cell(nSD,1);   % pooled per-(obs,illum) reduction, optimised white
% Loop over primary bandwidths.
for j=1:nSD
    % Display model and per-observer display response/inverse for this bandwidth.
    si=find(sdList==SDsel(j),1); Mdisp=Mdisp_set{si};
    C=pagemtimes(T_indiv,Pset{si}); invC=zeros(3,3,nObs);
    for o=1:nObs, invC(:,:,o)=inv(C(:,:,o)); end
    % Accumulators for the pooled per-(obs,illum) reductions.
    accD1=[]; accOpt=[];
    for k=okAll
        keepIdx=B(k).keepIdx; idxWhite=B(k).idxWhite; Wref=B(k).XYZ_white;
        % Colour signals and reference CAM16-UCS at D=0 and D=1 (fixed ref white).
        SPD=Ill(k).spd(:).*R_all(:,keepIdx); XYZw=T_xyz*SPD;
        refD0=cam16_s3(XYZw,Wref,0,M16,cam);
        refD1=cam16_s3(XYZw,Wref,1,M16,cam);
        % Per-observer metameric display matches, and each observer's display white.
        LMS=pagemtimes(T_indiv,SPD); XYZ=pagemtimes(Mdisp,pagemtimes(invC,LMS));
        Wown=squeeze(XYZ(:,idxWhite,:));           % 3 x nObs
        % d0/d1 = median CAM16-UCS distance to ref at D=0/D=1 (white-anchored);
        % dopt = same at D=1 but with the numerically optimised adapting white.
        d0=zeros(nObs,1); d1=zeros(nObs,1); dopt=zeros(nObs,1);
        parfor o=1:nObs
            % This observer's matched display XYZ and anchored white.
            Xo=XYZ(:,:,o); Wo=Wown(:,o);
            % Unadapted and fully-adapted metamerism (median over surfaces).
            s0=cam16_s3(Xo,Wo,0,M16,cam); d0(o)=median(vecnorm(s0-refD0),'omitnan');
            s1=cam16_s3(Xo,Wo,1,M16,cam); d1(o)=median(vecnorm(s1-refD1),'omitnan');
            % Search (in log space, so white stays positive) the D=1 adapting
            % white that minimises this observer's median metamerism.
            fun=@(lw) median(vecnorm(cam16_s3(Xo,exp(lw),1,M16,cam)-refD1),'omitnan');
            lwopt=fminsearch(fun,log(max(Wo,1e-6)),optset);
            dopt(o)=fun(lwopt);
        end
        % Per-observer % reduction vs D=0: white-anchored and optimised-white.
        accD1 =[accD1;  100*(1-d1./d0)];
        accOpt=[accOpt; 100*(1-dopt./d0)];
    end
    redD1{j}=accD1; redOpt{j}=accOpt;
end

% --- Summary table: median reduction, white-anchored vs optimised, per bandwidth ---
fprintf('\n=== Median metamerism reduction at D = 1 (pooled over all illuminant classes) ===\n');
fprintf('  SD(nm) | white-anchored | optimised white | gain\n');
for j=1:nSD
    % a = median white-anchored reduction; b = median optimised-white reduction.
    a=median(redD1{j},'omitnan'); b=median(redOpt{j},'omitnan');
    % gain (b-a) = extra reduction the best white leaves on the table.
    fprintf('   %3g   |   %6.1f %%     |   %6.1f %%      | +%.1f pts\n', SDsel(j), a, b, b-a);
end
% Compact restatement of the two rows (per SD 5/10/20).
fprintf('optimised-white values (D=1, per SD 5/10/20): %.1f / %.1f / %.1f %%\n', ...
    median(redOpt{1},'omitnan'), median(redOpt{2},'omitnan'), median(redOpt{3},'omitnan'));
fprintf('white-anchored values (D=1, per SD 5/10/20):  %.1f / %.1f / %.1f %%\n', ...
    median(redD1{1},'omitnan'), median(redD1{2},'omitnan'), median(redD1{3},'omitnan'));

% Save pooled reduction distributions and settings for the manuscript.
save(fullfile(projectRoot,'results','optimal_white_upperlimit.mat'),'redD1','redOpt','SDsel','nSub');

% ---- local CAM16-UCS forward (identical to the main pipeline) ----
% cam16_s3: CIECAM16 forward transform to CAM16-UCS J'a'b' (see the step-by-step
% comments in stats_ellipse_and_daylight.m for the same computation).
function s3=cam16_s3(XYZ,W,D,M16,cam)
% Viewing conditions; Yw = white luminance.
LA=cam.LA;Yb=cam.Yb;c=cam.c;Nc=cam.Nc;Yw=W(2,1,:);
% Luminance adaptation factor FL.
k=1/(5*LA+1);FL=0.2*k^4*(5*LA)+0.1*(1-k^4)^2*(5*LA)^(1/3);
% Background induction factors n, z, Nbb=Ncb.
n=max(Yb./Yw,1e-6);z=1.48+sqrt(n);Nbb=0.725.*(1./n).^0.2;Ncb=Nbb;
% CAT16 von Kries gain Dr and adapted/nonlinear white response.
RGBw=pagemtimes(M16,W);Dr=D.*Yw./RGBw+(1-D);RGBcw=Dr.*RGBw;
tw=(FL.*abs(RGBcw)/100).^0.42;RGBaw=sign(RGBcw).*400.*tw./(27.13+tw)+0.1;
% White achromatic response Aw.
Aw=(2*RGBaw(1,1,:)+RGBaw(2,1,:)+RGBaw(3,1,:)/20-0.305).*Nbb;
% Adapt + nonlinearity for the test stimuli.
RGB=pagemtimes(M16,XYZ);RGBc=Dr.*RGB;
t=(FL.*abs(RGBc)/100).^0.42;RGBa=sign(RGBc).*400.*t./(27.13+t)+0.1;
R1=RGBa(1,:,:);R2=RGBa(2,:,:);R3=RGBa(3,:,:);
% Opponent signals a,b and achromatic response Ac.
a=R1-12*R2/11+R3/11;b=(R1+R2-2*R3)/9;Ac=(2*R1+R2+R3/20-0.305).*Nbb;
% Lightness J.
J=100.*sign(Ac).*(abs(Ac)./Aw).^(c.*z);
% Hue hh and eccentricity et.
hh=atan2d(b,a);et=0.25*(cos(deg2rad(hh)+2)+3.8);
% Chroma C and colourfulness M.
denC=R1+R2+(21/20)*R3;tC=(50000/13*Nc.*Ncb).*et.*hypot(a,b)./max(denC,1e-9);
Cc=(max(tC,0)).^0.9.*sqrt(max(J,0)/100).*(1.64-0.29.^n).^0.73;Mc=Cc.*FL^0.25;
% UCS compression to J', M'; a'=M'cos h, b'=M'sin h.
Jp=(1+100*0.007).*J./(1+0.007.*J);Mp=(1/0.0228)*log(1+0.0228.*Mc);
s3=cat(1,Jp,Mp.*cosd(hh),Mp.*sind(hh));
end
