% ========================================================================
% adapting_white_variants.m
%
% Robustness of the metamerism reduction to the DEFINITION of the adapting
% white. Throughout the paper the adapting white is each observer's own
% reproduction of the highest-luminance (display-white) surface. Here the
% same CAM16-UCS pipeline is rerun at full adaptation (D = 1) with two other
% plausible scene-statistics anchors, each computed per observer from that
% observer's own reproduced stimuli (and, for the reference, from the
% standard-observer world spectra):
%
%   (0) brightest   : highest-luminance surface            (paper default)
%   (1) meanLMS     : mean cone excitation over the surface set ("grey
%                     world"). Because XYZ is linear in LMS this is the mean
%                     XYZ over surfaces, i.e. a LUMINANCE-WEIGHTED mean
%                     chromaticity (bright surfaces dominate).
%   (2) meanChrom   : mean CIE 1931 xy chromaticity over the surface set
%                     (every surface weighted equally) placed at the mean
%                     luminance. This differs from (1) because the mean of
%                     chromaticities is not the chromaticity of the mean.
%   (2b) meanChromUV: as (2) but averaging u'v' instead of xy (diagnostic
%                     only; shows the sensitivity to the chromaticity space).
%
% The unadapted baseline d0 (D = 0, brightest white, exactly as in the main
% analysis) is shared by all variants, so the reported reduction
%     100 * (1 - d1 / d0)   per (observer, illuminant)
% is directly comparable with the white-anchored D = 1 reduction in the
% Results and with the numerically optimised white (optimal_white_upperlimit).
% All 10,000 observers and all 36 illuminants are used (no subsampling).
%
% The heavy N = 10,000 run is NOT repeated: per-observer XYZ is rebuilt from
% the stored cone-fundamental family + display + world spectra.
%
% OUTPUT
%   results/adapting_white_variants.mat : per-variant pooled reductions and
%   D = 1 medians, per bandwidth and per illuminant class, plus a summary.
% ========================================================================
clearvars; clc;
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(projectRoot,'src')));
dataDir = fullfile(projectRoot,'data');
inMat   = fullfile(projectRoot,'results','results_AsanoCIE2006_metamerism_CALC_multiIlluminant.mat');
outMat  = fullfile(projectRoot,'results','adapting_white_variants.mat');

% CIECAM16 average-surround viewing conditions and CAT16 matrix (as main run).
cam = struct('LA',20,'Yb',20,'F',1.0,'c',0.69,'Nc',1.0);
M16 = [0.401288 0.650173 -0.051461; -0.250268 1.204414 0.045854; -0.002079 0.048952 0.953127];
variants = ["brightest","meanLMS","meanChrom","meanChromUV"];
nVar = numel(variants);

% ---- load the light parts of the main results (not byIll) -------------
S = load(inMat,'wls','params','display','cone_fund','illuminants','keepIdx_global');
wls = S.wls(:); sdList = S.params.sdList(:)'; nSD = numel(sdList);
Pset = S.display.Pset; Mdisp_set = S.display.Mdisp_set;
Ill = S.illuminants; nIll = numel(Ill);
keepIdx = S.keepIdx_global(:);
T_indiv = permute(S.cone_fund,[3 2 1]);            % 3 x nW x nObs
nObs = size(T_indiv,3);
clear S

% Reference CIE 1931 CMFs and the CIE 99 CFI reflectances on the working grid.
Axyz  = readmatrix(fullfile(dataDir,'CIE1931XYZ.txt'));
T_xyz = interp1(Axyz(:,1),Axyz(:,2:4),wls,'pchip',0)';
Rmat  = readmatrix(fullfile(dataDir,'CIE_srf_cfi_1nm.csv'));
R_all = min(max(interp1(Rmat(:,1),Rmat(:,2:end),wls,'pchip','extrap'),0),1);
R_keep = R_all(:,keepIdx); nSamples = numel(keepIdx);
catVec = string({Ill.category});
fprintf('observers = %d, illuminants = %d, surfaces = %d\n', nObs, nIll, nSamples);

% ---- per-bandwidth inverse display response per observer --------------
invC_set = cell(nSD,1);
for si = 1:nSD
    C = pagemtimes(T_indiv, Pset{si});
    invC = zeros(3,3,nObs);
    for o = 1:nObs, invC(:,:,o) = inv(C(:,:,o)); end
    invC_set{si} = invC;
end

% ---- main loop ---------------------------------------------------------
% d0  : nSD x nObs x nIll   unadapted (D=0, brightest white)
% d1  : nSD x nObs x nIll x nVar   fully adapted (D=1) per white variant
d0 = nan(nSD,nObs,nIll);
d1 = nan(nSD,nObs,nIll,nVar);
Wref_all = nan(3,nVar,nIll);
tStart = tic;
for k = 1:nIll
    SPD  = Ill(k).spd(:) .* R_keep;                 % nW x nSamples
    XYZw = T_xyz * SPD;                             % 3 x nSamples (reference)
    [~, idxWhite] = max(XYZw(2,:));
    Wref = white_variants(XYZw, idxWhite);          % 3 x 1 x nVar
    Wref_all(:,:,k) = squeeze(Wref);
    refD0 = cam16_s3(XYZw, Wref(:,1,1), 0, M16, cam);
    refD1 = cell(nVar,1);
    for v = 1:nVar, refD1{v} = cam16_s3(XYZw, Wref(:,1,v), 1, M16, cam); end
    LMS = pagemtimes(T_indiv, SPD);                 % 3 x nSamples x nObs
    for si = 1:nSD
        XYZ = pagemtimes(Mdisp_set{si}, pagemtimes(invC_set{si}, LMS));  % 3 x nS x nObs
        Wobs = white_variants(XYZ, idxWhite);       % 3 x nObs x nVar
        s0 = cam16_s3(XYZ, reshape(Wobs(:,:,1),3,1,nObs), 0, M16, cam);
        d0(si,:,k) = median(squeeze(vecnorm(s0 - refD0,2,1)),1,'omitnan');
        for v = 1:nVar
            s1 = cam16_s3(XYZ, reshape(Wobs(:,:,v),3,1,nObs), 1, M16, cam);
            d1(si,:,k,v) = median(squeeze(vecnorm(s1 - refD1{v},2,1)),1,'omitnan');
        end
    end
    fprintf('  [%2d/%2d] %-8s  sigma5 D1 median dE'' : %s   (%.0f s)\n', k, nIll, Ill(k).name, ...
        mat2str(round(squeeze(median(d1(1,:,k,:),2,'omitnan'))',2)), toc(tStart));
end

% ---- pooling: per class and over all classes ----------------------------
classes = ["daylight","smooth","peaky"];
groups  = [classes, "all"];
red_med = nan(nSD,nVar,numel(groups));   % median % reduction vs D=0
d1_med  = nan(nSD,nVar,numel(groups));   % median D=1 dE'
d0_med  = nan(nSD,numel(groups));
for g = 1:numel(groups)
    if groups(g) == "all", sel = true(1,nIll); else, sel = catVec == groups(g); end
    for si = 1:nSD
        a0 = reshape(d0(si,:,sel),[],1);
        d0_med(si,g) = median(a0,'omitnan');
        for v = 1:nVar
            a1 = reshape(d1(si,:,sel,v),[],1);
            red_med(si,v,g) = median(100*(1 - a1./a0),'omitnan');
            d1_med(si,v,g)  = median(a1,'omitnan');
        end
    end
end

fprintf('\n=== Median metamerism reduction at D = 1 (%% of unadapted), by adapting-white definition ===\n');
for g = 1:numel(groups)
    fprintf('\n  [%s]  D0 median dE'' (sigma 5/10/20): %.2f / %.2f / %.2f\n', groups(g), d0_med(:,g));
    fprintf('  %-12s | %s\n', 'white', strjoin(compose('sigma=%2g', sdList), ' | '));
    for v = 1:nVar
        fprintf('  %-12s | %s   (D1 dE'': %s)\n', variants(v), ...
            strjoin(compose('%5.1f %%', red_med(:,v,g)), '  | '), ...
            strjoin(compose('%.2f', d1_med(:,v,g)), ' / '));
    end
end

% how different are the whites themselves? (reference, per class, xy)
fprintf('\n=== Reference-white chromaticity (xy) per variant, averaged within class ===\n');
for c = 1:numel(classes)
    sel = find(catVec == classes(c));
    fprintf('  %-9s', classes(c));
    for v = 1:nVar
        W = squeeze(Wref_all(:,v,sel)); xy = W(1:2,:)./sum(W,1);
        fprintf('  %s x=%.4f y=%.4f Y=%.1f |', variants(v), mean(xy(1,:)), mean(xy(2,:)), mean(W(2,:)));
    end
    fprintf('\n');
end

save(outMat,'variants','classes','groups','sdList','red_med','d1_med','d0_med','d0','d1','Wref_all','-v7.3');
fprintf('\nSaved %s\n', outMat);

% ======================= local functions ================================
function W = white_variants(XYZ, idxWhite)
% XYZ : 3 x nSamples x nPages. Returns W : 3 x nPages x nVar
%   1 brightest surface, 2 mean XYZ (= mean LMS), 3 mean xy at mean Y,
%   4 mean u'v' at mean Y.
nP = size(XYZ,3);
W = zeros(3,nP,4);
W(:,:,1) = reshape(XYZ(:,idxWhite,:),3,nP);
W(:,:,2) = reshape(mean(XYZ,2,'omitnan'),3,nP);
Ym = reshape(mean(XYZ(2,:,:),2,'omitnan'),1,nP);
% mean xy
sXYZ = sum(XYZ,1);
x = XYZ(1,:,:)./sXYZ; y = XYZ(2,:,:)./sXYZ;
xm = reshape(mean(x,2,'omitnan'),1,nP); ym = reshape(mean(y,2,'omitnan'),1,nP);
W(:,:,3) = [xm./ym.*Ym; Ym; (1-xm-ym)./ym.*Ym];
% mean u'v'
den = XYZ(1,:,:) + 15*XYZ(2,:,:) + 3*XYZ(3,:,:);
u = 4*XYZ(1,:,:)./den; v = 9*XYZ(2,:,:)./den;
um = reshape(mean(u,2,'omitnan'),1,nP); vm = reshape(mean(v,2,'omitnan'),1,nP);
dd = 6*um - 16*vm + 12;
xu = 9*um./dd; yu = 4*vm./dd;
W(:,:,4) = [xu./yu.*Ym; Ym; (1-xu-yu)./yu.*Ym];
W = max(W, 1e-6);
end

function s3 = cam16_s3(XYZ, W, D, M16, cam)
% CIECAM16 forward transform to CAM16-UCS J'a'b' (identical to the main
% pipeline; see stats_ellipse_and_daylight.m for step-by-step comments).
LA=cam.LA;Yb=cam.Yb;c=cam.c;Nc=cam.Nc;Yw=W(2,1,:);
k=1/(5*LA+1);FL=0.2*k^4*(5*LA)+0.1*(1-k^4)^2*(5*LA)^(1/3);
n=max(Yb./Yw,1e-6);z=1.48+sqrt(n);Nbb=0.725.*(1./n).^0.2;Ncb=Nbb;
RGBw=pagemtimes(M16,W);Dr=D.*Yw./RGBw+(1-D);RGBcw=Dr.*RGBw;
tw=(FL.*abs(RGBcw)/100).^0.42;RGBaw=sign(RGBcw).*400.*tw./(27.13+tw)+0.1;
Aw=(2*RGBaw(1,1,:)+RGBaw(2,1,:)+RGBaw(3,1,:)/20-0.305).*Nbb;
RGB=pagemtimes(M16,XYZ);RGBc=Dr.*RGB;
t=(FL.*abs(RGBc)/100).^0.42;RGBa=sign(RGBc).*400.*t./(27.13+t)+0.1;
R1=RGBa(1,:,:);R2=RGBa(2,:,:);R3=RGBa(3,:,:);
a=R1-12*R2/11+R3/11;b=(R1+R2-2*R3)/9;Ac=(2*R1+R2+R3/20-0.305).*Nbb;
J=100.*sign(Ac).*(abs(Ac)./Aw).^(c.*z);
hh=atan2d(b,a);et=0.25*(cos(deg2rad(hh)+2)+3.8);
denC=R1+R2+(21/20)*R3;tC=(50000/13*Nc.*Ncb).*et.*hypot(a,b)./max(denC,1e-9);
Cc=(max(tC,0)).^0.9.*sqrt(max(J,0)/100).*(1.64-0.29.^n).^0.73;Mc=Cc.*FL^0.25;
Jp=(1+100*0.007).*J./(1+0.007.*J);Mp=(1/0.0228)*log(1+0.0228.*Mc);
s3=cat(1,Jp,Mp.*cosd(hh),Mp.*sind(hh));
end
