function Eout = ScalarDiffraction_ASR_AP(Beam,Scope,number_u,number_v)
%SCALARDIFFRACTION_ASR_AP calculate scalar diffraction by angular spectrum
%rearrangement (ASR) in Cartesian coordinate system.
%
% INPUT********************************************************************
% Beam.wavelength: scalar value, wavelength of light
% Beam.amp: LRy*LRx matrix, amplitude distribution on diffraction aperture
% Beam.phs: LRy*LRx matrix, phase distribution on diffraction aperture
% Beam.PixelSize: scalar value, pixel size of diffraction aperture after
% discretization
% Scope.us: 1*lu array, representing u axis of observation plane
% Scope.vs: 1*lv array, representing v axis of observation plane
% Scope.ws: 1*lw array, representing w axis of observation volume
% Scope.zs: scalar value, distance from source plane to observation plane
% number_u: number of target samples in the frequency domain after merging
% along the u-axis
% number_v: number of target samples in the frequency domain after merging
% along the v-axis
%
% OUTPUT*******************************************************************
% Eout: diffraction field on observation plane
%
% *************************************************************************
% LIU Xin
% liuxin24@hku.hk
% Apr.23, 2021
%
% updated by HU Yiwen
% huyw@zju.edu.cn

%% data initialization
n = 1;  % refractive index of medium

lambda = Beam.wavelength/n;  % wavelength in medium

%% diffraction aperture
% resolution of diffraction aperture
[LRy, LRx] = size(Beam.amp);

zs = Scope.zs;

E0 = Beam.amp.*exp(1i*Beam.phs);

% real size of diffraction aperture
LSx = (LRx-1)*Beam.PixelSize;
LSy = (LRy-1)*Beam.PixelSize;

% real coordinates of diffraction aperture
xd = linspace(-LSx/2,LSx/2,LRx);
yd = linspace(-LSy/2,LSy/2,LRy);

%% observation plane
[uus,vvs] = meshgrid(Scope.us,Scope.vs);

global theta2 phi2
xxs = cos(theta2)*cos(phi2)*uus - sin(phi2)*vvs + sin(theta2)*cos(phi2)*Scope.ws;
yys = cos(theta2)*sin(phi2)*uus + cos(phi2)*vvs + sin(theta2)*sin(phi2)*Scope.ws;
zzs = -sin(theta2)*uus + cos(theta2)*Scope.ws;

Lx = max(xxs(:)) - min(xxs(:));
Ly = max(yys(:)) - min(yys(:));

minZ = zs + min(zzs(:));

% bandwidth of the aperture plane
Lfx = 1./Beam.PixelSize;
Lfy = 1./Beam.PixelSize;

%% identify sampling interval in frequency domain
fmax_fft = 1/(2*Beam.PixelSize);

% maximum sampling interval limited by TF
dfMax1 = sqrt(1-(lambda*fmax_fft)^2)/(lambda*minZ*fmax_fft);

% maximum sampling interval limited by observation plane
dfxMax2 = 1/Lx;
dfyMax2 = 1/Ly;

% minimum requirements of sampling interval
dfx = min(dfxMax2,dfMax1);
dfy = min(dfyMax2,dfMax1);

s = 2;
LRfx = max(ceil(Lfx/dfx*s), LRx);
LRfy = max(ceil(Lfy/dfy*s), LRy);

% spatial frequency coordinate
fx = reshape(linspace(-Lfx/2,Lfx/2,LRfx),1,LRfx);
fy = reshape(linspace(-Lfy/2,Lfy/2,LRfy),LRfy,1);

% spatial frequency grid
[fxx,fyy] = meshgrid(fx,fy);

fzz = sqrt(1/lambda^2-(fxx.^2+fyy.^2));
fzz(lambda^2.*(fxx.^2+fyy.^2)>1) = nan;

%% calculation
F0 = mdft(E0, xd, yd, fx, fy);  % FT

% coordinate projection
fuu = cos(theta2)*cos(phi2)*fxx+cos(theta2)*sin(phi2)*fyy-sin(theta2)*fzz;
fvv = -sin(phi2)*fxx+cos(phi2)*fyy;
fww = sin(theta2)*cos(phi2)*fxx+sin(theta2)*sin(phi2)*fyy+cos(theta2)*fzz;

% merging
spa = lambda^2.*(fxx.^2+fyy.^2)<=1;
[fu_temp, idx_fu] = unique_tol(fuu(spa),number_u);
[fv_temp, idx_fv] = unique_tol(fvv(spa),number_v);

H = exp(1i*2*pi*zs*fzz);  % TF

Fd = F0.*H.*exp(1i*2*pi*Scope.ws*fww);

if (theta2 == 0)&&(phi2 == 0)
    fu_eff = fx;
    fv_eff = fy;
    Fw = Fd;
else
    Fw = sparse(idx_fv,idx_fu,Fd(spa));  % Place coordinates in xy at the counterparts in the tilted plane

    % weight calculation
    absFw = sparse(idx_fv,idx_fu,abs(Fd(spa)));

    fuFw = abs(Fd).*fuu;
    fvFw = abs(Fd).*fvv;

    fu_eff = full(sparse(ones(size(idx_fu)), idx_fu, fuFw))./full(sum(absFw));
    fv_eff = full(sparse(idx_fv, ones(size(idx_fv)), fvFw))./full(sum(absFw,2));
    if find(isnan(fu_eff))
        idx = isnan(fu_eff);
        fu_eff(idx) = fu_temp(idx);
    end
    if find(isnan(fv_eff))
        idx = isnan(fv_eff);
        fv_eff(idx) = fv_temp(idx);
    end
end

Eout = midft(Fw, Scope.us, Scope.vs, fu_eff, fv_eff);

disp('ScalarDiffraction_ASR_AP has completed!');

end