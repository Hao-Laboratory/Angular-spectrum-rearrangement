function Eout = ScalarDiffraction_GT_AP(Beam,Scope)
%SCALARDIFFRACTION_GT_AP calculate scalar diffraction by direct integration
%based angular spectrum method in Cartesian coordinate system.
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
lu = length(Scope.us); lv = length(Scope.vs); 

global theta2 phi2
[uus,vvs] = meshgrid(Scope.us,Scope.vs);
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
Eout = zeros(lv,lu);

totalNum = lv*lu;

F0 = mdft(E0, xd, yd, fx, fy);  % FT

for ii = 1:totalNum
    
    % TF
    H = exp(1i*2*pi*(zs+zzs(ii))*fzz);

    Fd = F0.*H;

    Eout(ii) = midft(Fd, xxs(ii), yys(ii), fx, fy);  % IFT
    
    textwaitbar(ii, totalNum, 'ScalarDiffraction_GT_AP in progress');
end

end