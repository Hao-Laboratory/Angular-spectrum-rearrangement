function Eout = ScalarDiffraction_CZT_AP(Beam,Scope)
%SCALARDIFFRACTION_CZT_AP calculate scalar diffraction by CZT based angular
%spectrum method in Cartesian coordinate system.
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
fzz(lambda^2.*(fxx.^2+fyy.^2)>1) = 0;

%% calculation
F0 = myczt(E0, LRfy, yd, fy);  % FT
F0 = F0.';
F0 = myczt(F0, LRfx, xd, fx);
F0 = F0.';

fuu = cos(theta2)*cos(phi2)*fxx+cos(theta2)*sin(phi2)*fyy-sin(theta2)*fzz;
fvv = -sin(phi2)*fxx+cos(phi2)*fyy;
fww = sin(theta2)*cos(phi2)*fxx+sin(theta2)*sin(phi2)*fyy+cos(theta2)*fzz;

minFu = min(fuu(:)); maxFu = max(fuu(:));
minFv = min(fvv(:)); maxFv = max(fvv(:));
fu = linspace(minFu,maxFu,LRfx);
fv = linspace(minFv,maxFv,LRfy);
[fuu2,fvv2] = meshgrid(fu,fv);
fww2_p = sqrt(abs(1/lambda^2-(fuu2.^2+fvv2.^2)));
fww2_n = -fww2_p;

% preparation for interpolation
fxx2_p = cos(theta2)*cos(phi2)*fuu2-sin(phi2)*fvv2+sin(theta2)*cos(phi2)*fww2_p;
fyy2_p = cos(theta2)*sin(phi2)*fuu2+cos(phi2)*fvv2+sin(theta2)*sin(phi2)*fww2_p;
fzz2_p = -sin(theta2)*fuu2+cos(theta2)*fww2_p;
fxx2_n = cos(theta2)*cos(phi2)*fuu2-sin(phi2)*fvv2+sin(theta2)*cos(phi2)*fww2_n;
fyy2_n = cos(theta2)*sin(phi2)*fuu2+cos(phi2)*fvv2+sin(theta2)*sin(phi2)*fww2_n;
fzz2_n = -sin(theta2)*fuu2+cos(theta2)*fww2_n;

% avoid numerical error
if any((fxx2_p(fxx2_p > Lfx/2) - Lfx/2) < 1e-15)
    fxx2_p((fxx2_p > Lfx/2)&(fxx2_p-Lfx/2<1e-15)) = Lfx/2;
end

if any((fxx2_p(fxx2_p < -Lfx/2) + Lfx/2) > -1e-15)
    fxx2_p((fxx2_p < -Lfx/2)&(fxx2_p+Lfx/2>-1e-15)) = -Lfx/2;
end

if any((fyy2_p(fyy2_p > Lfy/2) - Lfy/2) < 1e-15)
    fyy2_p((fyy2_p > Lfy/2)&(fyy2_p-Lfy/2<1e-15)) = Lfy/2;
end

if any((fyy2_p(fyy2_p < -Lfy/2) + Lfy/2) > -1e-15)
    fyy2_p((fyy2_p < -Lfy/2)&(fyy2_p+Lfy/2>-1e-15)) = -Lfy/2;
end

if any((fxx2_n(fxx2_n > Lfx/2) - Lfx/2) < 1e-15)
    fxx2_n((fxx2_n > Lfx/2)&(fxx2_n-Lfx/2<1e-15)) = Lfx/2;
end

if any((fxx2_n(fxx2_n < -Lfx/2) + Lfx/2) > -1e-15)
    fxx2_n((fxx2_n < -Lfx/2)&(fxx2_n + Lfx/2 > -1e-15)) = -Lfx/2;
end

if any((fyy2_n(fyy2_n > Lfy/2) - Lfy/2) < 1e-15)
    fyy2_n((fyy2_n > Lfy/2)&(fyy2_n-Lfy/2<1e-15)) = Lfy/2;
end

if any((fyy2_n(fyy2_n < -Lfy/2) + Lfy/2) > -1e-15)
    fyy2_n((fyy2_n < -Lfy/2)&(fyy2_n+Lfy/2>-1e-15)) = -Lfy/2;
end

mask1 = ((fxx2_p.^2+fyy2_p.^2)>1/lambda^2); % Exclude non-existent coordinate points
mask2 = ((fxx2_n.^2+fyy2_n.^2)>1/lambda^2);
mask3 = (lambda^2*(fuu2.^2+fvv2.^2)>1);
mask4 = (fww2_p == 0);
mask5 = (fzz2_p<0|fzz2_p>1/lambda);
mask6 = (fzz2_n<0|fzz2_n>1/lambda);
mask7 = ((fxx2_p<(-Lfx/2))|(fxx2_p>(Lfx/2))|(fyy2_p<(-Lfy/2))|(fyy2_p>(Lfy/2)));
mask8 = ((fxx2_n<(-Lfx/2))|(fxx2_n>(Lfx/2))|(fyy2_n<(-Lfy/2))|(fyy2_n>(Lfy/2)));

% Jacobian determinant
J_p = -sin(theta2)*fuu2./fww2_p + cos(theta2);
J_n = -sin(theta2)*fuu2./fww2_n + cos(theta2);

H = exp(1i*2*pi*zs*fzz);

Fd = F0.*H.*exp(1i*2*pi*Scope.ws*fww);

Fw_p = interp2(fxx,fyy,Fd,fxx2_p,fyy2_p,'linear',0);
Fw_n = interp2(fxx,fyy,Fd,fxx2_n,fyy2_n,'linear',0);

Fw_p = Fw_p.*abs(J_p);
Fw_n = Fw_n.*abs(J_n);

Fw_p(~isfinite(Fw_p)) = 0;
Fw_n(~isfinite(Fw_n)) = 0;

Fw_p(mask1 | mask3 | mask4 | mask5 | mask7) = 0;
Fw_n(mask2 | mask3 | mask6 | mask8) = 0;

Fw = Fw_p + Fw_n;

Eout = myiczt(Fw, lv, fv, Scope.vs);
Eout = Eout.';
Eout = myiczt(Eout, lu, fu, Scope.us);
Eout = Eout.';

disp('ScalarDiffraction_CZT_AP has completed!');

end