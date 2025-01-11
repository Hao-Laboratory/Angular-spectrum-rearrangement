function [Ex,Ey,Ez] = VectorialDiffraction_CZT_AP(Obj,Beam,Scope)
%VECTORIALDIFFRACTION_CZT_AP calculates vectorial diffraction by 2D chirp-z
%transform (CZT) in Cartesian coordinate system.
%
% INPUT********************************************************************
% Obj.NA: scalar value, numerical aperture of objective lens
% Obj.n: scalar value, refractive index in focal space
% Beam.wavelength: scalar value, Beam.wavelength of light
% Beam.amp: M*M matrix, amplitude distribution on pupil plane
% Beam.phs: M*M matrix, phase distribution on pupil plane
% Beam.plr: M*M*3 matrix, polarization distribution on pupil plane
% Scope.us: 1*N array, representing u axis of PSF
% Scope.vs: 1*N array, representing v axis of PSF
% Scope.ws: 1*N array, representing w axis of PSF
% Beam.PupilRes: scalar value, resolution of pupil
%
% OUTPUT*******************************************************************
% [Ex,Ey,Ez]: electric field of focused field
% 
% *************************************************************************
% LIU Xin
% liuxin24@hku.hk
% Apr.23, 2021
% 
% updated by HU Yiwen
% huyw@zju.edu.cn

%% pupil sampling
r0 = 1;  % radius of pupil
[xx,yy] = meshgrid(linspace(-r0,r0,Beam.PupilRes));  % cartesian coordinate of pupil plane
[phi,rho] = cart2pol(xx,yy);  % polar coordinate of pupil plane

thetaMax = asin(Obj.NA/Obj.n);  % maximum convergent angle of objective
sinTheta = sin(thetaMax).*rho;
sinTheta(rho>r0) = 0;
theta = asin(sinTheta);  % convergance angle of each ray (theta)

%% interpolation of pupil
if ~(size(Beam.amp) == Beam.PupilRes)
    Beam.amp = pupilInterp(Beam.amp,Beam.PupilRes);
end
if ~(size(Beam.phs) == Beam.PupilRes)
    Beam.phs = pupilInterp(Beam.phs,Beam.PupilRes);
end

px = Beam.plr(:,:,1);
py = Beam.plr(:,:,2);

if ~(size(px) == Beam.PupilRes)
    px = pupilInterp(px,Beam.PupilRes);
end

if ~(size(py) == Beam.PupilRes)
    py = pupilInterp(py,Beam.PupilRes);
end

%% pupil function
E0 = Beam.amp.*exp(1i.*Beam.phs);

% remove parts outside numerical aperture
E0(rho>r0) = 0;

%% check sampling condition
fc = Obj.NA/Beam.wavelength;  % cut-off frequency of objective
df = 2*fc/(Beam.PupilRes-1);  % sampling period
fs_xy = 1/df;  % sampling frequency

[uus,vvs] = meshgrid(Scope.us,Scope.vs);

global theta2 phi2
xxs = cos(theta2)*cos(phi2)*uus - sin(phi2)*vvs + sin(theta2)*cos(phi2)*Scope.ws;
yys = cos(theta2)*sin(phi2)*uus + cos(phi2)*vvs + sin(theta2)*sin(phi2)*Scope.ws;
zzs = -sin(theta2)*uus + cos(theta2)*Scope.ws;

fx_max = max(abs(xxs(:)));
fy_max = max(abs(yys(:)));
fz_max = max(abs(zzs(:)));

dfz_max = Obj.n/Beam.wavelength*...
    (sqrt(1-(sin(thetaMax)*(Beam.PupilRes-3)/(Beam.PupilRes-1))^2)-cos(thetaMax));
fs_z_min = 1/dfz_max;

% Nyquist–Shannon sampling theorem
sampling_condition_x = fs_xy > 2*fx_max;
sampling_condition_y = fs_xy > 2*fy_max;
sampling_condition_z1 = fs_z_min > 2*fz_max;

% sampling condition for kz
% the phase kz*z must not change by more than pi between neighboring
% sampling points in the pupil plane.
M0 = 2*Obj.NA^2/sqrt(Obj.n^2-Obj.NA^2)*(max(abs(zzs(:)))./Beam.wavelength);
M0 = round(M0);
sampling_condition_z2 = Beam.PupilRes > 2*M0;

if ~(sampling_condition_x && sampling_condition_y &&...
        sampling_condition_z1 && sampling_condition_z2)
    pupilResX = round((Beam.PupilRes-1)*2*fx_max/fs_xy) + 1;
    pupilResY = round((Beam.PupilRes-1)*2*fy_max/fs_xy) + 1;
    
    Chi = Beam.wavelength/Obj.n*(1/(2*fz_max))+cos(thetaMax);
    if fz_max == 0
        pupilResZ1 = 1;
    else
        pupilResZ1 = 2/(1-sqrt(1-Chi^2)/sin(thetaMax)) + 1;
    end
    pupilResZ2 = 2*M0;
    
    % the required least pupil resolution
    pupilResReq = max([pupilResX,pupilResY,pupilResZ1,pupilResZ2]);
    warning(['Pupil resolution should be larger than ',...
        num2str(pupilResReq), ' avoiding aliasing!']);
end

%% polarization in image space (exit pupil)
[Px,Py,Pz] = coortrans(px,py,0,theta,phi,'o2i');

%% amplitude apodization factor for energy conservation
% different objective may have different AF!!!
AAF = sqrt(cos(theta));  % only for objective obeying the sine condition

%% plane waves in image space
% amplitude projection factor from (theta, phi) to (kx, ky) coordinate for
% integral (FFT)
APF = 1./cos(theta);

%% amplitude transformation
E0 = E0.*AAF.*APF;

%% the coefficient in front of Debye integral
f = r0/sin(thetaMax);  % focal length
prefix = -1i*Obj.n*f/Beam.wavelength;
E0 = prefix*E0;

%% plane waves in image space (three polarization components)
Ex0 = E0.*Px;
Ey0 = E0.*Py;
Ez0 = E0.*Pz;
Ex0 = rot90(Ex0,2);
Ey0 = rot90(Ey0,2);
Ez0 = rot90(Ez0,2);

%% frequency domain
fn = Obj.n/Beam.wavelength;
fxy1 = linspace(-fc,fc,Beam.PupilRes);
[fxx1,fyy1] = meshgrid(fxy1);
fzz1_sq = fn.^2 - fxx1.^2 - fyy1.^2;
fzz1_sq(rho>r0) = 0;
fzz1 = sqrt(fzz1_sq);

fuu = cos(theta2)*cos(phi2)*fxx1+cos(theta2)*sin(phi2)*fyy1-sin(theta2)*fzz1;
fvv = -sin(phi2)*fxx1+cos(phi2)*fyy1;

minFu = min(fuu(:)); maxFu = max(fuu(:));
minFv = min(fvv(:)); maxFv = max(fvv(:));
fu = linspace(minFu,maxFu,Beam.PupilRes);
fv = linspace(minFv,maxFv,Beam.PupilRes);
[fuu2,fvv2] = meshgrid(fu,fv);
fww2_p = sqrt(abs(fn.^2-fuu2.^2-fvv2.^2));
fww2_n = -fww2_p;

% preparation for interpolation
fxx2_p = cos(theta2)*cos(phi2)*fuu2-sin(phi2)*fvv2+sin(theta2)*cos(phi2)*fww2_p;
fyy2_p = cos(theta2)*sin(phi2)*fuu2+cos(phi2)*fvv2+sin(theta2)*sin(phi2)*fww2_p;
fzz2_p = -sin(theta2)*fuu2+cos(theta2)*fww2_p;
fxx2_n = cos(theta2)*cos(phi2)*fuu2-sin(phi2)*fvv2+sin(theta2)*cos(phi2)*fww2_n;
fyy2_n = cos(theta2)*sin(phi2)*fuu2+cos(phi2)*fvv2+sin(theta2)*sin(phi2)*fww2_n;
fzz2_n = -sin(theta2)*fuu2+cos(theta2)*fww2_n;

mask1 = (fxx2_p.^2+fyy2_p.^2>fc.^2); % Exclude non-existent coordinate points
mask2 = (fxx2_n.^2+fyy2_n.^2>fc.^2);
mask3 = (fuu2.^2+fvv2.^2>fn.^2);
mask4 = (fzz2_p<fn*cos(thetaMax));
mask5 = (fzz2_n<fn*cos(thetaMax));
mask6 = (fww2_p == 0);
mask7 = (fzz2_p>fn|fzz2_p<0);
mask8 = (fzz2_n>fn|fzz2_n<0);

Ex0_p = interp2(fxx1,fyy1,Ex0,fxx2_p,fyy2_p,'linear',0);
Ey0_p = interp2(fxx1,fyy1,Ey0,fxx2_p,fyy2_p,'linear',0);
Ez0_p = interp2(fxx1,fyy1,Ez0,fxx2_p,fyy2_p,'linear',0);
Ex0_n = interp2(fxx1,fyy1,Ex0,fxx2_n,fyy2_n,'linear',0);
Ey0_n = interp2(fxx1,fyy1,Ey0,fxx2_n,fyy2_n,'linear',0);
Ez0_n = interp2(fxx1,fyy1,Ez0,fxx2_n,fyy2_n,'linear',0);
Ex0_p(~isfinite(Ex0_p)) = 0;
Ex0_n(~isfinite(Ex0_n)) = 0;
Ey0_p(~isfinite(Ey0_p)) = 0;
Ey0_n(~isfinite(Ey0_n)) = 0;
Ez0_p(~isfinite(Ez0_p)) = 0;
Ez0_n(~isfinite(Ez0_n)) = 0;

Ex0_p(mask1 | mask3 | mask4 | mask6 | mask7) = nan;
fww2_p(mask1| mask3 | mask4 | mask6 | mask7) = nan;
Ey0_p(mask1 | mask3 | mask4 | mask6 | mask7) = nan;
Ez0_p(mask1 | mask3 | mask4 | mask6 | mask7) = nan;
Ex0_n(mask2 | mask3 | mask5 | mask8) = nan;
fww2_n(mask2 | mask3 | mask5 | mask8) = nan;
Ey0_n(mask2 | mask3 | mask5 | mask8) = nan;
Ez0_n(mask2 | mask3 | mask5 | mask8) = nan;

% Jacobian determinant
J_p = -sin(theta2)*fuu2./fww2_p + cos(theta2);
J_n = -sin(theta2)*fuu2./fww2_n + cos(theta2);
J_p(isnan(J_p)) = 0;
J_n(isnan(J_n)) = 0;

Ex0_p = Ex0_p.*abs(J_p); Ex0_n = Ex0_n.*abs(J_n); 
Ey0_p = Ey0_p.*abs(J_p); Ey0_n = Ey0_n.*abs(J_n);
Ez0_p = Ez0_p.*abs(J_p); Ez0_n = Ez0_n.*abs(J_n);

Ex0_p(isnan(Ex0_p)) = 0; Ex0_n(isnan(Ex0_n)) = 0;
Ey0_p(isnan(Ey0_p)) = 0; Ey0_n(isnan(Ey0_n)) = 0;
Ez0_p(isnan(Ez0_p)) = 0; Ez0_n(isnan(Ez0_n)) = 0;

fww2_p(isnan(fww2_p)) = 0;
fww2_n(isnan(fww2_n)) = 0;

lu = length(Scope.us); lv = length(Scope.vs); 

%% calculation
defocusTerm_p = exp(1i*2*pi*fww2_p.*Scope.ws);
defocusTerm_n = exp(1i*2*pi*fww2_n.*Scope.ws);

Ewx_p = Ex0_p.*defocusTerm_p;
ExHold_p = myiczt(Ewx_p,lv,fv,Scope.vs);
ExHold_p = permute(ExHold_p,[2 1]);
ExHold_p = myiczt(ExHold_p,lu,fu,Scope.us);

Ewx_n = Ex0_n.*defocusTerm_n;
ExHold_n = myiczt(Ewx_n,lv,fv,Scope.vs);
ExHold_n = permute(ExHold_n,[2 1]);
ExHold_n = myiczt(ExHold_n,lu,fu,Scope.us);

Ex = permute(ExHold_p+ExHold_n,[2 1]);

Ewy_p = Ey0_p.*defocusTerm_p;
EyHold_p = myiczt(Ewy_p,lv,fv,Scope.vs);
EyHold_p = permute(EyHold_p,[2 1]);
EyHold_p = myiczt(EyHold_p,lu,fu,Scope.us);

Ewy_n = Ey0_n.*defocusTerm_n;
EyHold_n = myiczt(Ewy_n,lv,fv,Scope.vs);
EyHold_n = permute(EyHold_n,[2 1]);
EyHold_n = myiczt(EyHold_n,lu,fu,Scope.us);
Ey = permute(EyHold_p+EyHold_n,[2 1]);


Ewz_p = Ez0_p.*defocusTerm_p;
EzHold_p = myiczt(Ewz_p,lv,fv,Scope.vs);
EzHold_p = permute(EzHold_p,[2 1]);
EzHold_p = myiczt(EzHold_p,lu,fu,Scope.us);

Ewz_n = Ez0_n.*defocusTerm_n;
EzHold_n = myiczt(Ewz_n,lv,fv,Scope.vs);
EzHold_n = permute(EzHold_n,[2 1]);
EzHold_n = myiczt(EzHold_n,lu,fu,Scope.us);
Ez = permute(EzHold_p+EzHold_n,[2 1]);

disp('VectorialDiffraction_CZT_AP has completed!');

end