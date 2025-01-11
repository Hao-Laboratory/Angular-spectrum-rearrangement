function [Ex,Ey,Ez] = VectorialDiffraction_GT_AP(Obj,Beam,Scope)
%VECTORIALDIFFRACTION_GT_AP calculates vectorial diffraction by direct 
%integral in Cartesian coordinate system.
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
r0 = 1;  % radius of the pupil is unit
[xx,yy] = meshgrid(linspace(-r0,r0,Beam.PupilRes));  % cartesian coordinate of pupil plane
[phi,rho] = cart2pol(xx,yy);  % polar coordinate of pupil plane

thetaMax = asin(Obj.NA/Obj.n);  % maximum convergent angle of objective
sinTheta = sin(thetaMax).*rho;
sinTheta(rho>1) = 0;
theta = asin(sinTheta);  % convergance angle of each ray(theta)

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

%% wave vector
k0 = 2*pi/Beam.wavelength;  % wavenumber in vacuum
kc = k0*Obj.NA;  % cut-off frequency of objective (in k-space)
kx = linspace(-kc,kc,Beam.PupilRes);  % spatial frequency coordinates
ky = kx;
[kxx,kyy] = meshgrid(kx,ky);
kzz = k0*Obj.n*cos(theta);

lu = length(Scope.us);
lv = length(Scope.vs);
lw = length(Scope.ws);
totalNum = lu*lv*lw;

%% vectorial diffraction integral calculation
Ex = zeros(lv,lu,lw);
Ey = zeros(lv,lu,lw);
Ez = zeros(lv,lu,lw);
for ii = 1:totalNum
    Ex(ii) = sum(Ex0.*exp(1i*(kxx*xxs(ii) + kyy*yys(ii) + kzz*zzs(ii))),'all');
    Ey(ii) = sum(Ey0.*exp(1i*(kxx*xxs(ii) + kyy*yys(ii) + kzz*zzs(ii))),'all');
    Ez(ii) = sum(Ez0.*exp(1i*(kxx*xxs(ii) + kyy*yys(ii) + kzz*zzs(ii))),'all');
    textwaitbar(ii, totalNum, 'VectorialDiffraction_GT_AP in progress') 
end

end