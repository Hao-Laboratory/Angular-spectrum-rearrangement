function [Ex,Ey,Ez] = VectorialDiffraction_ASR_AP(Obj,Beam,Scope,number_u,number_v)
%VECTORIALDIFFRACTION_ASR_AP calculates vectorial diffraction by angular 
%spectrum rearrangement (ASR) in Cartesian coordinate system.
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
% number_u: number of target samples in the frequency domain after merging 
% along the u-axis
% number_v: number of target samples in the frequency domain after merging 
% along the v-axis
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
f = r0/sin(thetaMax);
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
fzz1_sq(rho>r0) = nan;
fzz1 = sqrt(fzz1_sq);

fuu = cos(theta2)*cos(phi2)*fxx1+cos(theta2)*sin(phi2)*fyy1-sin(theta2)*fzz1;
fvv = -sin(phi2)*fxx1+cos(phi2)*fyy1;
fww = sin(theta2)*cos(phi2)*fxx1+sin(theta2)*sin(phi2)*fyy1+cos(theta2)*fzz1;

if (theta2 == 0)&&(phi2 == 0)
    Ewx = Ex0;
    Ewy = Ey0;
    Ewz = Ez0;
    fux_eff = fxy1; fuy_eff = fux_eff; fuz_eff = fux_eff;
    fvx_eff = fux_eff.'; fvy_eff = fvx_eff; fvz_eff = fvx_eff;
    fww = fzz1;
else
    spa = fxx1.^2+fyy1.^2<=fc^2;
    Ewx = Ex0(spa);
    Ewy = Ey0(spa);
    Ewz = Ez0(spa);
    fww = fww(spa);

    [fu, idx_fu] = unique_tol(fuu(spa),number_u);
    [fv, idx_fv] = unique_tol(fvv(spa),number_v);
end

lu = length(Scope.us);
lv = length(Scope.vs);

%% calculation
Ewx_temp = Ewx.*exp(1i*2*pi*fww.*Scope.ws);
Ewy_temp = Ewy.*exp(1i*2*pi*fww.*Scope.ws);
Ewz_temp = Ewz.*exp(1i*2*pi*fww.*Scope.ws);
Ewx_temp(isnan(Ewx_temp)) = 0;
Ewy_temp(isnan(Ewy_temp)) = 0;
Ewz_temp(isnan(Ewz_temp)) = 0;

if (theta2 ~= 0 || phi2 ~= 0)
    % Place coordinates in xy plane at the counterparts in the tilted plane
    Ewx = sparse(idx_fv,idx_fu,Ewx_temp);
    Ewy = sparse(idx_fv,idx_fu,Ewy_temp);
    Ewz = sparse(idx_fv,idx_fu,Ewz_temp);

    % weight calculation
    absEwx_temp = sparse(idx_fv,idx_fu,abs(Ewx_temp));
    absEwy_temp = sparse(idx_fv,idx_fu,abs(Ewy_temp));
    absEwz_temp = sparse(idx_fv,idx_fu,abs(Ewz_temp));

    fuEwx = fuu(spa).*abs(Ewx_temp);
    fvEwx = fvv(spa).*abs(Ewx_temp);
    fuEwy = fuu(spa).*abs(Ewy_temp);
    fvEwy = fvv(spa).*abs(Ewy_temp);
    fuEwz = fuu(spa).*abs(Ewz_temp);
    fvEwz = fvv(spa).*abs(Ewz_temp);

    fux_eff = full(sparse(ones(size(idx_fu)), idx_fu, fuEwx))./full(sum(absEwx_temp));
    fvx_eff = full(sparse(idx_fv, ones(size(idx_fv)), fvEwx))./full(sum(absEwx_temp,2));
    if find(isnan(fux_eff))
        idx = isnan(fux_eff);
        fux_eff(idx) = fu(idx);
    end
    if find(isnan(fvx_eff))
        idx = isnan(fvx_eff);
        fvx_eff(idx) = fv(idx);
    end

    fuy_eff = full(sparse(ones(size(idx_fu)), idx_fu, fuEwy))./full(sum(absEwy_temp));
    fvy_eff = full(sparse(idx_fv, ones(size(idx_fv)), fvEwy))./full(sum(absEwy_temp,2));
    if find(isnan(fuy_eff))
        idx = isnan(fuy_eff);
        fuy_eff(idx) = fu(idx);
    end
    if find(isnan(fvy_eff))
        idx = isnan(fvy_eff);
        fvy_eff(idx) = fv(idx);
    end

    fuz_eff = full(sparse(ones(size(idx_fu)), idx_fu, fuEwz))./full(sum(absEwz_temp));
    fvz_eff = full(sparse(idx_fv, ones(size(idx_fv)), fvEwz))./full(sum(absEwz_temp,2));
    if find(isnan(fuz_eff))
        idx = isnan(fuz_eff);
        fuz_eff(idx) = fu(idx);
    end
    if find(isnan(fvz_eff))
        idx = isnan(fvz_eff);
        fvz_eff(idx) = fv(idx);
    end
end

Mux = exp(1i*2*pi*fux_eff.'*Scope.us);
Mvx = exp(1i*2*pi*Scope.vs.'*fvx_eff.');
Muy = exp(1i*2*pi*fuy_eff.'*Scope.us);
Mvy = exp(1i*2*pi*Scope.vs.'*fvy_eff.');
Muz = exp(1i*2*pi*fuz_eff.'*Scope.us);
Mvz = exp(1i*2*pi*Scope.vs.'*fvz_eff.');

M = length(fux_eff);
N = length(fvx_eff);
R = lv;
S = lu;

CostMM_Right = M*S*(N+R);
CostMM_Left = N*R*(M+S);

if CostMM_Left<=CostMM_Right
    Ex = Mvx*Ewx*Mux;
    Ey = Mvy*Ewy*Muy;
    Ez = Mvz*Ewz*Muz;
else
    Ex = Mvx*(Ewx*Mux);
    Ey = Mvy*(Ewy*Muy);
    Ez = Mvy*(Ewz*Muy);
end

disp('VectorialDiffraction_ASR_AP has completed!');

end