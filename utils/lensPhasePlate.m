function phs = lensPhasePlate(pupilRes,D,f,lambda,n)
%LENSPHASEPLATE generates lens phase
%
% LIU Xin
% liuxin24@hku.hk
% Mar.27, 2022

[xd,yd] = meshgrid(linspace(-D/2,D/2,pupilRes));
[~,rho] = cart2pol(xd,yd);
s = 1;
Nmin = ceil(s*4*(D/2)^2/(lambda*f/n));
if pupilRes < Nmin
    warning(['The pupil resolution must be larger than ',num2str(Nmin),' avoiding aliasing!']);
end
k = 2*pi*n/lambda;  % wavenumber
phs = -k./(2*f).*(xd.^2+yd.^2);  % quadratic phase factor
phs(rho>D/2) = 0;
phs = mod(phs,2*pi);
end

