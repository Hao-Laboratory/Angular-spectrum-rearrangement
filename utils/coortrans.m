function [Ex_out,Ey_out,Ez_out] = coortrans(Ex_in,Ey_in,Ez_in,theta,phi,option)
%coortrans is used to transform coordinate between the object space and the
% image space of single lens
%
% LIU Xin
% liuxin24@hku.hk
% Dec.11, 2020

switch option
    case 'o2i'  % from object space to image space
        M11 = 1+(cos(theta)-1).*cos(phi).^2;
        M12 = (cos(theta)-1).*cos(phi).*sin(phi);
        M21 = (cos(theta)-1).*cos(phi).*sin(phi);
        M22 = 1+(cos(theta)-1).*sin(phi).^2;
        M31 = sin(theta).*cos(phi);
        M32 = sin(theta).*sin(phi);
        M13 = 0;
        M23 = 0;
        M33 = 0;
    case 'i2o'  % from image space to object space
        M11 = 1+(cos(theta)-1).*cos(phi).^2;
        M21 = (cos(theta)-1).*cos(phi).*sin(phi);
        M12 = (cos(theta)-1).*cos(phi).*sin(phi);
        M22 = 1+(cos(theta)-1).*sin(phi).^2;
        M13 = sin(theta).*cos(phi);
        M23 = sin(theta).*sin(phi);
        M31 = 0;
        M32 = 0;
        M33 = 0;
    otherwise
        error('No such an option, please try again!');
end
Ex_out = (M11.*Ex_in + M12.*Ey_in + M13.*Ez_in);
Ey_out = (M21.*Ex_in + M22.*Ey_in + M23.*Ez_in);
Ez_out = (M31.*Ex_in + M32.*Ey_in + M33.*Ez_in);
end