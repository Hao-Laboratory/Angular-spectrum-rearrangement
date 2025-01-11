function Mout = pupilInterp(Min,Nout,interpMethod)
%PUPILINTERP: interpolation of the entrance pupil matrix.
%
%--------------------------------------------------------------------------
% Inputs:
% Min   -previous pupil matrix
% Nout  -pixel number of the pupil matrix after interpolation
% interpMethod -the method for interpolation
%
% Output:
% Mout  -pupil matrix after interpolation
%--------------------------------------------------------------------------
%
% Author: Xin Liu
% Email: liuxin24@hku.hk
% Aug.17, 2020

if nargin == 2
    interpMethod = 'nearest';
end

[M,N] = size(Min);
[xx,yy] = meshgrid(linspace(-1,1,N),...
    linspace(-1,1,M)); % create a Cartesian coordinate for Min

[xxi,yyi] = meshgrid(linspace(-1,1,Nout));
Mout = interp2(xx,yy,Min,xxi,yyi,interpMethod);
end