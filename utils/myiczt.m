function Xz = myiczt(xn, Mz, xin, xout)
%MYICZT perform inverse chirp-z transform in the first dimension of xn
%
% Xin Liu
% liuxin24@hku.hk
% Dec.19, 2020
% 
% updated by Yiwen Hu
% huyw@zju.edu.cn

%% check input
sz = size(xn);

if length(sz) > 3
    error('Invalid input dimensions!');
end

%% Preparation of czt
xout1 = xout(1);
pixelsize0 = xin(2) - xin(1);
fs = 1/pixelsize0;
pixelsize1 = xout(2) - xout(1);

% start point on unit circle
A = exp(-1i*2*pi*xout1/fs);

% step size on unit circle
W = exp( 1i*2*pi*pixelsize1/fs);

%% chirp z transform
% size of input signal
[Nx1, Nx2] = size(xn);

kk = ((-Nx1+1):max(Mz-1,Nx1-1)).';  % full length for FFT
ind_x = (0:(Nx1-1)).';  % index of the sequence of the first dimension of xn
WW = W.^((kk.^2)./2);
AW = A.^(-ind_x).*WW(Nx1+ind_x); % chose the indices 'm+nn' of ww
AW = repmat(AW, [1,Nx2]);

% two functions of n
fn = xn.*AW;
hn = WW(1:Nx1+Mz-1).^-1;

%% fast convolution via FFT
% length for power-of-two FFT
nfft = 2^nextpow2(Nx1+Mz-1);
Fr = fft(fn,nfft);
Hr = fft(hn,nfft);
Hr = repmat(Hr, [1,Nx2]);

Gr = Fr.*Hr;  % multiplication in frequency domain
Xz = ifft(Gr);  % inverse fourier transform
Xz = Xz(Nx1:(Nx1+Mz-1),:,:);  % extract effective data

%% multiply prefix
prefix = repmat(WW(Nx1:Nx1+Mz-1), [1,Nx2]);
Xz = prefix.*Xz;

%% phase shift
Mshiftx = min(xin(:))/(xin(2)-xin(1));

xout = reshape(xout, [length(xout) 1]);
xss = xout.*Mshiftx/fs;
xss = repmat(xss,[1 Nx2]);
Pshift = exp( 1i*2*pi*xss);

Xz = Xz.*Pshift;
end