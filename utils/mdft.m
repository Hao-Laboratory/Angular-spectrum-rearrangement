function Out = mdft(In,x,y,fx,fy)
%MDFT calculates 2D DFT using matrix triple product (MTP)
%
% LIU Xin
% liuxin24@hku.hk
% Dec.29, 2021

x = reshape(x,[],1);
y = reshape(y,1,[]);
fx = reshape(fx,1,[]);
fy = reshape(fy,[],1);
Mx = exp(-2*pi*1i*x*fx);
My = exp(-2*pi*1i*fy*y);
Out = My*In*Mx;
end