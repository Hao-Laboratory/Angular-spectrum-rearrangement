clear;clc;
close all;

%% calculation mode
DiffractionMethod = {'assmm','asczt','asap'};
n_method = length(DiffractionMethod);

% number of target samples in the frequency domain after merging
number_u = 1100;  % This value is chosen to demonstrate that our method is more accurate than the control method while maintaining comparable computational efficiency.
number_v = number_u;

LS = 6387.5;  % physical size of source plane (um)
LT = 255.5;  % physical size of observation plane

% colormap for visualization
cmap_E = 'hot';
cmap_P = 'jet';

%% set the angle of the tilted plane
global theta2 phi2
theta2 = 50/180*pi; phi2 = 30/180*pi;
disp([rad2deg(theta2), rad2deg(phi2)]);

%% beam parameters
Beam.Res = 256;  % number of samples on the source plane

Beam.PixelSize = LS/(Beam.Res-1);
[xx,yy] = meshgrid(-(Beam.Res-1)/2:(Beam.Res-1)/2,-(Beam.Res-1)/2:(Beam.Res-1)/2);
XX = xx.*Beam.PixelSize;
YY = yy.*Beam.PixelSize;
[phi,rho] = cart2pol(XX,YY);

% amplitude
load('amp.mat');
Beam.amp = amp;

Beam.wavelength = 785e-3;
n = 1;  % refractive index of the medium

% phase
load('phs1.mat');
Beam.phs = phs + lensPhasePlate(Beam.Res,LS,225e3,Beam.wavelength,n);

figure('Name','Field Distribution on the Source Plane');
tiledlayout(1,2);

nexttile;
pupilshow(Beam.amp);
colorbar('Ticks',[min(Beam.amp(:)),max(Beam.amp(:))]);
colormap(gca,cmap_E);
title('Amplitude');

nexttile;
pupilshow(wrapTo2Pi(phs));
axis image xy off;
colorbar('Ticks',[0 2*pi]);caxis([0 2*pi]);
colormap(gca,cmap_P);
title('Phase');

%% camera parameters
d = 225e3;  % wave propagation distance

% resolution of camera
c = 256;

% pixel size of camera
pixelSize = LT/(c-1);

Scope.us = -LT/2:pixelSize:LT/2;
Scope.vs = -LT/2:pixelSize:LT/2;
Scope.ws = 0;
Scope.zs = d;

%% phase induced by tilted plane
[uu,vv] = meshgrid(Scope.us,Scope.vs);
zz = -sin(theta2)*uu;
phs_tilt = 2*pi*zz/Beam.wavelength;

%% calcalation
E_gt = ScalarDiffraction(Beam,Scope,DiffractionMethod{n_method},number_u,number_v);  % ground truth
I_gt_tmp = abs(E_gt).^2;
I_gt = I_gt_tmp./max(I_gt_tmp(:));
I_gt_min = min(I_gt(:)); I_gt_max = max(I_gt(:));

I_all = cell(1,n_method);
phs_all = cell(1,n_method);
E_all = cell(1,n_method);

I_all{n_method} = I_gt;
phs_gt = wrapTo2Pi(angle(E_gt)-phs_tilt);
phs_all{n_method} = phs_gt;
E_all{n_method} = E_gt/max(abs(E_gt(:)));

fig_db = figure('Name','Simulation Results');
subplot(2,n_method,n_method);
imagesc(I_all{n_method});
axis image xy off;
caxis([0 I_gt_max]);
colorbar('Ticks',[0 I_gt_max],'TickLabels', ...
    {sprintf('%.1e',0),sprintf('%.1e',I_gt_max)});
colormap(gca,cmap_E);
title('Intensity_{GT}');

subplot(2,n_method,2*n_method);
imagesc(phs_all{n_method});
axis image xy off;
caxis([0 2*pi]);
colorbar('Ticks',[0 2*pi],'TickLabels', ...
    {sprintf('%.1e',0),sprintf('%.1e',2*pi)});
colormap(gca,cmap_P);
title('Phase_{GT}');

for ii = 1:(n_method-1)
    Eout = ScalarDiffraction(Beam,Scope,DiffractionMethod{ii},number_u,number_v);

    E_all{ii} = Eout/max(abs(Eout(:)));
    I = abs(Eout).^2;
    I = I./max(I(:));
    I_all{ii} = I;
    phs_all{ii} = wrapTo2Pi(angle(Eout)-phs_tilt);

    figure(fig_db);
    subplot(2,n_method,ii);
    imagesc(I_all{ii});
    axis image xy off;
    caxis([0 I_gt_max]);
    colorbar('Ticks',[0 I_gt_max],'TickLabels', ...
        {sprintf('%.1e',0),sprintf('%.1e',I_gt_max)});
    colormap(gca,cmap_E);
    if strcmp(DiffractionMethod{ii},'assmm')
        title('Intensity_{Ours}');
    elseif strcmp(DiffractionMethod{ii},'asczt')
        title('Intensity_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end

    subplot(2,n_method,ii+n_method);
    imagesc(phs_all{ii});
    axis image xy off;
    caxis([0 2*pi]);
    colorbar('Ticks',[0 2*pi],'TickLabels', ...
        {sprintf('%.1e',0),sprintf('%.1e',2*pi)});
    colormap(gca,cmap_P);
    if strcmp(DiffractionMethod{ii},'assmm')
        title('Phase_{Ours}');
    elseif strcmp(DiffractionMethod{ii},'asczt')
        title('Phase_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end
end

%% compare accuracy
fig_cmp = figure('Name','Comparison of Accuracy');

dev = cell(1,n_method-1);
for jj = 1:(n_method-1)
    figure(fig_cmp);
    subplot(1,n_method-1,jj);
    dev{jj} = abs(E_all{jj}-E_all{n_method});
    min_dev = min(dev{jj}(:));
    max_dev = max(dev{jj}(:));
    imagesc(dev{jj});
    caxis([min_dev,max_dev]);
    axis image xy off;
    colorbar('Ticks',[min_dev,max_dev],'TickLabels', ...
        {sprintf('%.1e',min_dev),sprintf('%.1e',max_dev)});
    colormap(gca,cmap_E);
    if strcmp(DiffractionMethod{jj},'assmm')
        title('Deviation_{Ours}');
    elseif strcmp(DiffractionMethod{jj},'asczt')
        title('Deviation_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end
end

%% evaluate errors and time consumption
sigma_mtp = sum(abs(E_all{1} - E_all{3}).^2,'all')/sum(abs(E_all{3}).^2,'all')
sigma_czt = sum(abs(E_all{2} - E_all{3}).^2,'all')/sum(abs(E_all{3}).^2,'all')

time_MTP = timeit(@() ScalarDiffraction_ASR_AP(Beam,Scope,number_u,number_v))
time_CZT = timeit(@() ScalarDiffraction_CZT_AP(Beam,Scope))