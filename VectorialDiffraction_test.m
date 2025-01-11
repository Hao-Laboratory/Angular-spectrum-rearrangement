clear;clc;
close all;

%% calculation mode
FocusingMethod = {'smm_ap','czt_ap','int_ap'};
n_method = length(FocusingMethod);

% number of target samples in the frequency domain after merging
number_u = 256; number_v = number_u;  % This value is chosen to demonstrate that our method is more accurate than the control method while maintaining comparable computational efficiency.

% colormap for visualization
cmap_E = 'hot';
cmap_P = 'jet';

%% set the angle of the tilted plane
global theta2 phi2
theta2 = 130/180*pi; phi2 = 30/180*pi;
disp([rad2deg(theta2), rad2deg(phi2)]);

%% beam parameters
Beam.PupilRes = 128;  % number of samples on the pupil plane

%% scope
psfRes = 100;
HalfScope = 1550e-3;  % um

pixelSize = 2*HalfScope/(psfRes-1)
Scope.us = linspace(-HalfScope,HalfScope,psfRes);
Scope.vs = Scope.us;
Scope.ws = 0;

%% objective
Obj.NA = 1.35; Obj.n = 1.406;

%% beam wavelength
Beam.wavelength = 785e-3;

%% amplitude
load('amp.mat');
Beam.amp = amp;

%% phase
load('phs2.mat');
Beam.phs = phs;

%% polarization
load plr.mat;
Beam.plr = plr;

%% show pupil function
fig_pupil = figure('Name','Field Distribution on the Pupil Plane');
tiledlayout(1,3);

nexttile;
dis_plr = imread('plr.png');
image(dis_plr);
axis image off;
title('Polarization');

nexttile;
pupilshow(Beam.amp);colormap(gca,cmap_E);
colorbar('Ticks',[min(Beam.amp(:)),max(Beam.amp(:))]);
title('Amplitude');

nexttile;
pupilshow(wrapTo2Pi(Beam.phs));colorbar('Ticks',[0 2*pi]);caxis([0 2*pi]);colormap(gca,cmap_P);
title('Phase');

%% phase induced by tilted plane
[uu,vv] = meshgrid(Scope.us,Scope.vs);
zz = -sin(theta2)*uu;
phs_tilt = 2*pi*zz/Beam.wavelength;

%% calcalation
fig_PSF = figure('Name','Simulation Results');
Ix_all = cell(1,n_method);
Phsx_all = cell(1,n_method);
Iy_all = cell(1,n_method);
Phsy_all = cell(1,n_method);
Iz_all = cell(1,n_method);
Phsz_all = cell(1,n_method);
I_all = cell(1,n_method);
Ex_all = cell(1,n_method);
Ey_all = cell(1,n_method);
Ez_all = cell(1,n_method);

[Ex_gt,Ey_gt,Ez_gt] = VectorialDiffraction(Obj,Beam,Scope,FocusingMethod{n_method},number_u,number_v);

Ex_all{n_method} = Ex_gt/max(abs(Ex_gt(:)));
Ey_all{n_method} = Ey_gt/max(abs(Ey_gt(:)));
Ez_all{n_method} = Ez_gt/max(abs(Ez_gt(:)));

phsEx_gt = wrapTo2Pi(angle(Ex_gt) - phs_tilt);
phsEy_gt = wrapTo2Pi(angle(Ey_gt) - phs_tilt);
phsEz_gt = wrapTo2Pi(angle(Ez_gt) - phs_tilt);

Ix_gt_tmp = abs(Ex_gt).^2; Iy_gt_tmp = abs(Ey_gt).^2; Iz_gt_tmp = abs(Ez_gt).^2;
I_gt_tmp = Ix_gt_tmp + Iy_gt_tmp + Iz_gt_tmp;

Ix_gt = Ix_gt_tmp./max(Ix_gt_tmp(:));
Iy_gt = Iy_gt_tmp./max(Iy_gt_tmp(:));
Iz_gt = Iz_gt_tmp./max(Iz_gt_tmp(:));
I_gt = I_gt_tmp./max(I_gt_tmp(:));

Ixmax = max(Ix_gt(:));
Iymax = max(Iy_gt(:));
Izmax = max(Iz_gt(:));
Imax = max(I_gt(:));

figure(fig_PSF);
% plot the intensity of X polarization
subplot(n_method,7,(n_method-1)*7+1);
imagesc(Scope.us,Scope.vs,Ix_gt);
axis image xy off;colormap(gca,cmap_E);title('IntensityX_{GT}');
Ix_all{n_method} = Ix_gt;
caxis([0,Ixmax]);
colorbar('Ticks',[0,Ixmax],'TickLabels',...
    {sprintf('%.2e',0),sprintf('%.2e',Ixmax)});

% plot the intensity of Y polarization
subplot(n_method,7,(n_method-1)*7+3);
imagesc(Scope.us,Scope.vs,Iy_gt);
axis image xy off;colormap(gca,cmap_E);title('IntensityY_{GT}');
Iy_all{n_method} = Iy_gt;
caxis([0,Iymax]);
colorbar('Ticks',[0,Iymax],'TickLabels',...
    {sprintf('%.2e',0),sprintf('%.2e',Iymax)});

% plot the intensity of Z polarization
subplot(n_method,7,(n_method-1)*7+5);
imagesc(Scope.us,Scope.vs,Iz_gt);
axis image xy off;colormap(gca,cmap_E);title('IntensityZ_{GT}');
Iz_all{n_method} = Iz_gt;
caxis([0, Izmax]);
colorbar('Ticks',[0,Izmax],'TickLabels',...
    {sprintf('%.2e',0),sprintf('%.2e',Izmax)});

% plot the intensity
subplot(n_method,7,n_method*7);
imagesc(Scope.us,Scope.vs,I_gt);
axis image xy off;colormap(gca,cmap_E);title('Intensity_{GT}');
I_all{n_method} = I_gt;
caxis([0, Imax]);
colorbar('Ticks',[0,Imax],'TickLabels',...
    {sprintf('%.2e',0),sprintf('%.2e',Imax)});

% plot the phase of X polarization
subplot(n_method,7,(n_method-1)*7+2);
imagesc(Scope.us,Scope.vs,wrapTo2Pi(phsEx_gt));
Phsx_all{n_method} = wrapTo2Pi(phsEx_gt);
axis image xy off;colormap(gca,cmap_P);caxis([0 2*pi]);
colorbar('Ticks',[0 2*pi],'TickLabels', ...
    {sprintf('%.2e',0),sprintf('%.2e',2*pi)});
title('PhaseX_{GT}');

% plot the phase of Y polarization
subplot(n_method,7,(n_method-1)*7+4);
imagesc(Scope.us,Scope.vs,wrapTo2Pi(phsEy_gt));
Phsy_all{n_method} = wrapTo2Pi(phsEy_gt);
axis image xy off;colormap(gca,cmap_P);caxis([0 2*pi]);
colorbar('Ticks',[0 2*pi],'TickLabels', ...
    {sprintf('%.2e',0),sprintf('%.2e',2*pi)});
title('PhaseY_{GT}');

% plot the phase of Z polarization
subplot(n_method,7,(n_method-1)*7+6);
imagesc(Scope.us,Scope.vs,wrapTo2Pi(phsEz_gt));
Phsz_all{n_method} = wrapTo2Pi(phsEz_gt);
axis image xy off;colormap(gca,cmap_P);caxis([0 2*pi]);
colorbar('Ticks',[0 2*pi],'TickLabels', ...
    {sprintf('%.2e',0),sprintf('%.2e',2*pi)});
title('PhaseZ_{GT}');

for ii = 1:(n_method-1)
    [Ex,Ey,Ez] = VectorialDiffraction(Obj,Beam,Scope,FocusingMethod{ii},number_u,number_v);

    Ex_all{ii} = Ex/max(abs(Ex(:)));
    Ey_all{ii} = Ey/max(abs(Ey(:)));
    Ez_all{ii} = Ez/max(abs(Ez(:)));

    phsEx = wrapTo2Pi(angle(Ex)-phs_tilt);
    phsEy = wrapTo2Pi(angle(Ey)-phs_tilt);
    phsEz = wrapTo2Pi(angle(Ez)-phs_tilt);

    Ix = abs(Ex).^2;
    Iy = abs(Ey).^2;
    Iz = abs(Ez).^2;
    I = Ix + Iy + Iz;

    Ix = Ix./max(I(:));
    Iy = Iy./max(I(:));
    Iz = Iz./max(I(:));
    I = I./max(I(:));

    figure(fig_PSF);
    % plot the intensity of X polarization
    subplot(n_method,7,(ii-1)*7+1);
    imagesc(Scope.us,Scope.vs,Ix);
    axis image xy off;colormap(gca,cmap_E);
    if strcmp(FocusingMethod{ii},'smm_ap')
        title('IntensityX_{Ours}');
    elseif strcmp(FocusingMethod{ii},'czt_ap')
        title('IntensityX_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end
    Ix_all{ii} = Ix;
    caxis([0,Ixmax]);
    colorbar('Ticks',[0,Ixmax],'TickLabels',...
        {sprintf('%.2e',0),sprintf('%.2e',Ixmax)});

    % plot the intensity of Y polarization
    subplot(n_method,7,(ii-1)*7+3);
    imagesc(Scope.us,Scope.vs,Iy);
    axis image xy off;colormap(gca,cmap_E);
    if strcmp(FocusingMethod{ii},'smm_ap')
        title('IntensityY_{Ours}');
    elseif strcmp(FocusingMethod{ii},'czt_ap')
        title('IntensityY_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end
    Iy_all{ii} = Iy;
    caxis([0,Iymax]);
    colorbar('Ticks',[0,Iymax],'TickLabels',...
        {sprintf('%.2e',0),sprintf('%.2e',Iymax)});

    % plot the intensity of Z polarization
    subplot(n_method,7,(ii-1)*7+5);
    imagesc(Scope.us,Scope.vs,Iz);
    axis image xy off;colormap(gca,cmap_E);
    if strcmp(FocusingMethod{ii},'smm_ap')
        title('IntensityZ_{Ours}');
    elseif strcmp(FocusingMethod{ii},'czt_ap')
        title('IntensityZ_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end
    Iz_all{ii} = Iz;
    caxis([0, Izmax]);
    colorbar('Ticks',[0,Izmax],'TickLabels',...
        {sprintf('%.2e',0),sprintf('%.2e',Izmax)});

    % plot the intensity
    subplot(n_method,7,ii*7);
    imagesc(Scope.us,Scope.vs,I);
    axis image xy off;colormap(gca,cmap_E);
    if strcmp(FocusingMethod{ii},'smm_ap')
        title('Intensity_{Ours}');
    elseif strcmp(FocusingMethod{ii},'czt_ap')
        title('Intensity_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end
    I_all{ii} = I;
    caxis([0, Imax]);
    colorbar('Ticks',[0,Imax],'TickLabels',...
        {sprintf('%.2e',0),sprintf('%.2e',Imax)});

    % plot the phase of X polarization
    subplot(n_method,7,(ii-1)*7+2);
    imagesc(Scope.us,Scope.vs,wrapTo2Pi(phsEx));
    Phsx_all{ii} = wrapTo2Pi(phsEx);
    axis image xy off;colormap(gca,cmap_P);caxis([0 2*pi]);
    colorbar('Ticks',[0 2*pi],'TickLabels', ...
        {sprintf('%.2e',0),sprintf('%.2e',2*pi)});
    if strcmp(FocusingMethod{ii},'smm_ap')
        title('PhaseX_{Ours}');
    elseif strcmp(FocusingMethod{ii},'czt_ap')
        title('PhaseX_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end

    % plot the phase of Y polarization
    subplot(n_method,7,(ii-1)*7+4);
    imagesc(Scope.us,Scope.vs,wrapTo2Pi(phsEy));
    Phsy_all{ii} = wrapTo2Pi(phsEy);
    axis image xy off;colormap(gca,cmap_P);caxis([0 2*pi]);
    colorbar('Ticks',[0 2*pi],'TickLabels', ...
        {sprintf('%.2e',0),sprintf('%.2e',2*pi)});
    if strcmp(FocusingMethod{ii},'smm_ap')
        title('PhaseY_{Ours}');
    elseif strcmp(FocusingMethod{ii},'czt_ap')
        title('PhaseY_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end

    % plot the phase of Z polarization
    subplot(n_method,7,(ii-1)*7+6);
    imagesc(Scope.us,Scope.vs,wrapTo2Pi(phsEz));
    Phsz_all{ii} = wrapTo2Pi(phsEz);
    axis image xy off;colormap(gca,cmap_P);caxis([0 2*pi]);
    colorbar('Ticks',[0 2*pi],'TickLabels', ...
        {sprintf('%.2e',0),sprintf('%.2e',2*pi)});
    if strcmp(FocusingMethod{ii},'smm_ap')
        title('PhaseZ_{Ours}');
    elseif strcmp(FocusingMethod{ii},'czt_ap')
        title('PhaseZ_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end
end

%% compare accuracy
fig_cmp = figure('Name','Comparison of Accuracy');

devx = cell(1,n_method-1);
devy = cell(1,n_method-1);
devz = cell(1,n_method-1);
dev = cell(1,n_method-1);

for jj = 1:(n_method-1)
    figure(fig_cmp);
    devx{jj} = abs(Ex_all{jj} - Ex_all{n_method});
    devy{jj} = abs(Ey_all{jj} - Ey_all{n_method});
    devz{jj} = abs(Ez_all{jj} - Ez_all{n_method});
    dev{jj} = abs(I_all{jj} - I_all{n_method});
    min_devx = min(devx{jj}(:)); max_devx = max(devx{jj}(:));
    min_devy = min(devy{jj}(:)); max_devy = max(devy{jj}(:));
    min_devz = min(devz{jj}(:)); max_devz = max(devz{jj}(:));
    max_dev = max(dev{jj}(:));

    subplot(n_method-1,4,(jj-1)*4+1);
    imagesc(devx{jj});
    caxis([min_devx, max_devx]);
    axis image xy off;
    colorbar('Ticks',[min_devx,max_devx],'TickLabels', ...
        {sprintf('%.1e',min_devx),sprintf('%.1e',max_devx)});
    colormap(gca,cmap_E);
    if strcmp(FocusingMethod{jj},'smm_ap')
        title('DeviationX_{Ours}');
    elseif strcmp(FocusingMethod{jj},'czt_ap')
        title('DeviationX_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end

    subplot(n_method-1,4,(jj-1)*4+2);
    imagesc(devy{jj});
    caxis([min_devy, max_devy]);
    axis image xy off;
    colorbar('Ticks',[min_devy,max_devy],'TickLabels', ...
        {sprintf('%.1e',min_devy),sprintf('%.1e',max_devy)});
    colormap(gca,cmap_E);
    if strcmp(FocusingMethod{jj},'smm_ap')
        title('DeviationY_{Ours}');
    elseif strcmp(FocusingMethod{jj},'czt_ap')
        title('DeviationY_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end

    subplot(n_method-1,4,(jj-1)*4+3);
    imagesc(devz{jj});
    caxis([min_devz, max_devz]);
    axis image xy off;
    colorbar('Ticks',[min_devz,max_devz],'TickLabels', ...
        {sprintf('%.1e',min_devz),sprintf('%.1e',max_devz)});
    colormap(gca,cmap_E);
    if strcmp(FocusingMethod{jj},'smm_ap')
        title('DeviationZ_{Ours}');
    elseif strcmp(FocusingMethod{jj},'czt_ap')
        title('DeviationZ_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end

    subplot(n_method-1,4,jj*4);
    imagesc(dev{jj});
    caxis([0, max_dev]);
    axis image xy off;
    colorbar('Ticks',[0,max_dev],'TickLabels', ...
        {sprintf('%.1e',0),sprintf('%.1e',max_dev)});
    colormap(gca,cmap_E);
    if strcmp(FocusingMethod{jj},'smm_ap')
        title('Deviation_{Ours}');
    elseif strcmp(FocusingMethod{jj},'czt_ap')
        title('Deviation_{Ctrl}');
    else
        error('Invalid diffraction method!');
    end
end

%% evaluate errors and time consumption
sigmax = cell(1,2);
sigmay = cell(1,2);
sigmaz = cell(1,2);

sigmax{1} = sum(abs(Ex_all{1} - Ex_all{3}).^2,'all')/sum(abs(Ex_all{3}).^2,'all');
sigmax{2} = sum(abs(Ex_all{2} - Ex_all{3}).^2,'all')/sum(abs(Ex_all{3}).^2,'all');
sigmay{1} = sum(abs(Ey_all{1} - Ey_all{3}).^2,'all')/sum(abs(Ey_all{3}).^2,'all');
sigmay{2} = sum(abs(Ey_all{2} - Ey_all{3}).^2,'all')/sum(abs(Ey_all{3}).^2,'all');
sigmaz{1} = sum(abs(Ez_all{1} - Ez_all{3}).^2,'all')/sum(abs(Ez_all{3}).^2,'all');
sigmaz{2} = sum(abs(Ez_all{2} - Ez_all{3}).^2,'all')/sum(abs(Ez_all{3}).^2,'all');

sigma_mtp = (sigmax{1}+sigmay{1}+sigmaz{1})/3
sigma_czt = (sigmax{2}+sigmay{2}+sigmaz{2})/3

time_MTP = timeit(@() VectorialDiffraction_ASR_AP(Obj,Beam,Scope,number_u,number_v))
time_CZT = timeit(@() VectorialDiffraction_CZT_AP(Obj,Beam,Scope))