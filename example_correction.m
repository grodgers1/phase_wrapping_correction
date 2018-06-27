%% Demo for phase wrapping correction algorithm
clear all

addpath ./functions/
addpath ./data/

%% Load sample data
wrap_sino = load('wrap_sino.mat'); wrap_sino = wrap_sino.wrap_sino;
abs_sino = load('abs_sino.mat'); abs_sino = abs_sino.abs_sino;

[sx, n_proj] = size(wrap_sino);
%% Measurement parameters
p2 = 7e-06; % [m]
d_t = 0.8; % [m]
pixsize = 2*2.3e-07; % pixel size of camera
E = 20000; % energy [eV]

angle_step = 0.15; % angle step between projections
angles = 0:angle_step:180; % angle steps of the projections

%% Reconstruct wrapped phase sinogram
wrap_slice = iradondpc(wrap_sino, angles, 'linear', 'Ram-Lak',sx);
wrap_slice = -(p2/(2*pi*d_t))*wrap_slice;
figure, imagesc(wrap_slice), axis equal, colormap gray

%% Perform correction with a range of model delta values
delta_range = (4.75:0.025:5.75)*1e-07;
buffer = 25;
edge_method = 1;
edge_bump = 1;
smoothedges = 1;
smoothfactor = 15;
usefits = 1;
bump = 20;

slice = zeros(sx,sx,length(delta_range));
for d = 1:length(delta_range)
    tic
    fprintf(['Working on delta value ' num2str(d) ' of ' num2str(length(delta_range)) '\n'])
    delta = delta_range(d);
	corr_sino = no_tank_correction(wrap_sino, abs_sino, d_t, p2, delta, buffer, bump, smoothedges, smoothfactor, usefits);
	slice(:,:,d) = iradondpc(corr_sino, angles, 'linear', 'Ram-Lak', sx);
	toc
end
slice = -(p2/(2*pi*d_t))*slice; % making units of delta

%% Find mean and std of regions of interest to determine best delta:
% % ROIs selected
% % Paraffin rois:
roi1_pa = [640 680 790 820]; % [x1, x2, y1, y2]
roi2_pa = [1050 1100 555 590]; % [x1, x2, y1, y2]
roi3_pa = [810 865 670 695]; % [x1, x2, y1, y2] 
roi4_pa = [1140 1195 475 525]; % [x1, x2, y1, y2] 
% % ML rois:
roi1_ml = [610 685 940 980];
roi2_ml = [500 560 770 810];
roi3_ml = [835 870 500 550];
% % White matter:
roi1_wm = [940 1060 1420 1460];
roi2_wm = [1290 1365 1260 1290];
% % GL rois:
roi1_gl = [625 685 530 590]; 
roi2_gl = [600 670 1170 1220];
% % visualize these:
figure, imagesc(slice(:,:,1), [4.5e-07 5.5e-07]), colormap gray, axis equal
hold on
rectangle('Position', [roi1_gl(1) roi1_gl(3) roi1_gl(2)-roi1_gl(1) roi1_gl(4)-roi1_gl(3)], 'EdgeColor', 'r')
rectangle('Position', [roi2_gl(1) roi2_gl(3) roi2_gl(2)-roi2_gl(1) roi2_gl(4)-roi2_gl(3)], 'EdgeColor', 'r')
rectangle('Position', [roi1_ml(1) roi1_ml(3) roi1_ml(2)-roi1_ml(1) roi1_ml(4)-roi1_ml(3)], 'EdgeColor', 'b')
rectangle('Position', [roi2_ml(1) roi2_ml(3) roi2_ml(2)-roi2_ml(1) roi2_ml(4)-roi2_ml(3)], 'EdgeColor', 'b')
rectangle('Position', [roi3_ml(1) roi3_ml(3) roi3_ml(2)-roi3_ml(1) roi3_ml(4)-roi3_ml(3)], 'EdgeColor', 'b')
rectangle('Position', [roi1_pa(1) roi1_pa(3) roi1_pa(2)-roi1_pa(1) roi1_pa(4)-roi1_pa(3)], 'EdgeColor', 'g')
rectangle('Position', [roi2_pa(1) roi2_pa(3) roi2_pa(2)-roi2_pa(1) roi2_pa(4)-roi2_pa(3)], 'EdgeColor', 'g')
rectangle('Position', [roi3_pa(1) roi3_pa(3) roi3_pa(2)-roi3_pa(1) roi3_pa(4)-roi3_pa(3)], 'EdgeColor', 'g')
rectangle('Position', [roi4_pa(1) roi4_pa(3) roi4_pa(2)-roi4_pa(1) roi4_pa(4)-roi4_pa(3)], 'EdgeColor', 'g')
rectangle('Position', [roi1_wm(1) roi1_wm(3) roi1_wm(2)-roi1_wm(1) roi1_wm(4)-roi1_wm(3)], 'EdgeColor', 'y')
rectangle('Position', [roi2_wm(1) roi2_wm(3) roi2_wm(2)-roi2_wm(1) roi2_wm(4)-roi2_wm(3)], 'EdgeColor', 'y')

m_gl_all = zeros(size(delta_range));
m_ml_all = zeros(size(delta_range));
m_pa_all = zeros(size(delta_range));
m_wm_all = zeros(size(delta_range));
std_gl_all = zeros(size(delta_range));
std_ml_all = zeros(size(delta_range));
std_pa_all = zeros(size(delta_range));
std_wm_all = zeros(size(delta_range));
for d = 1:length(delta_range)
    vol1_gl = slice(roi1_gl(3):roi1_gl(4), roi1_gl(1):roi1_gl(2), d);
    vol2_gl = slice(roi2_gl(3):roi2_gl(4), roi2_gl(1):roi2_gl(2), d);
    vol1_ml = slice(roi1_ml(3):roi1_ml(4), roi1_ml(1):roi1_ml(2), d);
    vol2_ml = slice(roi2_ml(3):roi2_ml(4), roi2_ml(1):roi2_ml(2), d);
    vol3_ml = slice(roi3_ml(3):roi3_ml(4), roi3_ml(1):roi3_ml(2), d);
    vol1_pa = slice(roi1_pa(3):roi1_pa(4), roi1_pa(1):roi1_pa(2), d);
    vol2_pa = slice(roi2_pa(3):roi2_pa(4), roi2_pa(1):roi2_pa(2), d);
    vol3_pa = slice(roi3_pa(3):roi3_pa(4), roi3_pa(1):roi3_pa(2), d);
    vol4_pa = slice(roi4_pa(3):roi4_pa(4), roi4_pa(1):roi4_pa(2), d);
    vol1_wm = slice(roi1_wm(3):roi1_wm(4), roi1_wm(1):roi1_wm(2), d);
    vol2_wm = slice(roi2_wm(3):roi2_wm(4), roi2_wm(1):roi2_wm(2), d);

    m_gl_all(d) = mean([vol1_gl(:)', vol2_gl(:)']);
    m_ml_all(d) = mean([vol1_ml(:)', vol2_ml(:)',vol3_ml(:)']);
    m_pa_all(d) = mean([vol1_pa(:)', vol2_pa(:)',vol3_pa(:)',vol4_pa(:)']);
    m_wm_all(d) = mean([vol1_wm(:)', vol2_wm(:)']);
    std_gl_all(d) = std([vol1_gl(:)', vol2_gl(:)']);
    std_ml_all(d) = std([vol1_ml(:)', vol2_ml(:)',vol3_ml(:)']);
    std_pa_all(d) = std([vol1_pa(:)', vol2_pa(:)',vol3_pa(:)',vol4_pa(:)']);
    std_wm_all(d) = std([vol1_wm(:)', vol2_wm(:)']);
end

% % Find the best delta for each material
d_gl_all = delta_range(std_gl_all == min(std_gl_all));
d_ml_all = delta_range(std_ml_all == min(std_ml_all));
d_pa_all = delta_range(std_pa_all == min(std_pa_all));
d_wm_all = delta_range(std_wm_all == min(std_wm_all));
d_vec = [d_gl_all, d_ml_all, d_pa_all, d_wm_all];

figure, plot(delta_range, std_gl_all, 'r')
hold on, plot(delta_range, std_ml_all, 'g')
plot(delta_range, std_pa_all, 'b')
plot(delta_range, std_wm_all, 'k')
plot(d_gl_all, std_gl_all((std_gl_all == min(std_gl_all))), 'r*')
plot(d_ml_all, std_ml_all((std_ml_all == min(std_ml_all))), 'g*')
plot(d_pa_all, std_pa_all((std_pa_all == min(std_pa_all))), 'b*')
plot(d_wm_all, std_wm_all((std_wm_all == min(std_wm_all))), 'k*')

% % choose a delta value
model_delta = mean(d_vec);

%% Perform correction with best model delta and reconstruct
corr_sino = no_tank_correction(wrap_sino, abs_sino, d_t, p2, model_delta, buffer, bump, smoothedges, smoothfactor, usefits);
corr_slice = iradondpc(corr_sino, angles, 'linear', 'Ram-Lak', sx);
corr_slice = -(p2/(2*pi*d_t))*corr_slice;

%% Compare corrected to uncorrected
figure, imagesc(corr_slice), colormap gray, axis equal
figure, imagesc(wrap_slice), colormap gray, axis equal


