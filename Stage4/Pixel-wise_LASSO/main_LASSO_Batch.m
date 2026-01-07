close all
clear
clc

%% Step 1: Load basis spectra

% Calibrate spectral window
% Find positions of two peaks (r1, r2) and frame number (f1, f2)
% Use standard chemicals for calibration.
% Users should change accordingly.

r1 = 1655;
r2 = 1745;
f1 = 20;
f2 = 58;
Nz = 80; % Set number of frames in the hyperspectral image

rpf =(r2-r1)/(f2-f1);
Raman_shift = linspace(r1-f1*rpf, r2+(Nz-f2+1)*rpf, Nz);

MatRefFile = 'pure_chemicals_1650.mat'; % Set basis spectra file name
S = load(MatRefFile);  % <-- load into struct instead of workspace

% Canonical component order (controls ref column order + output order)
components = {'BSA','TAG','CHL'};
k = numel(components);

% Mapping: canonical name -> variable name inside MAT file
varmap = struct( ...
    'BSA','BSA', ...
    'TAG','TAG', ...
    'CHL','CHL');

% Do not change the remaining lines in this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build normalized reference matrix ref (Nz x k)
ref = zeros(numel(Raman_shift), k);
for ii = 1:k
    cname = components{ii};
    vname = varmap.(cname);

    assert(isfield(S, vname), 'MAT file missing variable "%s" (mapped from %s).', vname, cname);

    x = S.(vname);
    x = x(:); % force column

    assert(numel(x) == numel(Raman_shift), ...
        'Variable "%s" has length %d but Raman_shift has %d.', vname, numel(x), numel(Raman_shift));

    ref(:,ii) = (x - min(x)) ./ (max(x) - min(x) + eps);
end

% Plot normalized references with offsets (same look as your original)
figure (1);
for ii = 1:k
    plot(Raman_shift, ref(:,ii) + (ii-1), 'LineWidth', 1); hold on
end
hold off
set(gca,'XMinorTick','on','YMinorTick','on');
xlabel('Raman shift (cm^{-1})'); ylabel('Int (a.u.)')
pause (4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 2: Set image batch directory and load file names

% Set input file type
%
% .txt --> hyperspectral text image in 2D montage
% .tif --> 3D tif image

Dir      = 'Data/';
filetype = '.tif';
Filelist = dir([Dir '*' filetype]); % Set file pattern

cells    = struct2cell(Filelist);
sortvals = cells(6,:);
mat      = cell2mat(sortvals);
[~,ix]   = sort(mat);

Filelist = Filelist(ix);

file_num = size(Filelist,1);
clear ix cells sortvals mat

%% Begin the loop for all images
for jj = 1:file_num

    %%%%%%%%%%% Loop Step 1: Load image %%%%%%%%%
    
    close all
    % Load specific image
    filename = Filelist(jj).name;

    % Do not change the lines in between 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if endsWith([filename], '.tif', 'IgnoreCase', true)
        info = imfinfo([Dir filename]);
        Nx = info(1).Height;
        Ny = info(1).Width;
        y =  zeros(info(1).Height, info(1).Width, Nz);
        for Nk = 1:Nz
            y(:,:,Nk) = imread([Dir filename], Nk);
        end
    elseif endsWith([filename], '.txt', 'IgnoreCase', true)
        y_montage = load([Dir filename]);
        Nx        = size(y_montage,1);
        Ny        = size(y_montage,2)/Nz;
        y         = permute(reshape(y_montage,[Nx,Ny,Nz]),[1,2,3]);
    else
        error('Unsupported file type: %s', [filename]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do not change the lines in between 

    
    %%%%%%%%%%% Loop Step 2: Subtract background %%%%%%%%%
    % Set the threshold value using the intensity histogram
    avg_threshold = -0.1;

    % Do not change the lines in between
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_sum = squeeze(mean(y,3));
    figure(2); histogram(y_sum);

    BGmask = zeros(size(y_sum));
    BG = find(y_sum < avg_threshold);
    BGmask(BG) = 1;
    figure(3); imagesc(BGmask); axis image off

    BG_spectrum = zeros(1,Nz);
    for i = 1:Nz
        y_temp = y(:,:,i);
        BG_spectrum(i) = mean(y_temp(BGmask == 1));
    end

    figure(4); plot(Raman_shift, BG_spectrum, 'Linewidth', 1.5);
    xlabel('Raman shift (cm^{-1})'); ylabel('Int (a.u.)')
    xlim([2820, 3030])
    set(gca,'FontSize',18,'LineWidth',2)

    BG_spectrum_n = (BG_spectrum - min(BG_spectrum)) ./ (max(BG_spectrum) - min(BG_spectrum) + eps);

    y_sub = zeros(size(y), 'like', y);
    for i = 1:Nz
        y_sub(:,:,i) = y(:,:,i) - BG_spectrum(i);
    end

    C    = zeros(Nx,Ny,k);  % Allocate empty concentration maps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do not change the lines in between

    
    %%%%%%%%%%% Loop Step 3: Run pixel-wise LASSO %%%%%%%%%

    % Set the sparsity level for all channels, same for all by default, change
    % if the initial outcome is unsatisfatory
    l = [5e-2, 5e-2, 5e-2];

    % Controls which channel to update, update all by default.
    update = true(1,k);

    % These two parameters do not need to change in most cases
    a    = 1;            % Controls convergence speed of ADMM
    iter = 5;            % Number of iterations for ADMM

    for ll = 1:k
        if update(ll)
            L  = l(ll);
            CC = nneg_lasso(y_sub, ref, L, a, iter);
            C(:,:,ll) = CC(:,:,ll);
        end
    end

    % Display concentration maps and examine image quality

    % Use the first channel (C(:,:,1)) to set up the display threshold.
    % Can use other channels if necessary
    disp_min = prctile(reshape(C(:,:,1), [Nx*Ny,1]), 0.3);
    disp_max = prctile(reshape(C(:,:,1), [Nx*Ny,1]), 99.7);

    figure(5);
    clims = [disp_min disp_max];

    for ii = 1:k
        subplot(2,3,ii); % Change subplot size if needed
        imagesc(C(:,:,ii), clims);
        colormap bone; axis off; axis square
    end
    drawnow
    pause(4);

    %%%%%%%%%%% Loop Step 4: Output image %%%%%%%%%

    % Set output filepath
    opt_filepath = 'Output_maps/';
    if ~exist(opt_filepath, 'dir')
        mkdir(opt_filepath);
    end

    % Do not change the remaining lines in this section
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    outExt = '.txt';
    for ii = 1:k
        cname = components{ii};
        outFile = fullfile(opt_filepath, [filename, '_', cname, outExt]);
        writematrix(C(:,:,ii), outFile, 'Delimiter', '\t');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end