function s02_voxel_degree_SD(subj_path, session)

%%% calculating degree centrality for scrubbing data

datetime('now')
restoredefaultpath;
addpath(genpath('/brain/guixue/Jintao/tools/spm12'));

% Iuput
base_dir = '/brain/guixue/Jintao/HCP_S1200';
Data_mask = fullfile(base_dir, 'GrayMask_02_2mm.nii');
Output_path = fullfile(base_dir,'VoxelDegreeFull', 'DifferentThresh_3dTproject', ...
    ['Scrubbing_' session]);
R_thr = 0.1;
Dis = 75;
LocDis = 20;
Nloop = 150;

subj_list = load(subj_path);

% mask
Vmask = spm_vol(Data_mask);
[Ymask, XYZ] = spm_read_vols(Vmask);
Ymask(isnan(Ymask)) = 0;

Index = find(Ymask);
XYZ = XYZ(:,Index)';
row_XYZ = size(XYZ,2);
[I, J, K] = ind2sub(size(Ymask),Index);

% runing: loop for subjects
for i = 1:length(subj_list)
    tic
    subj_id = num2str(subj_list(i));
    fprintf('Calculating Degree and SD for %s...\n', subj_id)
    
    Vin = spm_vol(fullfile(base_dir, 'RawData', subj_id, ...
        ['scr_bp_WGSR_3dTproject_' session '_LR_clean.nii']));
    Ydata = spm_get_data(Vin,[I J K]');
    Ydata = single(Ydata);
    
    Vin = spm_vol(fullfile(base_dir, 'RawData', subj_id, ...
        ['scr_hp_WGSR_3dTproject_' session '_LR_clean.nii']));
    Yff = spm_get_data(Vin,[I J K]');
    Yff = single(Yff);

    fprintf('Computing fractional SD...\n')
    lSD = nanstd(Ydata, 0, 1);
    fSD = lSD./nanstd(Yff, 0, 1);
    mfSD = fSD./nanmean(fSD(:));
    
    lSD_out = zeros(size(Ymask));
    lSD_out(Index) = lSD;
    fSD_out = zeros(size(Ymask));
    fSD_out(Index) = fSD;

    mfSD_out = zeros(size(Ymask));
    mfSD_out(Index) = mfSD;
    
    fprintf('Computing Degree ...\n')
    
    [m,n] = size(Ydata);
    Ydata = Ydata - repmat(mean(Ydata), m, 1);
    Ydata = Ydata./repmat(std(Ydata, 0, 1), m, 1);
    
    R_pos_wei_exloc = repmat({zeros(size(Ymask))}, 1, length(R_thr));
    R_pos_wei_long = repmat({zeros(size(Ymask))}, 1, length(R_thr));
    
    FisherZ_pos_wei_exloc = repmat({zeros(size(Ymask))}, 1, length(R_thr));
    FisherZ_pos_wei_long = repmat({zeros(size(Ymask))}, 1, length(R_thr));
    
    a = fix(n/Nloop);
    b = n - a*(Nloop-1);
    
    x = mat2cell(Ydata, m, [repmat(a, 1, Nloop-1), b]);
    xyz = mat2cell(XYZ, [repmat(a, 1, Nloop-1), b], row_XYZ);
    idx = mat2cell(Index, [repmat(a, 1, Nloop-1), b], 1);
    
    for ii = 1:Nloop
        fprintf('%d...', ii);
        D = pdist2(xyz{ii},XYZ);
        idx_exloc = D > LocDis;
        idx_long = D >= Dis;
        
        r = x{ii}' * Ydata / (m - 1);
        r(r>=1) = nan;
        
        for jj = 1:length(R_thr)
            r_tmp = r;
            r_tmp(r_tmp < R_thr(jj)) = NaN;
            
            r_exloc = r_tmp.*idx_exloc;
            r_long = r_tmp.*idx_long;
            
            z_exloc = atanh(r_exloc);
            z_long = atanh(r_long);
            
            R_pos_wei_exloc{jj}(idx{ii})       = nansum(r_exloc,2);
            R_pos_wei_long{jj}(idx{ii})        = nansum(r_long,2);
            
            FisherZ_pos_wei_exloc{jj}(idx{ii}) = nansum(z_exloc,2);
            FisherZ_pos_wei_long{jj}(idx{ii})  = nansum(z_long,2);
        end
    end
    fprintf('\n\t Writing voxel-wise degree and SD for %s...\n', subj_id);
    cd (Output_path)
    
    if ~exist(subj_id, 'dir')
        mkdir(subj_id)
    end
    
    cd(subj_id)
    Vout = Vin(1);
    Vout.dt(1) = 16;
    
    Vout.fname = [session '_rSD_WGSR_' subj_id '.nii'];
    Vout = spm_write_vol(Vout, lSD_out);
    
    Vout.fname = [session '_fSD_WGSR_' subj_id '.nii'];
    Vout = spm_write_vol(Vout, fSD_out);
    
    Vout.fname = [session '_mfSD_WGSR_' subj_id '.nii'];
    Vout = spm_write_vol(Vout, mfSD_out);
    
    for k = 1:length(R_thr)
        R_thr_name = R_thr(k)*100;
        
        Vout.fname = ['degree_pos_wei_exloc_' ...
            num2str(R_thr_name) '_' subj_id '.nii'];
        Vout = spm_write_vol(Vout, R_pos_wei_exloc{k});
        
        Vout.fname = ['degree_pos_wei_long_' ...
            num2str(R_thr_name) '_' subj_id '.nii'];
        Vout = spm_write_vol(Vout, R_pos_wei_long{k});
        
        Vout.fname = ['degree_pos_wei_exloc_' ...
            num2str(R_thr_name) '_' subj_id '_FisherZ.nii'];
        Vout = spm_write_vol(Vout, FisherZ_pos_wei_exloc{k});
        
        Vout.fname = ['degree_pos_wei_long_' ...
            num2str(R_thr_name) '_' subj_id '_FisherZ.nii'];
        Vout = spm_write_vol(Vout, FisherZ_pos_wei_long{k});
    end
    toc
    datetime('now')
end
return
