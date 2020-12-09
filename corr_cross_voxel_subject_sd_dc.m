%% correlation main results
clear
dirSD = '/brain/guixue/Jintao/HCP_S1200/VoxelSD/3dTproject_fSD';
dirDC = '/brain/guixue/Jintao/HCP_S1200/VoxelDegree/Based3TprojectMask';
dirResi = '/brain/guixue/Jintao/HCP_S1200/VoxelResi/Based3TprojectMask';
dirResult = '/brain/guixue/Jintao/HCP_S1200/results/Based3TprojectMask';

[Vmsk,Ymsk,~] = gretna_read_image('/brain/guixue/Jintao/HCP_S1200/mask_3dTproject_GM.nii');
Ymsk(isnan(Ymsk)) = 0;
idx = find(Ymsk);

subj = load('/brain/guixue/Jintao/HCP_S1200/HCP927.txt');

out_range = {'full','long'};
filter_cova = {'mean','r1','r2'};

n = 0;
for session = {'Mean12','REST1','REST2'}
    n = n+1;
    Cova = load(['/brain/guixue/Jintao/HCP_S1200/cova_' filter_cova{n} '_927.txt']);
    N_cova = size(Cova,2);
    
    load(fullfile(dirSD,[session{1} '_whole_brain_SD.mat']),'Y_Divided_fSD')
    zY_Divided_fSD = zscore(Y_Divided_fSD);
    Num_subj = size(zY_Divided_fSD,2);
    
    for thres = {'005','007','010','015'}
        m = 0;
        for Range = {'Exloc','Long'}
            m = m+1;
            load(fullfile(dirDC,[session{1} '_Mean_whole_brain_' ...
                Range{1} 'DC_' thres{1} '.mat']))
            zYout = zscore(Yout);
            
            % cross-subject correlation
            r = arrayfun(@(x) partialcorr(zY_Divided_fSD(x,:)',zYout(x,:)',...
                Cova,'Rows','pairwise'), 1:length(idx));
            DOF = Num_subj - 2 - N_cova;
            ztmp = gretna_fishertrans(r);
            R = zeros(size(Ymsk));
            R(idx) = r;
            z = zeros(size(Ymsk));
            
            z(idx) = ztmp;
            VoxelSize = repmat(Vmsk.mat(2,2), 1, 3);
            [dLh, ~, FWHM, ~] = y_Smoothest(z, Ymsk, DOF, VoxelSize);
            
            Vout = Vmsk;
            Vout.dt(1) = 16;
            
            Vout.descrip = sprintf('DPABI{R_[%d]}{dLh_%.6f}{FWHMx_%.6fFWHMy_%.6fFWHMz_%.6fmm}', ...
                DOF, dLh, FWHM(1), FWHM(2), FWHM(3));
            
            cd(fullfile(dirResult, session{1}, 'corr_fSD_DC'))
            
            Vout.fname = [session{1} '_partialcorr_Divided_fSD_' Range{1} ...
                'DC_' thres{1} '.nii'];
            spm_write_vol(Vout, R);
            
            [~,name,~] = fileparts(Vout.fname);
            Vout.fname = [name '_FsiherZ.nii'];
            spm_write_vol(Vout, z);
            
            % cross-voxel correlation
            [r,p] = corr(zY_Divided_fSD,zYout);
            R_all = diag(r); P_all = diag(p);
            
            lm = arrayfun(@(x) fitlm(zY_Divided_fSD(:,x),zYout(:,x)), 1:Num_subj,...
                'UniformOutput', false);
            
            for i = 1:Num_subj
                T_all(i,1) = lm{1,i}.Coefficients.tStat(2);
                
                Resid = zeros(size(Ymsk));
                Resid(idx) = lm{1,i}.Residuals.Raw;
                cd(fullfile(dirResi, session{1}, thres{1},...
                    out_range{m}, 'Divided_Zscore_fSD_DC'))
                
                if exist('ImageFiles', 'dir') == 0
                    mkdir ImageFiles
                end
                cd ImageFiles
                
                Vout = Vmsk;
                Vout.fname = ['Residuals_' session{1} '_' Range{1} 'DC_' ...
                    thres{1} '_' num2str(subj(i)) '.nii'];
                Vout.dt(1) = 16;
                spm_write_vol(Vout, Resid);
            end
            save(fullfile(dirResult, session{1}, 'corr_fSD_DC',[session{1} ...
                '_corr_Divided_fSD_' Range{1} 'DC_' thres{1} '.mat']),...
                'R_all','P_all','T_all')
        end
    end
end
