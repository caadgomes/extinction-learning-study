
clearvars;

warning('off', 'MATLAB:rmpath:DirNotFound');
% remove PsPM from path because it conflicts with SPM
rmpath(genpath('/usr/local/MATLAB/TOOLBOXES/PsPM_v6.0.0/'));

% Set-Up and Configuration

% Initialise SPM
%--------------------------------------------------------------------------
addpath('/usr/local/MATLAB/TOOLBOXES/spm12')
spm('Defaults','fMRI');


% Directory containing the data
%--------------------------------------------------------------------------
AG = 'A03';
exp = '3T';
base_dir = '/home/data/';
bids_dir = fullfile(base_dir, AG, exp);
rawdata_dir = fullfile(bids_dir, 'rawdata');
derivatives_dir = fullfile(bids_dir, 'derivatives');
fmriprep_dir = fullfile(derivatives_dir, 'fmriprep');
denoised_dir = fullfile(derivatives_dir, 'func_denoised');
ROIs_dir = fullfile(derivatives_dir, 'ROIs');
wf_dir = fullfile(bids_dir, 'workflows');
mri_space = 'T1w';
suit_template = fullfile(base_dir,'templates','atl-Anatom_space-SUIT_dseg.nii');
suit_lut = fullfile(base_dir,'templates','atl-Anatom.tsv');
idxs = 29:34;
temp_folder = '/usr/data/';

desc_atlas = '5rois';

opts = detectImportOptions(fullfile(rawdata_dir,'participants.tsv'), 'FileType','delimitedtext', 'Delimiter','\t');
opts = setvartype(opts, 'SFB1280_code', 'char');
subs_table = readtable(fullfile(rawdata_dir,'participants.tsv'), opts);
% subs_table = subs_table(1,:);

name_wf = 'suit';
wf_dir = fullfile(wf_dir, 'wf_cerebellum_ROIs', name_wf);
if ~exist(wf_dir,'dir'); mkdir(wf_dir); end

% Data Extraction and Aggregation

for sub_idx = 1:height(subs_table)

    ptc = subs_table(sub_idx,:).participant_id{:};

    fprintf('Processing subject: %s',ptc)

    argsnode = struct('p_wf_dir',fullfile(wf_dir, ptc), 'rerun',0, 'name_node','');

    p_cerebel_dir = fullfile(ROIs_dir, ptc, 'cerebellum');
    if ~exist(p_cerebel_dir,'dir'); mkdir(p_cerebel_dir); end

    p_wf_dir = fullfile(wf_dir, ptc);
    if ~exist(p_wf_dir,'dir'); mkdir(p_wf_dir); end

	% Data Cleaning and Preprocessing
	
    anat_file = glob(fullfile(fmriprep_dir, ptc, 'anat', sprintf('%s_desc-preproc_T1w.nii.gz',ptc)));

    acpc_temp_dir = fullfile(temp_folder, 'acpc_detect', AG, exp);

    args = {'-center-AC', '-no-tilt-correction','-output-orient RAS'};

    success = 0;
    n = 0;
    while success == 0 
        testf = sprintf('test%d',n);
        p_temp_dir = fullfile(acpc_temp_dir,ptc,testf);
        if ~exist(p_temp_dir,'dir'); mkdir(p_temp_dir); end        
        
        acpc_file = run_node('run_acpcdetect', argsnode, anat_file, p_temp_dir, p_cerebel_dir, args);
        if ~isempty(acpc_file)
            success = 1;
        else
            rmdir(p_temp_dir, 's');
            n = n+1;
        end 
    end

    suit_temp_nat = run_node('run_suit', argsnode, acpc_file, p_cerebel_dir, suit_template);

    sub_T1w = glob(fullfile(p_cerebel_dir,'*_desc-preproc_T1w.nii'));
    infovol = spm_vol(sub_T1w{:});

    p_ROIs = fullfile(ROIs_dir,ptc);
    optbids = struct('dir',p_ROIs, 'sub',ptc, 'space','T1w', 'desc','suit', 'desc_atlas',desc_atlas);

    resample_to_im_gz = glob(fullfile(fmriprep_dir, ptc, '**', 'func','*_space-T1w_boldref.nii.gz'));
    atlas_filegz = glob(fullfile(p_ROIs, 'func', sprintf('*_desc-%s_dseg.nii.gz',desc_atlas)));
    files_gz = {resample_to_im_gz(1), atlas_filegz};
    gunzip_files = run_node('run_gunzip', argsnode, files_gz, p_cerebel_dir);

    resample_to_im = gunzip_files{1};
    % argsnode.rerun = 1;
    [rois_anat, rois_func] = run_node('get_ROIs', argsnode, suit_temp_nat, suit_lut, idxs, optbids, infovol, resample_to_im);
    % argsnode.rerun = 0;

    if subs_table(sub_idx,:).rsfMRI
        atlas_json = glob(fullfile(p_ROIs, 'func', sprintf('*_desc-%s_dseg.json',desc_atlas)));
        rois = rois_func(contains(rois_func, '_roi-CerebellumNuclei_'));
        atlas_file = gunzip_files{2};
        [new_atlas, new_atlas_json] = run_node('update_atlas', argsnode, atlas_file, atlas_json, rois, optbids);
    end

end

