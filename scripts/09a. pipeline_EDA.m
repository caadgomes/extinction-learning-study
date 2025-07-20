%
clear all
% pspm_data_editor
AG = 'A02';
exp = 'Extinction_EEG_fMRI';

hiwis_dir = '/home/data';
rawdata_dir = fullfile(hiwis_dir, AG, exp, 'rawdata');
wf_dir = fullfile(hiwis_dir, AG, exp, 'workflows');
pspm_derivatives_dir = fullfile(hiwis_dir, AG, exp, 'derivatives', 'pspm');
gen_miss = 0;
filter_type = 'median'; %'butter', 'median', ''
T_struct = struct();
model_conds = struct('name','Context+CSs+USs', ...
    'values',{{'Context','CS+US','CS+noUS','CS-','US','noUS+'}});

% Different groups use different names for the events
AG_settings = check_json(AG, fullfile(hiwis_dir, 'PsPM_settings.json'));

% place where all the inputs and outputs will be saved
name_wf = 'miss';
wf_dir = fullfile(wf_dir, 'wf_pspm_miss', name_wf);
if ~exist(wf_dir,'dir'); mkdir(wf_dir); end

% initialise PsPM
pspm_init;

% get subjects
exclude = {
    };

subs_table = readtable(fullfile(rawdata_dir,'participants.tsv'), 'FileType','delimitedtext', 'Delimiter','\t');

% helper functions
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% IMPORT %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
for sub_idx = 1:height(subs_table)
    %     sub_idx = 1

    argsnode = struct('p_wf_dir','', 'rerun',0, 'name_node','');

    ptc = subs_table(sub_idx,:).participant_id{:};
    sfbcode = subs_table(sub_idx,:).SFB1280_code{:};
    %     if any(strcmp(ptc, exclude))
    %         continue
    %     end

    p_rawdata_dir = fullfile(rawdata_dir, ptc);
    p_derivatives_dir = fullfile(pspm_derivatives_dir, ptc);
    if ~exist(p_derivatives_dir,'dir'); mkdir(p_derivatives_dir); end

    p_wf_dir = fullfile(wf_dir, ptc);
    if ~exist(p_wf_dir,'dir'); mkdir(p_wf_dir); end

    task_files = glob(fullfile(p_rawdata_dir,'func','*task-*_physio.tsv.gz'));

    tasks = cellfun(@(x) extractBetween(x,'task-','_'), task_files);
    if isempty(tasks); continue; end

    tasks = intersect(tasks, {'acquisition'});

    for t = 1:numel(tasks)

        task = tasks{t};

        argsnode.p_wf_dir = fullfile(p_wf_dir, task);

        events_file = fullfile(p_rawdata_dir,'func',sprintf('%s_task-%s_events.tsv',ptc,task));
        onsets_dat = readtable(events_file, 'FileType','delimitedtext', 'Delimiter','\t');

        %argsnode.rerun = 1;
        [eda_fstruct, pars] = run_node('write_eda', argsnode, ptc, p_rawdata_dir,p_derivatives_dir,task);
        %argsnode.rerun = 0;

        info_markers = onsets_dat(startsWith(onsets_dat.trial_type, AG_settings.onsets.trial_type_names.(task)),:);
        [import_fstruct,completed] = run_node('import_eda', argsnode, eda_fstruct, pars, info_markers);

        if filter_type
            json_fname = fullfile(pspm_derivatives_dir, sprintf('task-%s_desc-Filtering_param.json',task));
            filt_pars_sub = check_json(sfbcode, json_fname, 0, AG_settings.filter);
            filt_opts = filt_pars_sub.(filter_type);
%                         argsnode.rerun=1;
            filtered_eda = run_node('filter_eda', argsnode, ptc, task, import_fstruct, p_derivatives_dir, filter_type, filt_opts);
%                         argsnode.rerun=0;
            json_fname = fullfile(pspm_derivatives_dir, sprintf('task-%s_desc-Trim_param.json',task));
            trim_pars_sub = check_json(sfbcode, json_fname, 0, AG_settings.trim);
            argsnode.rerun=1;
            [trimmed_fstruct, onsets_dat] = run_node('trim_eda', argsnode, filtered_eda, trim_pars_sub, p_derivatives_dir, ptc, task, onsets_dat);
            info_markers = onsets_dat(startsWith(onsets_dat.trial_type,AG_settings.onsets.trial_type_names.(task)),:); %onsets_dat has been updated
            argsnode.rerun=0;
        end

        missing_epochs = [];
        if gen_miss
            json_fname = fullfile(pspm_derivatives_dir, sprintf('task-%s_desc-MissingEpochs_param.json',task));
            default_pars = AG_settings.missing;
            miss_pars_sub = check_json(sfbcode, json_fname, 0, default_pars);
            miss_pars_sub.missing_epochs_filename = fullfile(p_derivatives_dir, sprintf('%s_task-%s_MissingEpochs.mat',ptc,task));

            if filter_type  
                file_to_miss = trimmed_fstruct;
                fd_name = fieldnames(filt_opts);
                fdfilt = struct(['afterfilt_',filter_type],1, ['afterfilt_',fd_name{:}],filt_opts.(fd_name{:}));
                miss_pars_sub = mergestructs(miss_pars_sub,fdfilt);
            else
                file_to_miss = import_fstruct;
            end
            missing_epochs = run_node('gen_missing_epochs', argsnode, ptc, task, file_to_miss, p_derivatives_dir, pars, miss_pars_sub, info_markers);
        end

    end
end
