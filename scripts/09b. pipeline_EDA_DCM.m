
clearvars;

warning('off', 'MATLAB:rmpath:DirNotFound');
% remove SPM from path because it conflicts with SPM
rmpath(genpath('/usr/local/MATLAB/TOOLBOXES/spm12/'));
addpath('/media/sf_G_80tb/installers/Toolboxes/PsPM6.1/');

% pspm_data_editor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% DEFINE GROUP, STUDY AND BIDS DIR %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AG = 'A09';
exp = 'Extinction_Generalization_II';
bids_dir = '/home/data';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rawdata_dir = fullfile(bids_dir, AG, exp, 'rawdata');
wf_dir = fullfile(bids_dir, AG, exp, 'workflows');
pspm_derivatives_dir = fullfile(bids_dir, AG, exp, 'derivatives', 'pspm');
group_dir = fullfile(pspm_derivatives_dir, 'group');
T_stats = table(); T_con = table(); T_diagn = table();
combine_tasks = false;
conds_dcm = 'all';
% Without a delay, the CS_onset and CS_interval impulses cannot be separated
% If there is a response to CS_onset but not during CS_interval, the latter may explain away the former (by assigning a very short latency),
% and the amplitudes of the two events become ambiguous
% So, add a delay after CS_onset (e.g., 2 seconds) in first column of 3rd event of the DCM events struct.
AG_settings = check_json(AG, fullfile(bids_dir, 'PsPM_settings.json'));
dcm_opt = struct('substhresh',4, 'constrained',1);

AG_settings.dcm.delay = struct('acquisition',[2.5,2.5], 'extinction',[2.5,2.5], 'renewal',[2.5,2.5], 'recall',[2.5,2.5]);
addFix = 0; % this is mostly useful when the CS-US interval is short and you either want to model the fixed response OR the flexible response, but not both

name_wf = 'DelayRen2p5Const';
wf_dir = fullfile(wf_dir, 'wf_pspm_DCM', name_wf);
if ~exist(wf_dir,'dir'); mkdir(wf_dir); end

name_DCM_model = name_wf;

% double curvy brackets ensure fields are treated as cells
cons = struct('acquisition',containers.Map(), 'extinction',containers.Map(), 'recall',containers.Map(), 'renewal',containers.Map(), 'all',containers.Map());

cons.acquisition('CS+US_VS_CS+noUS') = struct('con',{{'CS+US','CS+noUS'}}, 'type','simple', 'dir','pos');
cons.acquisition('CS+noUS_linear') = struct('con',{{'CS+noUS'}}, 'type','linear', 'dir','pos');
cons.acquisition('CS+_linear') = struct('con',{{'CS+US','CS+noUS'}}, 'type','linear', 'dir','pos');
cons.acquisition('CS+_exp') = struct('con',{{'CS+US','CS+noUS'}}, 'type','expcdf', 'dir','pos');
cons.acquisition('CS+noUS_exp') = struct('con',{{'CS+noUS'}}, 'type','expcdf', 'dir','pos');

cons.extinction('CS+noUS_linear') = struct('con',{{'CS+noUS'}}, 'type','linear', 'dir','pos');
cons.extinction('CS+noUS_exp') = struct('con',{{'CS+noUS'}}, 'type','exppdf', 'dir','pos');

if AG == 'A03'
    cons.renewal('CS+noUS_linear') = struct('con',{{'CS+noUS'}}, 'type','linear', 'dir','pos');
    cons.renewal('CS+noUS_exp') = struct('con',{{'CS+noUS'}}, 'type','exppdf', 'dir','pos');
elseif any(strcmp(AG, {'A05','A09'}))
    cons.recall('CS+noUS_linear') = struct('con',{{'CS+noUS'}}, 'type','linear', 'dir','pos');
    cons.recall('CS+noUS_exp') = struct('con',{{'CS+noUS'}}, 'type','exppdf', 'dir','pos');
end

cons.all('CS+US_VS_CS+noUS') = struct('con',{{'CS+US','CS+noUS'}}, 'type','simple', 'dir','pos');
cons.all('CS+noUS_linear') = struct('con',{{'CS+noUS'}}, 'type','linear', 'dir','pos');
cons.all('CS+_linear') = struct('con',{{'CS+US','CS+noUS'}}, 'type','linear', 'dir','pos');
cons.all('CS+_exp') = struct('con',{{'CS+US','CS+noUS'}}, 'type','expcdf', 'dir','pos');
cons.all('CS+noUS_exp') = struct('con',{{'CS+noUS'}}, 'type','expcdf', 'dir','pos');

% diagnostics
con_diag = struct('acquisition',containers.Map(), 'extinction',containers.Map(), 'recall',containers.Map(), 'renewal',containers.Map(), 'all',containers.Map());
con_diag.acquisition('CS+US_VS_CS+noUS') = struct('A02',{{'CS+US','CS+noUS'}}, 'A03',{{'CS+US','CS+noUS'}}, 'A05',{{'CS+US','CS+noUS'}}, 'A09',{{'CS+US_N','CS+noUS_N'}}, 'tail','right');

% initialise PsPM
pspm_init;

% get subjects
excludesubs = readtable(fullfile(pspm_derivatives_dir,'desc-excludeSubs_table.tsv'), 'FileType','delimitedtext', 'Delimiter','\t');
subs_table = readtable(fullfile(rawdata_dir,'participants.tsv'), 'FileType','delimitedtext', 'Delimiter','\t');

% helper functions
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

subs_not_run = array2table(cell(0,2), 'VariableNames',{'name','task'});

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% RUN %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

get_group_stats = 1;

choose_tasks = {'recall'};

parfor sub_idx = 1:height(subs_table)

    argsnode = struct('p_wf_dir','', 'rerun',0, 'name_node','');

    ptc = subs_table(sub_idx,:).participant_id{:};
    exc_sub = excludesubs(strcmp(excludesubs.participant_id,ptc),:);

    sfbcode = subs_table(sub_idx,:).SFB1280_code{:};

    p_rawdata_dir = fullfile(rawdata_dir, ptc);

    p_wf_dir = fullfile(wf_dir, ptc);
    if ~exist(p_wf_dir,'dir'); mkdir(p_wf_dir); end

    argsnode.p_wf_dir = p_wf_dir;

    task_files = glob(fullfile(p_rawdata_dir,'func','*task-*_physio.tsv.gz'));

    tasks = cellfun(@(x) extractBetween(x,'task-','_'), task_files);
    tasks = intersect(tasks, choose_tasks);
    
    if isempty(tasks); continue; end

    p_derivatives_dir = fullfile(pspm_derivatives_dir, ptc);

    valid_tasks = tasks;

    dcm_struct = struct('timing',{{}}, 'condition',{{}}, 'run_sub',false, 'options',struct());
    temp = struct('data_file',{{}}, 'missing_file',{{}});

    for t = 1:numel(tasks)

        task = tasks{t};

        if exc_sub.(task) % if sub needs to be excluded
            valid_tasks(strcmp(valid_tasks,task)) = []; % do not process this phase for this sub
            continue
        end

        fprintf('\n\nProcessing subject: %s\nProcessing task: %s\n\n',ptc,task);

        % create subject-specific folder only if data are valid
        if ~exist(p_derivatives_dir,'dir'); mkdir(p_derivatives_dir); end

        events_file = fullfile(p_derivatives_dir,sprintf('%s_task-%s_desc-trimmed_events.tsv',ptc,task));
        onsets_dat = readtable(events_file, 'FileType','delimitedtext', 'Delimiter','\t');
        % CS and US often co-terminate which is problematic if US
        % presentation is very long. Change the duration for CS from CS
        % onset till US onset
        onsets_dat = change_dur_cs(onsets_dat);
        sess_onsets_dat = {onsets_dat};
        % uncomment next line, if you wish to model all noUS as a single condition, i.e., combine noUS+ and noUS-
        %         onsets_dat.trial_type(startsWith(onsets_dat.trial_type,'noUS')) = {'noUS'};

        % the final processed eda datafile is the trimmed one
        sess_imported_file = glob(fullfile(p_derivatives_dir,sprintf('t*%s_task-%s_physio.mat',ptc,task)));
        sess_missing_epochs = glob(fullfile(p_derivatives_dir,sprintf('%s_task-%s_MissingEpochs.mat',ptc,task)));

        optsp = struct('task',task, 'ptc',ptc);

        if strcmp(AG,'A12') && strcmp(task, 'extinction')
            optsp.trim_at = {'Ratings_'}; optsp.rm_part = {{'US_','US_'}};
            sfiles = run_node('split_EDA', argsnode, sess_imported_file,sess_missing_epochs,onsets_dat,p_derivatives_dir,optsp);
            sess_imported_file = {sfiles(1:2).eda};
            sess_missing_epochs = {sfiles(1:2).miss};
            sess_onsets_dat = {sfiles(1:2).onsets};
        elseif strcmp(AG,'A09') && strcmp(task, 'recall')
            optsp.split = {struct('cond',["US"],'keep',0)};
            sfiles = run_node('split_EDA', argsnode, sess_imported_file,sess_missing_epochs,onsets_dat,p_derivatives_dir,optsp);
            sess_imported_file = {sfiles(1:2).eda};
            sess_missing_epochs = {sfiles(1:2).miss};
            sess_onsets_dat = {sfiles(1:2).onsets};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% DCM CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        further_cols = table(); %this is to add further info about the trials to the stats table
        for ses = 1:numel(sess_imported_file)

            imported_file = sess_imported_file(ses);
            missing_epochs = sess_missing_epochs(ses);
            onsets_dat = sess_onsets_dat{ses};

            dcm_delay = AG_settings.dcm.delay.(task);
            timing = {};
            if any(strcmp(AG, {'A02','A03','A05','A09','A12'}))
                if any(strcmp(AG, {'A02','A03'}))
                    context_onsets = onsets_dat.onset(matches(onsets_dat.trial_type,"Context","IgnoreCase",true));
                    timing(length(timing)+1) = {context_onsets};
                end
                where_CSs = startsWith(onsets_dat.trial_type,["CS+","CS-"]); %AG_settings.onsets.trial_type_names

                if AG=='A02'
                    further_cols = [further_cols; onsets_dat(where_CSs,{'CS_type'})]; 
                end

                CS_onsets = onsets_dat.onset(where_CSs);
                if sum(dcm_delay) > 0 || addFix
                    timing(length(timing)+1) = {CS_onsets};
                end

                if any(strcmp(AG, {'A03','A05','A09','A12'}))
                    % US is always presented "US_dur" before CS offset
                    US_onsets = CS_onsets + onsets_dat.duration1(where_CSs);
                elseif strcmp(AG, {'A02'})
                    % US is presented before both the CS and the interval, so
                    % slice based on actual US onsets
                    where_USs = matches(onsets_dat.trial_type,["noUS","US"]);
                    US_onsets = onsets_dat.onset(where_USs);
                end

                CS_interval = [CS_onsets + dcm_delay(1), US_onsets - dcm_delay(2)];
                if sum(dcm_delay) > 0 || ~addFix
                    timing(length(timing)+1) = {CS_interval};
                end

                timing(length(timing)+1) = {US_onsets};

                if any(strcmp(AG, {'A02','A03','A05'}))
                    CSpnoUS_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS");
                    CSm_idx = find(onsets_dat.trial_type(where_CSs)=="CS-");
                    if (strcmp(task, 'acquisition'))
                        CSpUS_idx = find(onsets_dat.trial_type(where_CSs)=="CS+US");
                        index_conds = struct('name',{'CS+US';'CS+noUS';'CS-'}, 'index', {CSpUS_idx; CSpnoUS_idx; CSm_idx});
                    else
                        if strcmp(AG,'A02')
                            CSpUS_idx = find(onsets_dat.trial_type(where_CSs)=="CS+US");
                            CSmUS_idx = find(onsets_dat.trial_type(where_CSs)=="CS-US");
                            index_conds = struct('name',{'CS+US';'CS+noUS';'CS-US';'CS-'}, 'index', {CSpUS_idx; CSpnoUS_idx; CSmUS_idx; CSm_idx});
                        else    
                            index_conds = struct('name',{'CS+noUS';'CS-'}, 'index', {CSpnoUS_idx; CSm_idx});
                        end
                    end

                elseif strcmp(AG, 'A09')
                    if strcmp(task, 'acquisition')
                        CSpnoUS_N_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS_N");
                        CSpnoUS_G_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS_G");
                        CSm_idx = find(onsets_dat.trial_type(where_CSs)=="CS-");
                        CSpUS_N_idx = find(onsets_dat.trial_type(where_CSs)=="CS+US_N");
                        CSpUS_G_idx = find(onsets_dat.trial_type(where_CSs)=="CS+US_G");
                        index_conds = struct('name',{'CS+US_N';'CS+US_G';'CS+noUS_N';'CS+noUS_G';'CS-'}, ...
                            'index', {CSpUS_N_idx; CSpUS_G_idx; CSpnoUS_N_idx; CSpnoUS_G_idx; CSm_idx});
                    elseif strcmp(task, 'extinction')
                        CSpnoUS_N_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS_N");
                        CSpnoUS_G_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS_G");
                        CSpnoUS_G2_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS_G2");
                        CSpnoUS_G3_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS_G3");
                        CSpnoUS_G4_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS_G4");
                        CSm_idx = find(onsets_dat.trial_type(where_CSs)=="CS-");
                        index_conds = struct('name',{'CS+noUS_N';'CS+noUS_G';'CS+noUS_G2';'CS+noUS_G3';'CS+noUS_G4';'CS-'}, ...
                            'index', {CSpnoUS_N_idx; CSpnoUS_G_idx; CSpnoUS_G2_idx; CSpnoUS_G3_idx; CSpnoUS_G4_idx; CSm_idx});
                    elseif strcmp(task, 'recall')
                        CSpnoUS_N_idx = find(startsWith(onsets_dat.trial_type(where_CSs),"CS+noUS_N_"));
                        CSpnoUS_G_idx = find(startsWith(onsets_dat.trial_type(where_CSs),"CS+noUS_G_"));
                        CSpnoUS_N0_idx = find(startsWith(onsets_dat.trial_type(where_CSs),"CS+noUS_N0_"));
                        CSpnoUS_G0_idx = find(startsWith(onsets_dat.trial_type(where_CSs),"CS+noUS_G0_"));
                        CSm_idx = find(startsWith(onsets_dat.trial_type(where_CSs),regexpPattern("CS-_.$")));
                        CSm_0_idx = find(startsWith(onsets_dat.trial_type(where_CSs),"CS-_0_"));

                        index_conds = struct('name',{'CS+noUS_N';'CS+noUS_G';'CS+noUS_N0';'CS+noUS_G0';'CS-';'CS-0'}, ...
                            'index', {CSpnoUS_N_idx; CSpnoUS_G_idx; CSpnoUS_N0_idx; CSpnoUS_G0_idx; CSm_idx; CSm_0_idx});
                    end

                elseif strcmp(AG, 'A12')
                    CSpnoUS_Vis_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS_Visceral");
                    CSpnoUS_Aud_idx = find(onsets_dat.trial_type(where_CSs)=="CS+noUS_Auditive");
                    CSm_idx = find(onsets_dat.trial_type(where_CSs)=="CS-");
                    if strcmp(task, 'acquisition')
                        CSpUS_Vis_idx = find(onsets_dat.trial_type(where_CSs)=="CS+US_Visceral");
                        CSpUS_Aud_idx = find(onsets_dat.trial_type(where_CSs)=="CS+US_Auditive");
                        index_conds = struct('name',{'CS+US_Visceral';'CS+US_Auditive';'CS+noUS_Visceral';'CS+noUS_Auditive';'CS-'}, ...
                            'index', {CSpUS_Vis_idx; CSpUS_Aud_idx; CSpnoUS_Vis_idx; CSpnoUS_Aud_idx; CSm_idx});
                    else
                        US_Vis_idx = find(onsets_dat.trial_type(where_CSs)=="US_Visceral");
                        US_Aud_idx = find(onsets_dat.trial_type(where_CSs)=="US_Auditive");
                        index_conds = struct('name',{'CS+noUS_Visceral';'CS+noUS_Auditive';'CS-';'US_Visceral';'US_Auditive'}, ...
                            'index', {CSpnoUS_Vis_idx; CSpnoUS_Aud_idx; CSm_idx; US_Vis_idx; US_Aud_idx});
                    end
                end
            end

            dcm_struct.ptc = ptc;
            % store the imported files and missing epochs in case I wish to
            % combine tasks. temp structure is needed when parallelising
            temp.data_file(end+1) = imported_file;
            temp.missing_file(end+1) = missing_epochs;
            data_file = temp.data_file;
            missing_file = temp.missing_file;

            dcm_struct.timing(end+1) = {timing};
            dcm_struct.condition(end+1) = {index_conds};

            dcm_struct.run_sub = true;
            dcm_struct.name_model = name_DCM_model;
            dcm_struct.task = task;
            dcm_struct.interpnan = 0;
            dcm_struct.modargs = dcm_opt;
            dcm_struct.modopts = struct('');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% DCM - SEPARATE TASKS %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_out_dir = fullfile(pspm_derivatives_dir, ptc, 'DCM');
        if ~exist(p_out_dir,'dir'); mkdir(p_out_dir); end
        if dcm_struct.run_sub && ~combine_tasks
            dcm_struct.p_out_dir = p_out_dir;

            % if running tasks individually, change argsnode.name_node to
            % the name of the task, so that results of the workflow are
            % stored in the task's subfolder
            argsnode.name_node = task;

%                         argsnode.rerun = 1;
            dcm = run_node('run_pspm_dcm', argsnode, data_file, missing_file, dcm_struct);
%                         argsnode.rerun = 0;

            if ~isempty(dcm)
                stat = 'dcm';
                info_betas = struct('sub',{ptc}, 'AG',{AG}, 'exp',{exp}, 'task',{task}, 'p_out_dir',p_out_dir, ...
                    'further_cols',further_cols, 'name_model',name_DCM_model);
                argsnode.name_node = fullfile(task, 'stats');

                argsnode.rerun = 1;
                t_stats = run_node('extract_betas', argsnode, stat, dcm, info_betas);
                argsnode.rerun = 0;

                if any(strcmp(task, {'renewal','recall'}))
                    if (AG=='A03')
                        if strcmp(exp,'7T')
                            eg = {'renewal'};
                        else
                            eg = subs_table{strcmp(subs_table.participant_id,ptc),'experimental_group'}; %if isnan(eg); eg = {'NaN'}; end
                        end
                        t_stats.task(:) = eg;
                    elseif AG=='A09'
                        t_stats(t_stats.session == 2, 'task') = {'reinstatement'};
                    end
                elseif any(strcmp(task, {'extinction'}))
                    if strcmp(AG,'A12')
                        first_US = find(startsWith(onsets_dat.trial_type,'US'),1,'first');
                        t_stats(first_US:end,'task') = {'reinstatement'};
                    end
                end

                tsk = unique([t_stats.task]);
                for tt = 1:numel(tsk)
                    tempt = t_stats(strcmp(t_stats.task,tsk{tt}),:);
                    writetable(tempt, fullfile(p_out_dir, sprintf('%s_stat-DCM_task-%s_desc-%s_df.tsv',ptc,tsk{tt},name_DCM_model)), 'FileType','Text', 'Delimiter','\t');
                end
                T_stats = [T_stats; t_stats];

            else
                subs_not_run = [subs_not_run;[ptc, {task}]];
            end
            % reset dcm_struct, otherwise PsPM will append next task to
            % dcm_struct
            dcm_struct = struct('datafile',{{}} ,'timing',{{}}, 'condition',{{}}, 'missing',{{}}, 'run_sub',false, 'options',struct());
            temp = struct('data_file',{{}}, 'missing_file',{{}});
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% DCM - COMBINE TASKS %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if dcm_struct.run_sub && combine_tasks

        % name node explicitly because some subjects have missing tasks and
        % node will be renamed whichever task comes last
        argsnode.p_wf_dir = p_wf_dir;
        argsnode.name_node = 'allTasks';
        dcm_struct.task = 'all';
        dcm = run_node('run_pspm_dcm', argsnode, ptc, p_out_dir, dcm_struct);

        if ~isempty(dcm)
            stat = 'dcm';
            argsnode.name_node = fullfile(argsnode.name_node, 'stats');
            info = struct('tasks',tasks);
            t_stats = run_node('extract_betas', argsnode, stat, dcm, info);
            T_stats = [T_stats; t_stats];

            modelfile = dcm.modelfile;
            argsnode.name_node = fullfile(argsnode.name_node,'cons');
            dcmcon = run_node('run_cons', argsnode, modelfile, cons, task);
            stat = 'con';
            t_con = run_node('extract_betas', argsnode, stat, dcmcon);
            T_con = [T_con; t_con];
            argsnode.name_node = '';
        else
            subs_not_run = [subs_not_run;[ptc, {task}]];
        end
        modelfile = dcm.modelfile;
        dcmcon = run_node('run_cons', argsnode, modelfile, cons);
    end

end

fprintf('\n\n%s\nThese subjects could not be run: \n\n', repmat('#',1,80));
disp(sortrows(subs_not_run))
fprintf('%s\n', repmat('#',1,80));

if get_group_stats
    if ~exist(group_dir,'dir'); mkdir(group_dir); end

    writetable(T_stats, fullfile(group_dir, sprintf('stat-DCM_type-betas_desc-%s_df.tsv',name_DCM_model)), 'Delimiter','\t', 'FileType','text');
    writetable(T_con, fullfile(group_dir, sprintf('stat-DCM_type-cons_desc-%s_df.tsv',name_DCM_model)), 'Delimiter','\t', 'FileType','text');

end

