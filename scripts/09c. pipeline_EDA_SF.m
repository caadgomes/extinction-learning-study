
clearvars;

warning('off', 'MATLAB:rmpath:DirNotFound');
% remove SPM from path because it conflicts with SPM
rmpath(genpath('/usr/local/MATLAB/TOOLBOXES/spm12UP/'));
% addpath('/usr/local/MATLAB/TOOLBOXES/PsPM_v6.0.0/');

% pspm_data_editor
AG = 'A12';
exp = 'study_2';
bids_dir = '/home/data/';
rawdata_dir = fullfile(bids_dir, AG, exp, 'rawdata');
wf_dir = fullfile(bids_dir, AG, exp, 'workflows');
pspm_derivatives_dir = fullfile(bids_dir, AG, exp, 'derivatives', 'pspm');
group_dir = fullfile(pspm_derivatives_dir, 'group');
T_stats = table(); T_con = table(); T_diagn = table();
model_conds = struct('name','Context+CSs+USs', ...
    'values',{{'Context','CS+US','CS+noUS','CS-','US','noUS+'}});
combine_tasks = false;
method_sf = 'dcm';


%'cond' uses the condition names (specified by user below) to create onsets
%and offsets. 'break' uses much longer segments, and the user needs to
%specify where to break the onsets using break_at
%if you want to use the entire data as a single segment, use break and
%leave break_at empty
epoch_type = 'break';
%missLongAcq means use missing values option (miss), a single break (Long) and
%acquisition (Acq)
descname = 'missLongAcq'; 
break_at = {}; %where to break onsets

choose_tasks = {'acquisition'};

% options for run_extract_sf
tsec = 0; % take this amount (in seconds) from CS offset

conductance_delay = 10;

T_sf = table();
AG_settings = check_json(AG, fullfile(bids_dir, 'PsPM_settings.json'));

name_wf = sprintf('%s',[epoch_type,descname]);
wf_dir = fullfile(wf_dir, 'wf_pspm_SF', name_wf);
if ~exist(wf_dir,'dir'); mkdir(wf_dir); end

name_SF_model = name_wf;

fun_sf = 'mean';

% double curvy brackets ensure fields are treated as cells
cons = struct('acquisition',containers.Map(), 'extinction',containers.Map(), 'recall',containers.Map(), 'renewal',containers.Map(), 'all',containers.Map());

% diagnostics
con_diag = struct('acquisition',containers.Map(), 'extinction',containers.Map(), 'recall',containers.Map(), 'renewal',containers.Map(), 'all',containers.Map());
con_diag.acquisition('CS+US_VS_CS+noUS') = struct('A02',{{'CS+US','CS+noUS'}}, 'A03',{{'CS+US','CS+noUS'}}, ...
    'A05',{{'CS+US','CS+noUS'}}, 'A09',{{'CS+US_N','CS+noUS_N'}}, 'A12',{{'CS+US_Visceral','CS+noUS_Visceral'}}, 'tail','right');

% initialise PsPM
pspm_init;

% get subjects
exclude_subs = readtable(fullfile(pspm_derivatives_dir,'desc-excludeSubs_table.tsv'), 'FileType','delimitedtext', 'Delimiter','\t');
subs_table = readtable(fullfile(rawdata_dir,'participants.tsv'), 'FileType','delimitedtext', 'Delimiter','\t');

% helper functions
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

subs_not_run = array2table(cell(0,2), 'VariableNames',{'name','task'});

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% RUN %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
get_group_stats = 1;

parfor sub_idx = 1:height(subs_table)

    argsnode = struct('p_wf_dir','', 'rerun',0, 'name_node','');

    ptc = subs_table(sub_idx,:).participant_id{:};
    exc_sub = exclude_subs(strcmp(exclude_subs.participant_id,ptc),:);

    sfbcode = subs_table(sub_idx,:).SFB1280_code{:};

    p_rawdata_dir = fullfile(rawdata_dir, ptc);

    p_wf_dir = fullfile(wf_dir, ptc);
    if ~exist(p_wf_dir,'dir'); mkdir(p_wf_dir); end

    argsnode.p_wf_dir = p_wf_dir;

    task_files = glob(fullfile(p_rawdata_dir,'func','*task-*_physio.tsv.gz'));

    tasks = cellfun(@(x) extractBetween(x,'task-','_'), task_files);
    if isempty(tasks); continue; end

    % If you want to run the script just for some of the tasks then
    % uncomment the following following lines
    tasks = intersect(tasks, choose_tasks);

    p_derivatives_dir = fullfile(pspm_derivatives_dir, ptc);

    valid_tasks = tasks;

    sf_struct = struct('timing','', 'timeunits','seconds', 'run_sub',false, 'options',struct());

    for t = 1:numel(tasks)

        task = tasks{t};

        if exc_sub.(task)
            valid_tasks(strcmp(valid_tasks,task)) = [];
            continue
        end

        fprintf('\n\nProcessing subject: %s\nProcessing task: %s\n\n',ptc,task);

        % create subject-specific folder only if data are valid
        if ~exist(p_derivatives_dir,'dir'); mkdir(p_derivatives_dir); end

        events_file = fullfile(p_derivatives_dir,sprintf('%s_task-%s_desc-trimmed_events.tsv',ptc,task));
        onsets_dat = readtable(events_file, 'FileType','delimitedtext', 'Delimiter','\t');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% SF CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(epoch_type, 'cond')
            if strcmp(task, 'acquisition')
                if strcmp(AG,'A12')
                    onsets_dat = change_dur_cs(onsets_dat, 'Merge',{{'CS+US','US'}}, 'AddTime',{{'CS+US_Visceral'},[-1]});
                elseif strcmp(AG,'A03')
                    onsets_dat = change_dur_cs(onsets_dat, 'Merge',{{'CS+US','US'}}, 'ExcludeCond',{'Fixation','noUS-'});
                elseif ismember(AG,{'A05','A09'})
                    onsets_dat = change_dur_cs(onsets_dat, 'Merge',{{'CS+US','US'}});
                elseif ismember(AG,{'A02'})
                    onsets_dat = change_dur_cs(onsets_dat, 'Merge',{{'CS+US','US'},{'CS+noUS','noUS'},{'CS-','noUS'}},'ExcludeConds',{'interval'});
                end
            else
                if strcmp(AG,'A03')
                    onsets_dat = change_dur_cs(onsets_dat, 'Merge',{{'CS+US','US'}}, 'ExcludeCond',{'Fixation','noUS-'});
                else
                    onsets_dat = change_dur_cs(onsets_dat);
                end
            end
            onsets = onsets_dat{:,'onset'};
            offsets = onsets_dat{:,'onset'} + onsets_dat{:,'duration2'}; %duration2 is the duration with US
        elseif strcmp(epoch_type, 'break')
            onsets_dat = change_dur_cs(onsets_dat, 'Merge',{{'CS+US','US'}}, 'AddTime',{{'CS+US_Visceral'},[-1]});
            if ~isempty(break_at)
                breaks = find(startsWith(onsets_dat.trial_type,break_at));
                from_idx = 1;
                onsets = [];  offsets = [];
                for b = 1:numel(breaks)
                    onsets(end+1,:) = onsets_dat{from_idx,'onset'};
                    offsets(end+1,:) = onsets_dat{breaks(b),'onset'};

                    onsets(end+1,:) = onsets_dat{breaks(b),'onset'};
                    if (breaks(b)+1)>height(onsets_dat)
                        offsets(end+1,:) = onsets_dat{breaks(b),'onset'} + onsets_dat{breaks(b),'duration1'};
                    else
                        offsets(end+1,:) = onsets_dat{breaks(b)+1,'onset'};
                    end
                    from_idx = breaks(b)+1;
                end
            else
                onsets = onsets_dat.onset(1);
                offsets = onsets_dat.onset(end) + onsets_dat.duration(end);
            end
        end

        if AG=='A12'
            d = diff([0; find(startsWith(onsets_dat.trial_type,'Ratings_'))]);
            sessions = repelem(1:numel(d), d)';
            onsets_dat.session = sessions;
        end

        % the final processed eda datafile is the trimmed one
        imported_file = glob(fullfile(p_derivatives_dir,sprintf('t*%s_task-%s_physio.mat',ptc,task)));
        missing_epochs_file = glob(fullfile(p_derivatives_dir,sprintf('%s_task-%s_MissingEpochs.mat',ptc,task)));

        p_out_dir = fullfile(pspm_derivatives_dir, ptc, 'SF');
        if ~exist(p_out_dir,'dir'); mkdir(p_out_dir); end

        sf_struct.ptc = ptc;
        sf_struct.p_out_dir = p_out_dir;
        sf_struct.timing = [onsets, offsets + conductance_delay];
        sf_struct.timeunits = 'seconds';
        sf_struct.run_sub = true;
        sf_struct.name_model = name_SF_model;
        sf_struct.task = task;
        sf_struct.modargs = struct('method',method_sf);
        sf_struct.interp = 0; %should missing data (nans) be interpolated

        argsnode.name_node = task;

        % pspm_sf returns a FILE, not a structure
        %         argsnode.rerun = 1;
        sf = run_node('run_pspm_sf', argsnode, imported_file, missing_epochs_file, sf_struct);
        %         argsnode.rerun = 0;

        if ~isempty(sf)
            stat = 'sf';
            sf_dat = load(sf{:});  sf_dat = sf_dat.sf;
            excl_conds_sf = {{'Ratings_'}}; %you can give partial matches to pass to startsWith
            info = struct('sub',{ptc}, 'AG',{AG}, 'exp',{exp}, 'task',{task}, 'name_model',name_SF_model, 'tsec',tsec, 'exclude_conds',excl_conds_sf, 'fun',fun_sf, ...
                'epoch_type',epoch_type, 'onsoff', [onsets,offsets]); %onsoff: to be able to disambiguate by onset time

            %             argsnode.rerun = 1;
            t_sf = run_node('run_extract_sf', argsnode, sf_dat, onsets_dat, info);
            %             argsnode.rerun = 0;

            if strcmp(AG,'A12')
                last_US = find(startsWith(t_sf.cond,'US_'),1,'first');
                t_sf(last_US:end,'task') = {'reinstatement'};
            end

            tsk = unique(t_sf.task);
            for tt = 1:numel(tsk)
                tempt = t_sf(strcmp(t_sf.task,tsk{tt}),:);
                writetable(tempt, fullfile(p_out_dir, sprintf('%s_stat-SF_task-%s_desc-%s_df.tsv',ptc,tsk{tt},[name_SF_model,fun_sf])), 'FileType','Text', 'Delimiter','\t');
            end
            T_sf = [T_sf; t_sf];

            opts = struct('ptc',ptc, 'AG',AG, 'exp',exp, 'task',task, 'tbl',t_sf, 'amp','amplitude', ...
                'cons',con_diag, 'model',sf_dat);
        else
            subs_not_run = [subs_not_run;[ptc, {task}]];
        end
    end

end

fprintf('\n\n%s\nThese subjects could not be run: \n\n', repmat('#',1,80));
disp(sortrows(subs_not_run))
fprintf('%s\n', repmat('#',1,80));

if get_group_stats
    if ~exist(group_dir,'dir'); mkdir(group_dir); end

    pspm_SFs = glob(fullfile(pspm_derivatives_dir,'**', sprintf('*_stat-SF_desc-%s_model*',[name_SF_model,fun_sf])));
    writetable(T_sf, fullfile(group_dir, sprintf('stat-SF_type-betas_desc-%s_df.tsv',[name_SF_model,fun_sf])), 'Delimiter','\t','FileType','text');

end

