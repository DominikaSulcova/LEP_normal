%% params
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave masterfiles
folder.input = uigetdir(pwd, 'Coose the input folder');             % raw data
folder.output = uigetdir(pwd, 'Choose the OneDrive folder');        % output folder --> figures, loutput file, exports 
cd(folder.output)

% input & output 
study = 'AperiodicPFC';
input_file = sprintf('%s\\NLEP_output.mat', folder.output);
output_file = sprintf('%s\\%s_output.mat', folder.output, study);

% dataset
params.subjects = 45;
params.area = {'hand' 'foot'};
params.side = {'right' 'left'}; 
params.block = {'b1' 'b2'};
params.LEP_comps = {'N1' 'N2' 'P2'}; 

% graphics
figure_counter = 1;

%% import existing measures and pre-processed data
% ----- section input -----
params.prestim_time = 'ready';
params.prefix = {'icfilt ica_all RS' 'ep reref ds notch bandpass dc'};
% -------------------------
% load NLEP output & AperiodicPFC output, if it already existis
load(input_file, 'NLEP_info', 'NLEP_data', 'NLEP_data_1to35', 'NLEP_measures')
if exist(output_file) == 2
    load(output_file, 'AperiodicPFC_data', 'AperiodicPFC_measures')
end

% extract LEP measures and pain ratings
for s = 1:params.subjects
    % subject info
    AperiodicPFC_measures(s).ID = NLEP_info.single_subject(s).ID;
    AperiodicPFC_measures(s).age = NLEP_info.single_subject(s).age;
    AperiodicPFC_measures(s).male = NLEP_info.single_subject(s).male;
    AperiodicPFC_measures(s).handedness_score = NLEP_info.single_subject(s).handedness;
    if NLEP_info.single_subject(s).handedness >= 0.3
        AperiodicPFC_measures(s).handedness = 'right';
    elseif NLEP_info.single_subject(s).handedness <= -0.3
        AperiodicPFC_measures(s).handedness = 'left';
    else
        AperiodicPFC_measures(s).handedness = 'bilateral';
    end

    % session info
    for c = 1:2
        condition = split(NLEP_info.single_subject(s).condition{c}, '_');
        AperiodicPFC_measures(s).conditions(c).condition = c; 
        AperiodicPFC_measures(s).conditions(c).area = condition{1}; 
        AperiodicPFC_measures(s).conditions(c).side = condition{2}; 
    end

    % LEP peak measures
    for c = 1:2
        if contains(NLEP_measures(s).LEP_ST.conditions{c}, AperiodicPFC_measures(s).conditions(c).area) && ...
                contains(NLEP_measures(s).LEP_ST.conditions{c}, AperiodicPFC_measures(s).conditions(c).side)
            AperiodicPFC_measures(s).LEP(c).condition = c;
            for a = 1:length(params.LEP_comps)
                % extract amplitude and latency
                AperiodicPFC_measures(s).LEP(c).(params.LEP_comps{a}).amplitude = squeeze(NLEP_measures(s).LEP_ST.amplitude(a, c, :))';
                AperiodicPFC_measures(s).LEP(c).(params.LEP_comps{a}).latency = squeeze(NLEP_measures(s).LEP_ST.latency(a, c, :))';
                
                % check for missing events
                missing{a} = find([AperiodicPFC_measures(s).LEP(c).(params.LEP_comps{a}).amplitude] == 0);
                AperiodicPFC_measures(s).LEP(c).(params.LEP_comps{a}).amplitude(missing{a}) = [];
                AperiodicPFC_measures(s).LEP(c).(params.LEP_comps{a}).latency(missing{a}) = [];                
            end

            % encode trial number
            if all(cellfun(@(x) isequal(x, missing{1}), missing))
                AperiodicPFC_measures(s).LEP(c).trials = length([AperiodicPFC_measures(s).LEP(c).(params.LEP_comps{a}).amplitude]);
            else
                error('ERROR: subject %d - LEPs - condition %d (%s %s) - the numbers of trials of LEP components do not match!',...
                    s, c, AperiodicPFC_measures(s).conditions(c).area, AperiodicPFC_measures(s).conditions(c).side)
            end
        else
            error('ERROR: subject %d - LEPs - the conditions do not match!', s)
        end
    end

    % pain ratings
    for c = 1:2
        % identify ratings associated with this condition
        idx = false(1, length(NLEP_measures(s).conditions));
        if size(NLEP_measures(s).pain, 1) == length(NLEP_measures(s).conditions)
            for b = 1:size(NLEP_measures(s).pain, 1)
                if contains(NLEP_measures(s).conditions{b}, AperiodicPFC_measures(s).conditions(c).area) && ...
                        contains(NLEP_measures(s).conditions{b}, AperiodicPFC_measures(s).conditions(c).side) 
                    idx(b) = true;
                end
            end
        end
        ratings = NLEP_measures(s).pain(idx, :);

        % encode
        if AperiodicPFC_measures(s).LEP(c).trials == 60
            AperiodicPFC_measures(s).pain(c).condition = c;
            AperiodicPFC_measures(s).pain(c).ratings = ratings(:)';
            AperiodicPFC_measures(s).pain(c).trials = AperiodicPFC_measures(s).LEP(c).trials;

        else
            % first check if the missing trials can be explained by removed
            % epochs, if not, provide missing trials manually
            fprintf('subject %d - %s %s: missing events (%d)\n', s, AperiodicPFC_measures(s).conditions(c).area, ...
                AperiodicPFC_measures(s).conditions(c).side, AperiodicPFC_measures(s).LEP(c).trials)
            removed = {};
            for b = 1:length(NLEP_info.single_subject(s).preprocessing.ERP(5).params)
                if contains(NLEP_info.single_subject(s).preprocessing.ERP(5).params(b).dataset, ...
                        sprintf('LEP %s %s', AperiodicPFC_measures(s).conditions(c).area, AperiodicPFC_measures(s).conditions(c).side))
                    block = str2double(regexp(NLEP_info.single_subject(s).preprocessing.ERP(5).params(b).dataset, '\d+', 'match', 'once'));
                    if ~isempty(NLEP_info.single_subject(s).preprocessing.ERP(5).params(b).discarded)                        
                        removed{block} = NLEP_info.single_subject(s).preprocessing.ERP(5).params(b).discarded;
                    else
                       removed{block} = []; 
                    end
                end
            end
            if ~isempty(removed)
                % identify indexes of removed trials
                removed_idx = [removed{1}, removed{2} + 30];

                % filter data if trial numbers match
                if AperiodicPFC_measures(s).LEP(c).trials == 60 - length(removed_idx)
                    AperiodicPFC_measures(s).pain(c).condition = c;
                    AperiodicPFC_measures(s).pain(c).ratings = ratings(:)';
                    AperiodicPFC_measures(s).pain(c).ratings(removed_idx) = [];
                    AperiodicPFC_measures(s).pain(c).trials = AperiodicPFC_measures(s).LEP(c).trials;
                    continue

                else                 
                    % ask for manual input
                    fprintf('LEP trials do not match rating trials. Please enter missing trials manually.\n')
                    prompt = {'b1:', 'b2:'};
                    dlgtitle = sprintf('subject %d - %s %s', s, AperiodicPFC_measures(s).conditions(c).area, ...
                        AperiodicPFC_measures(s).conditions(c).side);
                    dims = [1 50];
                    definput = {'' ''};
                    input = inputdlg(prompt,dlgtitle,dims,definput);

                end
            else
                % ask for manual input
                fprintf('LEP trials do not match rating trials. Please enter missing trials manually.\n')
                prompt = {'b1:', 'b2:'};
                dlgtitle = sprintf('subject %d - %s %s', s, AperiodicPFC_measures(s).conditions(c).area, ...
                    AperiodicPFC_measures(s).conditions(c).side);
                dims = [1 50];
                definput = {'' ''};
                input = inputdlg(prompt,dlgtitle,dims,definput);
            end

            % filter ratings according to manually input trial info
            removed_idx = [str2num(input{1}), str2num(input{2}) + 30];
            AperiodicPFC_measures(s).pain(c).condition = c;
            AperiodicPFC_measures(s).pain(c).ratings = ratings(:)';
            AperiodicPFC_measures(s).pain(c).ratings(removed_idx) = [];
            AperiodicPFC_measures(s).pain(c).trials = AperiodicPFC_measures(s).LEP(c).trials;            
        end
    end
end
fprintf('done.\n')
save(output_file, 'AperiodicPFC_measures', '-append');

% prepare PSD data for sensor-based analysis
for s = 1:params.subjects
    % select source dataset
    if s <= 35
        data = NLEP_data_1to35;
    else
        data = NLEP_data;
    end

    % subject & session info
    AperiodicPFC_data(s).ID = NLEP_info.single_subject(s).ID;
    AperiodicPFC_data(s).conditions = AperiodicPFC_measures(s).conditions;

    % extract PSD
    for c = 1:2
        % subset data associated with this condition
        psd.original = []; psd.oscillatory = []; psd.fractal = []; 
        idx = false(1, length(data.RSEEG(s).dataset));
        for b = 1:length(data.RSEEG(s).dataset)
            if contains(data.RSEEG(s).dataset{b}, AperiodicPFC_data(s).conditions(c).area) && ...
                    contains(data.RSEEG(s).dataset{b}, AperiodicPFC_data(s).conditions(c).side) && ...
                    contains(data.RSEEG(s).dataset{b}, params.prestim_time)
                idx(b) = true;
            end
        end
        for a = 1:length(data.RSEEG(s).PSD_st)
            if idx(a)
                for b = fieldnames(data.RSEEG(s).PSD_st)'
                    psd.(b{1})(end + 1:end + size(data.RSEEG(s).PSD_st(a).(b{1}), 1), :, :) = data.RSEEG(s).PSD_st(a).(b{1});
                end
            end
        end

        % check if trial numbers match, encode
        if size(psd.original, 1) == AperiodicPFC_measures(s).LEP(c).trials
            AperiodicPFC_data(s).PSD(c).condition = c;
            AperiodicPFC_data(s).PSD(c).trials = AperiodicPFC_measures(s).LEP(c).trials;
            AperiodicPFC_data(s).PSD(c).freq = data.RSEEG(s).freq;
            AperiodicPFC_data(s).PSD(c).original = psd.original;
            AperiodicPFC_data(s).PSD(c).oscillatory = psd.oscillatory;
            AperiodicPFC_data(s).PSD(c).fractal = psd.fractal;
        else
            error('ERROR: sibject %d - %s %s - the reial numbers do not match!', ...
                s, AperiodicPFC_data(s).conditions(c).area, AperiodicPFC_data(s).conditions(c).side)
        end
    end
end
fprintf('done.\n')
save(output_file, 'AperiodicPFC_data', '-append');

% prepare pre-processed data for source-based analysis
for s = 1:params.subjects
    for c = 1:2
        % identify and load datasets
        files2load = dir(sprintf('%s\\NLEP_%s\\%s %s %s %s LEP %s %s*', folder.input, AperiodicPFC_data(s).ID, params.prefix{1}, ...
            params.prestim_time, params.prefix{2}, AperiodicPFC_data(s).ID, AperiodicPFC_data(s).conditions(c).area, AperiodicPFC_data(s).conditions(c).side));
        files_idx = false(1, length(files2load));
        for f = 1:length(files2load)
            % identify block
            block = regexp(files2load(f).name, ' b(\d+)', 'tokens', 'once');
            block = str2num(block{1});

            % load the header
            if contains(files2load(f).name, '.lw6')
                load(sprintf('%s\\%s', files2load(f).folder, files2load(f).name), '-mat')
                dataset(block).header = header;
            end 

            % load the data
            if contains(files2load(f).name, '.mat')
                load(sprintf('%s\\%s', files2load(f).folder, files2load(f).name))
                dataset(block).data = data;
            end            
        end

        % save recording parameters
        AperiodicPFC_data(s).EEG_ready(c).condition = c;
        AperiodicPFC_data(s).EEG_ready(c).SR = 1 / dataset(1).header.xstep;
        AperiodicPFC_data(s).EEG_ready(c).times = -1:dataset(1).header.xstep:-dataset(1).header.xstep;

        % check if number of trials match, append
        data = cat(1, squeeze(dataset(1).data), squeeze(dataset(2).data));
        if size(data, 1) == AperiodicPFC_measures(s).pain(c).trials
            AperiodicPFC_data(s).EEG_ready(c).trials = size(data, 1);
            AperiodicPFC_data(s).EEG_ready(c).data = data;

        else
            % first check if the missing trials can be explained by removed
            % epochs, if not, provide missing trials manually
            fprintf('subject %d - %s %s: wrong number of events found (%d)\n', s, AperiodicPFC_data(s).conditions(c).area, ...
                AperiodicPFC_data(s).conditions(c).side, size(data, 1))
            removed = {};
            for b = 1:length(NLEP_info.single_subject(s).preprocessing.ERP(5).params)
                if contains(NLEP_info.single_subject(s).preprocessing.ERP(5).params(b).dataset, ...
                        sprintf('LEP %s %s', AperiodicPFC_measures(s).conditions(c).area, AperiodicPFC_measures(s).conditions(c).side))
                    block = str2double(regexp(NLEP_info.single_subject(s).preprocessing.ERP(5).params(b).dataset, '\d+', 'match', 'once'));
                    if ~isempty(NLEP_info.single_subject(s).preprocessing.ERP(5).params(b).discarded)                        
                        removed{block} = NLEP_info.single_subject(s).preprocessing.ERP(5).params(b).discarded;
                    else
                       removed{block} = []; 
                    end
                end
            end
            if ~isempty(removed)
                % identify indexes of removed trials
                removed_idx = [removed{1}, removed{2} + 30];

                % filter data if trial numbers match
                if size(data, 1) == 60 - length(removed_idx)
                    data(removed_idx, :, :) = [];
                    AperiodicPFC_data(s).EEG_ready(c).trials = size(data, 1);
                    AperiodicPFC_data(s).EEG_ready(c).data = data;               
                    continue

                else                 
                    % ask for manual input
                    fprintf('EEG trials do not match LEP trials. Please enter missing trials manually.\n')
                    prompt = {'b1:', 'b2:'};
                    dlgtitle = sprintf('subject %d - %s %s', s, AperiodicPFC_data(s).conditions(c).area, ...
                        AperiodicPFC_data(s).conditions(c).side);
                    dims = [1 50];
                    definput = {'' ''};
                    input = inputdlg(prompt,dlgtitle,dims,definput);

                end
            else
                % ask for manual input
                fprintf('EEG trials do not match LEP trials. Please enter missing trials manually.\n')
                prompt = {'b1:', 'b2:'};
                dlgtitle = sprintf('subject %d - %s %s', s, AperiodicPFC_data(s).conditions(c).area, ...
                    AperiodicPFC_data(s).conditions(c).side);
                dims = [1 50];
                definput = {'' ''};
                input = inputdlg(prompt,dlgtitle,dims,definput);
            end

            % filter ratings according to manually input trial info
            removed_idx = [str2num(input{1}), str2num(input{2}) + 30];
            data(removed_idx, :, :) = [];
            AperiodicPFC_data(s).EEG_ready(c).trials = size(data, 1);
            AperiodicPFC_data(s).EEG_ready(c).data = data;           
        end
    end
end
fprintf('done.\n')
save(output_file, 'AperiodicPFC_data', '-append');

% clean and continue
clear a b c f s condition missing idx ratings removed removed_idx prompt dlgtitle dims definput input ...
    block data psd trials files2load files_idx data header dataset ...
    NLEP_data_1to35 NLEP_data NLEP_info NLEP_measures
fprintf('section finished.\n\n')

%% sensor-based analysis: extract aperiodic exponent
% ----- section input -----
% -------------------------



% clean and continue
clear 
fprintf('section finished.\n\n')

%% sensor-based analysis: visualization 
% ----- section input -----
% -------------------------



% clean and continue
clear 
fprintf('section finished.\n\n')

%% sensor-based analysis: export for R
% ----- section input -----
% -------------------------



% clean and continue
clear 
fprintf('section finished.\n\n')


%%
% ----- section input -----
param.prefix_data = {'icfilt ica_all chunked' 'icfilt ica_all RS'};
param.cond_continuous = {'open' 'close'; 'pre' 'post'};
param.cond_ST = {'relaxed' 'ready'};
param.head_model = 'BEM';
% ------------------------- 

% set directories and load info structure
if ~exist("folder")
    folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave masterfiles
    folder.input = uigetdir(pwd, 'Coose the input folder');             % raw data
    folder.data = uigetdir(pwd, 'Choose the data folder');              % processed data
    folder.output = uigetdir(pwd, 'Choose the OneDrive folder');        % output folder --> figures, logfiles, output .mat file
    study = 'NLEP';
    output_file = sprintf('%s\\%s_output.mat', folder.output, study);
end
cd(folder.data)
load(output_file, 'NLEP_info');

% ask for subject number
if ~exist('subject_idx')
    prompt = {'subject number:'};
    dlgtitle = 'subject';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    subject_idx = str2num(input{1,1});
end
clear prompt dlgtitle dims definput input

% add FieldTrip to the top of search path
addpath(genpath([folder.toolbox '\FieldTrip']));

% load template electrode positions
chanlocs_template = ft_read_sens(sprintf('%s\\FieldTrip\\template\\electrode\\standard_1020.elc', folder.toolbox));
chanlocs_template = ft_determine_coordsys(chanlocs_template, 'interactive', 'yes');

% load template headmodel
headmodel = ft_read_headmodel(sprintf('%s\\FieldTrip\\template\\headmodel\\standard_bem.mat', folder.toolbox));
ft_plot_headmodel(headmodel, 'facealpha', 0.6)

% prepare template mri
mri = ft_read_mri(sprintf('%s\\FieldTrip\\template\\headmodel\\standard_mri.mat', folder.toolbox));

% prepare cortical atlas
atlas = ft_read_atlas(sprintf('%s\\FieldTrip\\template\\atlas\\aal\\ROI_MNI_V4.nii', folder.toolbox));
cortical_rois = {'Precentral_L', 'Precentral_R', 'Frontal_Sup_L', 'Frontal_Sup_R', ...
    'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R', 'Frontal_Mid_L', 'Frontal_Mid_R', ...
    'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', ...
    'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', ...
    'Rolandic_Oper_L', 'Rolandic_Oper_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', ...
    'Olfactory_L', 'Olfactory_R', 'Postcentral_L', 'Postcentral_R', 'Parietal_Sup_L', ...
    'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'Supramarginal_L', 'Supramarginal_R', ...
    'Angular_L', 'Angular_R', 'Precuneus_L', 'Precuneus_R', 'Temporal_Sup_L', 'Temporal_Sup_R', ...
    'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', ...
    'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Inf_L', 'Temporal_Inf_R', ...
    'Occipital_Sup_L', 'Occipital_Sup_R', 'Occipital_Mid_L', 'Occipital_Mid_R', ...
    'Occipital_Inf_L', 'Occipital_Inf_R', 'Cuneus_L', 'Cuneus_R', 'Calcarine_L', 'Calcarine_R', ...
    'Lingual_L', 'Lingual_R', 'Insula_L', 'Insula_R'};
roi_indices = find(ismember(atlas.tissuelabel, cortical_rois));
centroid_positions = zeros(length(roi_indices), 3);
for i = 1:length(roi_indices)
    roi_labels{i} = atlas.tissuelabel{i};
    roi_mask = atlas.tissue == roi_indices(i);
    [x, y, z] = ind2sub(size(atlas.tissue), find(roi_mask));
    centroid_positions(i, :) = mean([x, y, z], 1);
end

% create a custom source model grid
source_model = [];
source_model.pos = centroid_positions;  
source_model.inside = true(length(centroid_positions), 1); 
source_model.unit = 'mm';  


figure;
ft_plot_headmodel(headmodel, 'facealpha', 0.6)
hold on;
ft_plot_mesh(source_model.pos, 'vertexcolor', 'r', 'vertexsize', 10);  

rois = {'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Heschl_L', 'Heschl_R'};
cfg = [];
cfg.atlas = atlas;
cfg.roi = rois;
cfg.inputcoord = 'mni';  

leadfield = load(sprintf('%s\\FieldTrip\\template\\sourcemodel\\standard_sourcemodel3d10mm.mat', folder.toolbox));
% figure;
% ft_plot_mesh(leadfield.sourcemodel.pos(leadfield.sourcemodel.inside,:), 'vertexsize', 20);

% create a source model based on the atlas ROIs

mask = ft_volumelookup(cfg, mri);

% Visualize the ROIs
figure;
ft_plot_vol(atlas);  % Plot the atlas
hold on;
ft_plot_mesh(roi_grid.pos, 'vertexsize', 20);  % Plot the ROI positions
title('ROIs based on AAL atlas');

% identify RS-EEG data 
files2process = dir(sprintf('%s\\%s_%s\\%s*.mat', folder.input, study, NLEP_info.single_subject(subject_idx).ID, param.prefix_data{1}));
files2process = [files2process; dir(sprintf('%s\\%s_%s\\%s*.mat', folder.input, study, NLEP_info.single_subject(subject_idx).ID, param.prefix_data{2}))];

% concatenate datasets into a FieldTrip data structure 
for f = 1:length(files2process)
    % prepare data structure
    if ~exist('data_all')
        data_all = struct;
        data_all.trial = []; data_all.time = []; data_all.trialinfo = [];
    end

    % load letswave data and header
    load(sprintf('%s\\%s', files2process(f).folder, files2process(f).name))
    load(sprintf('%s\\%s.lw6', files2process(f).folder, files2process(f).name(1:end-4)), '-mat')
    data = squeeze(data);   

    % fill in common information
    if ~exist('data_all.label')
        data_all.label = {header.chanlocs.labels};
        data_all.fsample = 1/header.xstep;
        data_all.cfg = [];  % start with an empty config
    end

    % prepare data structure
    data.trial = cell(1, size(data, 1));
    data.time = cell(1, size(data, 1));
    data.trialinfo = ones(length(data.trial), 1)*f;
    
    % fill in data
    for i = 1:size(data, 1)
        data.trial{i} = squeeze(data(i, :, :));
        data.time{i} = (1:size(data, 3))/data_all.fsample;
    end

    % concatenate
    data_all.trial = [data_all.trial, data.trial];  
    data_all.time = [data_all.time, data.time]; 
    data_all.trialinfo = [data_all.trialinfo; data.trialinfo]; 
end
data = data_all;
ft_checkdata(data);

% extract channel locations
chanlocs = chanlocs_template;
for a = 1:length(chanlocs.label)
    if ismember(chanlocs.label{a}, data.label)
        chan_idx(a) = true;
    else
        chan_idx(a) = false;
    end
end
chanlocs.label = chanlocs.label(chan_idx);
chanlocs.chantype = chanlocs.chantype(chan_idx);
chanlocs.chanunit = chanlocs.chanunit(chan_idx);
chanlocs.chanpos = chanlocs.chanpos(chan_idx, :);
chanlocs.elecpos = chanlocs.elecpos(chan_idx, :);
[~,idx] = ismember(data.label', chanlocs.label);
chanlocs.label = chanlocs.label(idx);
chanlocs.chantype = chanlocs.chantype(idx);
chanlocs.chanunit = chanlocs.chanunit(idx);
chanlocs.chanpos = chanlocs.chanpos(idx, :);
chanlocs.elecpos = chanlocs.elecpos(idx, :);

% add fiducials
fiducials = {'LPA', 'RPA', 'Nz'};
for f = 1:length(fiducials)
    [~,idx] = ismember(fiducials(f), chanlocs_template.label);
    chanlocs.label(end+1) = chanlocs_template.label(idx);
    chanlocs.chantype(end+1) = chanlocs_template.chantype(idx);
    chanlocs.chanunit(end+1) = chanlocs_template.chanunit(idx);
    chanlocs.chanpos(end+1,:) = chanlocs_template.chanpos(idx, :);
    chanlocs.elecpos(end+1,:) = chanlocs_template.elecpos(idx, :);
end

% align the electrodes with the template MRI
cfg = [];
cfg.method = 'interactive';
cfg.elec = data.elec;
cfg.headshape = headmodel.bnd(1); 
data.elec = ft_electroderealign(cfg);

% check data
data.elec = chanlocs;
% ft_plot_sens(data.elec, 'style', 'k*', 'label', 'label');
ft_checkdata(data);

% align electrode positions
cfg = [];
cfg.method = 'interactive';
cfg.elec = data.elec;
cfg.headshape = headmodel.bnd(1);
cfg.template = chanlocs_template;
cfg.channel = 'all';
data.elec = ft_electroderealign(cfg);
ft_checkdata(data);

% loop through trials
source_trials = cell(1, length(data.trial));
for i = 1:length(data.trial)
    % extract the i-th trial data
    cfg = [];
    cfg.trials = i;
    single_trial = ft_selectdata(cfg, data);

    % compute the covariance matrix for the single trial
    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = 'all';  % Use the whole time window for covariance computation
    timelock = ft_timelockanalysis(cfg, single_trial);

    % perform source analysis for the single trial
    cfg = [];
    cfg.method = 'lcmv';
    cfg.grid = leadfield;            % Use the prepared lead field
    cfg.headmodel = headmodel;       % Use the head model (BEM)
    cfg.lcmv.keepfilter = 'yes';     % Keep the spatial filter
    cfg.lcmv.lambda = '5%';          % Regularization parameter
    cfg.lcmv.fixedori = 'yes';       % Fixed orientation of dipoles
    cfg.lcmv.projectnoise = 'yes';   % Project noise for visualization
    source = ft_sourceanalysis(cfg, timelock);

    % store the result for the i-th trial
    source_trials{i} = source;
end