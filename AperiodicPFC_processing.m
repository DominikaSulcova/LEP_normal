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

%% final data check
% check for extra/missing events --> generate .txt report
name = sprintf('%s\\NLEP_event_analysis.txt', folder.output);
fileID = fopen(name, 'a');
data = NLEP_data;
for a = 1:length(data.LEP)
    if ~isempty(data.LEP(a).ID)
        fprintf(fileID, sprintf('subject %d (%s) - LEP\r\n', a, data.LEP(a).ID));
        % unfiltered LEP data
        for b = fieldnames(data.LEP(a).unfiltered)'
            for c = {'cond1' 'cond2'}
                if size(data.LEP(a).unfiltered.(b{1}).(c{1}), 1) > 60
                    fprintf(fileID, sprintf('     --> unfiltered %s %s: - too many events - %d\r\n', ...
                        b{1}, c{1}, size(data.LEP(a).unfiltered.(b{1}).(c{1}), 1)));
                elseif size(data.LEP(a).unfiltered.(b{1}).(c{1}), 1) < 60
                    fprintf(fileID, sprintf('     --> unfiltered %s %s: - missing events - %d\r\n', ...
                        b{1}, c{1}, size(data.LEP(a).unfiltered.(b{1}).(c{1}), 1)));
                end
            end
        end
        % CWT-filtered LEP data
        for b = fieldnames(data.LEP(a).CWT_filtered)'
            for c = {'cond1' 'cond2'}
                if size(data.LEP(a).CWT_filtered.(b{1}).(c{1}), 1) > 60
                    fprintf(fileID, sprintf('     --> CWT_filtered %s %s: - too many events - %d\r\n', ...
                        b{1}, c{1}, size(data.LEP(a).CWT_filtered.(b{1}).(c{1}), 1)));
                elseif size(data.LEP(a).CWT_filtered.(b{1}).(c{1}), 1) < 60
                    fprintf(fileID, sprintf('     --> CWT_filtered %s %s: - missing events - %d\r\n', ...
                        b{1}, c{1}, size(data.LEP(a).CWT_filtered.(b{1}).(c{1}), 1)));
                end
            end
        end
        fprintf(fileID, sprintf('subject %d (%s) - RSEEG\r\n', a, data.LEP(a).ID));
        % PSD data
        for c = 1:length(data.RSEEG(a).PSD_st)
            if contains(data.RSEEG(a).dataset{c}, 'RS')
                if size(data.RSEEG(a).PSD_st(c).original, 1) > 59
                    fprintf(fileID, sprintf('     --> %s: - too many events - %d\r\n', ...
                        data.RSEEG(a).dataset{c}, size(data.RSEEG(a).PSD_st(c).original, 1)));
                elseif size(data.RSEEG(a).PSD_st(c).original, 1) < 59
                    fprintf(fileID, sprintf('     --> %s: - missing events - %d\r\n', ...
                        data.RSEEG(a).dataset{c}, size(data.RSEEG(a).PSD_st(c).original, 1)));
                end
            elseif contains(data.RSEEG(a).dataset{c}, 'LEP')
                if size(data.RSEEG(a).PSD_st(c).original, 1) > 30
                    fprintf(fileID, sprintf('     --> %s: - too many events - %d\r\n', ...
                        data.RSEEG(a).dataset{c}, size(data.RSEEG(a).PSD_st(c).original, 1)));
                elseif size(data.RSEEG(a).PSD_st(c).original, 1) < 30
                    fprintf(fileID, sprintf('     --> %s: - missing events - %d\r\n', ...
                        data.RSEEG(a).dataset{c}, size(data.RSEEG(a).PSD_st(c).original, 1)));
                end
            end
        end
        fprintf(fileID, '\r\n');
    end
end
fclose(fileID);

% remove specific events from specific datasets
a = 18; 
dataset = 'cond1';
% dataset = {'LEP hand left b2 - ready'};
trial = 31;
% if RSEEG
for d = 1:length(dataset)
    c = find(strcmp(NLEP_data_1to35.RSEEG(a).dataset, dataset{d}));
    for b = fieldnames(NLEP_data_1to35.RSEEG(a).PSD_st)'
        NLEP_data_1to35.RSEEG(a).PSD_st(c).(b{1})(trial, :, :) = [];
    end
end
% if LEP
for b = {'unfiltered' 'CWT_filtered'}   
    for c = fieldnames(NLEP_data_1to35.LEP(a).(b{1}))'
        NLEP_data_1to35.LEP(a).(b{1}).(c{1}).(dataset)(trial, :) = [];
    end
end
NLEP_data_1to35.LEP(a).blocks{str2double(regexp(dataset, '\d+', 'match', 'once'))}(trial) = [];
save(input_file, 'NLEP_info', 'NLEP_data_1to35', '-append')

% remove specific measures
s = 17; 
c = 1;
trial = 32;
for a = fieldnames(NLEP_measures(s).LEP_ST)'
    if ~strcmp(a{1}, 'conditions') 
        for b = 1:3
            NLEP_measures(s).LEP_ST.(a{1})(b, c, trial) = 0;
            idx = NLEP_measures(s).LEP_ST.(a{1})(b, c, :) ~= 0;
            data = NLEP_measures(s).LEP_ST.(a{1})(b, c, idx);
            data(1, 1, end+1:60) = 0;
            NLEP_measures(s).LEP_ST.(a{1})(b, c, :) = data;
        end
    end
end
save(input_file, 'NLEP_measures', '-append')

% check data in letswave
addpath(genpath([folder.toolbox '\letswave 7']));
letswave

%% import data and existing measures 

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