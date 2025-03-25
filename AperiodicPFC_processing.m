%% params
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave masterfiles
folder.input = uigetdir(pwd, 'Coose the input folder');             % raw data
folder.data = uigetdir(pwd, 'Coose the data folder');               % processed data
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

%% sensor-space analysis: import existing measures and pre-processed data
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

% clean and continue
clear a b c f s condition missing idx ratings removed removed_idx prompt dlgtitle dims definput input ...
    block data psd trials files2load files_idx data header dataset ...
    NLEP_data_1to35 NLEP_data NLEP_info NLEP_measures
fprintf('section finished.\n\n')

%% sensor-space analysis: extract aperiodic exponent
% ----- section input -----
params.eoi_target = {'AF3' 'AFz' 'AF4' 'F3' 'F1' 'F2' 'F4'};
params.eoi_ctrl = {'PO3' 'POz' 'PO4' 'P3' 'P1' 'P2' 'P4'};
params.foi = {[5, 50] [5 30] [30 50]};
params.foi_labels = {'broad' 'low' 'high'};
% -------------------------

% reload output structures if necessary
if exist('AperiodicPFC_data') ~= 1 || exist('AperiodicPFC_measures') ~= 1
    load(output_file, 'AperiodicPFC_data', 'AperiodicPFC_measures')
end

% prepare flip dictionary
load('dataset_default.lw6', '-mat')
params.labels = {header.chanlocs.labels};
params.chanlocs = header.chanlocs;
labels_flipped = params.labels;
for i = 1:length(params.labels)
    electrode_n = str2num(params.labels{i}(end));
    if isempty(electrode_n)
    else
        if electrode_n == 0
            label_new = params.labels{i}(1:end-2);
            label_new = [label_new num2str(9)];
        elseif mod(electrode_n,2) == 1              % odd number --> left hemisphere                    
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n + 1)];
        else                                    % even number --> right hemisphere 
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n - 1)];
        end
        a = find(strcmpi(params.labels, label_new));
        if isempty(a)
        else
            labels_flipped{i} = label_new;
        end
    end
end
labels_dict = cat(1, params.labels, labels_flipped)';

% homogenize the fractal PSD data --> flip as if all stimuli were delivered 
% to the right hand 
addpath(genpath([folder.toolbox '\letswave 6']));
row_counter = 1;
for s = 1:params.subjects
    for c = 1:2
        % select the data
        data = []; 
        data(:, :, 1, 1, 1, :) = AperiodicPFC_data(s).PSD(c).fractal;  

        % flip if left side stimulated
        flipped = false;
        if strcmp(AperiodicPFC_data(s).conditions(c).side, 'left')
           [header, data, ~] = RLW_flip_electrodes(header, data, labels_dict);
           flipped = true;
        end

        % save to the new dataset
        dataset(row_counter).subject = s;
        dataset(row_counter).ID = AperiodicPFC_data(s).ID;
        dataset(row_counter).freq = AperiodicPFC_data(s).PSD(c).freq;
        dataset(row_counter).area = AperiodicPFC_data(s).conditions(c).area;
        dataset(row_counter).flipped = flipped;
        dataset(row_counter).data = squeeze(data);

        % update row counter
        row_counter = row_counter + 1;
    end
end

% extract aperiodic measures for each subject/area/trial/electrode/foi
for d = 1:length(dataset)
    % provide update
    fprintf('extracting aperiodic measures: subject %d - %s\n', dataset(d).subject, dataset(d).area)

    % encode doi
    for f = 1:length(params.foi)
        dataset(d).foi(f).label = params.foi_labels{f};
        dataset(d).foi(f).limits = params.foi{f};
    end
    
    % fit PSD
    for a = 1:size(dataset(d).data, 1)
        for b = 1:size(dataset(d).data, 2)
            for c = 1:length(params.foi)
                % subset the data
                data = squeeze(dataset(d).data(a, b, :))';
                freq_idx = dataset(d).freq >= dataset(d).foi(c).limits(1) & dataset(d).freq <= dataset(d).foi(c).limits(2);
                roi.freq = dataset(d).freq(freq_idx);
                roi.psd = data(freq_idx);

                % log-transform 
                roi.freq_log = log10(roi.freq);
                roi.psd_log  = log10(roi.psd);

                % perform standard linear regression
                fit.polyfit = polyfit(roi.freq_log, roi.psd_log, 1);

                % perform robust linear regression
                [fit.robust, ~] = robustfit(roi.freq_log, roi.psd_log);
                fit.robust = fit.robust([2, 1])';

                % save to the dataset
                dataset(d).exponent_polyfit(a, b, c) = -fit.polyfit(1);
                dataset(d).offset_polyfit(a, b, c) = fit.polyfit(2);
                dataset(d).exponent_robust(a, b, c) = -fit.robust(1);
                dataset(d).offset_robust(a, b, c) = fit.robust(2);
            end
        end
    end

    % output mean differences between polyfit and rubust regression at each
    % electrode
    screen_size = get(0, 'ScreenSize');
    figure(figure_counter)
    set(gcf, 'Position', [1, screen_size(4)/4, screen_size(3), screen_size(4) / 2])
    for f = 1:length(params.foi)
        % calculate mean difference per electrode
        exp_diff.mean = [];
        exp_diff.sd = [];
        for b = 1:size(dataset(d).exponent_polyfit, 2)
            % extract differences at each trial
            diff_trial = [];
            for a = 1:size(dataset(d).exponent_polyfit, 1)
                diff_trial(a) = dataset(d).exponent_robust(a, b, f) - dataset(d).exponent_polyfit(a, b, f);
            end

            % calculate mean and SD
            exp_diff.mean(b) = mean(diff_trial);
            exp_diff.sd(b) = std(diff_trial);
        end    

        % plot
        subplot(1, 3, f)
        bar(exp_diff.mean);
        hold on;
        errorbar(1:length(exp_diff.mean), exp_diff.mean, exp_diff.sd, 'k.', 'LineWidth', 0.8);
        set(gca, 'XTick', 1:length(exp_diff.mean));
        set(gca, 'XTickLabel', params.labels);
        xlabel('electrodes');
        ylabel('mean robust - polyfit');
        title(sprintf('%s frequencies ([%d %d]Hz)', dataset(d).foi(f).label, ...
            dataset(d).foi(f).limits(1), dataset(d).foi(f).limits(2)));
    end

    % save figure and update counter
    saveas(gcf, sprintf('%s\\figures\\exp_diff_%s_%s.png', folder.output, dataset(d).ID, dataset(d).area))
    figure_counter = figure_counter + 1;
end
AperiodicPFC_APC = dataset;
save(output_file, 'AperiodicPFC_APC', '-append');

% extract aperiodic exponents and slopes from target/ctrl electrodes
target_idx = find(ismember(params.labels, params.eoi_target));
ctrl_idx = find(ismember(params.labels, params.eoi_ctrl));
for d = 1:length(dataset)
    for a = 1:size(dataset(d).data, 1)
        % subset data and average across rois
        data.target = squeeze(mean(dataset(d).data(a, target_idx, :), 2))'; 
        data.ctrl = squeeze(mean(dataset(d).data(a, ctrl_idx, :), 2))'; 

        % perform robust linear regression in all frequency bands
        exponent.target = []; offset.target = [];
        exponent.ctrl = []; offset.ctrl = [];
        for c = 1:length(params.foi)
            % identify frequencies
            freq_idx = dataset(d).freq >= dataset(d).foi(c).limits(1) & dataset(d).freq <= dataset(d).foi(c).limits(2);
            freq = dataset(d).freq(freq_idx);

            % log-transform
            roi.psd_target = log10(data.target(freq_idx));
            roi.psd_ctrl = log10(data.ctrl(freq_idx));
            roi.freq = log10(freq);

            % fit for target roi
            [fit, ~] = robustfit(roi.freq, roi.psd_target);
            exponent.target(c) = -fit(2);
            offset.target(c) = fit(1);

            % fit for ctrl roi
            [fit, ~] = robustfit(roi.freq, roi.psd_ctrl);
            exponent.ctrl(c) = -fit(2);
            offset.ctrl(c) = fit(1);           
        end

        % save to output structure
        dataset(d).exponent_target(a, :) = exponent.target;
        dataset(d).exponent_ctrl(a, :) = exponent.ctrl;
        dataset(d).offset_target(a, :) = offset.target;
        dataset(d).offset_ctrl(a, :) = offset.ctrl;
    end
end
AperiodicPFC_APC = dataset;
save(output_file, 'AperiodicPFC_APC', '-append');

% save to measures structure as well
for d = 1:length(dataset)
    % identify subject and condition
    s = dataset(d).subject;
    area = dataset(d).area;
    if dataset(d).flipped
        side = 'left';
    else
        side = 'right';
    end
    for a = 1:length(AperiodicPFC_measures(s).conditions)
        if strcmp(AperiodicPFC_measures(s).conditions(a).area, area) && ...
             strcmp(AperiodicPFC_measures(s).conditions(a).side, side)  
            c = AperiodicPFC_measures(s).conditions(a).condition;
        end
    end

    % encode
    AperiodicPFC_measures(s).APC(c).condition = c;
    AperiodicPFC_measures(s).APC(c).flipped = dataset(d).flipped;
    AperiodicPFC_measures(s).APC(c).foi = dataset(d).foi;
    AperiodicPFC_measures(s).APC(c).exponent_target = dataset(d).exponent_target;
    AperiodicPFC_measures(s).APC(c).exponent_ctrl = dataset(d).exponent_ctrl;
    AperiodicPFC_measures(s).APC(c).offset_target = dataset(d).offset_target;
    AperiodicPFC_measures(s).APC(c).offset_ctrl = dataset(d).offset_ctrl;
    AperiodicPFC_measures(s).APC(c).trials = size(dataset(d).exponent_target, 1);

    % throw error if trial numbers don't match
    if AperiodicPFC_measures(s).APC(c).trials ~= AperiodicPFC_measures(s).LEP(c).trials  
        error('ERROR: subject %d, condition %d - trial numbers do not match!', s, c)
    end
end
save(output_file, 'AperiodicPFC_measures', '-append');

% clean and continue
clear a b c d f i s header data labels_flipped electrode_n label_new labels_dict row_counter side area ...
    row_counter dataset freq_idx roi fit exp_diff diff_trial screen_size flipped freq exponent offset ctrl_idx target_idx
fprintf('section finished.\n\n')

%% sensor-space analysis: export for R
% ----- section input -----
params.foi = {[5, 50] [5 30] [30 50]};
params.foi_labels = {'broad' 'low' 'high'};
params.region = {'target' 'ctrl'};
% -------------------------

% reload output structures if necessary 
if exist('AperiodicPFC_data') ~= 1 || exist('AperiodicPFC_measures') ~= 1 || exist('AperiodicPFC_APC') ~= 1
    load(output_file, 'AperiodicPFC_data', 'AperiodicPFC_measures', 'AperiodicPFC_APC')
end

% export measured variables in a long-format table
for f = 1:length(params.foi)
    fprintf('exporting %s frequency data:\n', params.foi_labels{f})
    table_export = table;
    row_counter = 1;
    for s = 1:params.subjects
        for c = 1:length(AperiodicPFC_measures(s).conditions)
            for a = 1:AperiodicPFC_measures(s).LEP(c).trials
                for b = 1:length(params.LEP_comps)
                    for r = 1:length(params.region)
                        % subject info
                        table_export.subject(row_counter) = s;
                        table_export.ID{row_counter} = AperiodicPFC_measures(s).ID;
                        table_export.age(row_counter) = AperiodicPFC_measures(s).age;
                        table_export.male(row_counter) = AperiodicPFC_measures(s).male;
                        table_export.handedness{row_counter} = AperiodicPFC_measures(s).handedness;
        
                        % session info
                        table_export.condition(row_counter) = c;
                        table_export.area{row_counter} = AperiodicPFC_measures(s).conditions(c).area;
                        table_export.side{row_counter} = AperiodicPFC_measures(s).conditions(c).side;
                        if strcmp(table_export.handedness{row_counter},table_export.side{row_counter})
                            table_export.dominant(row_counter) = 1;
                        else
                            table_export.dominant(row_counter) = 0;
                        end
                        table_export.flipped{row_counter} = AperiodicPFC_measures(s).APC(c).flipped;
    
                        % dependent variable: LEP measures
                        table_export.component{row_counter} = params.LEP_comps{b};
                        table_export.amplitude(row_counter) = AperiodicPFC_measures(s).LEP(c).(params.LEP_comps{b}).amplitude(a);
                        table_export.latency(row_counter) = AperiodicPFC_measures(s).LEP(c).(params.LEP_comps{b}).latency(a);
    
                        % dependent variable: pain rating
                        table_export.rating(row_counter) = AperiodicPFC_measures(s).pain(c).ratings(a);
        
                        % independent variable: aperiodic measures
                        table_export.foi{row_counter} = AperiodicPFC_measures(s).APC(c).foi(f).label;  
                        table_export.region{row_counter} = params.region{r}; 
                        if strcmp(params.region{r}, 'target')
                            table_export.exponent(row_counter) = AperiodicPFC_measures(s).APC(c).exponent_target(a, f);
                            table_export.offset(row_counter) = AperiodicPFC_measures(s).APC(c).offset_target(a, f);  
                        elseif strcmp(params.region{r}, 'ctrl')
                            table_export.exponent(row_counter) = AperiodicPFC_measures(s).APC(c).exponent_ctrl(a, f);
                            table_export.offset(row_counter) = AperiodicPFC_measures(s).APC(c).offset_ctrl(a, f); 
                        end
        
                        % update row counter
                        row_counter = row_counter + 1;
                    end
                end
            end
        end
    end

    % save table to output structure and as .csv
    fprintf('saving...')
    statement = sprintf('AperiodicPFC_export_%s = table_export;', params.foi_labels{f});
    eval(statement)
    statement = sprintf('save(output_file, ''AperiodicPFC_export_%s'', ''-append'');', params.foi_labels{f});
    eval(statement)
    statement = sprintf('writetable(AperiodicPFC_export_%s, ''AperiodicPFC_export_%s.csv'');', params.foi_labels{f}, params.foi_labels{f});
    eval(statement)  
    fprintf('done.\n')
end

% clean and continue
clear a b c f r s row_counter table_export statement
fprintf('section finished.\n\n')

%% sensor-space analysis: visualization 
% ----- section input -----
params.plot_foi = 'high';
params.region = {'target' 'ctrl'};
% -------------------------
% grand average LEPs
%   - import LEP data
%   - C3 and C4 together --> flip?
%   - plot a figure with hand and foot together
% topographical distribution of aperiodic exponent 

% clean and continue
clear 
fprintf('section finished.\n\n')

%% source-space analysis: import pre-processed data
% ----- section input -----
params.prestim_time = 'ready';
params.prefix = {'icfilt ica_all RS' 'ep reref ds notch bandpass dc'};
% -------------------------

% import pre-processed pre-stimulus data 
fprintf('loading data:\n')
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
clear dataset
save(output_file, 'AperiodicPFC_data', '-append');

% prepare flip dictionary
load(sprintf('%s\\dataset_default.lw6', folder.output), '-mat')
params.labels = {header.chanlocs.labels};
params.chanlocs = header.chanlocs;
labels_flipped = params.labels;
for i = 1:length(params.labels)
    electrode_n = str2num(params.labels{i}(end));
    if isempty(electrode_n)
    else
        if electrode_n == 0
            label_new = params.labels{i}(1:end-2);
            label_new = [label_new num2str(9)];
        elseif mod(electrode_n,2) == 1              % odd number --> left hemisphere                    
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n + 1)];
        else                                    % even number --> right hemisphere 
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n - 1)];
        end
        a = find(strcmpi(params.labels, label_new));
        if isempty(a)
        else
            labels_flipped{i} = label_new;
        end
    end
end
labels_dict = cat(1, params.labels, labels_flipped)';

% homogenize the data --> flip as if all stimuli were delivered to the right hand 
addpath(genpath([folder.toolbox '\letswave 6']));
row_counter = 1;
fprintf('flipping data: ')
for s = 1:params.subjects
    fprintf('. ')
    for c = 1:2
        % select the data
        data = []; 
        data(:, :, 1, 1, 1, :) = AperiodicPFC_data(s).EEG_ready(c).data;  

        % flip if left side stimulated
        flipped = false;
        if strcmp(AperiodicPFC_data(s).conditions(c).side, 'left')
           [header, data, ~] = RLW_flip_electrodes(header, data, labels_dict);
           flipped = true;
        end

        % save to the new dataset
        dataset(row_counter).subject = s;
        dataset(row_counter).ID = AperiodicPFC_data(s).ID;
        dataset(row_counter).condition = c;
        dataset(row_counter).area = AperiodicPFC_data(s).conditions(c).area;
        dataset(row_counter).side = AperiodicPFC_data(s).conditions(c).side;
        dataset(row_counter).flipped = flipped;
        dataset(row_counter).SR = AperiodicPFC_data(s).EEG_ready(c).SR;
        dataset(row_counter).times = AperiodicPFC_data(s).EEG_ready(c).times;      
        dataset(row_counter).data = squeeze(data);

        % update row counter
        row_counter = row_counter + 1;
    end
end
AperiodicPFC_EEG_flipped = dataset;
save(output_file, 'AperiodicPFC_EEG_flipped', '-append', '-v7.3'); 
fprintf('done.\n')

% clean and continue
clear a c f i s block files2load files_idx data header labels_dict labels_flipped electrode_n ...
    label_new removed removed_idx prompt dlgtitle dims definput input flipped 
fprintf('section finished.\n\n')

%% source-space analysis: estimate source activity
% ----- section input -----
params.method = 'lcmv';
params.single_trial = true;
params.path_atlas = 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
params.path_headmodel = 'standard_bem.mat';
params.colours = [1.0000    0.0745    0.6510; 
    0.4745    0.3176    0.9098];
params.ratio_threshold = 0.001;
params.source_threshold = 60;
params.trial_threshold = 25;
% -------------------------

% load dataset if necessary
if exist('dataset') ~= 1
    load(output_file, 'AperiodicPFC_EEG_flipped')
    dataset = AperiodicPFC_EEG_flipped;
end

% prepare template-compatible electrode locations
addpath(genpath([folder.toolbox '\FieldTrip']));
load(sprintf('%s\\dataset_default.lw6', folder.output), '-mat')
params.labels = {header.chanlocs.labels};
params.elec_template = ft_read_sens('standard_1005.elc');
elec = create_elec(params, header);

% check electrode positions on the scalp
load(params.path_headmodel,'vol');
figure;
ft_plot_headmodel(vol,'facealpha',0.1,'facecolor',[0.1 0.1 0.1],'edgecolor',[1 1 1],'edgealpha',0.5);
hold on;
ft_plot_sens(elec,'style','r','label','label','elec','true','elecshape','disc','elecsize',5,'facecolor','r');
view(90,0);

% prepare source model --> centroid positons from Schaefer atlas
atlas = readtable(params.path_atlas);
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas.R, atlas.A, atlas.S];
cfg.unit = 'mm';
cfg.headmodel = params.path_headmodel;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

% compute source estimation for all datasets/trials
fprintf('estimating source activity:\n')
for d = 1:length(dataset)
    % provide update
    fprintf('subject %d - condition %d:\n', dataset(d).subject, dataset(d).condition)

    % load data into a FieldTrip-compatible data structure
    data = [];
    data.label = {header.chanlocs.labels}';
    data.fsample = dataset(d).SR;
    data.trial = cell(size(dataset(d).data, 1), 1);
    data.time = cell(size(dataset(d).data, 1), 1);
    for a = 1:size(dataset(d).data, 1)
        data.trial{a} = squeeze(dataset(d).data(a, :, :));
        data.time{a} = squeeze(dataset(d).times);
    end
    data.elec = elec;

    % estimate source activity
    if params.single_trial
        % prepare forward model = leadfield
        cfg = [];
        cfg.sourcemodel = sourcemodel_atlas;
        cfg.headmodel = params.path_headmodel;
        cfg.elec = data.elec; 
        cfg.normalize = 'yes';
        leadfield = ft_prepare_leadfield(cfg, data);

        % prepare spatial filter
        cfg = [];
        cfg.method = params.method;
        cfg.keeptrials = 'yes';
        cfg.rawtrial = 'yes';
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.lambda = '5%';
        cfg.lcmv.fixedori = 'yes';
        cfg.lcmv.projectnoise = 'yes';
        cfg.lcmv.weightnorm = 'arraygain';
        cfg.sourcemodel = leadfield;
        sources = ft_sourceanalysis(cfg, data);

    else
        % normalize time axis of the data 
        temptime = data.time{1};
        [data.time{:}] = deal(temptime);

        % compute the covaraciance matrix from the sensor data
        cfg = [];
        cfg.covariance = 'yes';
        cfg.keeptrials = 'no';
        cfg.removemean = 'yes';
        timelock = ft_timelockanalysis(cfg, data);

        % prepare forward model = leadfield
        cfg = [];
        cfg.sourcemodel = sourcemodel_atlas;
        cfg.headmodel = params.path_headmodel;
        cfg.normalize = 'yes';
        leadfield = ft_prepare_leadfield(cfg, data);

        % prepare spatial filter
        cfg = [];
        cfg.method = params.method;
        cfg.keeptrials = 'yes';
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.lambda = '5%';
        cfg.lcmv.fixedori = 'yes';
        cfg.lcmv.projectnoise = 'yes';
        cfg.lcmv.weightnorm = 'arraygain';
        cfg.sourcemodel = leadfield;
        sources = ft_sourceanalysis(cfg, data);
    end

    % save source activity to the output structure
    AperiodicPFC_sources(d).subject = dataset(d).subject;
    AperiodicPFC_sources(d).ID = dataset(d).ID;
    AperiodicPFC_sources(d).condition = dataset(d).condition;
    AperiodicPFC_sources(d).area = dataset(d).area;
    AperiodicPFC_sources(d).side = dataset(d).side;
    AperiodicPFC_sources(d).flipped = dataset(d).flipped;
    AperiodicPFC_sources(d).SR = dataset(d).SR;
    AperiodicPFC_sources(d).times = dataset(d).times;   
    for b = 1:length(params.labels)
        AperiodicPFC_sources(d).channels(b).label = params.labels{b};   
        AperiodicPFC_sources(d).channels(b).position = elec.chanpos(b, :); 
    end
    AperiodicPFC_sources(d).sources = atlas; 
    AperiodicPFC_sources(d).orientation = sources.trial(1).ori'; 
    AperiodicPFC_sources(d).filter = sources.trial(1).filter; 
    for c = 1:length(sources.trial)
        for s = 1:height(atlas)
            AperiodicPFC_sources(d).data(c, s, :) = sources.trial(c).mom{s}; 
        end
    end
    AperiodicPFC_sources(d).EEG = dataset(d).data;
end
save('AperiodicPFC_sources.mat', 'AperiodicPFC_sources', '-v7.3'); 

% quality check
for s = 1:params.subjects
    % select data of this subject
    data_idx = false(1, length(AperiodicPFC_sources));
    for a = 1:length(AperiodicPFC_sources)
        if AperiodicPFC_sources(a).subject == s
            data_idx(a) = true;
        end
    end
    data = AperiodicPFC_sources(data_idx);

    % plot real-part activity over all sources to check signal timecourse
    fprintf('subject %d:\nplotting the source activity...\n', s)
    fig = figure(figure_counter);
    set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1])
    for d = 1:length(data)
        for b = 1:size(data(d).data, 2)
            subplot(10, 10, b)
            for c = 1:size(data(d).data, 1)
                plot(data(d).times, real(squeeze(data(d).data(c, b, :))'), Color = params.colours(d, :));
                hold on
            end
            title(data(d).sources.ROIName{b})
        end
    end
    sgtitle(sprintf('source activity - subject %d (%s)', s, data(1).ID))
    saveas(fig, sprintf('%s\\figures\\source_activity_%s.png', folder.output, data(1).ID))
    figure_counter = figure_counter + 1;

    % check for the contribution of imaginary components - claculate
    % imaginary/real ratio per trial-source bin
    fprintf('checking the contribution of imaginary parts...\n')
    ratio = [];
    for d = 1:length(data)
        for b = 1:size(data(d).data, 2)
            for c = 1:size(data(d).data, 1)
                data_st = squeeze(data(d).data(c, b, :))';
                ratio(d, b, c) = max(abs(imag(data_st))) / max(abs(real(data_st)));
            end
        end
    end
    flagged = ratio > params.ratio_threshold;
    ratio =  permute(ratio, [1, 3, 2]);
    ratio = reshape(ratio, 2 * size(ratio, 2), size(ratio, 3));
    source_stats(s).subject = s;
    source_stats(s).ID = data(1).ID;
    source_stats(s).ratio.mean = mean(ratio, 'all');
    source_stats(s).ratio.threshold = params.ratio_threshold;
    source_stats(s).ratio.flagged = (sum(flagged, 'all') / numel(ratio));
    source_stats(s).flagged_data = flagged;
    fprintf('--> mean ratio %d\n', source_stats(s).ratio.mean)
    fprintf('--> %d%% of trial-source bins flagged (threshold %d)\n', ...
        source_stats(s).ratio.flagged * 100, params.ratio_threshold)
    
    % check the distribution of high ratios across sources
    flagged_source =  permute(flagged, [1, 3, 2]);
    flagged_source = reshape(flagged_source, 2 * size(flagged_source, 2), size(flagged_source, 3));
    flagged_source = sum(flagged_source, 1);
    flagged_source_num = find(flagged_source > params.source_threshold);
    source_stats(s).bad_sources.threshold = params.source_threshold;
    if ~isempty(flagged_source_num)
        flagged_source_label = data(1).sources.ROIName(flagged_source_num);
        source_stats(s).bad_sources.number = flagged_source_num;        
        source_stats(s).bad_sources.label = flagged_source_label;
        fprintf('--> %d bad sources identified.\n', length(flagged_source_num))
    else
        source_stats(s).bad_sources.number = [];        
        source_stats(s).bad_sources.label = {};
    end
    fig = figure(figure_counter);
    bar(flagged_source);
    xticks(flagged_source_num);
    xticklabels(flagged_source_label);
    xtickangle(40);
    xlabel('source');
    ylabel('number of flagged trials');
    title(sprintf('flagged trials per source - subject %d (%s)', s, data(1).ID));
    saveas(fig, sprintf('%s\\figures\\source_bad_%s.png', folder.output, data(1).ID))
    figure_counter = figure_counter + 1;

    % check variance across time for each source and trial
    fprintf('checking variance across timepoints...\n')
    variance = [];
    for d = 1:length(data)
        for b = 1:size(data(d).data, 2)
            for c = 1:size(data(d).data, 1)
                data_st = squeeze(data(d).data(c, b, :))';
                variance(d, b, c) = var(data_st);
            end
        end
    end
    var_threshold = mean(variance(:)) + 3*std(variance(:));
    flagged = variance > var_threshold;
    variance =  permute(variance, [1, 3, 2]);
    variance = reshape(variance, 2 * size(variance, 2), size(variance, 3));
    source_stats(s).variance.mean = mean(variance, 'all');
    source_stats(s).variance.threshold = var_threshold;
    source_stats(s).variance.flagged = (sum(flagged, 'all') / numel(variance));
    fprintf('--> mean variance %d\n', source_stats(s).variance.mean)
    fprintf('--> %d%% of trial-source bins flagged (threshold %d)\n', ...
        source_stats(s).variance.flagged * 100, var_threshold)
    fig = figure(figure_counter);
    set(fig, 'units', 'normalized', 'outerposition', [0.2 0.3 0.7 0.6])
    subplot(1, 2, 1)
    imagesc(variance);
    colorbar; xlabel('source'); ylabel('trial');
    title('variance heatmap');
    subplot(1, 2, 2)
    histogram(variance(:));
    xlabel('variance'); ylabel('frequency');
    title('variance distribution');
    sgtitle(sprintf('source activity variance - subject %d (%s)', s, data(1).ID));
    saveas(fig, sprintf('%s\\figures\\source_variance_%s.png', folder.output, data(1).ID))
    figure_counter = figure_counter + 1;

    % check the distribution of high variance across trials
    fprintf('checking for bad trials...\n')
    for d = 1:length(data) 
        flagged_trial = squeeze(flagged(d, :, :))';
        flagged_trial = sum(flagged_trial, 2);
        flagged_trial_num = find(flagged_trial > params.trial_threshold);
        source_stats(s).bad_trials(d).condition = data(d).condition; 
        source_stats(s).bad_trials(d).threshold = params.trial_threshold;
        source_stats(s).bad_trials(d).trials = flagged_trial_num; 
    end
    flagged_trial_num = [source_stats(s).bad_trials(1).trials; source_stats(s).bad_trials(2).trials];
    if ~isempty(flagged_trial_num)
        fprintf('--> %d bad trials identified across both conditions.\n', length(flagged_trial_num))
    end
    fprintf('\n')
end
AperiodicPFC_quality_check = source_stats;
save(output_file, 'AperiodicPFC_quality_check', '-append'); 

% clean and continue
clear a b c d s elec atlas cfg sourcemodel_atlas leadfield data header  ...
    dataset temptime timelock row_counter sources surface vol data_idx ...
    fig data_st ratio flagged flagged_source flagged_source_num flagged_source_label ...
    variance var_threshold flagged_trial flagged_trial_num source_stats
 fprintf('section finished.\n\n')

%% source-space analysis: extract aperiodic exponent
% ----- section input -----
params.foi = [5 80];
params.foi_APC = [30 50];
params.roi.lateral = {'7Networks_RH_Cont_PFCl_1', '7Networks_RH_Cont_PFCl_3', ...
    '7Networks_RH_Cont_PFCl_4'};
params.roi.medial = {'7Networks_LH_SalVentAttn_Med_1', '7Networks_LH_SalVentAttn_Med_2', '7Networks_LH_SalVentAttn_Med_3', ...
    '7Networks_RH_SalVentAttn_Med_1', '7Networks_RH_SalVentAttn_Med_2'};
params.roi.visual = {'7Networks_LH_Vis_2', '7Networks_LH_Vis_4', '7Networks_LH_Vis_5', ...
    '7Networks_RH_Vis_4', '7Networks_RH_Vis_5'};
% -------------------------

% load data if necessary
if exist('AperiodicPFC_sources') ~= 1 
    load('AperiodicPFC_sources.mat')
end
if exist('AperiodicPFC_quality_check') ~= 1
    load(output_file, 'AperiodicPFC_quality_check')
end
if exist('AperiodicPFC_measures') ~= 1
    load(output_file, 'AperiodicPFC_measures')
end
if exist('AperiodicPFC_data') ~= 1
    load(output_file, 'AperiodicPFC_data')
end

% calculate PSD for all datasets/trials/sources
addpath(genpath([folder.toolbox '\FieldTrip']));
fprintf('computing single-trial PSD...\n')
for a = 1:length(AperiodicPFC_sources)
    % select data, remove bad sources and trials 
    data = AperiodicPFC_sources(a).data;  
    s = AperiodicPFC_sources(a).subject;
    c = AperiodicPFC_sources(a).condition;
    data(:, AperiodicPFC_quality_check(s).bad_sources.number, :) = [];  
    idx = false(1, size(data, 1));
    for d = 1:size(data, 1)
        if ismember(d, AperiodicPFC_quality_check(s).bad_trials(c).trials')
            idx(d) = true;
        end
    end
    data(idx, :, :) = [];

    % encode info to the output structure
    AperiodicPFC_data(s).PSD_source(c).condition = c;
    AperiodicPFC_data(s).PSD_source(c).trials = size(data, 1);
    AperiodicPFC_data(s).PSD_source(c).trials_removed = AperiodicPFC_quality_check(s).bad_trials(c).trials';
    AperiodicPFC_data(s).PSD_source(c).sources = AperiodicPFC_sources(a).sources;
    AperiodicPFC_data(s).PSD_source(c).sources(AperiodicPFC_quality_check(s).bad_sources.number, :) = [];

    % compute frequency power spectrum for each trial and source    
    PSD = struct; 
    for d = 1:size(data, 1)
        % create a FieldTrip data structure
        cfg = [];
        cfg.trial = {squeeze(data(d, :, :))};       
        cfg.time = {AperiodicPFC_sources(a).times};      
        cfg.fsample = AperiodicPFC_sources(a).SR;
        cfg.label = AperiodicPFC_data(s).PSD_source(c).sources.ROIName;  
        cfg.sampleinfo = [1 1000];
        data_trial = ft_datatype_raw(cfg);
        ft_checkdata(data_trial);

        % extract original spectra
        cfg = [];
        cfg.output = 'pow';
        cfg.foilim = params.foi;  
        cfg.pad = 'nextpow2'; 
        cfg.method = 'irasa';  
        cfg.tapsmofrq = []; 
        cfg.output = 'original';
        PSD(d).original = ft_freqanalysis(cfg, data_trial);

        % extract fractal specra
        cfg.output = 'fractal';
        PSD(d).fractal = ft_freqanalysis(cfg, data_trial);

        % compute oscillatory spectra
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.operation     = 'x2-x1';
        PSD(d).oscillatory = ft_math(cfg, PSD(d).fractal, PSD(d).original);
    end

    % encode data to the output structure
    AperiodicPFC_data(s).PSD_source(c).freq = PSD(1).original.freq;  
    for d = 1:length(PSD)
        AperiodicPFC_data(s).PSD_source(c).original(d, :, :) = PSD(d).original.powspctrm;
        AperiodicPFC_data(s).PSD_source(c).fractal(d, :, :) = PSD(d).fractal.powspctrm;
        AperiodicPFC_data(s).PSD_source(c).oscillatory(d, :, :) = PSD(d).oscillatory.powspctrm;
    end
end
fprintf('done.\n')
save('AperiodicPFC_data2.mat', 'AperiodicPFC_data', '-v7.3'); 

% extract aperiodic measures for all ROIs
fprintf('extracting aperiodic measures..\n')
for r = fieldnames(params.roi)'
    % identify involved sources
    idx_sources = false(1, length(AperiodicPFC_data(1).PSD_source(1).sources.ROIName));
    for a = 1:length(AperiodicPFC_data(1).PSD_source(1).sources.ROIName)
        if any(strcmp(AperiodicPFC_data(1).PSD_source(1).sources.ROIName{a}, params.roi.(r{1})))
            idx_sources(a) = true;
        end
    end

    for s = 1:params.subjects
        for c = 1:2
            % subset data and average across rois
            data = AperiodicPFC_data(s).PSD_source(c).fractal(:, idx_sources, :); 
            data = squeeze(mean(data, 2)); 

            % identify frequencies
            idx_freq = AperiodicPFC_data(s).PSD_source(c).freq >= params.foi_APC(1) & ...
                AperiodicPFC_data(s).PSD_source(c).freq <= params.foi_APC(2);
            freq = AperiodicPFC_data(s).PSD_source(c).freq(idx_freq);

            % perform robust linear regression for all trials
            output.exponent = []; 
            output.offset = [];
            for b = 1:size(data, 1)
                % subset and log-transform 
                data_log = log10(data(b, idx_freq));
                freq_log = log10(freq);

                % fit for target roi
                [fit, ~] = robustfit(freq_log, data_log);
                output.exponent(b) = -fit(2);
                output.offset(b) = fit(1);
            end

            % encode to the output structure
            AperiodicPFC_measures(s).APC_source(c).condition = c;
            AperiodicPFC_measures(s).APC_source(c).flipped = AperiodicPFC_measures(s).APC(c).flipped;
            AperiodicPFC_measures(s).APC_source(c).trials = AperiodicPFC_data(s).PSD_source(c).trials; 
            AperiodicPFC_measures(s).APC_source(c).trials_removed = AperiodicPFC_data(s).PSD_source(c).trials_removed;
            AperiodicPFC_measures(s).APC_source(c).foi = params.foi_APC;
            AperiodicPFC_measures(s).APC_source(c).sources.(r{1}) = params.roi.(r{1});
            AperiodicPFC_measures(s).APC_source(c).exponent.(r{1}) = output.exponent;
            AperiodicPFC_measures(s).APC_source(c).offset.(r{1}) = output.offset;         
        end
    end
end
fprintf('done.\n')
save(output_file, 'AperiodicPFC_measures', '-append');

% clean and continue
clear a b c d r s idx data cfg data_trial PSD idx_sources idx_freq freq data_log freq_log fit output bad_sources bad_numbers
fprintf('section finished.\n\n')

%% source-space analysis: export for R  

% reload output structures if necessary 
if exist('AperiodicPFC_data') ~= 1 || exist('AperiodicPFC_measures') ~= 1
    load(output_file, 'AperiodicPFC_data', 'AperiodicPFC_measures')
end

% export measured variables in a long-format table
fprintf('exporting data:\n')
table_export = table;
row_counter = 1;
for s = 1:params.subjects
    for c = 1:length(AperiodicPFC_measures(s).conditions)
        % identify trials to be removed (based on trials discarded during source-analysis)
        trials2remove = AperiodicPFC_measures(s).APC_source(c).trials_removed;  

        % subset the data, filter out removed trials
        data = struct;
        data.exponent = AperiodicPFC_measures(s).APC_source(c).exponent;
        data.offset = AperiodicPFC_measures(s).APC_source(c).offset;
        for a = params.LEP_comps
            data.LEP.(a{1}).amplitude = AperiodicPFC_measures(s).LEP(c).(a{1}).amplitude;
            data.LEP.(a{1}).latency = AperiodicPFC_measures(s).LEP(c).(a{1}).latency;
            data.LEP.(a{1}).amplitude(trials2remove) = [];
            data.LEP.(a{1}).latency(trials2remove) = [];
        end
        data.rating = AperiodicPFC_measures(s).pain(c).ratings;
        data.rating(trials2remove) = [];
        
        % encode to the table
        for a = 1:length(data.rating)
            for b = fieldnames(data.LEP)'
                for r = fieldnames(data.exponent)'
                    % subject info
                    table_export.subject(row_counter) = s;
                    table_export.ID{row_counter} = AperiodicPFC_measures(s).ID;
                    table_export.age(row_counter) = AperiodicPFC_measures(s).age;
                    table_export.male(row_counter) = AperiodicPFC_measures(s).male;
                    table_export.handedness{row_counter} = AperiodicPFC_measures(s).handedness;
    
                    % session info
                    table_export.condition(row_counter) = c;
                    table_export.area{row_counter} = AperiodicPFC_measures(s).conditions(c).area;
                    table_export.side{row_counter} = AperiodicPFC_measures(s).conditions(c).side;
                    if strcmp(table_export.handedness{row_counter},table_export.side{row_counter})
                        table_export.dominant(row_counter) = 1;
                    else
                        table_export.dominant(row_counter) = 0;
                    end
                    table_export.flipped{row_counter} = AperiodicPFC_measures(s).APC_source(c).flipped;

                    % dependent variable: LEP measures
                    table_export.component{row_counter} = b{1};
                    table_export.amplitude(row_counter) = data.LEP.(b{1}).amplitude(a);
                    table_export.latency(row_counter) = data.LEP.(b{1}).latency(a);

                    % dependent variable: pain rating
                    table_export.rating(row_counter) = data.rating(a);
    
                    % independent variable: aperiodic measures
                    table_export.region{row_counter} = r{1}; 
                    table_export.exponent(row_counter) = data.exponent.(r{1})(a);
                    table_export.offset(row_counter) = data.offset.(r{1})(a);  
    
                    % update row counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end

% save table to output structure and as .csv
AperiodicPFC_export_sources = table_export;
save(output_file, 'AperiodicPFC_export_sources', '-append');
writetable(AperiodicPFC_export_sources, 'AperiodicPFC_export_sources.csv');
fprintf('done.\n')

% clean and continue
clear a b c r s row_counter table_export trials2remove data
fprintf('section finished.\n\n')

%% source-space analysis: visualization 
% ----- section input -----
params.foi_APC = [30 50];
params.path_headmodel = 'standard_bem.mat';
params.path_atlas = 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
params.roi.lateral = {'7Networks_RH_Cont_PFCl_1', '7Networks_RH_Cont_PFCl_3', ...
    '7Networks_RH_Cont_PFCl_4'};
params.roi.medial = {'7Networks_LH_SalVentAttn_Med_1', '7Networks_LH_SalVentAttn_Med_2', '7Networks_LH_SalVentAttn_Med_3', ...
    '7Networks_RH_SalVentAttn_Med_1', '7Networks_RH_SalVentAttn_Med_2'};
params.roi.visual = {'7Networks_LH_Vis_2', '7Networks_LH_Vis_4', '7Networks_LH_Vis_5', ...
    '7Networks_RH_Vis_4', '7Networks_RH_Vis_5'};
% -------------------------

% prepare average PSD  
data = []; data_counter = 1;
for s = 1:params.subjects
    for c = 1:2
        dataset = squeeze(mean(AperiodicPFC_data(s).PSD_source(c).fractal, 1));  
        data(data_counter, :, :) = dataset;
        data_counter = data_counter + 1;
    end
end
data = squeeze(mean(data, 1));

% identify frequencies
idx_freq = AperiodicPFC_data(1).PSD_source(1).freq >= params.foi_APC(1) & ...
    AperiodicPFC_data(1).PSD_source(1).freq <= params.foi_APC(2);
freq = AperiodicPFC_data(1).PSD_source(1).freq(idx_freq);

% calculate average aperiodic exponent for each source
exponent_avg = [];
for d = 1:size(data, 1)
    % subset and log-transform 
    data_log = log10(data(d, idx_freq));
    freq_log = log10(freq);

    % fit for target roi
    [fit, ~] = robustfit(freq_log, data_log);
    exponent_avg(d) = -fit(2);
end

% prepare a source model
addpath(genpath([folder.toolbox '\fieldtrip-20250318']));
load(sprintf('%s\\dataset_default.lw6', folder.output), '-mat')
params.labels = {header.chanlocs.labels};
params.elec_template = ft_read_sens('standard_1005.elc');
elec = create_elec(params, header);
load(params.path_headmodel,'vol');
atlas = readtable(params.path_atlas);
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas.R, atlas.A, atlas.S];
cfg.unit = 'mm';
cfg.headmodel = params.path_headmodel;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

% prepare plotting parameters 
sources2remove = AperiodicPFC_quality_check(1).bad_sources.number;
exponent_source = [];
exponent_source.pos = sourcemodel_atlas.pos;  
exponent_source.pos(sources2remove, :) = [];
exponent_source.inside = sourcemodel_atlas.inside;  
exponent_source.inside(sources2remove) = [];
exponent_source.pow = nan(size(exponent_source.pos, 1), 1);
exponent_source.pow(exponent_source.inside) = exponent_avg;
cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.projmethod = 'nearest';
cfg.surface = 'brain'; 
cfg.camlight = 'no';
cfg.camlight = 'no';
cfg.colorbar = 'yes';
cfg.funcolorlim = [mean(exponent_avg) - 2*std(exponent_avg), mean(exponent_avg) + 2*std(exponent_avg)];

% add target locations
params.regions = fieldnames(params.roi)';
for r = 1:3
    % identify source locations
    idx_sources = false(1, length(AperiodicPFC_data(1).PSD_source(1).sources.ROIName));
    for a = 1:length(AperiodicPFC_data(1).PSD_source(1).sources.ROIName)
        if any(strcmp(AperiodicPFC_data(1).PSD_source(1).sources.ROIName{a}, params.roi.(params.regions{r})))
            idx_sources(a) = true;
        end
    end
    centroids = table2array(AperiodicPFC_data(1).PSD_source(1).sources(idx_sources, [3:5]));  

    % plot average exponent value over all cortical sources
    ft_sourceplot(cfg, exponent_source);
    colormap(cool)
    alpha(0.4)
    hold on
    
    % add centroids to target locations
    scatter3(centroids(:,1), centroids(:,2), centroids(:,3), 75, 'k', 'filled');
    hold off
    title(sprintf('average aperiodic exponent values - %s areas', params.regions{r}))
end

% clean and continue
clear s c data dataset idx_freq freq data_log freq_log vol elec header ...
    atlas sourcemodel_atlas cfg exponent_source
fprintf('section finished.\n\n')

%% functions
function export_EEGLAB(lwdata, filename, subj)
% =========================================================================
% exports data from letswave to EEGLAB format
% =========================================================================  
% dataset
EEG.setname = filename;
EEG.filename = [];
EEG.filepath = [];
EEG.subject = subj; 
EEG.session = 1;
    
% time properties
EEG.nbchan = lwdata.header.datasize(2);
EEG.trials = lwdata.header.datasize(1);
EEG.pnts = lwdata.header.datasize(6);
EEG.srate = 1/lwdata.header.xstep;
EEG.times = lwdata.header.xstart + (0:EEG.pnts-1)*lwdata.header.xstep;
EEG.xmin = EEG.times(1);
EEG.xmax = EEG.times(end);
EEG.data = permute(single(lwdata.data),[2,6,1,3,4,5]);
EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'SEEG_enabled');
EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'topo_enabled');
    
% create events with appropriate latencies
EEG.event = lwdata.header.events;
if ~isempty(EEG.event)
    [EEG.event.type] = EEG.event.code;
    for e = 1:length(EEG.event)
        EEG.event(e).latency = (e-1)*EEG.pnts + EEG.xmin*(-1)*EEG.srate;
    end
    EEG.event = rmfield(EEG.event,'code');
end
    
% create required empty fields
EEG.icawinv = [];
EEG.icaweights = [];
EEG.icasphere = [];
EEG.icaweights = [];
EEG.icaweights = [];
EEG.icaweights = [];
EEG.icaweights = [];
save([filename,'.set'], 'EEG');
end
function elec = create_elec(params, header)
% =========================================================================
% creates elec structure from loaded electrode location file
% =========================================================================  
elec_idx = ismember(params.elec_template.label, {header.chanlocs.labels}');
if sum(elec_idx) ~= length({header.chanlocs.labels})
    error('ERROR: subject %d, condition %d - channel labels do not match!', s, c)
end
elec.chanpos = params.elec_template.chanpos(elec_idx,:);
elec.chantype = params.elec_template.chantype(elec_idx,:);
elec.chanunit = params.elec_template.chanunit(elec_idx,:);
elec.elecpos = params.elec_template.elecpos(elec_idx,:);
elec.label = params.elec_template.label(elec_idx,:);
[~, elec_orig] = sort({header.chanlocs.labels}');
[~, elec_new] = sort(elec_orig); 
[~, elec_temp] = sort(elec.label);
elec.chanpos = elec.chanpos(elec_temp(elec_new), :);
elec.chantype = elec.chantype(elec_temp(elec_new), :);
elec.chanunit = elec.chanunit(elec_temp(elec_new), :);
elec.elecpos = elec.elecpos(elec_temp(elec_new), :);
elec.label = elec.label(elec_temp(elec_new), :);
elec.type = 'custom';
elec.unit = params.elec_template.unit;
end