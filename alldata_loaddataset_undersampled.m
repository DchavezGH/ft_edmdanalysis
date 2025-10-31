% ======================================================
% Load all runs for all subjects, resample to 250 Hz,
% keep only midline channels & trials with >=5 letters
% ==================================================s====
%numSubjects = 1;
%runsPerSubject = 1;
numSubjects = 24;
runsPerSubject = 4;

% Root path where all subject folders live
rootPath = ['C:\Users\dchavezhuerta\OneDrive - Delft University of Technology\' ...
    'Documents\PHD\PHD project\Show rerun\Dataset Sternberg Woring Memory\' ...
    'DMD\ds004117-1.0.1'];

allData = cell(numSubjects, runsPerSubject);
allHeaders = cell(numSubjects, runsPerSubject);
for subj = 1:numSubjects
%for subj = 22:23

    subjID = sprintf('sub-%03d', subj); % sub-001, sub-002, etc.
    subjPath = fullfile(rootPath, subjID, 'ses-01', 'eeg');
    
    for run = 1:runsPerSubject
        runID = sprintf('run-%d', run);
        
        setFile = fullfile(subjPath, ...
            sprintf('%s_ses-01_task-WorkingMemory_%s_eeg.set', subjID, runID));
        fdtFile = strrep(setFile, '.set', '.fdt');
        
        if ~isfile(setFile) || ~isfile(fdtFile)
            warning('File missing: %s', setFile);
            continue;
        end
        
        fprintf('Processing %s\n', setFile);
        
        % ---- Load header (.set) ----
        headerInfo = load(setFile, '-mat');
        
        % ---- Load EEG time-series from .fdt ----
        fid = fopen(fdtFile, 'r', 'ieee-le');
        data = zeros(headerInfo.nbchan, headerInfo.pnts, headerInfo.trials);
        for trialIdx = 1:headerInfo.trials
            currentTrialData = fread(fid, ...
                [headerInfo.nbchan headerInfo.pnts], 'float32');
            data(:,:,trialIdx) = currentTrialData;
        end
        fclose(fid);

        % ---- Resample to 250 Hz ----
        origFs = headerInfo.srate;
        targetFs = 250;

        % Round sampling rate to nearest nominal value
        if abs(origFs - 250) < 1
            origFs = 250;
        elseif abs(origFs - 500) < 1
            origFs = 500;
        elseif abs(origFs - 1000) < 1
            origFs = 1000;
        else
            warning('Unusual sampling rate: %.4f Hz (subject %d, run %d)', ...
                origFs, subj, run);
            origFs = round(origFs);
        end

        if origFs ~= targetFs
            fprintf('Resampling from %.4f Hz to %d Hz\n', headerInfo.srate, targetFs);

            resampledTrials = cell(1, size(data,3));
            for tr = 1:size(data,3)
                resampledTrials{tr} = resample(data(:,:,tr)', targetFs, origFs)';
            end

            % Match all trial lengths by truncating to shortest
            newPnts = min(cellfun(@(x) size(x,2), resampledTrials));
            resampledData = zeros(size(data,1), newPnts, size(data,3));
            for tr = 1:size(data,3)
                resampledData(:,:,tr) = resampledTrials{tr}(:,1:newPnts);
            end

            data = resampledData;

            % Update header info
            headerInfo.srate = targetFs;
            headerInfo.pnts = size(data,2);

            % Rescale event latencies (samples)
            scaleFactor = targetFs / origFs;
            for e = 1:length(headerInfo.event)
                headerInfo.event(e).latency = ...
                    round(headerInfo.event(e).latency * scaleFactor);
            end
        end
        
        % ---- Store results ----
        allData{subj, run} = data;
        allHeaders{subj, run} = headerInfo;
    end
end

% Preallocate for processed subsets
fmidData = cell(numSubjects, runsPerSubject);
fmidHeaders = cell(numSubjects, runsPerSubject);

% Filtering, segmentation, channel selection 
for subj = 1:numSubjects
    for run = 1:runsPerSubject
        if ~isempty(allData{subj, run})
            disp(['Processing subject ' num2str(subj) ', run ' num2str(run)])
            currentTrialData = allData{subj, run};
            headerInfo = allHeaders{subj, run};

            latency    = round([headerInfo.event.latency]');
            letters    = string({headerInfo.event.letter})';
            event_card = string({headerInfo.event.memory_cond})';
            trial      = string({headerInfo.event.trial})';
            task_role  = string({headerInfo.event.task_role})';

            eventArray = [num2cell(latency), cellstr(letters), ...
                          num2cell(event_card), num2cell(trial), ...
                          cellstr(task_role)];

            X_0 = filloutliers(currentTrialData, "linear");

            [pSegmentsA, eventInfo] = segmentByEventFiltered(X_0, eventArray);
            pSegmentstotal = pSegmentsA';

            % Regularize segment lengths
            cols = cellfun(@(x) size(x,2), pSegmentstotal);
            cols = cols(cols > 0);
            if isempty(cols), continue; end
            n = min(cols);
            pSegmentstotal = cellfun(@(M) M(:,1:n), ...
                                     pSegmentstotal, 'UniformOutput', false);

            % Keep frontal midline channels (hard-coded for now)
            rows = [14 25 35 45 6 13 15];
            pSegmentstotal = cellfun(@(M) M(rows,:), ...
                                     pSegmentstotal, 'UniformOutput', false);

            fmidData{subj, run} = pSegmentstotal;
            fmidHeaders{subj, run} = eventInfo;
        end
    end
end

clearvars -except fmidData fmidHeaders numSubjects runsPerSubject

% ======================================================
% Helper functions
% ======================================================
function [eventSegments, eventInfo] = segmentByEventFiltered(x, event)
    snapshotCol  = str2double(string(event(:,1)));
    eventCardCol = str2double(string(event(:,3)));
    trialCol     = str2double(string(event(:,4)));
    roleCol      = string(event(:,5));

    mask = (eventCardCol >= 5) & (roleCol == "to_remember");
    idx = find(mask);

    eventSegments = {};
    eventInfo = {};

    uniqueTrials = unique(trialCol(idx));
    for t = 1:numel(uniqueTrials)
        thisTrial = uniqueTrials(t);
        trialIdx = idx(trialCol(idx) == thisTrial);

        % Keep only first 5
        trialIdx = trialIdx(1:min(5, numel(trialIdx)));

        for k = 1:numel(trialIdx)
            startSnap = snapshotCol(trialIdx(k));
            if trialIdx(k) < size(event,1)
                endSnap = snapshotCol(trialIdx(k)+1) - 1;
            else
                endSnap = size(x, 2);
            end
            eventSegments{end+1} = x(:, startSnap:endSnap);
            eventInfo(end+1,:) = event(trialIdx(k), :);
        end
    end
end
