function [spiketimes, cids, Srise, Sfall, cpos] = extractspikes(KSPath, recPath, TTLPath)
% Kilosort: post-curation unit extraction
% INPUT - paths to sorted data (cluster_info, table), spike times (vector)
% and matched unit ids (vector), recording time stamps (vector)
% OUTPUT - spike times of each single unit
% based on postcuration in PostCuration_ABW.m

% check input arguments (Fs now initiated in 2 functions)

% check if paths contain needed files
if ~isfile([KSPath,'cluster_info.tsv']) || ~isfile([KSPath,'spike_clusters.npy']) || ~isfile([KSPath,'spike_times.npy'])
    error('Files not found in KSPath.');
elseif ~isfile([TTLPath,'sample_numbers.npy']) ||  ~isfile([TTLPath,'states.npy'])
    error('Files not found in TTLPath.')
elseif ~isfile([recPath, 'sample_numbers.npy'])
    error('Files not found in RecPath.')
end

disp('Automatically extracting units labelled ''good'' in Phy.')

cluster_info = readtable([KSPath,'cluster_info.tsv'],'FileType','text'); % info on clusters
cids = cluster_info.cluster_id(strcmp(cluster_info.group,'good'))';
fprintf('%i good units for analysis\n', length(cids));% show good units
[~, idx] = ismember(cids, cluster_info.cluster_id);
cpos(:,1) = cluster_info.cluster_id(idx);
cpos(:,2) = cluster_info.ch(idx);
cpos(:,3) = cluster_info.depth(idx);

% extract spike times of 'single' units
timestamps = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording
spike_times = readNPY([KSPath 'spike_times.npy']); % spike_times contains spike time indexing, not time/samplenr itself
spike_clusters = readNPY([KSPath 'spike_clusters.npy']); % matched cluster ids

% load TTLs
TTL_samples = readNPY([TTLPath 'sample_numbers.npy']);
TTL_states = readNPY([TTLPath 'states.npy']);
Fs = 30000; % Sampling frequency (in Hz)

% select on and offset stimuli (5=SOM, 2=AUD)
Srise = TTL_samples((TTL_states == 2) | (TTL_states == 5)); % Srise = column vector of length NStim where the spikes should be aligned
Sfall = TTL_samples((TTL_states == -2) | (TTL_states == -5)); % Sfall = column vector of length nStm indicating end of stimulus

% remove artefacts where Srise == Sfall
idx = find(ismember(Sfall,Srise));
Srise(idx) = [];
Sfall(idx) = [];

% get all spiketimes from each good unit
spiketimes = cell(length(cids), 1);

% % exclusion criteria: firingrate > 0.1Hz
total_rec = (timestamps(length(timestamps)) - timestamps(1))/Fs;
minimal_freq = total_rec * 0.1; % min amount of spike

for cluster = 1:length(cids)
    idx = (spike_clusters == cids(cluster));
    spiketimes{cluster} = timestamps(spike_times(idx));
    % if length(spiketimes{cluster}) < minimal_freq % to do: actually remove unit from all relevant variables
    %     spiketimes{cluster} = NaN;
    % end
end

disp('unit extraction done');

end
