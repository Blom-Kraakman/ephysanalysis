function [spiketimes, cids, Srise, Sfall] = extractspikes(BehaviorPath, KSPath, TTLPath, relevant_sessions, skip_sessions, rec_samples, sessions_TTLs, Fs, OutPath)
% Kilosort: post-curation unit extraction
% INPUT - paths to sorted data (cluster_info, table), spike times (vector)
% and matched unit ids (vector), recording time stamps (vector)
% OUTPUT - spike times of each single unit
% based on postcuration in PostCuration_ABW.m

%check if paths contain needed files
if ~isfile([KSPath,'cluster_info.tsv']) || ~isfile([KSPath,'spike_clusters.npy']) || ~isfile([KSPath,'spike_times.npy'])
    error('Files not found in KSPath.');
elseif ~isfile([TTLPath,'sample_numbers.npy']) ||  ~isfile([TTLPath,'states.npy'])
    error('Files not found in TTLPath.')
end

cluster_info = readtable([KSPath,'cluster_info.tsv'],'FileType','text'); % info on clusters
cids = cluster_info.cluster_id(strcmp(cluster_info.group,'good'))';
[~, idx] = ismember(cids, cluster_info.cluster_id);
% to do: improve name of table cpos
cpos(:,1) = cluster_info.cluster_id(idx);
cpos(:,2) = cluster_info.ch(idx);
cpos(:,3) = cluster_info.depth(idx);
cpos(:,4) = cluster_info.fr(idx);
cpos(:,5) = cluster_info.n_spikes(idx);
cpos = array2table(cpos, 'VariableNames', {'id' , 'channel', 'depth', 'firing_rate', 'nr_spikes'});

fprintf('Found %i good units for analysis\n', length(cids));

% load files
spike_times = readNPY([KSPath 'spike_times.npy']); % spike_times contains spike time indexing, not time/samplenr itself
spike_clusters = readNPY([KSPath 'spike_clusters.npy']); % matched cluster ids
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));

% get behaviour hile
Nr_sessions = (relevant_sessions(1):relevant_sessions(2))';
for file = 1:length(Nr_sessions)
    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);
end

% remove unexecuted sessions
idx_remove = find(sum(sessions_TTLs(:,1) == skip_sessions, 2));
sessions_TTLs(idx_remove,:) = [];

%keep only TTLs recorded during specific session
[Srise, Sfall] = TTLsToUse(sessions_TTLs, TTLPath, rec_samples);

% get all spiketimes from each good unit
spiketimes = cell(length(cids), 1);

% exclusion criteria: firingrate > 0.1Hz
total_rec = (rec_samples(length(rec_samples)) - rec_samples(1))/Fs;
minimal_freq = total_rec * 0.1; % min amount of spike

for cluster = 1:length(cids)
    idx = (spike_clusters == cids(cluster));
    spiketimes{cluster} = rec_samples(spike_times(idx));
    % if length(spiketimes{cluster}) < minimal_freq % to do: actually remove unit from all relevant variables
    %     spiketimes{cluster} = NaN;
    % end
end

filename = sprintf('M%.2i_S%02d-%02d_InfoGoodUnits', str2double(stimuli_parameters.Par.MouseNum), relevant_sessions(1), relevant_sessions(2));
save(fullfile(OutPath, filename), "cpos") %cpos variables: unit id, channel, depth, avg firing rate, nr spikes;

fprintf('unit extraction done\n');

end


