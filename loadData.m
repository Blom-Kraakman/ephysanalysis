function [cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, onsetDelay, StimResponseFiring] = loadData(OutPath, session, BehaviorPath)
% load several data files, if they exist

% load cluster unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']); %.name;
if ~isempty(cpos_file)
    cpos = load([OutPath '\' cpos_file.name]);
    cids = cpos.clusterinfo.id';
else
    cids = [];
    warning('no cluster details found')
end

% load behaviour files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

% load aligned spike times
aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
if ~isempty(aligned_spikes_files)
    aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
    if isfield(aligned_spikes,"Srise")
        disp("Srise/Sfall loaded from data file.")
        Srise = aligned_spikes.Srise;
        Sfall = aligned_spikes.Sfall;
    else
        Srise = [];
        Sfall = [];
        warning("Extract Srise before continuing")
    end
    aligned_spikes = aligned_spikes.SpkT;
else
    Srise = [];
    Sfall = [];
    aligned_spikes = [];
    warning('no aligned spikes file found')
end

% load sessions details
TTLs_file = dir([OutPath '\*_OE_TTLs.mat']);
if ~isempty(TTLs_file)
    sessions_TTLs = load([OutPath '\' TTLs_file.name]);
    sessions_TTLs = sessions_TTLs.sessions_TTLs;
else
    sessions_TTLs = [];
    warning('no session TTL file found')
end

% load response properties
if isfolder([OutPath '\ResponseProperties\'])
    files = ['\*_S' num2str(session, '%.2d') '_*.mat'];
    resp_files = dir(fullfile([OutPath '\ResponseProperties'], files));
    if ~isempty(resp_files)
        StimResponseFiring = load([resp_files.folder '\' resp_files.name]);
        StimResponseFiring = StimResponseFiring.StimResponseFiring;
        disp('analysed data loaded')
    else
        StimResponseFiring = [];
    end
else
    StimResponseFiring = [];
end

% add delay to "SO"
if strcmp(stimuli_parameters.Par.Rec, 'SxA')
    onsetDelay = stimuli_parameters.Stm.SomAudSOA./1000;
    onsetDelay(isnan(stimuli_parameters.Stm.SomAudSOA)) = 0;
else
    onsetDelay = zeros(size(stimuli_parameters.Stm, 1), 1);
end

end