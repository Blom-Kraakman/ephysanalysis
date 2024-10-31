function [uAmp, uFreq, conditions] = selectparameters(stimuli_parameters)
% INPUT: stimuli parameters with the relevant variable
% OUTPUT: vectors with unique variables

% select correct amplitude and frequency parameters
if strcmp(stimuli_parameters.Par.Rec, "AMn") % noise session
    uAmp = unique(stimuli_parameters.Stm.Intensity); 
    uFreq = unique(stimuli_parameters.Stm.Mf);
elseif strcmp(stimuli_parameters.Par.Rec, "FRA")
    uAmp = unique(stimuli_parameters.Stm.Intensity);
    uFreq = unique(stimuli_parameters.Stm.Freq);
elseif strcmp(stimuli_parameters.Par.Rec, "SxA") && strcmp(stimuli_parameters.Par.SomatosensoryWaveform, "Square") % SxA pressure session
    uAmp = unique(stimuli_parameters.Stm.Amplitude);
    uFreq = unique(stimuli_parameters.Stm.AudIntensity); % sound variable
else % SxA vibrotac and SOM sessions
    uAmp = unique(stimuli_parameters.Stm.Amplitude);
    uFreq = unique(stimuli_parameters.Stm.SomFreq);
end


% select multiple stimuli types in one session
if strcmp(stimuli_parameters.Par.Rec, "SxA") %dual modal session
    conditions = {'OO', 'OA', 'SO', 'SA'};
else % unimodal sessions
    conditions = 1;
end
end