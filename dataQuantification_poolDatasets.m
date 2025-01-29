function dataQuantification_poolDatasets(stimOrder, fn)
%% combines datasets across mice
% INPUT: .csv containing matched session type + session number of all animals to be pooled
% OUTPUT: *_UnitResponses per session type containing 
%   data_all: delta FR for all units
%   unitResponses_all: table with metadata 

for condition = 2:size(stimOrder, 2)

    %initiate variables
    unitResponses_all = table;
    data_all = [];
    mousenum_all =[];
    session_all = [];
    amplitudes_all = [];
    frequencies_all = [];
    PreT_all = [];
    PostT_all = [];

    disp(['Concatinating ' stimOrder.Properties.VariableNames{condition}])

    % get data from session to analyse
    for animal = 1:size(stimOrder,1)

        disp(['Mouse ' num2str(stimOrder{animal,1})])


        if isnan(stimOrder{animal,condition})
            continue
        end

        OutPath = ['D:\DATA\Processed\M' num2str(stimOrder{animal,1}, '%d') '\ICX\ResponseProperties\'];
        sessionFile = ['M' num2str(stimOrder{animal,1}, '%.2d') '_S' num2str(stimOrder{animal,condition}, '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
        stim_files = dir(fullfile(OutPath, sessionFile));
        dataS = load([stim_files.folder '\' stim_files.name]);

        %disp(['post stim time: ' num2str(dataS.StimResponseFiring.PostT) 'sec'])
        %disp(['pre stim time: ' num2str(dataS.StimResponseFiring.PreT) 'sec'])

        % concatinate data from all units
        % select common amp and freq data points
        % if strcmp(dataS.StimResponseFiring.type, 'SxA')
        %
        %     if animal(animal) == 9 && session(animal) == 4
        %         data = dataS.StimResponseFiring.firing_mean(:, [1,3,5], :, :);
        %     elseif animal(animal) == 8 || animal(animal) == 10 || animal(animal) == 11
        %         data = dataS.StimResponseFiring.firing_mean(1:8, :, :, :);
        %     else
        %         data = dataS.StimResponseFiring.firing_mean; % no changes needed
        %     end
        % else
        %     data = dataS.StimResponseFiring.firing_mean; % no changes needed
        % end

        data = dataS.StimResponseFiring.firing_mean;
        data_all = cat(4, data_all, data);

        % combine unitResponses tables
        % if strcmp(dataS.StimResponseFiring.type, 'SxA')
        %     sessionFile = ['M' num2str(animal(file), '%.2d') '_S' num2str(session(file), '%.2d') '*UnitResponses.mat']; % select based on stimulus type '_*.mat'
        %     stim_files = dir(fullfile(OutPath, sessionFile));
        %     unitResponses = load([stim_files.folder '\' stim_files.name]);
        %     unitResponses_all = vertcat(unitResponses_all, unitResponses.unitResponses);
        % else
        unitResponses_all = vertcat(unitResponses_all, dataS.StimResponseFiring.unitResponses);

        % add details
        mousenum_all = cat(2, mousenum_all, str2double(dataS.StimResponseFiring.MouseNum));
        session_all = cat(2, session_all, str2double(dataS.StimResponseFiring.session));
        amplitudes_all = cat(2, amplitudes_all, dataS.StimResponseFiring.amplitudes);
        frequencies_all = cat(2, frequencies_all, dataS.StimResponseFiring.frequencies);
        PreT_all = cat(2, PreT_all, dataS.StimResponseFiring.PreT);
        PostT_all = cat(2, PostT_all, dataS.StimResponseFiring.PostT);

    end

    StimResponseFiring_all.unitResponses = unitResponses_all;
    StimResponseFiring_all.firing_mean = data_all;
    StimResponseFiring_all.MouseNum = mousenum_all;
    StimResponseFiring_all.session = session_all;
    StimResponseFiring_all.amplitudes = amplitudes_all;
    StimResponseFiring_all.frequencies = frequencies_all;
    StimResponseFiring_all.PreT = PreT_all;
    StimResponseFiring_all.PostT = PostT_all;

    % save
    %filename = sprintf('M12-16_%s_UnitResponses',stimOrder.Properties.VariableNames{condition});
    save(fullfile('D:\DATA\Processed\', fn), "StimResponseFiring_all", "data_all")

    clear data_all

end

%% 5.5 save if needed
%filename = sprintf('M09-10-11_SxA_earplug_UnitResponses');
%save(fullfile(OutPath, filename), "unitResponses_all", "data_all")

end