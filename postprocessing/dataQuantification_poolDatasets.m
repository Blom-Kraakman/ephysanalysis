function dataQuantification_poolDatasets(stimOrder, fn)
%% combines datasets across mice
% INPUT: .csv containing matched session type + session number of all animals to be pooled
% OUTPUT: *_UnitResponses per session type containing 
%   data_all: delta FR for all units
%   unitResponses_all: table with metadata 

for condition = 2:(size(stimOrder, 2)-1)

    % skip analysis if already done
    StimType = stimOrder.Properties.VariableNames{condition};
    if ~isempty(dir(['\\store\department\neuw\shared\Aaron Wong\Data\ProcessedData\Blom\Processed\M' '*_' StimType '_UnitResponses.mat']))
        continue
    end

    disp(['Concatinating ' StimType])

    %initiate variables
    unitResponses_all = table;
    firing_mean = [];
    FSLmed = [];
    FSLiqr = [];
    hvalue = [];
    mousenum_all =[];
    session_all = [];
    amplitudes_all = [];
    frequencies_all = [];
    PreT_all = [];
    PostT_all = [];

    % get data from session to analyse
    for animal = 1:size(stimOrder,1)

        disp(['Mouse ' num2str(stimOrder{animal,1})])

        if isnan(stimOrder{animal,condition})
            continue
        end

        if stimOrder{animal,1} == 10 || stimOrder{animal,1} == 11 ||stimOrder{animal,1} == 19 || stimOrder{animal,1} == 20
            OutPath = ['\\store\department\neuw\shared\Aaron Wong\Data\ProcessedData\Blom\Processed\M' num2str(stimOrder{animal,1}, '%d') '\ICX\ResponseProperties\'];
        else
            OutPath = ['\\store\department\neuw\shared\Aaron Wong\Data\ProcessedData\Blom\Processed\M' num2str(stimOrder{animal,1}, '%d') '\ResponseProperties\'];
        end

        sessionFile = ['M' num2str(stimOrder{animal,1}, '%.2d') '_S' num2str(stimOrder{animal,condition}, '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
        stim_files = dir(fullfile(OutPath, sessionFile));
        dataS = load([stim_files.folder '\' stim_files.name]);

        % firing mean per stimulus combination
        firing_data = dataS.StimResponseFiring.firing_mean;
        firing_mean = cat(4, firing_mean, firing_data);

        % FSL per stimulus combination
        FSLmed_data = dataS.StimResponseFiring.FSLmed;
        FSLmed = cat(4, FSLmed, FSLmed_data);

        FSLiqr_data = dataS.StimResponseFiring.FSLiqr;
        FSLiqr = cat(4, FSLiqr, FSLiqr_data);

        % h-val per stimulus combination
        hval_data = dataS.StimResponseFiring.hvalue;
        hvalue = cat(3, hvalue, hval_data);

        % add units responses table
        unitResponses_all = vertcat(unitResponses_all, dataS.StimResponseFiring.unitResponses);

        % add details
        mousenum_all = cat(2, mousenum_all, str2double(dataS.StimResponseFiring.MouseNum));
        session_all = cat(2, session_all, str2double(dataS.StimResponseFiring.session));
        amplitudes_all = cat(2, amplitudes_all, dataS.StimResponseFiring.amplitudes);
        frequencies_all = cat(2, frequencies_all, dataS.StimResponseFiring.frequencies);
        PreT_all = cat(2, PreT_all, dataS.StimResponseFiring.PreT);
        PostT_all = cat(2, PostT_all, dataS.StimResponseFiring.PostT);

    end

    % add mousenum as prefix to unit nr to make unique
    cids = unitResponses_all.Cluster';
    mice = unitResponses_all.MouseNum';
    newcids = sprintf( '%1d%1d ', [mice;cids]);
    unitResponses_all.Cluster = str2num(newcids)';

    StimResponseFiring_all.StimType = StimType;
    StimResponseFiring_all.MouseNum = mousenum_all;
    StimResponseFiring_all.session = session_all;
    StimResponseFiring_all.amplitudes = amplitudes_all;
    StimResponseFiring_all.frequencies = frequencies_all;
    StimResponseFiring_all.PreT = PreT_all;
    StimResponseFiring_all.PostT = PostT_all;
    StimResponseFiring_all.unitResponses = unitResponses_all;
    StimResponseFiring_all.firing_mean = firing_mean;
    StimResponseFiring_all.FSLmed = FSLmed;
    StimResponseFiring_all.FSLiqr = FSLiqr;
    StimResponseFiring_all.hvalue = hvalue;

    % onsetdelay


    % save
    filename = sprintf('M10-11-19-20-29-31-33_%s_UnitResponses', StimType);
    save(fullfile('\\store\department\neuw\shared\Aaron Wong\Data\ProcessedData\Blom\Processed\', filename), "StimResponseFiring_all", "firing_mean")

    clear firing_mean

end
end