function dataQuantification_poolDatasets(stimOrder, fn)
%% combines datasets across mice
% INPUT: .csv containing matched session type + session number of all animals to be pooled
% OUTPUT: *_UnitResponses per session type containing 
%   data_all: delta FR for all units
%   unitResponses_all: table with metadata 

for condition = 2:(size(stimOrder, 2))

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

        % firing mean, FSL and h-val per stimulus combination
        if strcmp(stimOrder.Properties.VariableNames{condition}, 'AM')
            firing_data = dataS.StimResponseFiring.firing_mean;
            FSLmed_data = dataS.StimResponseFiring.FSLmed;
            FSLiqr_data = dataS.StimResponseFiring.FSLiqr;
            hval_data = dataS.StimResponseFiring.hvalue;
            
            tPre_T = dataS.StimResponseFiring.PreT;
            tPost_T = dataS.StimResponseFiring.PostT;
            ufreqs = dataS.StimResponseFiring.frequencies;
            uamps = dataS.StimResponseFiring.amplitudes;

        else

            % adjust data matrix per animal to accomodate difference in dimensions
            if str2num(dataS.StimResponseFiring.MouseNum) == 27
                firing_data = dataS.StimResponseFiring.firing_mean(1:8, 2:3,:,:);
                firing_data(1,:,:,:) = repmat(dataS.StimResponseFiring.firing_mean(1,1,:,:), 1,2);
                FSLmed_data = dataS.StimResponseFiring.FSLmed(1:8, 2:3,:,:);
                FSLmed_data(1,:,:,:) = repmat(dataS.StimResponseFiring.FSLmed(1,1,:,:),1,2);
                FSLiqr_data = dataS.StimResponseFiring.FSLiqr(1:8, 2:3,:,:);
                FSLiqr_data(1,:,:,:) = repmat(dataS.StimResponseFiring.FSLiqr(1,1,:,:),1,2);
                hval_data = dataS.StimResponseFiring.hvalue(1:8, 2:3,:,:);
                hval_data(1,:,:,:) = repmat(dataS.StimResponseFiring.hvalue(1,1,:,:),1,2);
                uamps = dataS.StimResponseFiring.amplitudes(dataS.StimResponseFiring.amplitudes > 0);
            else
                firing_data = dataS.StimResponseFiring.firing_mean(1:8, :,:,:);
                FSLmed_data = dataS.StimResponseFiring.FSLmed(1:8, :,:,:);
                FSLiqr_data = dataS.StimResponseFiring.FSLiqr(1:8, :,:,:);
                hval_data = dataS.StimResponseFiring.hvalue(1:8, :,:,:);
                uamps = dataS.StimResponseFiring.amplitudes;
            end

            ufreqs = dataS.StimResponseFiring.frequencies(1:8);
            tPre_T = dataS.StimResponseFiring.PreT(1);
            tPost_T = dataS.StimResponseFiring.PostT(1);

        end

        % add to matrix
        firing_mean = cat(4, firing_mean, firing_data);
        FSLmed = cat(4, FSLmed, FSLmed_data);
        FSLiqr = cat(4, FSLiqr, FSLiqr_data);
        hvalue = cat(3, hvalue, hval_data);

        unitResponses_all = vertcat(unitResponses_all, dataS.StimResponseFiring.unitResponses);

        % add units responses table
        mousenum_all = cat(2, mousenum_all, str2double(dataS.StimResponseFiring.MouseNum));
        session_all = cat(2, session_all, str2double(dataS.StimResponseFiring.session));
        amplitudes_all = cat(2, amplitudes_all, uamps);
        frequencies_all = cat(2, frequencies_all, ufreqs);
        PreT_all = cat(2, PreT_all, tPre_T);
        PostT_all = cat(2, PostT_all, tPost_T);
    end

    % add mousenum as prefix to unit nr to make unique
    cids = unitResponses_all.Cluster';
    mice = unitResponses_all.MouseNum';
    newcids = sprintf( '%1d%1d ', [mice;cids]);
    unitResponses_all.Cluster = str2num(newcids)';
    unitResponses_all.NewCluster = str2num(newcids)';

    StimResponseFiring_all.StimType = stimOrder.Properties.VariableNames{condition};
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

    % save
    filename = sprintf('M10-11-19-20_%s_UnitResponses', stimOrder.Properties.VariableNames{condition});
    save(fullfile('D:\DATA\Processed\', filename), "StimResponseFiring_all", "firing_mean")

    clear firing_mean firing_data

end
end