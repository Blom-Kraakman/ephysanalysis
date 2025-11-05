%% Archived plotting
% keep track of old plotting code, currently not in used / not functional

%% venn diagram: response ratio groups
%sets = ["sound" "vibrotactile" "pressure"];
sets = ["sound" "tactile"];
total = 64; % %size(data_all,4);
sound_resp = 50;
tactile_sound = 9; % 27 for all tac resp
tactile_resp = 17; % 35 for all tac resp
%tactile_only = 8;

%labels = [round((sum(sound_resp)/total)*100), round((sum(tactile_resp)/total)*100), round((sum(tactile_sound)/total)*100)]; % order: A, B, C, A&B, A&C, B&C and A&B&C
labels = [sound_resp-tactile_sound, tactile_resp-tactile_sound, tactile_sound]; % order: A, B, C, A&B, A&C, B&C and A&B&C

venn(2,'sets',sets,'labels',labels, 'colors', [0.808 0.529 0.666;0.247 0.635 0.831], 'alpha',0.75);

%% plotting different multisensory indexes

% % preference index
% subplot(1,3,1);
% imagesc(preference_index(:,:, cluster),'AlphaData', NaNindex_preference(:,:,cluster))
% cb = colorbar(gca, 'eastoutside');
% cb.Label.String = 'modality preference';
% clim([-1, 1]);
% cb.Ticks = [-1, 0, 1];
%cb.TickLabels = {'tactile preference' 'no preference' 'sound preference'};
% set(gca, 'YDir', 'normal', 'XTick',1:size(preference_index,2),'XTickLabel',umN,'YTick',1:size(preference_index,1),'YTicklabel',StimResponseFiring.frequencies(:,1))
% xlabel('stimulus intensity (mN)')
% ylabel('bbn intensity (dB SPL)')
% % modulation index
% subplot(2,2,2);
% imagesc(modulation_index(:,:, cluster),'AlphaData', NaNindex_modulation(:,:,cluster))
% cb = colorbar(gca, 'eastoutside');
% cb.Label.String = 'multisensory modulation';
% set(gca, 'YDir', 'normal', 'XTick',1:size(modulation_index,2),'XTickLabel',umN,'YTick',1:size(modulation_index,1),'YTicklabel',StimResponseFiring.frequencies(:,1))
% xlabel('stimulus intensity (mN)')
% ylabel('bbn intensity (dB SPL)')
% %% multisensory index tests
% StimResponseFiring = StimResponseFiring_all;
% data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids
% unit_idx = 57;
% data2plot = data(:,:,unit_idx);
% dFRmulti = data2plot;
% dFRsound = repmat(dFRmulti(:,1),[1,9]);
% dFRvibrotac = repmat(dFRmulti(1,:),[5,1]);
% dFRmulti_sign = sign(dFRmulti);
% dFRmax = max(dFRsound,dFRvibrotac);
% dFRmax = max(dFRsound, dFRvibrotac,"ComparisonMethod","abs");
% dFRmono_sign = sign(dFRmax);
% 
% msi = (data2plot - dFRmax)./dFRmax * 100; % msi = ((((dFRmulti_sign * dFRmono_sign(max_index)) * abs(dFRmulti)) - abs(dFRmax)) / abs(dFRmax))*100;
% enhancement_magnitude = ((((dFRmulti_sign .* dFRmono_sign) .* abs(dFRmulti)) - abs(dFRmax)) ./ abs(dFRmax))*100; % multisensory enhancement index = (dFRmulti - dFRmax) / dFRmax
% 
% ax1 = subplot(2,3,1); imagesc(data2plot);ax1.YDir='normal';colorbar
% subplot(2,3,2); plot(data2plot(1,:));title('somato')
% subplot(2,3,3); plot(data2plot(:,1));title('aud')
% ax2 = subplot(2,3,4); imagesc(enhancement_magnitude);ax2.YDir='normal';colorbar

%% select TTLs of session 4
nTrials = size(aligned_spikes, 1);
%nClusters = size(aligned_spikes, 2);
nClusters = 1;
NStim = 1;

% get TTL on and off for a session
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));
for file = 1:9

    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);

    disp(NStim);

    % session to plot
    if ismember(str2double(stimuli_parameters.Par.Set), session)
        % keep Srise and Sfall withing boundaries
        tSrise = Srise(NStim: (NStim + size(stimuli_parameters.Stm, 1)-1));
        tSfall = Sfall(NStim: (NStim + size(stimuli_parameters.Stm, 1)-1));
    end

    % update cummulative stimuli
    NStim = NStim + size(stimuli_parameters.Stm, 1);

end

%% barplot: average response strength per condition
% TO DO: check if works
% data format:
%   firing_mean(freq, amp, condition, :)

savefigs = 0;
conditions = ["OO", "OA", "SO", "SA"];
unitResponses = StimResponseFiring_all.unitResponses;

data_all = StimResponseFiring_all.firing_mean;

% select and format data
%sound_vibro_idx = max(all_stim, sound_vibrotac);
sound_vibro_idx = ismember(StimResponseFiring_all.unitResponses.Cluster, resp_cids);
%data = (mean(data, 1, "omitnan"));
%data = squeeze(mean(mean(data, 1, "omitnan"), 2, "omitnan"));

% make data matrix to plot
data_OO = squeeze(data_all(1,1,1,sound_vibro_idx));
data_OA = squeeze(data_all(1,1,2,sound_vibro_idx)); % if vibrotac
%data_OA = squeeze(data_all(1,1,2,sound_vibro_idx)); % if pressure TO DO
data_SO = squeeze(data_all(3,3,3,sound_vibro_idx)); % 20Hz, 0.3V
data_SA = squeeze(data_all(3,3,4,sound_vibro_idx));

% [maxValue, maxValueIdx] = max(max(data_all(:,:,3,sound_vibro_idx),[],2));
% maxValueIdx = squeeze(maxValueIdx);
% sound_vibro_units = find(sound_vibro_idx);
% data_SO = zeros(length(maxValueIdx),1);
% data_SA = zeros(length(maxValueIdx),1);
% for i = 1:length(maxValueIdx)
%     data_SO(i,1) = squeeze(data_all(maxValueIdx(i),3,3,sound_vibro_units(i))); % max Hz, 0.3V
%     data_SA(i,1) = squeeze(data_all(maxValueIdx(i),3,4,sound_vibro_units(i)));
% end

data = [data_OO'; data_OA'; data_SO'; data_SA']; % max vibrotac
data(5,:) = data(2,:) + data(3,:);
%data = log(data);
%errors = std(data(:, :), 0, 2); % / sqrt(size(data, 2));

ymin = -5;
ymax = 70;
%%plot bar graph
figure;
set(gcf,'position',[500,150,900,700])
hold on
fig = bar(mean(data,2)); % avg FR of all responsive units
x = repmat((1:5)',1,17);
for i = 1:5
    swarmchart(x(i,:), data(i,:), 20, 'k', 'filled','XJitterWidth',0.4);
end

% bar colors
set(fig, 'FaceColor', 'Flat')
fig.CData = [0 0 0; 0.808 0.529 0.666; 0.247 0.635 0.831; 0.580 0.455 0.651; 0.5 0.5 0.5];

% axis
set(gca,'fontsize',18)
%xlabel('conditions')
xticks(1:5)
xticklabels(["control", "sound", "vibrotactile", "multimodal", "sound + vibrotactile"])
ylabel('\Delta Firing rate (spikes/s)')
ylim([ymin ymax])
title(['Condition: ' StimResponseFiring_all.StimType])

hold off

% line graph version
% data format: condition (oo oa so sa oa+so) x units

% TODO: color code responsive units

figure;
subplot(1,2,1)
plot([data_OO,data_OA,data_SA,data_SO]','-o','MarkerSize', 8, 'MarkerFaceColor','k','Color',[.5,.5,.5], 'LineWidth', 1.25)
ylabel('\Delta Firing rate (spikes/s)')
xticks(1:4)
xticklabels(["control", "sound", "multimodal", "vibrotactile"])
ylim([ymin ymax])

if savefigs
    figname = sprintf('M10-11-19-20_stimrespchanges_%s_', StimResponseFiring_all.StimType);
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
end

%figure;
subplot(1,2,2)
plot([data_SA,data_OA+data_SO]','-o','MarkerSize', 8, 'MarkerFaceColor','k','Color',[.5,.5,.5], 'LineWidth', 1.25)
%ylabel('\Delta Firing rate (spikes/s)')
xticks(1:2)
xticklabels(["multimodal", "sound + vibrotactile"])
ylim([ymin ymax])

if savefigs
    figname = sprintf('M10-11-19-20_stimrespchanges_additivity_%s_', StimResponseFiring_all.StimType);
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
end


%% compare firing rate unimodal to multimodal stimulus presentation
% VERSION 1
% makes scatter plot of firing rates per responsive unit between two conditions
% data format:
%   firing_mean(freq, amp, condition, resp_cids)
%   conditions = {"OO", "OA", "SO", "SA"};

% select all units responding to sound and vibration
unitResponses = unitResponses_all;
sound_vibrotac = table2array(((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO) |unitResponses.SA);
%sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
all_stim = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO);
%all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);

% select and format data
sound_vibro_idx = max(all_stim, sound_vibrotac);
data = data_all(:,:,:,sound_vibro_idx);

%data = data_all;
%data = squeeze(mean(mean(data, 1, "omitnan"), 2, "omitnan"));
aud_data = abs(squeeze(data(1,1,2,:)));
som_data = abs(squeeze(mean(data(2:7,3,3,:), 1)));

% for unit = 1:size(data,4)
% 
%     % stat test
%     [p,h,stats] = signrank(aud_data, som_data, 'alpha', 0.01); % different from control trials?
% 
%     if h
%         if abs(aud_trials(unit)) > abs(som_trials(unit))
%             audPref = [audPref; unitResponses.responsive(unit)];
%         elseif abs(aud_trials(unit)) < abs(som_trials(unit))
%             somPref = [somPref; unitResponses.responsive(unit)];
%         elseif abs(aud_trials(unit)) == abs(som_trials(unit))
%             noPref = [noPref; unitResponses.responsive(unit)];
%         end
%     else
%         noPref = [noPref; unitResponses.responsive(unit)];
%     end
% end

% best modality per unit
figure;
hold on
scatter(aud_data, som_data, 'k', "filled")
% add 45 degree angle line
plot([-1, 15],[ -1, 15], 'k')

% format
xlabel('\Delta Firing rate sound trials (spikes/s)')
ylabel('\Delta Firing rate vibrotactile trials (spikes/s)')
axis([-1 15 -1 15])

set(gca,'fontsize',18)

%% compare firing rate unimodal to multimodal stimulus presentation
% VERSION 2
% makes scatter plot of firing rates per responsive unit between two conditions
% data format:
%   firing_mean(freq, amp, condition, :)
%   conditions = {"OO", "OA", "SO", "SA"};

data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids

% compare multi to sound

uSOM = StimResponseFiring.amplitudes;
uAUD = StimResponseFiring.frequencies;

figure
hold on
pos = 1;
for aud = 2:length(uAUD)
    for som = 2:length(uSOM)

        subplot((length(uAUD)-1),(length(uSOM)-1), pos);
        scatter(squeeze(data(aud,1,:)), squeeze(data(aud,som,:)), 'filled', 'k')
        hold on
        xlim(gca, [-10 100])
        ylim(gca, [-10 100])

        xlabel(['\Delta FR sound (' num2str(uAUD(aud)) 'dbSPL)'])
        ylabel(['\Delta FR (' num2str(uAUD(aud)) 'dbSPL *' num2str(uSOM(som)) 'V)'])

        % add 45 degree angle line
        plot([-10, 100],[ -10, 100], 'k')

        pos = pos+1;
    end
end

% compare multi to som pressure

figure
hold on
pos = 1;
for aud = 2:length(uAUD)
    for som = 2:length(uSOM)

        subplot((length(uAUD)-1),(length(uSOM)-1), pos);
        scatter(squeeze(data(1,som,:)), squeeze(data(aud,som,:)), 'filled', 'k')
        hold on
        xlim(gca, [-10 100])
        ylim(gca, [-10 100])

        xlabel(['\Delta FR pressure (' num2str(uSOM(som)) 'V)'])
        ylabel(['\Delta FR (' num2str(uAUD(aud)) 'dbSPL *' num2str(uSOM(som)) 'V)'])

        % add 45 degree angle line
        plot([-10, 100],[ -10, 100], 'k')

        pos = pos+1;
    end
end

%%  compare firing rate unimodal to unimodal stimuli presentation
% for all units plot dFR for each som and aud stimulus combination
% does not give useful plot

uSOM = StimResponseFiring.amplitudes;
uAUD = StimResponseFiring.frequencies;
data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids

figure
hold on
pos = 1;
for aud = 2:length(uAUD)
    for som = 2:length(uSOM)

        subplot((length(uAUD)-1),(length(uSOM)-1), pos);
        scatter(squeeze(data(aud,1,:)), squeeze(data(1,som,:)), 'filled', 'k')
        hold on
        xlim(gca, [-10 100])
        ylim(gca, [-10 100])
        xlabel([num2str(uSOM(som)) ' (V)']);
        ylabel([num2str(uAUD(aud)) ' dB SPL']);

        %axis([-10 90 -10 90])
        plot([-10, 100],[ -10, 100], 'k')

        pos = pos+1;
    end
end

%% quantify reactive units
% 3. firing in relation to phase stimulus (try cycle histogram)
% get and plot analog signal

% read and save into variable
fid = fopen('D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B\continuous.dat');
analog_trace = fread(fid, '*int16');
analog_samples = readNPY('D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B\sample_numbers.npy');

% select which session to analyse
session = 4;

% get start + end of session
idx = sessions_TTLs(:,1) == session;
session_time = sessions_TTLs(idx, 3);

% get first and last sample from session
session_start = analog_samples(analog_samples(:, 1) == session_time(1));
session_end = analog_samples(analog_samples(:, 1) == session_time(end));

figure;
plot(analog_samples(session_start:session_end), analog_trace(session_start:session_end)); % analog trace during session 4

%% ABW --- Start

Fs = 30000;
all_freqs = unique(stimuli_parameters.Stm.SomFreq);


% for cluster = 1:nClusters
cluster = 1;
    %for freq = 1:length(all_freqs)
freq = 3;
        figure;

        index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.1);
        SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
        Var = SOM_Hz;
        yaxistext = [num2str(all_freqs(freq)) ' Hz'];

        % make rasterplot
        [f, YTick, ~, ~, ~, YTickLim] = plotraster(gca, aligned_spikes(index, cluster), Var, [0, 0, 0], [5 5], 1);
        yrange = [min(YTick{end}) - 50, max(YTick{end}) + 50];

        % get analogue signal
        PreT = 0.2; PostT = 1.5; %s
        PreT_samp = -round(PreT*Fs);
        PostT_samp = round(PostT*Fs);
        ttSrise = double(tSrise(index));
        ttSfall = tSfall(index);
        [trialon,~] = find(analog_samples == ttSrise');
        [trialoff,~] = find(analog_samples == ttSfall');

        tt = (PreT_samp:PostT_samp) ./ Fs;
        motorSignal = double(analog_trace(trialon + (PreT_samp:PostT_samp)));
        hold(f,'on');
        plot(f,tt,motorSignal','k');
    %end
% end

% [CycT] = CycTimes(aligned_spikes,StimDur, Mf,SkipVal, SkipMethod);


% ABW --- End

%% ISI histogram

%StimResponseFiring = StimResponseFiring_all;

conditions = StimResponseFiring.conditions; % OO OA SO SA
uparamA = StimResponseFiring.amplitudes;
uparamB = StimResponseFiring.frequencies;
nparamA = length(uparamA);
nparamB = length(uparamB);
ISICVs = nan(nparamA, nparamB, length(conditions), length(cids)); % most freq ISI %amp,freq,con,cluster

binsize = 0.05; %50ms
%cluster = 1;
for cluster = 1%close all:length(StimResponseFiring.cids)
    for con = 1:length(conditions)
        for amp = 1:nparamA

            figure;
            hold on
            pos = 1;

            for freq = 1:nparamB         

                %spikeISI 4D cell array (amp, freq, con, cluster)
                toplot = spikeISI{amp, freq, con, cluster};
                [N, edges] = histcounts(toplot, -7:binsize:0); %-7: 0.001 sec ISI
                subplot(3, 3, pos)
                for i = 1:(length(edges)-1)
                    binCenters(:,i) = mean(edges(:,i:i+1));
                    % binCenters(:,i) = sqrt(edges(:,i)* edges(:,i+1)); % geometric mean
                end
              
                % plot ISI histograms
                bar(binCenters, N, 1, 'FaceColor', 'k')
                xline(log(1/uparamB(freq)),'r--');
                %xline(log(1/uFreq(freq)/2),'b--');
                %xline(log(0.75*1/uFreq(freq)),'k--');
                %xline(log(0.25*1/uFreq(freq)),'g--');
                % stairs(edges(1:end-1), N, 'Color', [0.5 0.5 0.5])
                xticks([log(0.001), log(0.01), log(0.1), log(1)])
                set(gca,'TickDir','out')
                xLab = [0.001, 0.01, 0.1, 1];
                xticklabels(xLab)
                xlabel('ISI (sec)')
                ylabel('# spikes')
                %ylim([0,max(N)+10])
                ylim([0,30])
                %legend;
                title(figTitle)
                clear binCenters

                pos = pos+1;

            end

           sgtitle(figSGTitle)

        end

    end

    % figname = sprintf('M%.2i_cluster%i_VxApre', animal, cids(cluster));
    % saveas(gcf, fullfile(OutPath, figname));
    % saveas(gcf, fullfile(OutPath, [figname '.jpg']));

end

% plot CVs
figure

for cluster = 1:length(cids)
    hold on
    scatter(0, squeeze(ISICVs(1,1,1,cluster)), 'r')
    scatter(-10, squeeze(ISICVs(1,1,2,cluster)), 'r')
    scatter(uparamB(2:8),squeeze(ISICVs(3,2:8,4,cluster)), 'c') %SA, 0.3V, all freqs, all units
    scatter(uparamB(2:8),squeeze(ISICVs(3,2:8,3,cluster)), 'k') %SO, 0.3V, all freqs, all units

end
ylabel('CV')
legend('cyan: multimodal, black: tactile')
%set(gca,'Xscale','log')

%% -------------------------- 8. Motor signal  -------------------------- %%
% plot analog motor signal per condition
for freq = 1:length(all_freqs)
    index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.3);

    % figure;
    ttSrise = tSrise(index);
    ttSfall = tSfall(index);

    figure;
    hold on

    for i = 1:length(ttSrise)
        trialon = find(analog_samples == ttSrise(i));
        trialoff = find(analog_samples == ttSfall(i));
        tanalog_samples = trialon:trialoff;

        plot(double(analog_samples(tanalog_samples)-analog_samples(trialon))./Fs, analog_trace(tanalog_samples));

    end

    xlabel('Samples (normalized)')

    
    title([num2str(all_freqs(freq)) 'Hz (session '  stimuli_parameters.Par.Set ')']);

end

%% plot difference in FR between no sound & multimodal condition
% data selection

load('\\store\department\neuw\shared\Aaron Wong\Data\ProcessedData\Blom\Processed\M10-11-19-20\M10-11-19-20_PxA_P_pre_UnitResponses.mat');

cids = StimResponseFiring_all.unitResponses.Cluster;
%resp_cids = [19153,19265,19287,19296,20277,20290,20303,20306];
resp_cids = [19153,19265,19296,20277,20290,20303,20306];

index = ismember(cids, resp_cids); %max(max(sound_resp_units, vibrotac_resp_units), multi_resp_units);

% variables
uAmp = unique(StimResponseFiring_all.amplitudes)';
umN = round((uAmp * 0.158)*1000); %0.158N/V callibation
uInt = StimResponseFiring_all.frequencies(:,end);

% select data
tdata = StimResponseFiring_all.firing_mean;
data = NaN(size(tdata, 1), size(tdata, 2), sum(index)); % initiate data and fill
data(1,1,:) = tdata(1,1,1,index);
data(2:end, 1, :) = tdata(2:end,1,2,index);
data(1, 2:end, :) = tdata(1,2:end,3,index);
data(2:end, 2:end, :) = tdata(2:end,2:end,4,index);

% Remove non responsive units
data = data(:,:,~(all(isnan(data) | isinf(data), [1 2])));

fontsize = 14;

% plot difference in FR between no sound & multimodal condition

% median FSL
%fr_median = median(data, 3); % size: [length(uInt) x length(umN)]
fr_median = mean(data,3);

% Compute difference from sound-only condition
fr_diff = fr_median(:, 2:end) - fr_median(:, 1); % subtract sound-only from each pressure condition

% Plot each BBN intensity as a separate line
colors = ([0 0 0;...
    0.6549019607843137 0.615686274509804 0.9254901960784314;...
    0.592156862745098 0.43137254901960786 0.8196078431372549; ...
    0.42745098039215684	0.27450980392156865	0.6313725490196078; ...
    0.24705882352941178	0.1568627450980392	0.3058823529411765]); % Distinct colors for each line

% plot
figure;
hold on

for i = 1:length(uInt)
    % Plot mean line
    plot(umN(2:end), fr_diff(i,:), '-o', ...
        'LineWidth', 3, ...
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), ...
        'DisplayName', sprintf('%d dB SPL', uInt(i)));
end

% Customize axes and labels
title('relative FR change')
xlabel('Pressure intensity (mN)');
ylabel('relative \Delta frirng rate (spikes/s)');
xticks(umN);
set(gca, 'FontSize', 18);
hold off;

%% heatmap: calculate multisensory integration indexes

% Only works with PxA (not VxA) stimuli

% basic formulas:
% preference_index = (dFRsom - dFRaud) / (dFRsom + dFRaud)
% additivity_index = (dFRmulti - dFRsom - dFRaud) / (dFRsom + dFRaud)
% enhancement_magnitude = (dFRmulti - dFRmax) / dFRmax
% modulation_index = dFRmulti - (dFRsound + dFRsom)

% TO DO:
% remove data if FR of th modalities < 1Hz
% add FR heatmap

%StimResponseFiring = StimResponseFiring_all;

% add monosensory FR to same matrix
%data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids

signvalue = sign(data); % get sign of dFR
signvalue(signvalue(:,:,:) == 0) = 1; % make sign 1 if dFR=0

for c = 1:size(data,3)
    for i = 2:size(data,1)
        for j = 2:size(data,2)

            dFRsound = data(i,1,c);
            dFRvibrotac = data(1,j,c);
            dFRmulti = data(i,j,c);
            [dFRmax, max_index] = max([dFRsound, dFRvibrotac],[],"ComparisonMethod","abs");

            % NaN if FR of both modalities < 1Hz
            if abs(dFRsound) <= 1 && abs(dFRvibrotac) <= 1
                dFRsound = NaN;
                dFRvibrotac = NaN;
                dFRmulti = NaN;
                dFRmax = NaN;
            end

            dFRmono_sign = [signvalue(i,1,c), signvalue(1,j,c)];
            dFRmulti_sign = signvalue(i,j,c);
            %dFRmax = max(abs(data(1,j,c)), abs(data(i,1,c)));
            %dFRmax = max(data(1,j,c),data(i,1,c),"ComparisonMethod","abs");

            % calculate multimodal enhancement index for each stim combination
            %enhancement_magnitude(i,j,c) = (dFRmulti - dFRmax) / dFRmax; % multisensory enhancement index
            enhancement_magnitude(i,j,c) = ((((dFRmulti_sign * dFRmono_sign(max_index)) * abs(dFRmulti)) - abs(dFRmax)) / abs(dFRmax))*100; % multisensory enhancement index = (dFRmulti - dFRmax) / dFRmax
            modulation_index(i,j,c) = dFRmulti - (dFRsound + dFRvibrotac); % multisensory modulation index, sign indicated direction multi
            additivity_index(i,j,c) = ((dFRmulti - dFRvibrotac - dFRsound) / abs(dFRvibrotac + dFRsound))*100;
            preference_index(i,j,c) = (dFRvibrotac - dFRsound) / (abs(dFRvibrotac) + abs(dFRsound)); % - pref sound; + pref vibrotac

        end
    end
end

%% add NaN to irrelevant combinations (unisensory except control)
enhancement_magnitude(1,2:end,:) = NaN;
enhancement_magnitude(2:end,1,:) = NaN;

modulation_index(1,2:end,:) = NaN;
modulation_index(2:end,1,:) = NaN;

additivity_index(1,2:end,:) = NaN;
additivity_index(2:end,1,:) = NaN;

preference_index(1,2:end,:) = NaN;
preference_index(2:end,1,:) = NaN;

NaNindex_enhancement = ~isnan(enhancement_magnitude) & ~isinf(enhancement_magnitude);
NaNindex_modulation = ~isnan(modulation_index) & ~isinf(modulation_index);
NaNindex_additivity = ~isnan(additivity_index) & ~isinf(additivity_index);
NaNindex_preference = ~isnan(preference_index) & ~isinf(preference_index);

clear i j c

%% Plotting multisensory integration measures
% heatmaps

for cluster = 1:size(data, 3) % loop through clusters
    figure('Position',[100,100,1800,500]);
    colormap('parula'); % Choose a color map
    sgtitle(['unit ' num2str(resp_cids(cluster))])

    % delta FR
    subplot(1,3,1);
    imagesc(data(:,:, cluster))
    cb = colorbar(gca, 'eastoutside');
    cb.Label.String = '	\Delta spike rate (spikes/sec)';
    cb.FontSize = fontsize;
    set(gca, 'YDir', 'normal', 'XTick',1:length(umN),'XTickLabel',umN','YTick',1:length(uInt),'YTicklabel',uInt)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')
    set(gca,'fontsize',fontsize)
    title('Average firing rate')

    % enhancement index
    subplot(1,3,2);
    imagesc(enhancement_magnitude(:,:, cluster),'AlphaData', NaNindex_enhancement(:,:,cluster))
    cb = colorbar(gca, 'eastoutside');
    %clim([minV, maxV]);
    cb.Label.String = 'multisensory index (%)'; %cb.Label.String = 'multisensory enhancement index';
    cb.FontSize = fontsize;
    set(gca, 'YDir', 'normal', 'XTick',1:length(umN),'XTickLabel',umN,'YTick',1:length(uInt),'YTicklabel',uInt)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')
    set(gca,'fontsize',fontsize)
    title('Multisensory modulation')

    % additivity index
    subplot(1,3,3);
    imagesc(additivity_index(:,:, cluster),'AlphaData', NaNindex_additivity(:,:,cluster))
    cb = colorbar(gca, 'eastoutside');
    cb.Label.String = 'additivity index (%)';
    cb.FontSize = fontsize;
    set(gca, 'YDir', 'normal', 'XTick',1:length(umN),'XTickLabel',umN,'YTick',1:length(uInt),'YTicklabel',uInt)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')
    set(gca,'fontsize',fontsize)
    title('Multisensory additivity')
 
    % % save responsive units only
    %  if ismember(StimResponseFiring.unitResponses.Cluster(cluster), resp_cids)
    %      sgtitle(['Responsive unit (' num2str(StimResponseFiring.unitResponses.Cluster(cluster)) ')'])
    %      figname = sprintf('M10-11-19-20_unit-%i_multisensory-indexes_%s', StimResponseFiring.unitResponses.Cluster(cluster), StimResponseFiring.StimType);
    %      saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
    %      saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
    %  else
    %      sgtitle(['Non responsive unit (' num2str(StimResponseFiring.unitResponses.Cluster(cluster)) ')'])
    %      close(gcf)
    %  end

end


%% Plotting multisensory integration measures
% line graphs
close all
for cluster = 1:size(data, 3) % loop through clusters
    figure;
    sgtitle(['unit ' num2str(resp_cids(cluster))])
    ymin = min(enhancement_magnitude(2:end,:, cluster), [], 'all') * 1.1;
    ymax = max(enhancement_magnitude(2:end,:, cluster), [], 'all') * 1.1;


    for soundint = 2:size(data, 1)
        subplot(length(uInt)-1,1,soundint-1)

        bar(enhancement_magnitude(soundint,:,cluster))

        xlabel('pressure (mN)')
        xticklabels(umN)
        title([num2str((uInt(soundint))) 'dB SPL'])
        ylabel('MSI')

        %ylim([0 max(enhancement_magnitude(soundint,:,1)) * 1.1]); % Optional: consistent y-axis
        ylim([ymin ymax])

    end

end

%% plot average modulation
%close all

% Compute average across clusters
%mean_enhancement = mean(enhancement_magnitude, 3); % size: [soundint x pressure]
mean_enhancement = median(enhancement_magnitude, 3); % size: [soundint x pressure]
nClusters = size(enhancement_magnitude, 3);

figure;
sgtitle('median MSI');

for soundint = 2:size(data, 1)

    % Determine y-axis limits per sound intensity
    ymin = floor(min(enhancement_magnitude(soundint,:,:), [], 'all'));
    ymax = ceil(max(enhancement_magnitude(soundint,:,:), [], 'all'));

    subplot(2,2,soundint-1)
    title([num2str(uInt(soundint)) ' dB SPL'])

    % Plot average bar
    bar(mean_enhancement(soundint,:), 'FaceColor', [0.3 0.6 0.8]);
    hold on;

    % Overlay individual cluster data points
    for cluster = 1:nClusters
        scatter(1:length(umN), enhancement_magnitude(soundint,:,cluster), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

    xlabel('Pressure (mN)');
    xticklabels(umN);
    title([num2str(uInt(soundint)) ' dB SPL']);
    ylabel('MSI');
    ylim([ymin ymax]);
    set(gca, 'FontSize', 12);
    hold off;
end


