%% Single unit plotting
% plotting of single cluster data such as:
% - FRA analysis & heatmaps
% - raster plots
% - example traces
% - tuning curves / rate level functions

%% ---------------------- Set directories  ---------------------- %%
% aligned spike data

close all
clearvars


%% ----------------------- FRA analysis & Plotting ----------------------- %%
% output: FRA & MedFSL 4D: intensity, frequency, set number, cluster

% load FRA session(s) aligned spikes
aligned_spikes_files = dir(fullfile(OutPath, '*FRA_AlignedSpikes.mat'));

for file = 1:size(aligned_spikes_files, 1)

    % load each FRA aligned spikes file
    %aligned_spikes = load([aligned_spikes_files(file).folder '\' aligned_spikes_files(file).name]);
    % load corresponding stimuli file
    %session = aligned_spikes_files(file).name(6:7);
    %stim_files = dir(fullfile(BehaviorPath, ['*_S' session '_*.mat']));
    %stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

    session = str2double(aligned_spikes_files(file).name(6:7));
    [cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

    % FRA analysis saves heatmap figures
    FSL = 0;
    heatmap = 1;
    FRAanalysis(stimuli_parameters, aligned_spikes, cids, OutPath, heatmap, FSL);

end

%% ---------------- OVERVIEW RASTER PLOT ---------------- %%
% plotting single sessions - in use
% raster plot of single units
% saves figures to OutPath

close all

session = 2;
[cids, stimuli_parameters, aligned_spikes, ~, ~, sessions_TTLs, onsetDelay, ~, clusterinfo] = loadData(OutPath, session, BehaviorPath);
plotResponses(stimuli_parameters, aligned_spikes, cids, OutPath);

%% add spacing where needed
% never needed anymore, automatically done with loadData.m

idx = find(ismember(stimuli_parameters.Stm.MMType,["SO","OO"]));
for cluster = 1:length(cids)
    for ii = idx'
        aligned_spikes{ii,cluster} = aligned_spikes{ii,cluster} + 0.25;
    end
end

%% plot raster + PSTH of single SOM stimulus type - keep
% plots stimulus shape + raster + psth
% per cluster pr vibration frequency seperatly

close all
condition = 'SO';
saveplots = 0; %0 don't save, 1 save plots in OutPath
cidstoplot = [121];
SOMplotting(stimuli_parameters, aligned_spikes, cidstoplot, cids, OutPath, saveplots, condition);

%% plot single unit - in use --> single unit
% raster plot + PSTH

% define shared x-lim parameters
preT  = -0.2;
postT = 1.2;
xrange = [preT, postT];
binsize = 0.01;
start_aud = 0;
start_som = 0.25;
end_som = 0.75;
end_aud = 1;
xlinerange = [start_aud start_som end_som end_aud];

index = (stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq ~= 600);
SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
SOM_Amp = stimuli_parameters.Stm.Amplitude(index);
StimType = stimuli_parameters.Stm.Var25(index);
Var = [SOM_Hz, StimType];
%Var =  [stimuli_parameters.Stm.Var25, stimuli_parameters.Stm.SomFreq, stimuli_parameters.Stm.Amplitude];
raster_yinc = [5,10,10];
raster_color = [0, 0, 0];

for cluster = 1:length(cids)
%fig = figure; % rasterplot 
figure;
fig = subplot(3,1,1:2); % rasterplot
ax = gca;

[f, YTick, ~, ~, ~, YTickLim] = plotraster(ax, aligned_spikes(index, cluster), Var, raster_color, raster_yinc, 1);

yticks(YTick{1});
yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];
ylim(f,yrange);
xlim(f,xrange);
xline(xlinerange) % on/off set

for i = 1:size(YTickLim,1) % delimit groups
    yline(ax,YTickLim(i,1)-3,':k');
    yline(ax,YTickLim(i,2)+3,':k');
end

%xlabel('Time (s)')
all_freqs = unique(stimuli_parameters.Stm.SomFreq);
yaxislabels = {'Control',all_freqs(2:8)};
yticklabels(yaxislabels)
ylabel('Vibrotactile frequency (Hz)')
ax.FontSize = 16;

% make histogram / PSTH
binsize = 0.02;

fig = subplot(3,1,3); 
hold on
% select groups for hist
OO = stimuli_parameters.Stm.Var25 == 1;
OA = stimuli_parameters.Stm.Var25 == 4;
SO = (stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq ~= 600) & stimuli_parameters.Stm.Var25 == 3;
SA = (stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq ~= 600) & stimuli_parameters.Stm.Var25 == 2;

[N, edges] = histcounts(vertcat(aligned_spikes{OO, cluster}), preT:binsize:postT); % OO
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(OO))/binsize)) % spike/s
%plot(edges(1:end-1), ((N/sum(OO))/binsize), 'Color', '#D95319', 'LineWidth',1.5) % spike/s
plot(edges(1:end-1), ((N/sum(OO))/binsize), 'k','LineWidth',2.5) % spike/s

[N, edges] = histcounts(vertcat(aligned_spikes{OA, cluster}), preT:binsize:postT); % OA
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(OA))/binsize)) % spike/s
plot(edges(1:end-1), ((N/sum(OA))/binsize), 'Color', "#c21069",'LineWidth',2.5) % spike/s

[N, edges] = histcounts(vertcat(aligned_spikes{SO, cluster}), preT:binsize:postT); % SO
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(SO))/binsize)) % spike/s
plot(edges(1:end-1), ((N/sum(SO))/binsize), 'Color', "#3eb5f0",'LineWidth',2.5) % spike/s

[N, edges] = histcounts(vertcat(aligned_spikes{SA, cluster}), preT:binsize:postT); % SA
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(SA))/binsize)) % spike/s
plot(edges(1:end-1), ((N/sum(SA))/binsize), 'Color', "#9475CB",'LineWidth',2.5) % spike/s

%format axis
xlabel('Time (s)')
ylabel('Spike rate (spikes/s)')
xline(xlinerange) % on/off set
legend('control', 'sound only', 'vibrotactile only', 'multimodal', 'Location', 'northeast')

ax2 = gca;
ax2.FontSize = 16;
xlim(ax2,xrange);
%ylim(ax2, [0 60])

title(cids(cluster))


end
