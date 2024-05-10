clearvars
% load behavioural stimuli
filename = 'M9_S03_SxA';
%load([filename,'.mat']);
load(['D:\DATA\Behavioral Stimuli\M9\' filename,'.mat'],'Stm','Aud', 'Par');

%% load microphone data
filename = [filename '_Sound'];
load(['D:\DATA\Behavioral Stimuli\M9\Sound recording\' filename,'.mat'], 'Aud', 'Fs');
%load([filename,'.mat'],'Aud','Stm','Par');

%% Format Aud, if necessary

%Select correct trials from Aud
AudTemp = Aud;
Aud = AudTemp(Stm.AudDur==0);

%% Format Aud (cell array --> matrix)
AudTemp = [];
for i = 1:size(Aud,2)
    AudTemp(i,:) = Aud{i};
end
Aud = AudTemp;

%%
nSamp = size(Aud,2);
sound = Aud(1,:);
%
window_dur = 0.02 ;% s
window = round(window_dur*Fs);
noverlap = round(0.9*window);
% spectrogram(sound,window,noverlap,[],Fs,'yaxis');

%
numTrials = size(Aud,1);
tic;
[s,f,t,ps_all] = spectrogram(Aud(1,:),window,noverlap,[],Fs,'yaxis');
fprintf('%4d/%-4d\r',1,numTrials);
for ii = 2:numTrials
    [~,~,~,ps] = spectrogram(Aud(ii,:),window,noverlap,[],Fs,'yaxis');
    ps_all = cat(3,ps_all,ps);
    fprintf('\b\b\b\b\b\b\b\b\b\b');
    fprintf('%4d/%-4d\r',ii,numTrials);
end
toc;

[nF,nT,~] = size(ps_all);

% adapt to select no sound only trials (for SxA sessions)
UStim = unique(Stm(:,{'SomFreq','Amplitude'}),'rows');
nUStim = size(UStim,1);

ps_mean = nan(nF,nT,nUStim);
for ii = 1:nUStim
    idx = Stm.SomFreq == UStim.SomFreq(ii) & Stm.Amplitude == UStim.Amplitude(ii);
    ps_mean(:,:,ii) = mean(ps_all(:,:,idx),3);
end

%% Calculate instantaneous power

fIdx = f > 10;
fLow = f > 500 & f <= 2000;
fMid = f > 2000 & f <= 10000;
fHigh = f > 10000 ;

instPower_raw = reshape(sum(ps_mean(fIdx,:,:)),[nT,nUStim]);
instPower_raw_dB = 10*log10(instPower_raw);
instPower_raw_dBSPL = dbv2spl(instPower_raw_dB);

instPower_low = reshape(sum(ps_mean(fLow,:,:)),[nT,nUStim]);
instPower_low_dB = 10*log10(instPower_low);
instPower_low_dBSPL = dbv2spl(instPower_low_dB);

instPower_mid = reshape(sum(ps_mean(fMid,:,:)),[nT,nUStim]);
instPower_mid_dB = 10*log10(instPower_mid);
instPower_mid_dBSPL = dbv2spl(instPower_mid_dB);

instPower_high = reshape(sum(ps_mean(fHigh,:,:)),[nT,nUStim]);
instPower_high_dB = 10*log10(instPower_high);
instPower_high_dBSPL = dbv2spl(instPower_high_dB);

%% plot average
figure('Position',[10,10,1000,800]);
clim = [-110,-60];
freqRange = [0,40];%kHz
dBRange_pow = [-Inf,Inf];
dBRange_pow = [-2,60];
tRange = [0,0.6];

nRows = floor(sqrt(nUStim));
nCols = ceil(nUStim / nRows);

for ii = 1:nUStim
    ax1 = subplot(nRows,nCols,ii);
    yyaxis(ax1,'left')
    imagesc(t,f./1000,squeeze(10*log10(ps_mean(:,:,ii))),clim)
    ax1.YDir = 'normal';
    ylim(freqRange);
    ylabel(ax1,'Frequency (kHz)')
    xlabel(ax1,'Time (s)')
    title([num2str(UStim.SomFreq(ii),'%d Hz'),'   ',num2str(UStim.Amplitude(ii),'%.3f V  ')])
    yyaxis(ax1,'right')
    plot(ax1,t,instPower_raw_dBSPL(:,ii),'w-'); hold(ax1,'on'); 
    plot(ax1,t,instPower_low_dBSPL(:,ii),'-','Color',[1,.5,.5])
    plot(ax1,t,instPower_mid_dBSPL(:,ii),'-','Color',[1,0.1,0.1])
    plot(ax1,t,instPower_high_dBSPL(:,ii),'-','Color',[.5,0,0]); hold(ax1,'off');
   % legend(ax1,{'all (>10Hz)','low (500-2000Hz)','mid (2-10kHz)','high (>10kHz)'},'Location','best');
    ylabel(ax1,'Instantaneous power (dB SPL)')
    ylim(dBRange_pow)
end

% make space for legend in subplot
ax = subplot(nRows,nCols,ii+1,'Visible','off');
axPos = ax.Position;
hL = legend(ax1,{'all (>10Hz)','low (500-2000Hz)','mid (2-10kHz)','high (>10kHz)'});
hL.Position(1:2) = axPos(1:2); % move legend to position of extra axes

sgtitle(['Ramp: ',num2str(Stm.Ramp(1)),'ms','   ',...
         'Actuator: ',Stm.Actuator(1,:),    '   ', ...
         'Waveform: ',Stm.Waveform(1,:)])

saveas(gcf,[filename,'_spectogram.png'])
%% plot individual
figure('Position',[10,10,1400,900]);
stimToPlot = [1,4:4:numTrials];
nToPlot = length(stimToPlot)+1;
clim = [-110,-60];
freqRange = [0,10];%kHz
dBRange_pow = [-Inf,Inf];
tRange = [0,0.6];

nRows = floor(sqrt(nToPlot));
nCols = ceil(nToPlot / nRows);

for ii = 1:nToPlot
    ax1 = subplot(nRows,nCols,ii);
    yyaxis(ax1,'left')
    if ii < nToPlot
    idx = stimToPlot(ii);
    imagesc(t,f./1000,squeeze(10*log10(ps_all(:,:,idx))),clim)
    title([num2str(Stm.Rep(idx),'Rep: %d ')])
    else
    imagesc(t,f./1000,squeeze(10*log10(ps_mean(:,:,1))),clim)
    title(['mean'])
    end
    ax1.YDir = 'normal';
    ylim(freqRange);
    ylabel(ax1,'Frequency (kHz)')
    xlabel(ax1,'Time (s)')
end

sgtitle(['Ramp: ',num2str(Stm.Ramp(1)),'ms','   ',...
         'Actuator: ',Stm.Actuator(1,:),    '   ', ...
         'Waveform: ',Stm.Waveform(1,:)])
saveas(gcf,[filename,'_spectogram_trials.png'])

%% local functions
function dbspl = dbv2spl(dbv)
mic = '7012';
switch mic
    case '7012'
        micSens = 16.4e-3;  % V / 1 Pa; 16.4mV/Pa ~ -35.7 dBV @ 1Pa
    case '7016'
        micSens = 3.89e-3;  % V / 1 Pa; 3.89mV/Pa ~ -48.2 dBV @ 1Pa
    otherwise % 7016
        micSens = 3.89e-3;  % V / 1 Pa; 3.89mV/Pa ~ -48.2 dBV @ 1Pa
end
    micDBV  = 20*log10(micSens); 
    refdB   = 94;       % dB SPL = 1Pa
    postAmp = 40;       % dB

    dbspl = dbv - micDBV + refdB - postAmp;
end