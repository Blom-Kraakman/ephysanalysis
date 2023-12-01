function Wav = ChannelWaveForms(F,KSPath,cid,cpos,nWav)
% a function to load and plot raw waveforms.

% --- default parameter(s) ---
if isempty(nWav);    nWav    =   200; end
% ----------------------------

% --- Channel map ---
% ChMap       =   1+ [26    24    22    19    21    30    28    17    29    31    14    20    18    23    25    27     1     3     4 6     8    10    12    16    15    13    11     9     7     5     2     0    63    61    58    56    54    52 50    46    49    51    53    55    57    59    60    62    36    38    40    45    43    32    34    47    35 33    48    42    44    41    39    37];
% load('W:\KS\configfiles\CambridgeH3chronic.mat','chanMap');
load('C:\Users\TDT\Documents\GitHub\KilosortSettings\Kilosort3Config', 'OE_Intan_cambridgeneurotech_ASSY-77-H10_67ch_disconnected');
ChMap = chanMap;
% -------------------

% --- Read info of individual (sorted) spikes ---
% all spikes
ids     =   readNPY([KSPath '\spike_clusters.npy']);
ST      =   readNPY([KSPath '\spike_times.npy']);
amp     =   readNPY([KSPath '\amplitudes.npy']);

% select spike with specified cluster ID "cid"
sel     =   ids == cid;
ST      =   ST(sel);
ids     =   ids(sel);
amp     =   amp(sel);
% -----------------------------------------------

% --- Parameters ---
samp    =   -30:31; % number of samples to plot (time);
NAdjCh  =   2; % number of adjacent channels (each side) to plot
sampRate = 30000;
t       =   (samp/sampRate)*1000;
nSamp   = length(samp);

gwfparams = struct;

gwfparams.dataDir = 'KSPath';    % KiloSort/Phy output folder
%gwfparams.fileName = [Mouse '.dat'];         % .dat file containing the raw
%  gwfparams.dataDir = 'W:\temp';    % KiloSort/Phy output folder
gwfparams.fileName = 'continuos.dat';   % 'temp_wh.dat';         % .dat file containing the raw
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 64;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [samp(1) samp(end)];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = nWav;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    ST; % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = ids; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

wf = getWaveForms(gwfparams,ChMap);

mWav = wf.waveFormsMean; mWav = reshape(mWav,[64 62]); mWav = 0.195*mWav;
Wav = wf.waveForms; Wav = reshape(Wav,[nWav,64,62]); Wav = 0.195*Wav;

%plots

%make subplot handles
%  w1 = subplot(4,4,1);
%  w2 = subplot(4,4,5);
%  w3 = subplot(4,4,9);
w1 = subplot(4,4,[1,9]);

a = subplot(4,4,13:16);
plot(a,(ST/30000),amp,'k.'); ylabel(a,'amplitude (a.u)'); xlabel(a,'time re: experiment');
set(a,'FontName','Arial','FontSize',12);
ylim(a,[0,Inf]);

W = subplot(4,4,[2.5:4 6.5:8 10.5:12]);

ch2plot = cpos + (-NAdjCh:NAdjCh);
ch2plot = ch2plot(ch2plot>0 & ch2plot <= 64);
mn = min(mWav(ch2plot,:),[],'all');
mx = max(mWav(ch2plot,:),[],'all');
offset = 1.2*(mx-mn);
plot(w1,t,mWav(ch2plot,:)' + ch2plot*offset,'k-')

%  if cpos == 1
%      plot(w1,t,mWav(cpos+1,:),'k-'); title(w1,['E' num2str(cpos+1)]);
%      plot(w2,t,mWav(cpos,:),'k-'); title(w2,['E' num2str(cpos)]);
%  elseif cpos == 64
%      plot(w2,t,mWav(cpos,:),'k-'); title(w2,['E' num2str(cpos)]);
%      plot(w3,t,mWav(cpos-1,:),'k-');title(w3,['E' num2str(cpos-1)]);
%  else
%      plot(w1,t,mWav(cpos+1,:),'k-'); title(w1,['E' num2str(cpos+1)]);
%      plot(w2,t,mWav(cpos,:),'k-'); title(w2,['E' num2str(cpos)]);
%      plot(w3,t,mWav(cpos-1,:),'k-');title(w3,['E' num2str(cpos-1)]);
%  end

%  ylabel(w2,'Vpp (\muV)');
% yticks(w1,ch2plot*offset);yticklabels(ch2plot);
set(w1,'YTick',ch2plot*offset,'YTickLabel',ch2plot)
xlabel(w1,'time (ms)');
set(w1,'FontName','Arial','FontSize',12);
%  set(w2,'FontName','Arial','FontSize',12);
%  set(w3,'FontName','Arial','FontSize',12);
%  linkaxes([w1 w2 w3],'y');
%  mn = min(mWav(cpos,:));
%  mx = max(mWav(cpos,:));
%  df = diff([mn mx]);
%  ylim(w2,[mn-(df/10) mx+(df/10)]);

Wo = reshape(Wav(:,cpos,:),[nWav, nSamp]);
for k=1:nWav
    w=Wo(k,:) - (mean(Wo(k,1:20)));
    plot(t,w,'k-','LineWidth',0.5);
    hold on;
end
xlim(W,[min(t) max(t)]);
set(W,'FontName','Arial','FontSize',12);
ylabel('Vpp (\muV)');
xlabel('time (ms)');

sgtitle([Mouse ' - unit ' num2str(cid)],'FontName','Arial','FontSize',25);
end