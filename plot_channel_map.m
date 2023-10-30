% 27/09/2023

channel_map = readNPY('D:\DATA\EphysRecordingsSorted\M01_S07\kilosort3\channel_map.npy');
channel_positions = readNPY('D:\DATA\EphysRecordingsSorted\M01_S07\kilosort3\channel_positions.npy');

fig = scatter(channel_positions(:,1), channel_positions(:,2), ".", 'k');
fig.SizeData = 100;
hold on;

for i = 1:64
    text((channel_positions(i,1) +1) ,(channel_positions(i,2) +1),num2str(channel_map(i)))
end


title('Channel map')
xlim([-10 200])
ylabel('Relative depth')

