function fig = plot_in_channel_map(KSPath, cpos)

channel_map = readNPY([KSPath 'channel_map.npy']);
channel_positions = readNPY([KSPath 'channel_positions.npy']);

xcoords = channel_positions(:,1);
ycoords = channel_positions(:,2);
margin = 50;
xMin = min(xcoords); xMax = max(xcoords);
xSpan = xMax-xMin; xMid = 0.5*(xMax+xMin);
yMin = min(ycoords); yMax = max(ycoords);
ySpan = yMax-yMin; yMid = 0.5*(yMax+yMin);
maxSpan = max(xSpan,ySpan);

fig = scatter(xcoords, ycoords, ".", 'k');
hold on;

for i = 1:length(channel_map)
    text((channel_positions(i,1)+1), (channel_positions(i,2)+1),num2str(channel_map(i)))
end

% add position of units included in analysis
hold on
units = table2array(cpos(:,2));
idx = ismember(channel_map, units);
scatter(channel_positions(idx,1), channel_positions(idx,2), 'o', 'XJitter', 'density'); 

% to do: reflect number of clusters on 1 channel. apply jitter in circle
% marking to achieve this
% to do: add unit number with corresponding channel
%text((channel_positions(unique(idx),1)+1), (channel_positions(unique(idx),1)+1),num2str(channel_map()))

% format figure
title('Channel map')
ylabel('Relative depth')
fig.SizeData = 100;
fontsize(13,"points")
axis square
xlim([xMid - 0.5*maxSpan - margin, xMid + 0.5*maxSpan + margin]);
ylim([yMid - 0.5*maxSpan - margin, yMid + 0.5*maxSpan + margin]);
xticks(unique(xcoords));
yticks(unique(ycoords));

hold off;
end