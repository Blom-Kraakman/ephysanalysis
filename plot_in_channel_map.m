function [figa, figb, channelpos] = plot_in_channel_map(KSPath, cpos)

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

figa = scatter(xcoords, ycoords, ".", 'k');
hold on;

for i = 1:length(channel_map)
    text((channel_positions(i,1)+1), (channel_positions(i,2)+1),num2str(channel_map(i)))
end

% add position of units included in analysis
hold on
unitspos = table2array(cpos(:,2));
idx = ismember(channel_map, unitspos);
scatter(channel_positions(idx,1), channel_positions(idx,2), 'o', 'XJitter', 'density'); 

% format figure
title('Channel map')
ylabel('Relative depth')
figa.SizeData = 100;
fontsize(13,"points")
axis square
xlim([xMid - 0.5*maxSpan - margin, xMid + 0.5*maxSpan + margin]);
ylim([yMid - 0.5*maxSpan - margin, yMid + 0.5*maxSpan + margin]);
xticks(unique(xcoords));
yticks(unique(ycoords));

hold off;

% to do: reflect number of clusters on 1 channel. apply jitter in circle
% marking to achieve this
% to do: add unit number with corresponding channel
%text((channel_positions(unique(idx),1)+1), (channel_positions(unique(idx),1)+1),num2str(channel_map()))

%fig = scatter(1:length(cids)/10, cpos.depth', ".", 'k');
%xticks(unique(xcoords));

figure;
figb = scatter(channel_positions(idx,1), channel_positions(idx,2), '.k');
%units = cpos.id';
channelpos = [channel_positions(idx,1), channel_positions(idx,2), channel_map(idx)];
for unit_to_plot = 1:length(channelpos)
    index = cpos.channel == channelpos(unit_to_plot, 3);
    text(channelpos(index,1)+1, channelpos(index,2),num2str(cpos.id(index)))
end

title('Channel map')
ylabel('Relative depth')
figb.SizeData = 100;
fontsize(13,"points")
axis square
xlim([-10 30])

end