clear 
close all

opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = [4, Inf];

opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["coords", "x", "y", "likelihood", "x1", "y1", "likelihood1", "x2", "y2", "likelihood2"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ImportErrorRule = "omitrow";
opts.MissingRule = "omitrow";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

filename = uigetfile({'*.csv'}, 'Pick the Coordinate file');
tbl = readtable(filename, opts);

% Convert to output type
frame = tbl.coords;


% Clear temporary variables
%clear opts 
fps = 24; %define frames per second

filtx = tbl.x;
filty = tbl.y;
del=find(tbl.likelihood < 0.90 );
filtx(del) = NaN; filty(del) = NaN;
filtx = fillmissing(filtx, 'nearest');
filty = fillmissing(filty, 'nearest');
clear del

filtx1 = tbl.x1;
filty1 = tbl.y1;
del=find(tbl.likelihood1 < 0.99);
filtx1(del) = NaN; filty1(del) = NaN;
filtx1 = fillmissing(filtx1, 'nearest');
filty1 = fillmissing(filty1, 'nearest');
clear del

filtx2 = tbl.x2;
filty2 = tbl.y2;
del=find(tbl.likelihood2 < 0.99);
filtx2(del) = NaN; filty2(del) = NaN;
filtx2 = fillmissing(filtx2, 'nearest');
filty2 = fillmissing(filty2, 'nearest');
clear del


%using a lowpass filter to remove jittering and noise
cutoff_frequency = 3; 
[b, a] = butter(4, cutoff_frequency / (fps / 2), 'low');
X = filtfilt(b, a, filtx);
Y = filtfilt(b, a, filty);

distance = zeros (length(X), 1);

for i = 2:length(X)
    dist =  sqrt((X(i)-X(i-1))^2 + (Y(i)-Y(i-1))^2);
    distance (i, 1)= distance (i-1, 1) + dist;
end

speed = diff(distance);
sspeed = [NaN  smoothdata(speed, 'gaussian', 60)']';
accel = [NaN diff(sspeed)']';
jerk = [NaN diff(accel)']';

figure

colormap turbo

subplot (1,2,1)
hp = patch([X' NaN], [Y' NaN], 0);
set(hp,'cdata', [sspeed' NaN], 'edgecolor','interp','facecolor','none');
axis off
title(filename, 'HorizontalAlignment', 'left')

subplot (1,2,2)
scatter3 (X, Y, frame/30/60,[], sspeed,"filled", 'MarkerFaceAlpha', 0.5)
%set(gca, 'ZDir', 'reverse')
zlabel('Time (min)')


% just another way to color the plot 
%figure
%c = sspeed/max(sspeed); 
%for i = 1:length(X)-1
%    plot3(X(i:i+1), Y(i:i+1), frame(i:i+1), ...
%        'LineWidth', 2, 'Color', [c(i), 0, 1-c(i)]);
%    hold on;
%end
%% open video frame

vidquery = input("Want to open a video file? (y/n) ", 's');
if strcmp (vidquery, 'y') == 1
vidObj = VideoReader(uigetfile({'*.MP4'}, 'Pick the associetes video file'));
frame = read(vidObj, 1000);
figure
subplot(2,1,1)
image(frame);hold on
set(gca, 'YDir', 'normal')
end

%% bin data

bin = 1*fps*60; %bin in minutes
j=1;
binned_distance(j)=distance(bin);
j=2;
for i = 2*bin:bin:length(distance)
    binned_distance(j)=distance(i)-distance(i-bin);
    j= j+1;
    end
binned_distance=binned_distance';
figure
plot (binned_distance, 'o-', 'LineWidth', 2)
ylabel('Distance')
xlabel('Time bin')
%% time in zones

%define Zone!!
Xleft = 523; Xright = 570; %(start of both sides of the open arms)
%you can skip the Y coordinate if the + is aligned
% & = AND operator, | = OR operator
% filtx = nose
% filtx1 = center
% filtx2 = tail base

% in this example, nose, center and tail need to be in the open arms
%inZone = find((filtx < Xleft & filtx1 < Xleft & filtx2 < Xleft) |  (filtx > Xright & filtx1 > Xright & filtx2 > Xright));
inZone = find(X < Xleft  |  X > Xright);

distance_inzone = zeros (length(inZone),1);

for i = 2:length(inZone)
   zdist = sqrt((X(inZone(i)) - X(inZone(i-1)))^2 + (Y(inZone(i)) - Y(inZone(i-1)))^2);
    distance_inzone (i)= distance_inzone (i-1) + zdist;
end


if strcmp (vidquery, 'y') == 0
figure
subplot(2,1,1)
hold on;
end

plot (X, Y); plot (X(inZone),Y(inZone),'o');
%plot ([Xleft Xright Xright Xleft Xleft], [Ymin Ymin Ymax Ymax Ymin], '-k');
plot ([Xright Xright], [0 1080], '-k')
axis off
title(filename)
timeInZone = length(inZone)/fps;
totalTime = length(X)/fps;
hold on

totalDistance = distance(end);
disp(['Animal spent ' num2str(timeInZone) 's (' num2str(timeInZone/totalTime*100) ' %) in the defined Zone, of a total of ' num2str(totalTime) 's.'])
disp(['Animal traveled ' num2str(totalDistance) 'px in TOTAL (remember to convert de metric units)'])
disp(['Animal traveled ' num2str(distance_inzone(end)) 'px in the selected ZONE'])
%% plot a Heatmap from the coordinates


gridResolution = 1;  
xGrid = min(X)-10:gridResolution:max(X)+10;
yGrid = min(Y)-10:gridResolution:max(Y)+10;

density = zeros(length(yGrid), length(xGrid));

for i = 1:length(xGrid)
    for j = 1:length(yGrid)
        radius = 8; % aprox volumen of the mouse
        inCircle = (X - xGrid(i)).^2 + (Y - yGrid(j)).^2 <= radius^2;
        density(j, i) = sum(inCircle)/30;
    end
end
density(density == 0) = NaN;
% Plot heatmap
subplot(2,1,2)
imagesc(xGrid, yGrid, density)
colorbar
colormap parula
alpha('color') % Enable transparency
axis off
c = colorbar;
set(gca, 'YDir', 'normal')
ylabel(c, 'Seconds')
clim([0 20])
%% save

name = strsplit (filename, '.');
name = string(name(1));
save (name, 'name', 'distance', 'speed',... 
    'accel', 'sspeed', 'timeInZone', 'totalTime', 'binned_distance', 'bin', ...
    'Xright', 'Xleft')