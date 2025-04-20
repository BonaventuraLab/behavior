clear 
close all

opts = delimitedTextImportOptions("NumVariables", 16);

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["coords", "x", "y", "likelihood", "x1", "y1", "likelihood1", "x2", "y2", "likelihood2", "x3", "y3", "likelihood3", "x4", "y4", "likelihood4"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
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
fps = 30; %define frames per second

filtx = tbl.x;
filty = tbl.y;
del=find(tbl.likelihood < 0.99);
filtx(del) = NaN; filty(del) = NaN;
filtx = fillmissing(filtx, 'nearest');
filty = fillmissing(filty, 'nearest');
clear del

filtx1 = tbl.x1;
filty1 = tbl.y1;
del=find(tbl.likelihood1 < 0.95);
filtx1(del) = NaN; filty1(del) = NaN;
filtx1 = fillmissing(filtx1, 'nearest');
filty1 = fillmissing(filty1, 'nearest');
clear del

filtx2 = tbl.x2;
filty2 = tbl.y2;
del=find(tbl.likelihood2 < 0.95);
filtx2(del) = NaN; filty2(del) = NaN;
filtx2 = fillmissing(filtx2, 'nearest');
filty2 = fillmissing(filty2, 'nearest');
clear del

filtx3 = tbl.x3;
filty3 = tbl.y3;
del=find(tbl.likelihood3 < 0.95 | tbl.y3<0);
filtx3(del) = NaN; filty3(del) = NaN;
filtx3 = fillmissing(filtx3, 'nearest');
filty3 = fillmissing(filty3, 'nearest');
clear del

filtx4 = tbl.x4;
filty4 = tbl.y4;
del=find(tbl.likelihood4 < 0.95);
filtx4(del) = NaN; filty4(del) = NaN;
filtx4 = fillmissing(filtx4, 'nearest');
filty4 = fillmissing(filty4, 'nearest');
clear del



%using a lowpass filter to remove jittering and noise
%MIGHT NEED TO BE ADJUSTED
cutoff_frequency = 3; 
[b, a] = butter(4, cutoff_frequency / (fps / 2), 'low');
X = filtfilt(b, a, filtx1);
Y = filtfilt(b, a, filty1);

%% trim file
% this will overwrite X and Y, if you want to take another section
% of the video, you need to re-run the previous section before.
% skip/silence this section if you don't want to trim it

START = input("Enter the START time (in seconds): ");
DURATION = input("Enter the DURATION (in minutes): " );

X = X(START*fps:(START+DURATION*60)*fps);
Y = Y(START*fps:(START+DURATION*60)*fps);
frame = frame(START*fps:(START+DURATION*60)*fps);

%% calculate distance and plots

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
vidObj = VideoReader(uigetfile({'*.mp4'}, 'Pick the associetes video file'));
frame = read(vidObj, 100);
figure
subplot(2,1,1)
image(frame);hold on
end


%% time in zones

%define Zone!!
Xmin = 925; Xmax = Xmin+160;
Ymin = 280; Ymax = Ymin+90; 

inZone = find(X >= Xmin & X <= Xmax & Y >= Ymin & Y <= Ymax);

if strcmp (vidquery, 'y') == 0
figure
subplot(2,1,1)
hold on;
end

plot (X, Y); plot (X(inZone),Y(inZone));
plot ([Xmin Xmax Xmax Xmin Xmin], [Ymin Ymin Ymax Ymax Ymin], '-k');
axis off
title(filename)
timeInZone = length(inZone)/fps;
totalTime = length(X)/fps;
hold on
%subplot (2,1,2)
%labels = {'Time in Zone', 'Other'};
%pie ([timeInZone, totalTime-timeInZone],[1 0], labels)
totalDistance = distance(end);
disp(['Animal spent ' num2str(timeInZone) 's (' num2str(timeInZone/totalTime*100) ' %) in the defined Zone, of a total of ' num2str(totalTime) 's.'])
disp(['Animal traveled ' num2str(totalDistance) 'px (remember to convert de metric units'])
%% animation

%figure
%for i = 1:length(X)
%    plot (x(i)-X(i),y(i)-Y(i), 'o')
%    hold on
%    plot (x1(i)-X(i),y1(i)-Y(i), 'o')
%    plot (x2(i)-X(i),y2(i)-Y(i), 'o')
%    plot (x3(i)-X(i),y3(i)-Y(i), 'o')
%    plot (x4(i)-X(i),y4(i)-Y(i), 'o')

%    plot (x(1:i)-X(1:i),y(1:i)-Y(1:i), 'b-')
    
%    ylim([-400 400]);
%    xlim([-400 400]);
%    pause (1/fps)
%    cla
%end

%calculate angle between tail end and centroid
%
nose = cat(2, filtx, filty);
r_ear = cat(2, filtx1, filty1);
l_ear = cat(2, filtx2, filty2);
centroid = cat(2, filtx3, filty3);
tailend = cat (2, filtx4, filty4);

centered_tailend = tailend - centroid;
centered_nose = nose - centroid;
centered_r_ear = r_ear - centroid;
centered_l_ear = l_ear - centroid;
centered_centroid = centroid -centroid;

theta = atan2(centered_tailend(:, 2), centered_tailend(:, 1));
N = size(centered_nose, 1);
R = zeros(N, 2, 2);
for i = 1:N
    angle = theta(i);
    R(i, :, :) = [cos(angle), -sin(angle); sin(angle), cos(angle)];
end

% rotate_tailend
for i = 1:size(centered_tailend, 1)
    coord_row = centered_tailend(i, :);
    rotation_matrix = squeeze(R(i, :, :));
    rotated_coord_row = coord_row * rotation_matrix;
    rotated_tailend(i, :) = rotated_coord_row;
end

% rotate_nose
for i = 1:size(centered_nose, 1)
    coord_row = centered_nose(i, :);
    rotation_matrix = squeeze(R(i, :, :));
    rotated_coord_row = coord_row * rotation_matrix;
    rotated_nose(i, :) = rotated_coord_row;
end

% rotate_r_ear
for i = 1:size(centered_r_ear, 1)
    coord_row = centered_r_ear(i, :);
    rotation_matrix = squeeze(R(i, :, :));
    rotated_coord_row = coord_row * rotation_matrix;
    rotated_r_ear(i, :) = rotated_coord_row;
end

% rotate_l_ear
for i = 1:size(centered_l_ear, 1)
    coord_row = centered_l_ear(i, :);
    rotation_matrix = squeeze(R(i, :, :));
    rotated_coord_row = coord_row * rotation_matrix;
    rotated_l_ear(i, :) = rotated_coord_row;
end

%for i = 1:length(rotated_nose)
%scatter(rotated_nose(i,1), rotated_nose(i,2));hold on
%scatter(rotated_tailend(i,1), rotated_tailend(i,2));
%ylim([-400 400]);
%xlim([-400 400]);
%pause(1/(3*fps));
%cla;
%end


nose_angle = sin(atan2(rotated_nose(:, 2), rotated_nose(:, 1)));
nose_to_center = sqrt(rotated_nose(:,1).^2 + rotated_nose(:,2).^2);
tail_to_center = sqrt(rotated_tailend(:,1).^2 + rotated_tailend(:,2).^2);
nose_to_tail = nose_to_center + tail_to_center;
centroid_to_centerarena = (sqrt((centroid(:,1)-(max(X)-min(X))).^2 + (centroid(:,2)-(max(Y)-min(Y))).^2))/(max(X)-min(X));


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
subplot (2,1,2)
plot (binned_distance, 'o-', 'LineWidth', 2)
ylabel('Distance')
xlabel('Time bin')

%% export for PCA analysis
%currently in another script: DLC_PCA.m

A = [centroid_to_centerarena, sspeed, accel, nose_angle, nose_to_tail,...
    nose_to_center];
%for i = 2:length(frame)-1;
%    A(i,7)= sspeed(i-1);
%    A(i,8)= sspeed(i+1);
%end

%% save

name = strsplit (filename, '.');
name = string(name(1));
%save (name, 'name', 'A', 'distance', 'speed',... 
%    'accel', 'sspeed', 'timeInZone', 'totalTime', 'binned_distance', 'bin', ...
%    'Xmax', 'Xmin', 'Ymax', 'Ymin')