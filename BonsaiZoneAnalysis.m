% Load data (assuming data is in a table format)
data = readtable(uigetfile({'*.csv'})); 

% process time info
data.Timestamp = datetime(data.Timestamp, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSSSXXX", ...
     'TimeZone', 'UTC');
data.Timestamp.Format = "yyyy-MM-dd HH:mm:ss.SSSSSSS";
data.elapsedtime = seconds(data.Timestamp - data.Timestamp(1));
fps = length(data.elapsedtime)/data.elapsedtime(end);
disp(['Frames per second: ', num2str(fps)]);
disp(['Total time: ', num2str(data.elapsedtime(end)), ' seconds, (', num2str(data.elapsedtime(end)/60), ' min)']);
%plot(data.elapsedtime, 'o');hold on;
% Extract the zone column
zones = data.Value; 

% Identify transitions (when the zone changes from the previous frame)
zone_changes = [true; ~strcmp(zones(1:end-1), zones(2:end))];

% Count the number of times the animal enters each zone
unique_zones = unique(zones);
entry_counts = arrayfun(@(z) sum(strcmp(zones(zone_changes), z)), unique_zones);

% Compute the duration (frames per stay)
stay_lengths = diff([find(zone_changes); length(zones) + 1]);
zone_stays = zones(zone_changes); % Zones corresponding to detected stays

% Compute average duration per zone
average_durations = arrayfun(@(z) mean(stay_lengths(strcmp(zone_stays, z))), unique_zones);

% Cunt total zone occurences
counts = cellfun(@(x) sum(strcmp(zones, x)), unique_zones);
total_frames = length(zones);

% Compute percentage time spent in each zone
percentage_time = (counts / total_frames) * 100;

% Display results
results = table(unique_zones,percentage_time, entry_counts, average_durations, ...
    'VariableNames', {'Zone', '%time', 'EntryCount', 'AvgDuration'});
disp(results);


%% plot Trajectory


plot (data.Value_Item1_X, data.Value_Item1_Y)
ylim ([0 480])
xlim ([0 640])
axis off
