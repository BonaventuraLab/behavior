%% Set directory to folder containing CSV files

close all; clear;

folder = uigetdir('','Select a folder');
name = strsplit (folder, '\');
name = (string(name(length(name))));
% Get a list of all CSV files in the folder
file_list = dir(fullfile(folder, '*.csv'));

%% plot time course
% Loop through each file and plot the data
figure;
hold on;
for i = 1:length(file_list)
    % Load data from current file
    warning_state = warning('off', 'all');
    data = readtable(fullfile(folder, file_list(i).name));
    warning(warning_state);
    data.Properties.VariableNames{1} = 'Time';
    data.InSessionTime  = data.Time(:,1)-data.Time(1:1);
  
    % Plot data
   plot (data.InSessionTime, data.Pellet_Count)
end


% Add legend to the plot using the file names

legend(strrep({file_list.name},'.csv',''));

%% raster plots
figure();

for i = 1:length(file_list)
    warning_state = warning('off', 'all');
    data = readtable(fullfile(folder, file_list(i).name));
    warning(warning_state);
    data.Properties.VariableNames{1} = 'Time';
    data.InSessionTime  = data.Time(:,1)-data.Time(1:1);
time = data.InSessionTime;
event = data.Event;

hold on;
  for j = 1:numel(event)
    if strcmp(event{j}, 'Pellet')
        plot([time(j), time(j)], [i-1, i], 'k-', 'LineWidth', 2);
    else
       if strcmp(event{j}, 'Left')
        plot([time(j), time(j)], [i-1, i], 'b-', 'LineWidth', 0.5);
       else
          if strcmp(event{j}, 'Right')
        plot([time(j), time(j)], [i-1, i], 'r-', 'LineWidth', 0.5);
         end
       end
    end
  end
end

xlim([min(time), max(time)]);
%ylim([0, 1]);
xlabel('Time');
ylabel('Trials');
clearvars time event

%% do histograms

summary = table(strrep({file_list.name},'.csv','')');
summary.Properties.VariableNames{1} = 'File';
bin_edges = 0:10:60;

for i = 1:length(file_list)
    warning_state = warning('off', 'all');
    data = readtable(fullfile(folder, file_list(i).name));
    warning(warning_state);
    data.Properties.VariableNames{1} = 'Time';
    data.InSessionTime  = minutes(data.Time(:,1)-data.Time(1:1));

    summary.Session_type(i) = data.Session_type(2);
    Pellet = data(strcmp(data.Event, 'Pellet'), :);
    pellet_counts (:,i) = histcounts (Pellet.InSessionTime, bin_edges);

    
    LeftPokes = data(strncmpi(data.Event, 'Left', 4), :); 
        %using strncmpi instead of strcmp compares the first n characters
        %this is to include LeftWithPellet and LeftTimeout
    Left_counts (:,i) = histcounts (LeftPokes.InSessionTime, bin_edges);

    RightPokes = data(strncmpi(data.Event, 'Right', 4), :);
    Right_counts (:,i) = histcounts (RightPokes.InSessionTime, bin_edges);

   retrieval{i} = data.Retrieval_Time(~isnan(data.Retrieval_Time));

     
end 
clearvars  LeftPokes RightPokes

percent_correct = Left_counts./(Left_counts+Right_counts)*100;
summary.Pellets = sum(pellet_counts, 1)';
summary.Left_total = sum(Left_counts, 1)';
summary.Right_total = sum(Right_counts, 1)';


figure()
subplot (1,2,1)
plot (summary.Pellets, 'o-', 'LineWidth', 2)
ylim([0, max(summary.Pellets)+1])
ylabel('Pellets')
xlabel('Trials')

subplot(1,2,2)
plot (summary.Left_total, 'o-', 'LineWidth', 2);hold on;
plot (summary.Right_total, 'o-', 'LineWidth', 2)
ylabel('Pokes')
xlabel('Trials')
legend("Left Pokes", "Right Pokes")

figure()
subplot (1, 2, 1)
pellet_mean = mean (pellet_counts');
pellet_sem = std(pellet_counts') / sqrt(length(pellet_counts'));
bar (bin_edges(2:end),pellet_mean);hold on;
errorbar(bin_edges(2:end),pellet_mean, pellet_sem, '.', 'color', 'black');
title('Pellets')

subplot (1, 2, 2)
Left_mean = mean (Left_counts');
Left_sem = std(Left_counts') / sqrt(length(Left_counts'));



Right_mean = mean (Right_counts');
Right_sem = std(Right_counts') / sqrt(length(Right_counts'));

bar (bin_edges(2:end),Left_mean);hold on;
plot (bin_edges(2:end),Right_mean,'color', 'r', 'LineWidth', 2);hold on;

errorbar(bin_edges(2:end),Left_mean, Left_sem, '.', 'color', 'k');
errorbar(bin_edges(2:end),Right_mean, Right_sem, '.', 'color', 'r');
title('Nose Pokes')
legend ('Left (Active)', 'Right (Inactive)')

%% calc/average cummulative counts
cumulative_pellets = pellet_counts;
for ii = 2:length(cumulative_pellets)
    cumulative_pellets(ii,:) = cumulative_pellets(ii,:) +cumulative_pellets(ii-1,:);
end
mean_cumul = mean(cumulative_pellets, 2);
SEM_cumul = (std(cumulative_pellets')/(size(cumulative_pellets,2)-1)')';
STD_cumul = std(cumulative_pellets')';

cumulative_leftpokes = Left_counts;
for ii = 2:length(cumulative_leftpokes)
    cumulative_leftpokes(ii,:) = cumulative_leftpokes(ii,:) +cumulative_leftpokes(ii-1,:);
end

cumulative_rightpokes = Right_counts;
for ii = 2:length(cumulative_rightpokes)
    cumulative_rightpokes(ii,:) = cumulative_rightpokes(ii,:) +cumulative_rightpokes(ii-1,:);
end
%figure;
plot (bin_edges(2:end), mean_cumul, 'k', 'LineWidth', 2);hold on
plot (bin_edges(2:end), mean_cumul+SEM_cumul, 'k--')
plot (bin_edges(2:end), mean_cumul-SEM_cumul, 'k--')
fill ([1:length(SEM_cumul) fliplr(1:length(SEM_cumul))], ...
[(mean_cumul+SEM_cumul)' flip(mean_cumul-SEM_cumul)'], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');


ylabel('Rewards obtained')
xlabel('Time')
title('Average Cumulative Rewards')
%% SAVE
clearvars i j data ans
save (name)
writetable (summary,strcat(name, "_summary_",".xls"))