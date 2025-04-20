% Load data from the file (replace 'filename.csv' with your file name)



filename = uigetfile({'*.csv'}, 'Pick the coordinate file');

name = input("Enter a name: ", 's');
if isempty(name)
    name = strsplit(filename,".");
    name = string(name(1));
end


opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Coordinate_ 1", "Coordinate_2", "Time"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Coordinates = readtable(filename, opts);
data = table2array(Coordinates);

time = data(:, 3);
x = data(:, 1);
y = data(:, 2);
subplot(2,1,1)
plot (x, y)
axis off

% Set the number of bins for x and y coordinates
num_bins_x = 30;
num_bins_y = 30;

% Create a 2D histogram
subplot (2,1,2)
histogram2d = hist3([x, y], [num_bins_x, num_bins_y]);
axis off


% Create the heatmap using imagesc
imagesc(histogram2d);
%colorbar; % Add color bar for reference
title('Heatmap from Coordinate File');
