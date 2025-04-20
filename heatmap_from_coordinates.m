%% plot a Heatmap from the coordinates


gridResolution = 5;  
xGrid = min(X):gridResolution:max(X);
yGrid = min(Y):gridResolution:max(Y);

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
imagesc(xGrid, yGrid, density)
colorbar
colormap parula
alpha('color') % Enable transparency
axis off
c = colorbar;
ylabel(c, 'Seconds')
title(name)
clim([0 25])