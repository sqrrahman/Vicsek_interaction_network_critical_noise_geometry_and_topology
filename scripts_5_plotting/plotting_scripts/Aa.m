clearvars; close all; clc;

cc = parula(6); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:); colors = {cc1, cc2, cc3};

groups = {1:3, 4:6, 7:9};

figure('Units','pixels','Position',[200 200 400 400]);
ax = axes('Units','pixels','Position',[80 80 300 300]);

for i = 2
    hold on;
    axis square;
    for j = groups{i}
        data = readmatrix(sprintf('../post_binned_data/N512/perimeter_binned/perimeter_set_%d.txt', j), 'NumHeaderLines', 1);
        x = data(:, 1); % Cluster size
        n = data(:, 2); % Number of clusters in bins
        y = data(:, 17); % Perimeter area
        std = data(:, 18); % Population standard deviation
    
        width = std;
        
        % Only consider valid data (positive values)
        valid = x > 0 & y > 0 & std > 0 & n > 10;
        x = x(valid); y = y(valid);
        width = width(valid);
    
        % Shade for std or 95% confidence interval
        xx = [x; flipud(x)];
        yy = [y-width; flipud(y+width)];
    
        plot(x, y, 'Color', colors{mod(j-1,3)+1}, 'LineWidth', 1.5);
        fill(xx, yy, colors{mod(j-1,3)+1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

        % if j == 4
        %     draw_slope(x,y,4,length(y),1.6);
        % end
    end
    
    plot_styling()
end

filename = sprintf('Aa.png');
print(gcf, filename, '-dpng', '-r600');


function plot_styling()

    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16);
    
    % Ticks
    xticks([0.1 1 10 100 1000 10000]);
    xticklabels({'$10^-1$', '$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$'});
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    yticks([0.1 1 10 100 1000 10000 100000]);
    yticklabels({'$10^{-1}$', '$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'});
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    xlim([1 10000]);
    ylim([0.1 10000]);
    grid on; box on;
    
    % label, legend, text, title
    xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$A$', 'Interpreter', 'latex', 'FontSize', 16);
    legend({'$s = 0.1$' '$s = 0.25$' '$s = 1.0$'}, 'Interpreter','latex','FontSize',16,'Location','nw');
    legend off;
end


function draw_slope(x,y,a,b,shift)
    x = x(a:b);
    y = y(a:b);
    mdl = fitlm(log(x), log(y));
    [yp, yci] = predict(mdl, log(x));
    yfit = exp(yp)*shift;
    slope_linear = (log(yfit(end)) - log(yfit(1))) / (log(x(end)) - log(x(1)));
    % plot(x/shift, yfit*shift, 'k-', 'LineWidth', 1.2, 'HandleVisibility','off');
    
    % Draw right triangle to show slope
    % Choose a position for the triangle (around 60% along the line)
    start_idx = 5;
    a_x = x(start_idx) / shift; % Adjust for the scaling factor
    a_y = yfit(start_idx) * shift;


    end_idx = 10;
    b_x = x(end_idx) / shift; % Adjust for the scaling factor
    b_y = yfit(end_idx) * shift;

    c_x = a_x;
    c_y = b_y;
    
    
    % Triangle vertices
    x_tri = [a_x, c_x, b_x, a_x];
    y_tri = [a_y, c_y, b_y, a_y];
    
    % Draw the triangle
    cc = parula(6); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:); colors = {cc1, cc2, cc3};
    plot(x_tri, y_tri, '-',  'Color', cc1, 'LineWidth', 1.2, 'HandleVisibility','off');
    

    at_x = a_x - 14;
    at_y = a_y + 15;

    bt_x = b_x - 90;
    bt_y = c_y + 40;
    % Add annotations
    % Horizontal side annotation (1 in log scale = 10x)
    text(at_x, at_y, sprintf('$%.2f$', slope_linear), ...
         'Interpreter', 'latex', 'FontSize', 16, 'HorizontalAlignment', 'center');
    
    % Vertical side annotation (slope value)
    text(bt_x, bt_y, ...
         '$1$', 'Interpreter', 'latex', 'FontSize', 16, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

