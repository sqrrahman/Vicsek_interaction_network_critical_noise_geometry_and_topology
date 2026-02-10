clearvars; close all; clc;

cc = parula(6); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:); colors = {cc1, cc2, cc3};

groups = {1:3, 4:6, 7:9};

figure('Units','pixels','Position',[200 200 400 400]);
ax = axes('Units','pixels','Position',[80 80 300 300]);

for i = 2
    hold on;
    axis square;
    for j = groups{i}
        data = readmatrix(sprintf('../post_binned_data/N512/cc_binned/cc_set_%d.txt', j), 'NumHeaderLines', 1);
        x = data(:, 1); % Cluster size
        n = data(:, 2); % Number of clusters in bins
        y = data(:, 3); % clustering coefficient
        std = data(:, 4); % Population standard deviation
    
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
    end
    plot_styling()
end

filename = sprintf('cc.png');
print(gcf, filename, '-dpng', '-r600');


function plot_styling()

    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16);
    % set(gca, 'XScale', 'log', 'FontSize', 16);
    
    % Ticks
    xticks([0.1 1 10 100 1000 10000]);
    xticklabels({'$10^-1$', '$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$'});
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    % yticks([0.4 0.5 0.6 0.7]);
    % yticklabels({'$0.4$', '$0.5$', '$0.6$', '$0.7$'});
    yticks([0.45 0.60 0.75]);
    yticklabels({'$0.45$', '$0.60$', '$0.75$'});
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    xlim([1 10000]);
    ylim([0.45 0.75]);
    grid on; box on;
    
    % label, legend, text, title
    xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$$\langle c \rangle$$', 'Interpreter', 'latex', 'FontSize', 16);
    legend({'$s = 0.1$' '$s = 0.25$' '$s = 1.0$'}, 'Interpreter','latex','FontSize',16,'Location','se');
    % legend off;
end

