clearvars; close all; clc;

cc = parula(6); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:); colors = {cc1, cc2, cc3};

folders = {'N512', 'N4096', 'N32768'};

styles = {'-','-',':'};

groups = {1:3, 4:6, 7:9};

figure('Units','pixels','Position',[100 100 500 500]); 
ax = axes('Units','pixels','Position',[80 80 400 400]); 

    hold on;
    axis square;

        data = readmatrix(sprintf('../post_binned_data/N512/perimeter_binned/perimeter_set_4.txt'), 'NumHeaderLines', 1);
        x = data(:, 1); % Cluster size
        n = data(:, 2); % Number of clusters in bins
        y1 = data(:, 11); % k on the boundary
        std1 = data(:, 12); % Population standard deviation
        y2 = data(:, 13); % k in the interior
        std2 = data(:, 14); % Population standard deviation
        y3 = data(:, 15); % k of leaf points
        std3 = data(:, 16); % Population standard deviation
        
        % Only consider valid data (positive values)
        valid = n > 10 & y1 >1 & y2 >1 & y3 > 1;
        x = x(valid);
        y1 = y1(valid); std1 = std1(valid);
        y2 = y2(valid); std2 = std2(valid);
        y3 = y3(valid); std3 = std3(valid);
    
        % Shade for std
        xx = [x; flipud(x)];
        yy1 = [y1-std1; flipud(y1+std1)];
        yy2 = [y2-std2; flipud(y2+std2)];
        yy3 = [y3-std3; flipud(y3+std3)];
    
        plot(x, y1, 'Color', colors{1}, 'LineStyle', '-', 'LineWidth', 1.5);
        fill(xx, yy1, colors{1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');
        plot(x, y2, 'Color', colors{2}, 'LineStyle', '-', 'LineWidth', 1.5);
        fill(xx, yy2, colors{2}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');
        plot(x, y3, 'Color', colors{3}, 'LineStyle', '-', 'LineWidth', 1.5);
        fill(xx, yy3, colors{3}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');
        
    plot_styling()


filename = sprintf('kBIE.png');
print(gcf, filename, '-dpng', '-r600');


function plot_styling()
    
    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16);
    
    % Ticks
    xticks([1 10 100 1000 10000]);
    xticklabels({'$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$'});
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    yticks([1 10 100 1000 10000]);
    yticklabels({'$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$'});
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    xlim([1 10000]);
    ylim([1 100]);
    grid on; box on;
    
    % label, legend, text, title
    xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\langle k \rangle$', 'Interpreter', 'latex', 'FontSize', 16);
    legend({'$\langle k \rangle_b$' '$\langle k \rangle_i$' '$\langle k \rangle_e$'}, 'Interpreter','latex','FontSize',16,'Location','nw');
end