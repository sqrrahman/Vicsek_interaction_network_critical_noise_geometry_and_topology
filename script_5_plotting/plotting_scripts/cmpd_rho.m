clearvars; close all; clc;

cc = parula(6); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:); colors = {cc1, cc2, cc3};

figure('Units','pixels','Position',[200 200 400 400]);
ax = axes('Units','pixels','Position',[86 86 300 300]);

hold on; axis square;
for j = [1,4,7]
    data = readmatrix(sprintf('../post_binned_data/N512/cmpd_binned/cmpd_set_%d.txt', j), 'NumHeaderLines', 1);
    x = data(:, 1); % Cluster size
    num_clusters = data(:,2); %number of clusters in that bin
    y = data(:, 3); % Probability
    std   = data(:,4);         % std
    
    % Only consider valid data (positive values)
    valid = x > 0 & y > 0;
    x = x(valid); y = y(valid); std = std(valid);

    % Shade for std
    xx = [x; flipud(x)];                         % [x; reversed x]
    yy = [y-std; flipud(y+std)];                 % [lower; reversed upper]

    plot(x, y, 'Color', colors{(j-1)/3 + 1}, 'LineStyle', '-.', 'LineWidth', 1.5);
    fill(xx, yy, colors{(j-1)/3 + 1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

    % if j == 4
    %     draw_slope(x,y,1,14,1.6,colors{(j-1)/3 + 1});
    % end
end
plot_styling_b()

filename = sprintf('cmpd_rho.png');
print(gcf, filename, '-dpng', '-r600');


function plot_styling_b()
    
    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16);
    
    % Ticks
    xticks([1 10 100 1000 10000 100000]);
    xticklabels({'$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'});
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    yticks([1e-12 1e-9 1e-6 1e-3 1]);
    yticklabels({'$10^{-12}$', '$10^{-9}$', '$10^{-6}$', '$10^{-3}$', '$10^{0}$'});
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    xlim([1 100000]);
    ylim([10e-13 1]);
    grid on; box on;
    
    % label, legend, text, title
    xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$P(m)$', 'Interpreter', 'latex', 'FontSize', 16);
    
    
    legend_labels = {'$\rho = 0.1$' '$\rho = 0.25$' '$\rho = 1.0$'};
    legend(legend_labels, 'Interpreter','latex','FontSize',16,'Location','sw');

    % text(0.98, 0.98, '(b)', 'Units','normalized','FontSize',16,'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top');
end


function draw_slope(x,y,a,b,shift, cc)
    x = x(a:b);
    y = y(a:b);
    mdl = fitlm(log(x), log(y));
    [yp, yci] = predict(mdl, log(x));
    yfit = exp(yp)*shift;
    % slope_linear = (log(yfit(end)) - log(yfit(1))) / (log(x(end)) - log(x(1)));
    slope_linear = mdl.Coefficients.Estimate(2);
    % plot(x/shift, yfit*shift, 'k-', 'LineWidth', 1.2, 'HandleVisibility','off');
    
    % Draw right triangle to show slope
    % Choose a position for the triangle (around 60% along the line)
    start_idx = 7;
    a_x = x(start_idx) * shift; % Adjust for the scaling factor
    a_y = yfit(start_idx) * shift;


    end_idx = 12;
    b_x = x(end_idx) * shift; % Adjust for the scaling factor
    b_y = yfit(end_idx) * shift;

    c_x = b_x;
    c_y = a_y;
    
    
    % Triangle vertices
    x_tri = [a_x, c_x, b_x, a_x];
    y_tri = [a_y, c_y, b_y, a_y];
    
    % Draw the triangle
    plot(x_tri, y_tri, '-',  'Color', cc, 'LineWidth', 1.2, 'HandleVisibility','off');
    

    at_x = a_x + 300;
    at_y = a_y - 0.016;

    bt_x = b_x - 50;
    bt_y = c_y + 0.05;
    % Add annotations
    % Horizontal side annotation (1 in log scale = 10x)
    text(at_x, at_y, sprintf('$%.2f$', abs(slope_linear)), ...
         'Interpreter', 'latex', 'FontSize', 16, 'HorizontalAlignment', 'center');
    
    % Vertical side annotation (slope value)
    text(bt_x, bt_y, ...
         '$1$', 'Interpreter', 'latex', 'FontSize', 16, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end


