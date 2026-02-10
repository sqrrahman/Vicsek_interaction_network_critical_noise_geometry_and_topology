clearvars; close all; clc;

cc = parula(6); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:); colors = {cc1, cc2, cc3};

groups = {1:3, 4:6, 7:9};

figure('Units','pixels','Position',[200 200 400 400]);
ax = axes('Units','pixels','Position',[80 80 300 300]);

for i = 2
    hold on;
    axis square;
    for j = groups{i}
        data = readmatrix(sprintf('../post_binned_data/N512/mkl_binned/mkl_set_%d.txt', j), 'NumHeaderLines', 1);
        x = data(:, 1); % Cluster size
        n = data(:, 2); % Number of clusters in bins
        y = data(:, 3); % Average Degree
        std = data(:, 4); % Population standard deviation
    
        width = std;
        
        % Only consider valid data (positive values)
        valid = x > 0 & y >= 1 & std >= 0 & n > 10;
        x = x(valid); y = y(valid);
        width = width(valid);
    
        % Shade for std or 95% confidence interval
        xx = [x; flipud(x)];
        yy = [y-width; flipud(y+width)];
    
        plot(x, y, 'Color', colors{mod(j-1,3)+1}, 'LineWidth', 1.5);
        fill(xx, yy, colors{mod(j-1,3)+1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

        % if j == 4
        %     draw_slope(x,y,1,length(y),0.95,0.4-j/20, 1.6);
        % end

    end
    plot_styling();
end

filename = sprintf('k.png');
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
    legend({'$s = 0.1$' '$s = 0.25$' '$s = 1.0$'}, 'Interpreter','latex','FontSize',16,'Location','nw');
    legend off;
end

function draw_slope(x,y,a,b,txt_pos_x, txt_pos_y, shift)
    % restrict data to the given range
    x = x(a:b);
    y = y(a:b);

    % --- fit model: log(y) ~ log(log(x))
    mdl = fitlm(log(log(x)), log(y));

    % predictions in log-space
    [yp, yci] = predict(mdl, log(log(x)));
    yfit = exp(yp);

    % plot fitted curve
    plot(x/shift, yfit*shift, 'k-', 'LineWidth', 1.4, 'HandleVisibility','off');

    % % confidence band
    % fill([x; flipud(x)], ...
    %      [exp(yci(:,1)); flipud(exp(yci(:,2)))], ...
    %      'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
    %      'HandleVisibility','off');
    % 
    % % slope = exponent b, intercept = log(a)
    % ci = coefCI(mdl);
    % slope = mdl.Coefficients.Estimate(2);   % b
    % intercept = mdl.Coefficients.Estimate(1); % log(a)
    % a_est = exp(intercept);
    % r2 = mdl.Rsquared.Ordinary;
    % 
    % % equation text
    % txt = sprintf(['$y \\approx %.3f (\\log x)^{%.3f}$ \\\\ ' ...
    %            '\n 95\\%% CI for $b$: [%.3f, %.3f],  $R^2 = %.3f$'], ...
    %            a_est, slope, ci(2,1), ci(2,2), r2);
    % 
    % text(txt_pos_x, txt_pos_y, txt, ...
    %      'Units','normalized', ...
    %      'Interpreter','latex', ...
    %      'FontSize',14, ...
    %      'HorizontalAlignment','right');
end
