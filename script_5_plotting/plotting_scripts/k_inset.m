clearvars; close all; clc;

cc = parula(6); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:);
colors  = {cc1, cc2, cc3};
markers = {'*', '^', 'o'};          % for main pkk plot
styles  = {'-','-',':'};            % for avg <k> plot
groups  = {1:3, 4:6, 7:9};          % for avg <k> plot

% --- Figure & main axes
fig     = figure('Units','pixels','Position',[100 100 400 400]);
ax_main = axes('Parent',fig,'Units','pixels','Position',[80 80 300 300]);

% --- Plot pkk
plot_main_pkk(ax_main, colors, markers);

% --- Inset axes
% ax_inset = axes('Parent',fig,'Units','pixels','Position',[170 110 155 155]);
ax_inset = axes('Parent',fig,'Units','pixels','Position',[180 180 180 180]);

% --- Plot avg <k>
plot_inset_avgk(ax_inset, colors, styles, groups);


print(fig, 'k_inset.png', '-dpng', '-r600');



function plot_main_pkk(ax, colors, markers)

    bin_min = 50;
    bin_max = 68;

    % bin_min = 919;
    % bin_max = 1271;

    % bin_min = 1272;
    % bin_max = 1759;
    
    % bin_min = 1760;
    % bin_max = 2435;

    % bin_min = 2436;
    % bin_max = 3370;

    setIdx  = 4:6;

    % legend_entries = { ...
    %     '$\rho = 0.1, s = 0.1$', '$\rho = 0.1, s = 0.25$', '$\rho = 0.1, s = 1.0$', ...
    %     '$\rho = 0.25, s = 0.1$', '$\rho = 0.25, s = 0.25$', '$\rho = 0.25, s = 1.0$', ...
    %     '$\rho = 1.0, s = 0.1$', '$\rho = 1.0, s = 0.25$', '$\rho = 1.0, s = 1.0$' ...

    axes(ax); hold(ax,'on'); axis(ax,'square');

    for i = 1:length(setIdx)
        set_num = setIdx(i);

        % Determine color group and marker
        color_group  = mod(set_num - 1, 3) + 1;     % 1,4,7 -> 1; 2,5,8 -> 2; 3,6,9 -> 3
        marker_group = ceil(set_num / 3);           % 1-3 -> 1, 4-6 -> 2, 7-9 -> 3
        current_color  = colors{color_group};
        current_marker = markers{marker_group};

        filename = sprintf('../post_binned_data/N512/dd_binned/dd_set_%d_bin%d_%d.txt', set_num, bin_min, bin_max);
        data = readtable(filename, 'Delimiter', '\t', 'NumHeaderLines', 1);

        k_over_k = data{:, 1};
        mean_prob = data{:, 2};
        std_prob  = data{:, 3};

        [k_over_k, sortIdx] = sort(k_over_k);
        mean_prob = mean_prob(sortIdx);
        std_prob  = std_prob(sortIdx);

        idx = mean_prob > 0;
        k_over_k = k_over_k(idx);
        mean_prob = mean_prob(idx);
        std_prob  = std_prob(idx);


        % Scatter
        % plot(ax, k_over_k, mean_prob, current_marker, ...
        %     'MarkerSize', 6, 'MarkerEdgeColor', current_color, ...
        %     'MarkerFaceColor', current_color, 'LineStyle', '-');

        % fill(ax, [k_over_k; flipud(k_over_k)], ...
        %      [(mean_prob+std_prob); flipud((mean_prob-std_prob))], ...
        %      current_color, 'FaceAlpha',0.2,'EdgeColor','none', 'HandleVisibility','off');
        plot(ax, k_over_k, mean_prob, 'Color', current_color, ...
             'LineWidth',1.5,'Marker','none');
    end

    set(ax, 'XScale', 'lin', 'YScale', 'lin', 'FontSize', 16);

    % Ticks
    xticks(ax, 0:5:10);
    xticklabels(ax, {'$0$', '$5$','$10$'});
    yticks([0 0.5 1]);
    yticklabels({'$0$', '$0.5$', '$1.0$'});

    ax.TickLabelInterpreter = 'latex';

    xlim(ax, [0 10]);
    ylim(ax, [0 1]);
    grid(ax,'on'); box(ax,'on');

    ax.TickLabelInterpreter = 'latex';
    grid on; grid minor;

    xlabel(ax, '$k/\langle k \rangle _p$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel(ax, '$P(k/\langle k \rangle _p)$', 'FontSize', 16, 'Interpreter', 'latex');

    set(ax, 'FontSize', 16);

    % legend(ax, legend_entries(setIdx), 'Location', 'ne', 'FontSize', 16, 'Interpreter', 'latex');
end



function plot_inset_avgk(ax, colors, styles, groups)
    axes(ax); hold(ax,'on'); axis(ax,'square');

    i = 2;

    for j = groups{i}
        data = readmatrix(sprintf('../post_binned_data/N512/mkl_binned/mkl_set_%d.txt', j), 'NumHeaderLines', 1);
        x   = data(:, 1); % Cluster size
        n   = data(:, 2); % Number of clusters in bins
        y   = data(:, 3); % Average Degree
        std  = data(:, 4); % Population standard deviation

        valid = x > 0 & y >= 1 & std >= 0 & n > 10;
        x = x(valid); y = y(valid); std = std(valid);

        xx = [x; flipud(x)];
        yy = [y - std; flipud(y + std)];

        plot(ax, x, y, 'Color', colors{mod(j-1,3)+1}, 'LineStyle', styles{i}, 'LineWidth', 1.1);
        fill(xx, yy, colors{mod(j-1,3)+1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');
    end


    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 15);

    % Ticks
    xticks(ax, [1 10000]);
    xticklabels(ax, {'$10^0$', '$10^4$'});
    yticks(ax, [1 100 10000]);
    yticklabels(ax, {'$10^0$', '$10^2$', '$10^4$'});

    ax.TickLabelInterpreter = 'latex';

    xlim(ax, [1 10000]);
    ylim(ax, [1 100]);
    grid(ax,'on'); box(ax,'on');

    xlabel(ax, '$m$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(ax, '$\langle k \rangle$', 'Interpreter', 'latex', 'FontSize', 14);
    ax.XLabel.Units = 'normalized';
    ax.YLabel.Units = 'normalized';
    
    ax.XLabel.Position(2) = -0.04;
    ax.YLabel.Position(1) = -0.04;

    % legend({'$s = 0.1$' '$s = 0.25$' '$s = 1.0$'}, 'Interpreter','latex','FontSize',15,'Location','nw');
end
