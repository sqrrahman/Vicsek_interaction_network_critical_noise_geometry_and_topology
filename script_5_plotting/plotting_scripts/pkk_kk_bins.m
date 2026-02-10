clearvars; close all; clc;

set_num = 4;
% bin_ranges = {[50, 68], [69, 94], [347, 479], [1272, 1759], [3371, 4664], [12364, 17109]};
bin_ranges = {[36, 49], [50, 68], [69, 94]};

% [36, 49], [50, 68], [69, 94], [95, 130], [182, 250], [251, 346],
% [347, 479], [480, 663], [664, 918], [919, 1271], [1760, 2435],
% [2436, 3370], [3371, 4664], [4665, 6455], [6456, 8933], [8934, 12363],
% [12364, 17109], [17110, 23677], [23678, 32768]

cc = parula(12); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:);
cc4 = cc(7,:); cc5 = cc(9,:); cc6 = cc(11,:);
colors = {cc1, cc2, cc3, cc4, cc5, cc6};

figure('Units','pixels','Position',[100 100 500 500]); 
ax = axes('Units','pixels','Position',[80 80 400 400]);  

markers = {'o', 's', '^', 'd', 'v', 'p'};

title_entries = { ...
    '$\rho = 0.1, s = 0.1$', '$\rho = 0.1, s = 0.25$', '$\rho = 0.1, s = 1.0$', ...
    '$\rho = 0.25, s = 0.1$', '$\rho = 0.25, s = 0.25$', '$\rho = 0.25, s = 1.0$', ...
    '$\rho = 1.0, s = 0.1$', '$\rho = 1.0, s = 0.25$', '$\rho = 1.0, s = 1.0$' ...
};

legend_entries = {};
legend_handles = [];

hold on;

for i = 1:length(bin_ranges)
    bin_min = bin_ranges{i}(1);
    bin_max = bin_ranges{i}(2);
    
    filename = sprintf('../post_binned_data/N512/dd_binned/dd_set_%d_bin%d_%d.txt', set_num, bin_min, bin_max);
    
    data = readtable(filename, 'Delimiter', '\t', 'NumHeaderLines', 1);
    k_over_k = data{:, 1};
    mean_prob = data{:, 2};
    
    % h = scatter(k_over_k, mean_prob*10, 64, colors{i}, 'filled', ...
    %     'Marker', markers{i});

    h = plot(k_over_k, mean_prob, markers{i}, ...
        'MarkerSize', 6, 'MarkerEdgeColor', colors{i}, ...
        'MarkerFaceColor', 'none', 'LineStyle', 'none');
    
    legend_entries{end+1} = sprintf('$%d \\leq m \\leq %d$', bin_min, bin_max);
    legend_handles(end+1) = h;
    
    fprintf('Successfully plotted bin %d-%d for set %d\n', bin_min, bin_max, set_num);

end



xticks([0 1 2 3 4 5 6]);
xticklabels({'$0$', '$1$', '$2$', '$3$', '$4$', '$5$', '$6$'});
ax = gca;
ax.TickLabelInterpreter = 'latex';

yticks([0 0.2 0.4 0.6 0.8 1.0]);
yticklabels({'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'});
ax = gca;
ax.TickLabelInterpreter = 'latex';

grid on; box on;

xlabel('$k/\langle k \rangle_p$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$P(k/\langle k \rangle_p)$', 'FontSize', 16, 'Interpreter', 'latex');

set(gca, 'FontSize', 16);

legend(legend_handles, legend_entries, 'Location', 'ne', 'FontSize', 16, 'Interpreter', 'latex');



xlim([0 6]);
ylim([0 1]);
grid on; box on;

filename = sprintf('pkk_kk_bins.png');
print(gcf, filename, '-dpng', '-r600');

