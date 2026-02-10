clearvars; close all; clc;

% bin_min = 919;
% bin_max = 1271;

bin_min = 36;
bin_max = 49;

setIdx = 1:9;

cc = parula(6); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:); colors = {cc1, cc2, cc3};

markers = {'*', '^', 'o'};  % star, triangle, circle

figure('Units','pixels','Position',[100 100 500 500]); 
ax = axes('Units','pixels','Position',[80 80 400 400]);

legend_entries = { ...
    '$\rho = 0.1, s = 0.1$', '$\rho = 0.1, s = 0.25$', '$\rho = 0.1, s = 1.0$', ...
    '$\rho = 0.25, s = 0.1$', '$\rho = 0.25, s = 0.25$', '$\rho = 0.25, s = 1.0$', ...
    '$\rho = 1.0, s = 0.1$', '$\rho = 1.0, s = 0.25$', '$\rho = 1.0, s = 1.0$' ...
};

hold on; axis square;

for i = 1:length(setIdx)
    set_num = setIdx(i);
    
    % color group and marker
    color_group = ceil(set_num / 3);  % 1-3 -> 1, 4-6 -> 2, 7-9 -> 3
    marker_group = mod(set_num - 1, 3) + 1;  % 1,4,7 -> 1; 2,5,8 -> 2; 3,6,9 -> 3
    current_color = colors{color_group};
    current_marker = markers{marker_group};
    
    filename = sprintf('../post_binned_data/N512/dd_binned/dd_set_%d_bin%d_%d.txt', set_num, bin_min, bin_max);
    data = readtable(filename, 'Delimiter', '\t', 'NumHeaderLines', 1);
    k_over_k = data{:, 1};
    mean_prob = data{:, 2};
    std_prob = data{:, 3};
    
    % Plot scatter with error bars
    % errorbar(k_over_k, mean_prob, std_prob, 'LineStyle', 'none', ...
    %            'Marker', current_marker, 'MarkerSize', 8, ...
    %            'Color', current_color, 'MarkerFaceColor', current_color, ...
    %            'MarkerEdgeColor', current_color, 'LineWidth', 1.5);

    % plot(k_over_k, mean_prob, current_marker, ...
    %     'MarkerSize', 6, 'MarkerEdgeColor', current_color, ...
    %     'MarkerFaceColor', 'none', 'LineStyle', 'none');
    plot(k_over_k, mean_prob, current_marker, ...
        'MarkerSize', 6, 'MarkerEdgeColor', current_color, ...
        'MarkerFaceColor', current_color, 'LineStyle', 'none');
    fprintf('Successfully plotted: %s\n', filename);
end

% set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16);
xticks([0 1 2 3 4 5 6]);
xticklabels({'$0$', '$1$', '$2$', '$3$', '$4$', '$5$', '$6$'});
ax = gca;
ax.TickLabelInterpreter = 'latex';

yticks([0 0.2 0.4 0.6 0.8 1.0]);
yticklabels({'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'});
ax = gca;
ax.TickLabelInterpreter = 'latex';

xlim([0 6]);
ylim([0 1]);
grid on; box on;

xlabel('$k/\langle k \rangle_p$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$P(k/\langle k \rangle_p)$', 'FontSize', 16, 'Interpreter', 'latex');

set(gca, 'FontSize', 16);

% legend(legend_entries(setIdx), 'Location', 'ne', 'FontSize', 15, 'Interpreter', 'latex');


headers = {'', '$\rho = 0.1$', '$\rho = 0.25$', '$\rho = 1.0$'};
s_values = {'$s = 0.1$', '$s = 0.25$', '$s = 1.0$'};
create_custom_legend(colors, markers, headers, s_values)

filename = sprintf('pkk_kk_sets.png');
print(gcf, filename, '-dpng', '-r600');

function create_custom_legend(colors, markers, headers, s_values)
    % Create a table-style legend for linear plots
    % Define table position and size
    table_x = 1.5; % left edge
    table_y = 0.58; % bottom edge
    cell_width = 1.1; % width of each cell
    cell_height = 0.1; % height of each cell (linear spacing)
    bg_width = 4 * cell_width;
    bg_height = 4 * cell_height; % 4 rows
    
    rectangle('Position', [table_x, table_y-0.002, bg_width, bg_height], ...
        'FaceColor', 'white', 'EdgeColor', 'k');
    
    % Draw internal grid lines
    % Vertical lines (3 internal lines)
    for i = 1:3
        x_line = table_x + i * cell_width;
        line([x_line, x_line], [table_y-0.002, table_y-0.002+bg_height], ...
             'Color', 'k', 'LineWidth', 0.5);
    end
    
    % Horizontal lines (3 internal lines)  
    for i = 1:3
        y_line = table_y + i * cell_height;
        line([table_x, table_x+bg_width], [y_line, y_line], ...
             'Color', 'k', 'LineWidth', 0.5);
    end
    
    % Draw table grid and fill content
    for row = 1:4
        for col = 1:4
            x_pos = table_x + (col-1) * cell_width;
            y_pos = table_y + (4-row) * cell_height; % linear spacing, top to bottom
            
            % Add content
            if row == 1 % Header row
                if col > 1
                    text(x_pos + cell_width/2, y_pos + cell_height/2, headers{col}, ...
                        'FontSize', 14, 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end
            else % Data rows
                if col == 1 % First column: s values
                    text(x_pos + cell_width/2, y_pos + cell_height/2, s_values{row-1}, ...
                        'FontSize', 14, 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                else % Data cells: markers
                    set_num = (col-2)*3 + (row-1); % Calculate set number (1-9)
                    if set_num <= 9
                        color_group = ceil(set_num / 3);
                        marker_group = mod(set_num - 1, 3) + 1;
                        current_color = colors{color_group};
                        current_marker = markers{marker_group};
                        
                        % Plot marker in cell center
                        plot(x_pos + cell_width/2, y_pos + cell_height/2, current_marker, ...
                            'MarkerSize', 8, 'MarkerEdgeColor', current_color, ...
                            'MarkerFaceColor', current_color);
                    end
                end
            end
        end
    end
end


