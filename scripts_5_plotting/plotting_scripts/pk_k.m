clearvars; close all; clc;

setIdx = 1:9;

cc = parula(6); cc1 = cc(1,:); cc2 = cc(3,:); cc3 = cc(5,:); colors = {cc1, cc2, cc3};

markers = {'*', '^', 'o'}; % star, triangle, circle

figure('Units','pixels','Position',[100 100 500 500]); 
ax = axes('Units','pixels','Position',[80 80 400 400]); 

legend_entries = { ...
    '$\rho = 0.1, s = 0.1$', '$\rho = 0.1, s = 0.25$', '$\rho = 0.1, s = 1.0$', ...
    '$\rho = 0.25, s = 0.1$', '$\rho = 0.25, s = 0.25$', '$\rho = 0.25, s = 1.0$', ...
    '$\rho = 1.0, s = 0.1$', '$\rho = 1.0, s = 0.25$', '$\rho = 1.0, s = 1.0$' ...
};


hold on; axis square;

% Loop through each set
for i = 1:length(setIdx)
    set_num = setIdx(i);
    
    % Determine color group and marker
    color_group = ceil(set_num / 3); % 1-3 -> 1, 4-6 -> 2, 7-9 -> 3
    marker_group = mod(set_num - 1, 3) + 1; % 1,4,7 -> 1; 2,5,8 -> 2; 3,6,9 -> 3
    current_color = colors{color_group};
    current_marker = markers{marker_group};
    
    % Construct filename
    filename = sprintf('../post_binned_data/N512/pk_calc/set_%d_pk_calc.txt', set_num);
    
    % Read the file (two columns: mean_prob, std_prob)
    data = readmatrix(filename);
    probabilities = data(:, 1);  % First column: mean probabilities
    std_devs = data(:, 2);       % Second column: standard deviations
    
    % Define bin edges
    bin_size = 1; % Change this to adjust bin size
    max_degree = length(probabilities) - 1;
    bin_edges = 0:bin_size:(max_degree + bin_size);
    
    % Bin the probabilities and standard deviations
    binned_prob = [];
    binned_std = [];
    bin_centers = [];
    
    for b = 1:length(bin_edges)-1
        bin_start = bin_edges(b);
        bin_end = bin_edges(b+1) - 1;
        
        % Sum probabilities and combine standard deviations in this bin
        if bin_end < length(probabilities)
            bin_prob = sum(probabilities(bin_start+1:bin_end+1));
            % For combining standard deviations: sqrt(sum of variances)
            bin_std = sqrt(sum(std_devs(bin_start+1:bin_end+1).^2));
        else
            bin_prob = sum(probabilities(bin_start+1:end));
            bin_std = sqrt(sum(std_devs(bin_start+1:end).^2));
        end
        
        if bin_prob > 0 % Only keep non-zero bins
            binned_prob(end+1) = bin_prob;
            binned_std(end+1) = bin_std;
            bin_centers(end+1) = (bin_start + min(bin_end, max_degree)) / 2;
        end
    end
    
    k_plot = bin_centers';
    prob_plot = binned_prob';
    std_prob = binned_std';
    
    % Plot
    h = plot(k_plot, prob_plot, current_marker, ...
        'MarkerSize', 5, 'MarkerEdgeColor', current_color, ...
        'MarkerFaceColor', current_color, 'LineStyle', 'none');

    % h = errorbar(k_plot, prob_plot, std_prob, current_marker, ...
    %     'MarkerSize', 4, 'MarkerEdgeColor', current_color, ...
    %     'MarkerFaceColor', current_color, 'LineStyle', 'none');
    fprintf('Successfully plotted: %s (max degree: %d)\n', filename, length(probabilities)-1);
end

xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 16); 
ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 16); 

set(gca, 'YScale', 'log');

yticks([1e-8 1e-6 1e-4 1e-2 1]); 
yticklabels({'$10^{-8}$', '$10^{-6}$', '$10^{-4}$', '$10^{-2}$', '$10^{0}$'});
ax = gca; 
ax.TickLabelInterpreter = 'latex'; 

set(gca, 'FontSize', 16);
grid on; box on;

% legend(legend_entries(setIdx),'Location', 'ne', 'FontSize', 16, 'Interpreter', 'latex');

% headers = {'', '$\rho = 0.1$', '$\rho = 0.25$', '$\rho = 1.0$'};
% s_values = {'$s = 0.1$', '$s = 0.25$', '$s = 1.0$'};
% create_custom_legend(colors, markers, headers, s_values);

set(gcf, 'Color', 'white');
filename = sprintf('pk_k.png');
print(gcf, filename, '-dpng', '-r600');



function create_custom_legend(colors, markers, headers, s_values)
    % Create a table-style legend
    
    % Define table position and size
    table_x = 45;      % left edge
    table_y = 0.015;     % bottom edge
    cell_width = 50;    % width of each cell
    cell_height = 0.01; % height of each cell (in log scale)


    bg_width = 4 * cell_width;
    bg_height = 80 * cell_height;
    rectangle('Position', [table_x, table_y-0.002, bg_width, bg_height], ...
              'FaceColor', 'white', 'EdgeColor', 'k');
    
    
    
    % Draw table grid and fill content
    for row = 1:4
        for col = 1:4
            x_pos = table_x + (col-1) * cell_width;
            y_pos = table_y * (10^((row-1) * 0.5)); % logarithmic spacing
            
            
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

