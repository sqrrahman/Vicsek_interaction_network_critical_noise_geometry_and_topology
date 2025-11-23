% MATLAB script to plot clusters from C++ processed data for SET_PLACEHOLDER

clear; clc;

set_idx = SET_PLACEHOLDER;
cluster_data_dir = ['./cluster_data_set_' num2str(set_idx)];
output_dir = ['./set_' num2str(set_idx)];

% Create output directory
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Define colors
cc = parula(6); 
cc1 = cc(1,:); 
cc2 = cc(3,:); 

fprintf('Plotting clusters for set_%d...\n', set_idx);

% Get all cluster data files
cluster_files = dir([cluster_data_dir '/cluster_set_*.txt']);

for i = 1:length(cluster_files)
    filename = cluster_files(i).name;
    filepath = fullfile(cluster_data_dir, filename);
    
    try
        % Read cluster data
        fid = fopen(filepath, 'r');
        
        % Read header: cluster_size avg_degree centroid_x centroid_y box_size
        header = fscanf(fid, '%f', 5);
        cluster_size = header(1);
        avg_degree = header(2);
        centroid_x = header(3);
        centroid_y = header(4);
        box_size = header(5);
        
        % Read agent coordinates
        agent_coords = fscanf(fid, '%f %f', [2, cluster_size]);
        agent_x = agent_coords(1, :)';
        agent_y = agent_coords(2, :)';
        
        % Read number of edges
        num_edges = fscanf(fid, '%d', 1);
        
        % Read edge coordinates
        if num_edges > 0
            edge_coords = fscanf(fid, '%f %f %f %f', [4, num_edges]);
            edge_x1 = edge_coords(1, :);
            edge_y1 = edge_coords(2, :);
            edge_x2 = edge_coords(3, :);
            edge_y2 = edge_coords(4, :);
        else
            edge_x1 = [];
            edge_y1 = [];
            edge_x2 = [];
            edge_y2 = [];
        end
        
        fclose(fid);
        
        % Create figure
        fig = figure('Visible', 'off', 'Position', [100, 100, 800, 800]);
        hold on;
        
        % Plot nodes only (no edges)
        scatter(agent_x, agent_y, 20, cc1, 'filled');
        
        % Set limits and formatting
        box_size = 100; % Fixed box size
        half_box = box_size / 2;
        xlim([centroid_x - half_box, centroid_x + half_box]);
        ylim([centroid_y - half_box, centroid_y + half_box]);
        axis equal;
        grid on;
        
        % Add axis labels
        xlabel('X Coordinate', 'FontSize', 12);
        ylabel('Y Coordinate', 'FontSize', 12);
        
        % Add text annotation
        text_str = sprintf('N = %d, avg k = %.2f', cluster_size, avg_degree);
        text(centroid_x - half_box + 0.05*box_size, centroid_y + half_box - 0.1*box_size, ...
            text_str, 'FontSize', 12, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
        
        % Save figure
        output_filename = strrep(filename, 'cluster_', '');
        output_filename = strrep(output_filename, '.txt', '.png');
        output_file = fullfile(output_dir, output_filename);
        
        print(fig, output_file, '-dpng', '-r150');
        close(fig);
        
        fprintf('Plotted: %s (N=%d, k=%.2f)\n', output_filename, cluster_size, avg_degree);
        
    catch ME
        fprintf('Error plotting %s: %s\n', filename, ME.message);
    end
end

fprintf('Plotting complete for set_%d!\n', set_idx);