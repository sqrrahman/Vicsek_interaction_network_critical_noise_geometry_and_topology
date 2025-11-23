clearvars; close all; clc;

% Parameters
alpha = 0.5;

% Simulation parameters
sim_values = 1:30;
t_values = 0:10000:2000000;

output_folder = 'perimeter_data';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Function to check if cluster is valid (min_x > 1 and min_y > 1) also, ignoring 2-agent clusters
function valid = isValidCluster(cluster_agents, agents)
    if length(cluster_agents) < 3
        valid = false;
        return
    end
    
    min_x = min(agents(cluster_agents, 1));
    min_y = min(agents(cluster_agents, 2));
    
    valid = min_x > 1.0 && min_y > 1.0;
end

% Function to find matching p_vector and sim_data files
function [p_filename, sim_filename] = findMatchingFiles(sim, t)
    p_filename = fullfile('p_data', sprintf('sim_%02d', sim), ...
                          sprintf('p_vector_t_%07d.txt', t));
    sim_filename = fullfile('sim_data', sprintf('sim_%02d', sim), ...
                            sprintf('sim_data_t_%07d.txt', t));
end

% Main processing loop
for sim = sim_values
    for t = t_values
        % Find matching input files
        [p_filename, sim_filename] = findMatchingFiles(sim, t);
        
        if ~isfile(p_filename) || ~isfile(sim_filename)
            fprintf('Could not find matching files for sim %d t %d\n', sim, t);
            continue;
        end
        
        % Read parent vector and agent data
        try
            p_vector = readmatrix(p_filename, 'OutputType', 'double') + 1;
            agents = readmatrix(sim_filename, 'OutputType', 'double');
        catch
            fprintf('Error reading files for sim %d t %d\n', sim, t);
            continue;
        end
        
        if isempty(p_vector) || isempty(agents) || length(p_vector) ~= size(agents, 1)
            fprintf('Error: size mismatch for sim %d t %d\n', sim, t);
            continue;
        end
        
        % Group agents by parent
        unique_parents = unique(p_vector);
        clusters = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        
        for parent = unique_parents'
            clusters(parent) = find(p_vector == parent);
        end
        
        results = [];
        
        parents = keys(clusters);
        for i = 1:length(parents)
            parent = parents{i};
            cluster_agents = clusters(parent);
        
            % Check if cluster is valid
            if isValidCluster(cluster_agents, agents)
                x = agents(cluster_agents, 1);
                y = agents(cluster_agents, 2);
        
                try
                    shp = alphaShape(x, y, alpha);
                    L = perimeter(shp);
                    A = area(shp);
        
                    % Get boundary points with unique coordinates
                    [~, boundaryPoints] = boundaryFacets(shp);
                    [uniqueBoundaryCoords, unique_idx] = unique(boundaryPoints, 'rows');
                    numBoundaryPoints = numel(unique_idx);
                    perimeterLogical = ismember([x, y], uniqueBoundaryCoords, 'rows');
                    inShapeLogical = inShape(shp, x, y);
                    leafLogical = ~inShapeLogical;
                    insideLogical = inShapeLogical & ~perimeterLogical;
                    
                    % Count nodes of each type
                    numPerimeterPoints = sum(perimeterLogical);
                    numInsidePoints = sum(insideLogical);
                    numLeafPoints = sum(leafLogical);
                    
                    % Vectorized degree calculations
                    n = length(x);
                    dx_matrix = abs(x - x');
                    dy_matrix = abs(y - y');
                    dist_matrix = sqrt(dx_matrix.^2 + dy_matrix.^2);
                    neighbor_matrix = (dx_matrix <= 1) & (dy_matrix <= 1) & (dist_matrix <= 1) & (dist_matrix > 0);
                    
                    if numPerimeterPoints > 0
                        perimeter_degrees = sum(neighbor_matrix(perimeterLogical, :), 2);
                        avg_degree_perimeter = mean(perimeter_degrees);
                    else
                        avg_degree_perimeter = NaN;
                    end
                    
                    if numInsidePoints > 0
                        inside_degrees = sum(neighbor_matrix(insideLogical, :), 2);
                        avg_degree_inside = mean(inside_degrees);
                    else
                        avg_degree_inside = NaN;
                    end
                    
                    if numLeafPoints > 0
                        leaf_degrees = sum(neighbor_matrix(leafLogical, :), 2);
                        avg_degree_leaf = mean(leaf_degrees);
                    else
                        avg_degree_leaf = NaN;
                    end
        
                    % Store results
                    cluster_size = length(cluster_agents);
                    results = [results; cluster_size, L, numPerimeterPoints, numInsidePoints, numLeafPoints, ...
                        avg_degree_perimeter, avg_degree_inside, avg_degree_leaf, A, double(parent)];
        
                end
            end
        end
        
        % Write results
        outdir = fullfile(output_folder, sprintf('sim_%02d', sim));
        if ~exist(outdir, 'dir')
            mkdir(outdir);
        end
        output_filename = fullfile(outdir, sprintf('perimeter_data_t_%07d.txt', t));
        
        if ~isempty(results)
            writematrix(results, output_filename, 'Delimiter', ' ');
            fprintf('Processed: %s with %d clusters\n', output_filename, size(results, 1));
        else
            fid = fopen(output_filename, 'w'); fclose(fid);
            fprintf('Processed: %s with 0 clusters\n', output_filename);
        end
    end
end

fprintf('All processing complete.\n');
