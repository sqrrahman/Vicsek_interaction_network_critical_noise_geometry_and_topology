clear; clc; close all;
N = 512;

% Get all the text files in the directory
files = dir('sim_data/sim_data_*.txt');

% Initialize arrays to store noise values and time steps
noiseValues = zeros(1, numel(files));
timeSteps = zeros(1, numel(files));

% Loop through files to extract noise and time step numbers
for i = 1:numel(files)
    noiseValues(i) = sscanf(files(i).name, 'sim_data_noise_%f_sim_01_t_%*d.txt');
    timeSteps(i) = sscanf(files(i).name, 'sim_data_noise_%*f_sim_01_t_%d.txt');
end

% Sort files based on the extracted time step
[~, idx] = sort(timeSteps);
simDataFiles = files(idx);

% Number of frames (files)
numFrames = numel(simDataFiles);

% Set up parameters
rho = {RHO_VALUE};  % Change as needed
L = sqrt(N/ rho);
video_name = sprintf('AniMotion_set_{setIdx}');
video_writer = VideoWriter(video_name, 'Motion JPEG AVI');
video_writer.FrameRate = 60;  % Adjust frame rate as needed
open(video_writer);
figure_handle = figure('Position', [100, 100, 1200, 1200]);  % Set resolution
set(figure_handle, 'Toolbar', 'none');  % Disable toolbar

% Loop over the sorted files
for i = 1:numFrames
    fname = fullfile('sim_data', simDataFiles(i).name);
    data = readmatrix(fname);  % [x, y, θ]
    x = data(:, 1);
    y = data(:, 2);
    theta = data(:, 3);

    % Plot and write the frame to the video
    plt_motion(x, y, theta, i, L, figure_handle);
    frame = getframe(gcf);  % Capture current figure as a frame
    writeVideo(video_writer, frame);  % Write frame to video
end

% Close the video writer
close(video_writer);
fprintf('Saved: %s\n', video_name);  % Display saved video name

% Function for plotting
function plt_motion(xx, yy, theta, t, L, figure_handle)
    figure(figure_handle);  % Use provided figure handle
    set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');

    u = cos(theta); v = sin(theta);

    quiver(xx, yy, u, v, 0.1, 'b', 'HandleVisibility', 'off'); hold on;
    xlim([0 L]);
    ylim([0 L]);
    xlabel('x', 'Interpreter', 'latex');
    ylabel('y', 'Interpreter', 'latex');
    title(sprintf('t = %d', t), 'FontSize', 14, 'Interpreter', 'latex');
    hold off;
end
