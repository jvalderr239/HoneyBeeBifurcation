%{
Authors: Jose, Hemant, Barin
Description: Describe the system dynamics of honeybee nest 
decision-making
%}

close all; clear all; clc
rng(239, 'twister')
% time steps to take for simulation
max_iter = 250;
dt = 0.02;
% inertia 
ui = 1.5;
% social effort
us = 3;
%%
% Network Description
BA = 5;
BB = -5;
figure_to_plot = 2; % figure from honeybee paper
[D, A, B, N] = generate_network(figure_to_plot, BA, BB);

%% System dynamics

% in order to visualize the behavior of the system
% plot trajectories of x as a function of the control parameter u
u_vals = .25:0.05:2;
% number of points to plot for bifurcation diagram
steady_state_iterations = 0.1*max_iter;
% Initialize vectors
x = zeros(N,max_iter);
avg_opinion = zeros(1, max_iter);
%initial vals
avg_opinion(:,1) = (1/N) * sum(x(:,1));

% Initialize plot
if figure_to_plot == 4
    model_redux_fig = figure();
end
bifurcation_fig = figure(); hold on;
for u = 1:size(u_vals,2) 
    % Initialize zero velocities
    ux = zeros(N,1); 
    
    % Agents' Initial Positions
    x0 = BB + (BA - BB)*rand(N,1);
    x(:,1) = x0;
    %simulate for a certain number of steps to find steady-state opinion
    for sim = 1:max_iter-1
        % Compute dx/ds = -Dx + u * AS(x) + B
        ux = -D*x(:,sim) + u_vals(:,u) * A * arrayfun(@tanh,x(:,sim)) + B;
        % Integration step
        x(:,sim+1) = x(:,sim) + ux.*ui*dt;
        avg_opinion(:,sim+1) = (1/N) * sum(x(:,sim+1));
    end
    % Update plot if
    if u == 2 && figure_to_plot == 4
        for node = 1:N
            color = '-.b';
            if node < 4
                color = '-.g';
            elseif node > 5
                color = '-.m';
            end
            figure(model_redux_fig); hold on; plot(dt:dt:max_iter*dt, x(node,:), color)
        end
        xlabel('Time (s)'); ylabel('Agent and average opinion ($x_{i}$,y)',...
            'Interpreter','latex');
        hold off
    end
    %plot(1:max_iter, avg_opinion)
    figure(bifurcation_fig), plot(u_vals(:,u)*ones(steady_state_iterations,1)', ...
        avg_opinion(:,max_iter - steady_state_iterations + 1:end),...
        '-.or','markersize', 2);
    
end
figure(bifurcation_fig)
title('Bifurcation diagram');
xlabel('Bifurcation paramter (u)');ylabel('Average Opinion (y)');
hold off;
