%{
Authors: Jose, Hemant, Barin
Description: Describe the system dynamics of honeybee nest 
decision-making
%}

close all; clear all; clc
rng(239, 'twister')
% time steps to take for simulation
max_iter = 600;
dt = 0.02;
% inertia 
ui = 2;
% social effort
us = 1.2;
%%
% Network Description
BA = 1;
BB = -1;
figure_to_plot = 2; % figure from honeybee paper
[D, A, B, N] = adaptive_network(figure_to_plot, BA, BB);

%% System dynamics

% in order to visualize the behavior of the system
% plot trajectories of x as a function of the control parameter u
%u_vals = .25:0.05:2;
% number of points to plot for bifurcation diagram
%steady_state_iterations = 0.1*max_iter;
% Initialize vectors
x = zeros(N,max_iter);
avg_opinion = zeros(1, max_iter);
ubar_avg = zeros(1,max_iter)
U_0 = zeros(1,N)
%initial vals
avg_opinion(:,1) = (1/N) * sum(x(:,1));
avg_opinion_threshold = 2
epsilon = 0.07

% Initialize zero velocities
    ux = zeros(N,1); 
    
    % Agents' Initial Positions
    x0 = BB + (BA - BB)*rand(N,1);
    x(:,1) = x0;
% Initialize plot
%if figure_to_plot == 4
%    model_redux_fig = figure();
%end
bifurcation_fig = figure(); hold on;

for sim = 1:max_iter-1
    ubar_dot = (avg_opinion_threshold^2 - avg_opinion(:,sim)^2)*epsilon;
    ubar_avg(1,sim+1) = ubar_avg(1,sim) + ubar_dot*dt;
    U_0 = (ubar_avg(1,sim+1));
    U_matrix  = diag(U_0);
    ux = -D*x(:,sim) + U_matrix * A * arrayfun(@tanh,x(:,sim)) + B
    x(:,sim+1) = x(:,sim) + ux.*ui*dt;
    x(:,sim+1);
    avg_opinion(:,sim+1) = (1/N) * sum(x(:,sim+1))
    a = avg_opinion(:,sim+1);
    figure(bifurcation_fig), plot(ubar_avg(:,sim)', ...
        avg_opinion(1,sim),...
        'or','markersize', 2, 'DisplayName','cos(2x)');
    %legend('adaptive control')
end


    
    %simulate for a certain number of steps to find steady-state opinion
%     for sim = 1:max_iter-1
%         % Compute dx/ds = -Dx + u * AS(x) + B
%         ux = -D*x(:,sim) + u_vals(:,u) * A * arrayfun(@tanh,x(:,sim)) + B;
%         % Integration step
%         x(:,sim+1) = x(:,sim) + ux.*ui*dt;
%         avg_opinion(:,sim+1) = (1/N) * sum(x(:,sim+1));
%     end
    % Update plot if
%     if u == 2 && figure_to_plot == 4
%         for node = 1:N
%             color = '-.b';
%             if node < 4
%                 color = '-.g';
%             elseif node > 5
%                 color = '-.m';
%             end
%             figure(model_redux_fig); hold on; plot(dt:dt:max_iter*dt, x(node,:), color)
%         end
%         xlabel('Time (s)'); ylabel('Agent and average opinion ($x_{i}$,y)',...
%             'Interpreter','latex');
%         hold off
%     end
%     %plot(1:max_iter, avg_opinion)

% in order to visualize the behavior of the system
% plot trajectories of x as a function of the control parameter u
u_vals = 0:0.005:2.5;
% number of points to plot for bifurcation diagram
steady_state_iterations = 0.1*max_iter;
% Initialize vectors
x = zeros(N,max_iter);
avg_opinion = zeros(1, max_iter);
thresh_avg_opinion = avg_opinion_threshold* ones(1, max_iter);
%initial vals
avg_opinion(:,1) = (1/N) * sum(x(:,1));

% Initialize plot
if figure_to_plot == 4
    model_redux_fig = figure();
end
%bifurcation_fig = figure(); hold on;
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
    figure(bifurcation_fig), plot(u_vals(:,u), ...
        avg_opinion(:,sim+1),...
        'b--o','markersize', 1);
    figure(bifurcation_fig), plot(u_vals(:,u), ...
        thresh_avg_opinion(:,sim+1),...
        'ok','markersize', 1);
    %legend('adaptive', 'bifurcation')
    hold on;

h = zeros(3, 1);
h(1) = plot(0,0,'or', 'visible', 'on');
h(2) = plot(0,0,'ob', 'visible', 'on');
h(3) = plot(0,0,'ok', 'visible', 'on');
legend(h,'adaptive','bifurcation', 'y_{avg} threshold');
    
end
    
    

figure(bifurcation_fig)
title('Bifurcation diagram');
xlabel('Bifurcation paramter (u)');ylabel('Average Opinion (y)');
hold off;
