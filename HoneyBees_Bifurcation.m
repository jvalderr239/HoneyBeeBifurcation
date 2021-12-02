%{

%}
% close all;
 clear all; clc
rng(239, 'twister')

% Number of Time Steps
max_iter = 250;
% Time delta
dt = 0.02;
% Agent (honeybee) Inertia 
ui = 1.5;
% Agent (honeybee) Social Effort or Attention paid to opinions of other
% agents
us = 3;

%% Network Description
BA = -5;
BB = -4.5;
fig2plot = 2; % figure from honeybee paper
[D, A, B, N] = generate_network_corrected(fig2plot, BA, BB);

%% System dynamics

% Control Parameter (stop signalling cross-inhibition in honeybees)
u = .25:0.005:2;

% number of points to plot for bifurcation diagram
steady_state_iterations = floor(max_iter/10);

% Initialization
x = zeros(N,max_iter);
y = zeros(1, max_iter);
% Initial Condition for Average Opinion
y(:,1) = (1/N) * sum(x(:,1));

% Initialize plot
if fig2plot == 4
    model_redux_fig = figure();
end
bifurcation_fig = figure(); hold on;
for uidx = 1:size(u,2) 
    % Initialize zero velocities
    dxds = zeros(N,1); 
    
    % Agents' Initial Opinions
    x0 = BB + (BA - BB)*rand(N,1);
    x(:,1) = x0;
    %simulate for a certain number of steps to find steady-state opinion
    for sim = 1:max_iter-1
        % Compute dx/ds = -Dx + u * AS(x) + B
        dxds = -D*x(:,sim) + u(:,uidx) * A * arrayfun(@tanh,x(:,sim)) + B;
        % Integration step
        x(:,sim+1) = x(:,sim) + dxds*dt;
        y(:,sim+1) = (1/N) * sum(x(:,sim+1));
    end
    % Update plot if
    if u(uidx) == 2 && fig2plot == 4
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
    figure(bifurcation_fig), plot(u(:,uidx)*ones(steady_state_iterations,1)', ...
        y(:,max_iter - steady_state_iterations + 1:end),...
        '-.or','markersize', 2);
    
end
figure(bifurcation_fig)
title('Bifurcation diagram');
xlabel('Bifurcation paramter (u)'); ylabel('Average Opinion (y)');
hold off;
