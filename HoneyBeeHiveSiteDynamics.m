%{
HoneyBeeHiveSiteDynamics.m
Script that computes bifurcation for hive site selection dynamics with two
alternatives. A plot of the bifurcation diagram is generated in a figure. 
Steady state solutions (of the average opinions) are computed for each
bifurcation parameter value of interest using fsolve.m

Dependencies: generate_network.m

Authors: J. Valderrama, H. Kumawat, B.A. Moghimi
Date: 2021 Dec 01
%}

close all; 
clear all; clc
rng(239, 'twister')

% Number of Time Steps
max_iter = 250;
% Time delta
dt = 0.02;

%% Network Description
BA = 1;
BB = -1;
% Select figure to reproduce from "Multiagent Decision-Making Dynamics Inspired by Honeybees"
fig2plot = 2; 
[D, A, B, N] = generate_network_corrected(fig2plot, BA, BB);

%% System dynamics
% Control Parameter (stop signalling cross-inhibition in honeybees)
u = .25:0.005:2;

% Initialization of Opinions of Honeybees and Average Opinion
% x = zeros(N,max_iter);
y = zeros(1, max_iter);

% Initial Condition for Average Opinion
% y(:,1) = (1/N) * sum(x(:,1));

% Initialize plot
if fig2plot == 4
    model_redux_fig = figure();
end
bifurcation_fig = figure(); hold on;
% Loop over 
for uidx = 1:length(u) 
    
    % Agents' Initial Opinions
    x0 = BB + (BA - BB)*rand(N,1);
%     x(:,1) = x0;
    funcHive = @(x) hiveSiteConsensus(x,D,A,u(uidx),B);
    xidx = fsolve(funcHive,x0);
    yidx = mean(xidx);
    
%     %simulate for a certain number of steps to find steady-state opinion
%     for sim = 1:max_iter-1
%         % Compute dx/ds = -Dx + u * AS(x) + B
%         dxds = -D*x(:,sim) + u(uidx) * A * arrayfun(@tanh,x(:,sim)) + B;
%         % Integration step
%         x(:,sim+1) = x(:,sim) + dxds*dt;
%         y(:,sim+1) = (1/N) * sum(x(:,sim+1));
%     end
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
     figure(bifurcation_fig), plot(u(uidx),yidx,'-.or','markersize', 2);
%     figure(bifurcation_fig), plot(u(uidx)*ones(steady_state_iterations,1)', ...
%         y(:,max_iter - steady_state_iterations + 1:end),...
%         '-.or','markersize', 2);
    
end
figure(bifurcation_fig)
title('Bifurcation diagram');
xlabel('Bifurcation paramter (u)'); ylabel('Average Opinion (y)');
hold off;
