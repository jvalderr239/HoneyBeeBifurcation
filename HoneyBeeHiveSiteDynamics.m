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
% Time delta (s indep. variable)
ds = 0.02;

% Singularity value for bifurcation parameter u
singu = 1;

%% Network Description
BA = 1;
BB = -1;
% Select figure to reproduce from "Multiagent Decision-Making Dynamics Inspired by Honeybees"
fig2plot = 2; 
[D, A, B, N] = generate_network_corrected(fig2plot, BA, BB);

% N=3 network (strongly connected) for case of grid search
D = 2*eye(3);
A = ones(3);
A = A - eye(3);
B = [1; -1; 0];
N = 3; 

% N=4 network (strongly connected) for case of grid search
D = eye(4);
D(1,1) = 2; D(2,2) = 2; D(3,3) = 3; D(4,4) = 1;
A = [0 1 1 0 
     1 0 1 0 
     1 1 0 1 
     0 0 1 0];
B = [1; -1; 0; 0];
N = 4;

% N=6 network (strongly connected) for case of grid search
% D = eye(6);
% D(1,1) = 2; D(2,2) = 2; D(3,3) = 3;
% D(4,4) = 3; D(5,5) = 2; D(6,6) = 2;
% A = [0 1 1 0 0 0
%      1 0 1 0 0 0
%      1 1 0 1 0 0
%      0 0 1 0 1 1
%      0 0 0 1 0 1
%      0 0 0 1 1 0];
% B = [ones(3,1); -ones(3,1)];
% N = 6;
     
%% System dynamics
% Control Parameter (stop signalling cross-inhibition in honeybees)
u = .25:0.05:2;

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
    % Grid Search (brute force) N = 3
    if N == 3
        if u(uidx) > singu
            grid = BB:.2:BA;
            yidx_arr =[];
            for g1 = 1:length(grid)
                for g2 = 1:length(grid)
                    for g3 = 1:length(grid)
                        x0 = [grid(g1) grid(g2) grid(g3)]';
                        xidx = fsolve(funcHive,x0);
                        yidx = mean(xidx);
                        yidx_arr = [yidx_arr yidx];
                    end
                end
            end
            yidx = round(yidx_arr,4); % Round to 4 decimal places
            yidx = unique(yidx); % Select only for unique values
        end
    end
    % Grid Search (brute force) N = 4
    if N == 4
        if u(uidx) > singu
            grid = BB:.2:BA;
            yidx_arr =[];
            for g1 = 1:length(grid)
                for g2 = 1:length(grid)
                    for g3 = 1:length(grid)
                        for g4 = 1:length(grid)
                            x0 = [grid(g1) grid(g2) grid(g3) grid(g4)]';
                            xidx = fsolve(funcHive,x0);
                            yidx = mean(xidx);
                            yidx_arr = [yidx_arr yidx];
                        end
                    end
                end
            end
            yidx = round(yidx_arr,4); % Round to 4 decimal places
            yidx = unique(yidx); % Select only for unique values
        end
    end
    % Grid Search (brute force) N = 6
    if N == 6
        if u(uidx) > singu
            grid = BB:.2:BA;
            yidx_arr =[];
            for g1 = 1:length(grid)
                for g2 = 1:length(grid)
                    for g3 = 1:length(grid)
                          for g4 = 1:length(grid)
                                for g5 = 1:length(grid)
                                    for g6 = 1:length(grid)
                                        x0 = [grid(g1) grid(g2) grid(g3) grid(g4) grid(g5) grid(g6)]';
                                        xidx = fsolve(funcHive,x0);
                                        yidx = mean(xidx);
                                        yidx_arr = [yidx_arr yidx];
                                    end
                                end
                          end
                    end
                end
            end
            yidx = round(yidx_arr,4); % Round to 4 decimal places
            yidx = unique(yidx); % Select only for unique values
        end
    end
     
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
            figure(model_redux_fig); hold on; plot(ds:ds:max_iter*ds, x(node,:), color)
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
