%{
Authors: Jose, Hemant, Barin
Description: Describe the system dynamics of honeybee nest 
decision-making
%}

close all; clear all; clc
rng(239, 'twister')
% time steps to take for simulation
max_iter = 250;
% number of agents
N = 8;  
dt = 0.02;

% inertia 
ui = .5;
% social effort
us = .5;

%x0 = [-10 -5 0 -1 2]';
%%
% Network Description
% From Figure 4 of the paper

% degree matrix for strongly connected network
D = 4 * eye(N);
D(4,4) = 7;
D(5,5) = 7;

% adjacency matrix for fully connected network
A = zeros(N);
A(1,2:5) = 1;                           % node 1
A(2,1) = 1;A(2,3:5) = 1;                % node 2
A(3,1:2) = 1; A(3,4:5) = 1;             % node 3
A(4,1:3) = 1; A(4,5) = 1; A(4,6:8) = 1; % node 4
A(5,1:3) = 1; A(5,4) = 1; A(5,6:8) = 1; % node 5
A(6:8,4:5) = 1; A(6,7:8) = 1;           % node 6
A(7,6) = 1; A(7,8) = 1;                 % node 7
A(8,6:7) = 1;                           % node 8
%% System dynamics

BA = 1; BB = -1;

% in order to visualize the behavior of the system
% plot trajectories of x as a function of the control parameter u
u_vals = .25:0.01:2;
steady_state_iterations = 0.1*max_iter;
% Initialize vectors
x = zeros(N,max_iter);
xc = zeros(2,max_iter);
B = [1 1 1 0 0 -1 -1 -1]';
avg_opinion = zeros(1, steady_state_iterations);
%initial vals
avg_opinion(:,1) = (1/N) * sum(x(:,1));

% Initialize plot
figure; hold on
for u = 1:size(u_vals,2) 
    % Initialize zero velocities
    ux = zeros(N,1); 
    
    % Agents' Initial Positions
    x0 = BB + (BA - BB)*rand(N,1);
    x(:,1) = x0;
    %simulate for a certain number of steps to find steady-state opinion
    for sim = 1:max_iter-1
        % Compute dx/ds = -Dx + u * AS(x) + B
        ux = -D*x(:,sim) + u_vals(:,u) * A * arrayfun(@tanh,x(:,sim)) +B;
        % Integration step
        x(:,sim+1) = x(:,sim) + ux.*ui*dt;
        avg_opinion(:,sim+1) = (1/N) * sum(x(:,sim+1));
    end
    % plot average opinions per control parameter
    %avg_opinion = zeros(1, steady_state_iterations);
    %avg_opinion(:,1) = (1/N)*sum(x(:,1));
    %for iter = 2:steady_state_iterations
    %    avg_opinion(:,iter) = (1/N) * sum(x(:,end));
        %set(agenPlot,'xdata',u_vals(:,u),'ydata',avg_opinion(:,u))
        %set(xBplot,'xdata',xc(1,k),'ydata',xc(2,k))
    %end
    % Update plot
    %if u == 2
    %    for node = 1:N
    %        color = '-.b';
    %        if node < 4
    %            color = '-.g';
    %        elseif node > 5
    %            color = '-.m';
    %        end
    %        plot(1:max_iter, x(node,:), color)
    %    end
    %end
    %plot(1:max_iter, avg_opinion)
    plot(u_vals(:,u)*ones(steady_state_iterations,1)', ...
        avg_opinion(:,max_iter - steady_state_iterations),...
        '-.or','markersize',2);
    
end
%axis([0,3,-100,100])
title('Bifurcation diagram');
xlabel('Bifurcation paramter (u)');ylabel('Average Opinion (y)');
hold off;
