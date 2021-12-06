%{
HoneyBeeHiveSiteDynamicsN4.m
Script that computes bifurcation for hive site selection dynamics with two
alternatives. A plot of the bifurcation diagram is generated in a figure. 
Steady state solutions (of the average opinions) are computed for each
bifurcation parameter value of interest using fsolve.m

Dependencies: None

Authors: B.A. Moghimi, J. Valderrama, H. Kumawat
Date: 2021 Dec 04
%}

close all; 
clear all; clc
rng(239, 'twister')

% Singularity value for bifurcation parameter u (Beta = 0 case)
singu = 1;

%% Network Description & System dynamics
% Grid Search limits
BA = 1;
BB = -1;

% Select figure to reproduce from "Multiagent Decision-Making Dynamics Inspired by Honeybees"
% fig2plot = 2; 
% [D, A, B, N] = generate_network_corrected(fig2plot, BA, BB);

% N=4 network (strongly connected) for case of grid search
D = eye(4);
D(1,1) = 2; D(2,2) = 2; D(3,3) = 3; D(4,4) = 1;
A = [0 1 1 0 
     1 0 1 0 
     1 1 0 1 
     0 0 1 0];
% B = [1; -1; 0; 0]; % Pitchfork case
B = [1; 1; -.5; 0]; % Unfolding case 2
N = 4;

% Control Parameter (stop signalling cross-inhibition in honeybees)
u = .25:0.05:2;

% Initialize plot
% if fig2plot == 4
%     model_redux_fig = figure();
% end

tic
y = {}; % Average opinions 
% Loop over 
for uidx = 1:length(u) %1:length(u) 
    
    % Agents' Initial Opinions
%     x0 = BB + (BA - BB)*rand(N,1);
%     x(:,1) = x0;
    funcHive = @(x) hiveSiteConsensus(x,D,A,u(uidx),B);
%     xidx = fsolve(funcHive,x0);
%     yidx = mean(xidx);

    % Grid Search (brute force) N = 4
    if N == 4
%         if u(uidx) > singu
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
            y{uidx} = yidx;
%         end
    end
     
    % Update plot if
%     if u(uidx) == 2 && fig2plot == 4
%         for node = 1:N
%             color = '-.b';
%             if node < 4
%                 color = '-.g';
%             elseif node > 5
%                 color = '-.m';
%             end
%             figure(model_redux_fig); hold on; plot(ds:ds:max_iter*ds, x(node,:), color)
%         end
%         xlabel('Time (s)'); ylabel('Agent and average opinion ($x_{i}$,y)',...
%             'Interpreter','latex');
%         hold off
%     end
  
end

%% Visualization
bifurcation_fig = figure(); hold on;
figure(bifurcation_fig)
for ii=1:length(u)
    z = y{ii};
    for jj=1:length(z)
        plot(u(ii),z(jj),'rx') % Need to include a marker type such as 'x' or 'o' to see points on figure
    end
end
title({'Bifurcation diagram','N=4'});
xlabel('Bifurcation paramter (u)','FontSize',16,'FontWeight','bold');
ylabel('Average Opinion (y)','FontSize',16,'FontWeight','bold');
hold off;

toc

savefig(bifurcation_fig,'N4_Unfold2') % Save figure to .fig file
save('N4_Unfold2','y') % Save data to .mat file

%% ***************************************************************************************
% ******************************* UNUSED CODE *****************************
%     %simulate for a certain number of steps to find steady-state opinion
%     for sim = 1:max_iter-1
%         % Compute dx/ds = -Dx + u * AS(x) + B
%         dxds = -D*x(:,sim) + u(uidx) * A * arrayfun(@tanh,x(:,sim)) + B;
%         % Integration step
%         x(:,sim+1) = x(:,sim) + dxds*dt;
%         y(:,sim+1) = (1/N) * sum(x(:,sim+1));
%     end

%     figure(bifurcation_fig), plot(u(uidx)*ones(steady_state_iterations,1)', ...
%         y(:,max_iter - steady_state_iterations + 1:end),...
%         '-.or','markersize', 2);

% Initialization of Opinions of Honeybees and Average Opinion
% x = zeros(N,max_iter);
% y = zeros(1, max_iter);

% Initial Condition for Average Opinion
% y(:,1) = (1/N) * sum(x(:,1));

% c = @(x)[];
% ceq = @(x) -D*x(:,1) + u*A*sigmoidFunc(x(:,1)) + B;
% nonlinfcn = @(x) deal(c(x),ceq(x));
% problem = createOptimProblem('fmincon','x0',x0,'objective',objfun,'nonlcon',nonlinfcn);

% Global Optimization Solver Method N = 3
% if N == 3
%     if u(uidx) > singu
%         objfun = @(x) norm(hiveSiteConsensus(x,D,A,u(uidx),B),1);
%         problem = createOptimProblem('fmincon','x0',x0,'objective',objfun,'lb',0,'ub',0);
%         ms = MultiStart('XTolerance',1e-4,'MaxTime',500);
%         stpoints = RandomStartPointSet('NumStartPoints',100,'ArtificialBound',1);
%         [xidx,fval,exitflag,output,solutions] = run(ms,problem,stpoints);
%         yidx = mean(xidx);
%     end
% end