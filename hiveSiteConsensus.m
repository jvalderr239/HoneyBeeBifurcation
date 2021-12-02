function dxds = hiveSiteConsensus(x,D,A,u,B)
% hiveSiteConsensus.m 
% Computes the rate of change of honeybee preferences for two alternatives
% to hive site selection. 
% Inputs: 
%   x - input vector containing the preferences (opinions) for N honeybees
%   D - degree matrix corresponding to network graph representing the swarm
%   A - adjacency matrix for network graph
%   sigmoidFunc - nonlinear function applied to input x to saturate
%                 opinions that exceed a certain threshold.
%   u - stop-signalling cross-inhibition parameter (bifurcation parameter)
%   B - external preference values of each agent (honeybee) (e.g. external
%       preference for hive site A or B)
% Author: B.A. Moghimi
% Date: 2021 Dec 01
sigmoidFunc = @tanh;
dxds = -D*x + u*A*sigmoidFunc(x) + B;