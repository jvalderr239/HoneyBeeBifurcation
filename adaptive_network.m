function [Degree, Adjacency, Beta, N] = adaptive_network(figureNum, BA, BB)
%{
This function generates the fully connected network graphs from
Gray, Rebecca & Franci, Alessio & Srivastava, Vaibhav & Leonard, Naomi. 
(2017). Multiagent Decision-Making Dynamics Inspired by Honeybees. 
IEEE Transactions on Control of Network Systems. PP. 10.1109/TCNS.2018.2796301. 

@param figureNum is the network corresponding to the figure in the paper
@param BA is the value of the agent preference for nest A
@param BB is the value of the agent preference for nest B
%}
if figureNum == 2
    N = 12;
    % degree matrix
    D = eye(N);
    D(1:5) = 8 * D(1:5);
    D(6:9) = 6 * D(6:9);
    D(10:12) = 7 * D(10:12);
    d = [8,8,8,8,8,6,6,6,6,7,7,7]
    D = diag(d)
    %adjacency matrix
    A = ones(N);
    for row = 1:N
        for col = 1:N
            if row == col
                A(row,col) = 0;
            end
        end
    end
    A(1:5,6:9) = 0
    A(6:9,10:12) = 0
    A(10:12,1:5) =0
    B = [BA * ones(1,5) BB * ones(1,4) zeros(1,3)]';
elseif figureNum == 3
    N = 12;
    % degree matrix
    D = 7 * eye(N);
    D(5:8,:) = D(5:8,:) * 11/7;
    % adjacency matrix
    A = zeros(N);
    A(1,2:4) = 1; A(1:4,5:8) = 1;  %group 1
    A(2,1) = 1; A(2,3:4) = 1;
    A(3,1:2) = 1; A(3,4) = 1;
    A(4,1:3) = 1;
    A(5:8,1:4) = 1;A(5:8,6:12) = 1; %group 2
    A(5,6:8) = 1;
    A(6,5) = 1; A(6,7:8) = 1;
    A(7,5:6) = 1; A(7,8) = 1;
    A(8,5:7) = 1;
    A(9:12,5:8) = 1; % group 3
    A(9,10:12) = 1;
    A(10,9) = 1; A(10,11:12) = 1;
    A(11,9:10) = 1; A(11,12) = 1;
    A(12,9:11) = 1;
    
    B = [BA * ones(1,4) zeros(1,4) BB * ones(1,4)]';
elseif figureNum == 4
    N = 8; 
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
    % Agent preferences
    B = [BA * ones(1,3) zeros(1,2) BB * ones(1,3)]';
else
    error('Error. \n Only figures 2,3 or 4 can be generated from Honeybee Paper.');
end

Degree = D; Adjacency = A; Beta = B;
end

