% First Assignment Linear Programming for the SC42055 Otimization in Systems and
% Control course
% By Bart de Jong [4367146] and Mike Pesselse [4300564] from Group 19
% This script was written in MATLAB R2019b
clc
close all
clearvars

%% Initialisation
StudentNumber1 = [4 3 6 7 1 4 6];               % Student number Bart de Jong
StudentNumber2 = [4 3 0 0 5 6 4];               % Student number Mike Pesselse
E1 = StudentNumber1(5) + StudentNumber2(5);     % Calculate variable E1
E2 = StudentNumber1(6) + StudentNumber2(6);     % Calculate variable E2
E3 = StudentNumber1(7) + StudentNumber2(7);     % Calculate variable E3



%% 2) First Optimization problem in standard form
% Matrices for the optimization-algorithm
c = [-80 -150 2000+15*E3];                      % Calculate C-matrix
A = [4 6 0; 4 7 -160+2*E2; -1 0 0; 0 0 1];      % Calculate A-matrix
b = [12000+30*E1 0 -1000 50]';                  % Calculate B-matrix

lb = zeros(3, 1);                               % Lower bound is zero for every variable
ub = [];                                        % No upper bound for every variable

options = optimoptions('linprog','Algorithm', 'dual-simplex', 'display', 'off');

%% 2) Calculate optimum for first problem
% Optimize variables x1, x2 and x3 and calculate maximum profit
[x_Q2, ~, exitflag2, ~] = linprog(c, A, b, [], [], lb, ub, [], options);    % Find optimal values for x1, x2 and w

profit_Q2 = -c*x_Q2;                                                        % Calculate maximum profit
w_Q2 = x_Q2(3);                                                             % Obtain number of workers from x-vector

if exitflag2 ~= 1                                                           % Check if the exitflag is anything other than 1
    fprintf('ERROR Q2: No viable solution\n\n');
else
    str = ['For Question 2:\nMaximum profit      ', char(8364), num2str(profit_Q2, '%.2f') '\n' 'Number of Phone-Xs  ', num2str(x_Q2(1)), '\nNumber of Phone-Ss  ', num2str(x_Q2(2)), '\nNumber of workers   ', num2str(x_Q2(3)), '\n\n'];
    fprintf(str)
    clear str
end

%% 2) Plot constraint figure
x2 = 0:3000;                                    % Define linearly spaced vector for number of Phone-Ss

c1x1 = (b(1)-A(1,2)*x2)/A(1,1);                	% Constraint memory cells
c2x1 = (w_Q2*-A(2,3)-A(2,2)*x2)/A(2,1);         % Constraint work hours
c3x1 = ones(1,length(x2))*-b(3);                % Constraint minimum supply
p1x1 = (profit_Q2+c(2)*x2+w_Q2*c(3))/-c(1);     % Profit line

figure                                          % Plot constraint lines, profit lines and optimum
plot(x2, c1x1, x2, c2x1, x2, c3x1, x_Q2(2), x_Q2(1), 'r*', x2, p1x1, 'k:')
hold on

patch([0, 0, 428],[1750, 1000, 1000],'red')     % Visualise feasible solution set
alpha(0.2)

legend('Constraint: memory cells', 'Constraint: maximum work hours', 'Constraint: minimum supply Phone-Xs', 'Optimal solution', 'Profit lines', 'Feasible set')

for profit = -50000:50000:5000000                                                   % Draw profit lines and disable legend labels
    p1x1 = (profit+c(2)*x2+w_Q2*c(3))/-c(1);
    p = plot(x2, p1x1, 'k:');
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     % Turn of legend labels
end

ylim([0 max(max([c1x1', c2x1', c3x1']))*1.1])
xlim([0 max(x2)])
xlabel('Number of Phone-Ss (x_2)')
ylabel('Number of Phone-Xs (x_1)')
title(['Q2) Plot showing constraint lines and optimum (number of workers = ', num2str(w_Q2), ')'])



%% 3) Second optimization problem in standard form
minworkers = 50;                                % Define minimum number of workers
maxworkers = 70;                                % Define maximum number of workers

intcon = [1, 2];                                % Variables 1 and 2 are integers
c = [-80 -150];                                 % Re-define C-matrix
lb = zeros(2, 1);                               % Re-define lower bounds
lengthloop = maxworkers-minworkers+1;           % Calculate the length of the for-loop
xloop = zeros(2,lengthloop);                    % Store optimal variables for every number of workers
profitloop = zeros(1, lengthloop);              % Store profit for every number of workers
exitflag3 = zeros(1, lengthloop);               % Initialise vector for exitflags
options = optimoptions('intlinprog', 'display', 'off');

%% 3) Calculate optimum for second problem
% Loop over every number of workers and find corresponding optimal values
% for the variables and the profit
for w = minworkers:maxworkers                   % Find optima for x1 and x2 for various values of w
    A = [4 6; 4-(w-50)/20 7-(w-50)/20; -1 0];   % Define A-matrix
    b = [12000+30*E1 w*(160-2*E2) -1000]';      % Define B-matrix
    [xloop(:,w-minworkers+1), ~, exitflag3(w-minworkers+1), ~] = intlinprog(c, intcon, A, b, [], [], lb, ub, [], options);
    profitloop(w-minworkers+1) = -c*xloop(:,w-minworkers+1) - w*(2000+15*E3);   % Calculate maximum profit
end

[profit_Q3, index] = max(profitloop);          	% Find the maximum profit and corresponding index in the vector
w_Q3 = index+minworkers-1;                     	% Find the optimal number of workers
x_Q3 = xloop(:, index);                         % Find the optimal values of the variables x1 and x2

if any(exitflag3-ones(1, length(exitflag3)))    % Check for values other than 1 in the exitflags
    fprintf('ERROR Q3: No viable solution\n\n');
else
    str = ['For Question 3:\nMaximum profit      ', char(8364), num2str(profit_Q3, '%.2f') '\n' 'Number of Phone-Xs  ', num2str(x_Q3(1)), '\nNumber of Phone-Ss  ', num2str(x_Q3(2)), '\nNumber of workers   ', num2str(w_Q3), '\n\n'];
    fprintf(str)
    clear str
end

figure                                          % Plot profit as a function of number of workers
plot(minworkers:maxworkers, profitloop, 'k', w_Q3, profit_Q3, 'r*')
annotation('textarrow', [0.8 0.9], [0.85 0.88], 'String', ['Maximum profit: ', char(8364), num2str(profit_Q3, '%.2f')], 'Color', 'red')
grid on
xlim([minworkers, maxworkers])
xlabel('Number of workers')
ylabel(['Profit in euros [', char(8364), ']'])
title('Q3) Maximum profit as a function of number of workers')
line([minworkers, maxworkers],[profit_Q3 profit_Q3],'Color','red','LineStyle','--')

clear w xloop profitloop

%% 3) Plot constraint figure
c1x1 = (b(1)-A(1,2)*x2)/A(1,1);                                     % Constraint memory cells
c2x1 = (w_Q3*(160-2*E2)-(7-(w_Q3-50)/20)*x2)/(4-(w_Q3-50)/20);      % Constraint work hours
c3x1 = ones(1,length(x2))*-b(3);                                    % Constraint minimum supply
p1x1 = (profit_Q3+c(2)*x2+w_Q3*(2000+15*E3))/-c(1);                 % Profit line

figure                                                              % Plot constraint lines, profit lines and optimum
plot(x2, c1x1, x2, c2x1, x2, c3x1, x_Q3(2), x_Q3(1), 'r*', x2, p1x1, 'k:')
hold on

patch([0, 0, 1133, 443],[3045, 1000, 1000, 2380],'red')             % Visualise feasible solution set
alpha(0.2)

legend('Constraint: memory cells', 'Constraint: maximum work hours', 'Constraint: minimum supply Phone-Xs', 'Optimal solution', 'Profit lines', 'Feasible set')

for profit = -50000:50000:5000000                                   % Draw profit lines and disable legend labels
    p1x1 = (profit+c(2)*x2+w_Q3*-(2000+15*E3))/-c(1);
    p = plot(x2, p1x1, 'k:');
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

ylim([0 max(max([c1x1', c2x1', c3x1']))*1.1])
xlim([0 max(x2)])
xlabel('Number of Phone-Ss (x_2)')
ylabel('Number of Phone-Xs (x_1)')
title(['Q3) Plot showing constraint lines and optimum (number of workers = ', num2str(w_Q3), ')'])



%% 4) Calculate optimum for problem with no upper bound on workers
minworkers = 50;                                % Define minimum number of workers
maxworkers = 120;                               % Define maximum number of workers

lengthloop = maxworkers-minworkers+1;           % Calculate length of for-loop
xloop = zeros(2,lengthloop);                    % Store optimal variables for every number of workers
profitloop = zeros(1, lengthloop);              % Store profit for every number of workers
exitflag4 = zeros(1, lengthloop);               % Intialise vector for storing exitflags

for w = minworkers:maxworkers                   % Find optima for x1 and x2 for various values of w
    A = [4 6; 4-(w-50)/20 7-(w-50)/20; -1 0];   % Define A-matrix
    b = [12000+30*E1 w*(160-2*E2) -1000]';      % Define B-matrix
    [xloop(:,w-minworkers+1), ~, exitflag4(w-minworkers+1), ~] = intlinprog(c, intcon, A, b, [], [], lb, ub, [], options);
    profitloop(w-minworkers+1) = -c*xloop(:,w-minworkers+1) - w*(2000+15*E3);   % Calculate maximum profit
end

[profit_Q4, index] = max(profitloop);         	% Find the maximum profit and corresponding index in the vector
w_Q4 = index+minworkers-1;                    	% Find the optimal number of workers
x_Q4 = xloop(:, index);                         % Find the optimal values of the variables x1 and x2

if any(exitflag4-ones(1, length(exitflag4)))    % Check for values other than 1 in the exitflags
    fprintf('ERROR Q4: No viable solution\n\n');
else
    str = ['For Question 4:\nMaximum profit      ', char(8364), num2str(profit_Q4, '%.2f') '\n' 'Number of Phone-Xs  ', num2str(x_Q4(1)), '\nNumber of Phone-Ss  ', num2str(x_Q4(2)), '\nNumber of workers   ', num2str(w_Q4), '\n\n'];
    fprintf(str)
    clear str
end

figure                                          % Plot profit as a function of number of workers
plot(minworkers:maxworkers, profitloop, 'k', w_Q4, profit_Q4, 'r*')
annotation('textarrow', [0.45 0.41], [0.55 0.79], 'String', ['Maximum profit: ', char(8364), num2str(profit_Q4, '%.2f')], 'Color', 'red')
grid on
xlim([minworkers, maxworkers])
xlabel('Number of workers')
ylabel(['Profit in euros [', char(8364), ']'])
title('Q4) Maximum profit as a function of number of workers')
line([minworkers, maxworkers],[profit_Q4 profit_Q4],'Color','red','LineStyle','--')

clear w xloop profitloop

%% 4) Plot constraint figure
c1x1 = (b(1)-A(1,2)*x2)/A(1,1);                                     % Constraint memory cells
c2x1 = (w_Q4*(160-2*E2)-(7-(w_Q4-50)/20)*x2)/(4-(w_Q4-50)/20);      % Constraint work hours
c3x1 = ones(1,length(x2))*-b(3);                                    % Constraint minimum supply
p1x1 = (profit_Q4+c(2)*x2+w_Q4*(2000+15*E3))/-c(1);                 % Profit line

figure                                                              % Plot constraint lines, profit lines and optimum
plot(x2, c1x1, x2, c2x1, x2, c3x1, x_Q4(2), x_Q4(1), 'r*', x2, p1x1, 'k:')
hold on

patch([0, 0, 1348, 1308],[3045, 1000, 1000, 1083],'red')            % Visualise feasible solution set
alpha(0.2)

legend('Constraint: memory cells', 'Constraint: maximum work hours', 'Constraint: minimum supply Phone-Xs', 'Optimal solution', 'Profit lines', 'Feasible set')

for profit = -50000:50000:5000000                                   % Draw profit lines and disable legend labels
    p1x1 = (profit+c(2)*x2+w_Q4*-(2000+15*E3))/-c(1);
    p = plot(x2, p1x1, 'k:');
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

ylim([0 max(max([c1x1', c2x1', c3x1']))*1.1])
xlim([0 max(x2)])
xlabel('Number of Phone-Ss (x_2)')
ylabel('Number of Phone-Xs (x_1)')
title(['Q4) Plot showing constraint lines and optimum (number of workers = ', num2str(w_Q4), ')'])



%% 5) Calculate optimum for changed contract
minworkers = 64;                                    % Define minimum number of workers
maxworkers = 120;                                   % Define maximum number of workers

lengthloop = maxworkers-minworkers+1;               % Calculate length of for-loop
xloop = zeros(2,lengthloop);                        % Store optimal variables for every number of workers
profitloop = zeros(1, lengthloop);                  % Store profit for every number of workers
exitflag5 = zeros(1, lengthloop);                   % Initialise vector for storing exitflags

for w = minworkers:maxworkers                       % Find optima for x1 and x2 for various values of w
    A = [4 6; 4-(w-50)/20 7-(w-50)/20; -1 0; 0 -1]; % Define A-matrix
    b = [12000+30*E1 w*(160-2*E2) -800 -1000]';     % Define B-matrix
    [xloop(:,w-minworkers+1), ~, exitflag5(w-minworkers+1), ~] = intlinprog(c, intcon, A, b, [], [], lb, ub, [], options);
    profitloop(w-minworkers+1) = -c*xloop(:,w-minworkers+1) - w*(2000+15*E3);   % Calculate maximum profit
end

[profit_Q5, index] = max(profitloop);               % Find the maximum profit and corresponding index in the vector
w_Q5 = index+minworkers-1;                          % Find the optimal number of workers
x_Q5 = xloop(:, index);                             % Find the optimal values of the variables x1 and x2

if any(exitflag5-ones(1, length(exitflag5)))        % Check for values other than 1 in the exitflags
    fprintf('ERROR Q5: No viable solution\n\n');
else
    str = ['For Question 5:\nMaximum profit      ', char(8364), num2str(profit_Q5, '%.2f') '\n' 'Number of Phone-Xs  ', num2str(x_Q5(1)), '\nNumber of Phone-Ss  ', num2str(x_Q5(2)), '\nNumber of workers   ', num2str(w_Q5), '\n\n'];
    fprintf(str)
    clear str
end

figure                                              % Plot profit as a function of number of workers
plot(minworkers:maxworkers, profitloop, 'k', w_Q5, profit_Q5, 'r*')
annotation('textarrow', [0.35 0.30], [0.6 0.86], 'String', ['Maximum profit: ', char(8364), num2str(profit_Q5, '%.2f')], 'Color', 'red')
grid on
xlim([minworkers, maxworkers])
xlabel('Number of workers')
ylabel(['Profit in euros [', char(8364), ']'])
title('Q5) Maximum profit as a function of number of workers')
line([minworkers, maxworkers],[profit_Q5 profit_Q5],'Color','red','LineStyle','--')

clear w xloop profitloop

%% Plot constraint figure
x1 = 0:4400;                                                        % Define linearly spaced vector for number of Phone-Xs
c1x1 = (b(1)-A(1,2)*x2)/A(1,1);                                     % Constraint memory cells
c2x1 = (w_Q5*(160-2*E2)-(7-(w_Q5-50)/20)*x2)/(4-(w_Q5-50)/20);      % Constraint work hours
c3x1 = ones(1,length(x2))*-b(3);                                    % Constraint minimum supply Phone-Xs
c4x2 = ones(1,length(x1))*-b(4);                                    % Constraint minimum supply Phone-Ss
p1x1 = (profit_Q5+c(2)*x2+w_Q5*(2000+15*E3))/-c(1);                 % Profit line

figure                                                              % Plot constraint lines, profit lines and optimum
plot(x2, c1x1, x2, c2x1, x2, c3x1, c4x2, x1, x_Q5(2), x_Q5(1), 'r*', x2, p1x1, 'k:')
hold on

patch([1000, 1000, 1488, 1466],[1545, 800, 800, 846],'red')         % Visualise feasible solution set
alpha(0.2)

legend('Constraint: memory cells', 'Constraint: maximum work hours', 'Constraint: minimum supply Phone-Xs', 'Constraint: minimum supply Phone-Ss', 'Optimal solution', 'Profit lines', 'Feasible set')

for profit = -50000:50000:5000000                                   % Draw profit lines and disable legend labels
    p1x1 = (profit+c(2)*x2+w_Q5*-(2000+15*E3))/-c(1);
    p = plot(x2, p1x1, 'k:');
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

ylim([0 max(max([c1x1', c2x1', c3x1']))*1.1])
xlim([0 max(x2)])
xlabel('Number of Phone-Ss (x_2)')
ylabel('Number of Phone-Xs (x_1)')
title(['Q5) Plot showing constraint lines and optimum (number of workers = ', num2str(w_Q5), ')'])
