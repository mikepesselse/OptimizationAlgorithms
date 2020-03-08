% Third Assignment Nonlinear Programming for the SC42055 Otimization in Systems and
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

clear StudentNumber1 StudentNumber2

%% Define variables
N = 120;        % number of intervals of length T

% The following section defines the number of times various algorithms are
% run and has a large impact on the computation time of the script.
% The runtime of the script with the recommended values [10, 2, 2, 10, 0]
% is around 370 seconds on a 2019 TU Delft laptop.
n_Q5_ms = 10;   % Number of times SQP algorithm is run for unconstrained optimisation problem
% For reference: 10 runs take 53 seconds
n_Q5_sa = 2;    % Number of times simulated annealing algorithm is run for unconstrained optimisation problem
% For reference: 10 runs take 401 seconds
n_Q5_ga = 2;    % Number of times genetic algorithm is run for unconstrained optimisation problem
% For reference: 10 runs take 624 seconds
n_Q6_ms = 10;   % Number of times SQP algorithm is run for unconstrained optimisation problem
% For reference: 10 runs take 105 seconds
n_Q6_ga = 0;    % Number of times genetic algorithm is run for unconstrained optimisation problem
% For reference: 10 runs take 2269 seconds


%% Q4: Run algorithm for various starting points
% Define options for SQP algorithm, starting values and bounds on r(k)
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1e6, 'Display', 'off');
x0_1 = 0*ones(1, N);
x0_2 = 0.99*ones(1, N);
lb = zeros(1, N);
ub = ones(1, N);

% Run the algorithm for the two different starting values
[r_Q4_1, TTS_Q4_1, exitflag_Q4_1] = fmincon(@TTS,x0_1,[],[],[],[],lb,ub,[],options);
[r_Q4_2, TTS_Q4_2, exitflag_Q4_2] = fmincon(@TTS,x0_2,[],[],[],[],lb,ub,[],options);

% Output answers
str = ['For Question 4:\n', 'The optimal value of the TTS is ',num2str(TTS_Q4_1, '%.2f'), ' when using an SQP algorithm with a starting value of ',num2str(x0_1(1), '%.2f'), ' for r(k)\n',...
    'The optimal value of the TTS is ',num2str(TTS_Q4_2, '%.2f'), ' when using an SQP algorithm with a starting value of ',num2str(x0_2(1), '%.2f'), ' for r(k)', '\n\n'];
fprintf(str)
clear str

%% Q5: Multistart sqp algorithm
% Define options for the multistart SQP algorithm
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1e6, 'Display', 'off');
problem = createOptimProblem('fmincon','objective', @TTS, 'x0', rand(1,N),'lb',lb,'ub',ub,'options',options);
ms = MultiStart('Display', 'off');

% Run the multistart SQP algorithm
[r_Q5_ms, TTS_Q5_ms, exitflag_Q5_ms] = run(ms, problem, RandomStartPointSet('NumStartPoints',n_Q5_ms));

% Output answers
str = ['For Question 5:\n', 'The optimal value of the TTS is ',num2str(TTS_Q5_ms, '%.2f'), ' when using a multistart SQP algorithm with ',num2str(n_Q5_ms),' uniformly distributed starting values', '\n'];
fprintf(str)
clear str

%% Q5: Run simulated annealing algorithm
% Define options for the simulated annealing algorithm
options = optimoptions('simulannealbnd','PlotFcns',{@saplotbestx,@saplotbestf,@saplotx,@saplotf}, 'MaxIterations', 500, 'Display', 'off', 'ReannealInterval', 100);

% Initialise vectors
r_Q5_sa_loop = zeros(n_Q5_sa, N);
TTS_Q5_sa_loop = zeros(1, n_Q5_sa);
exitflag_Q5_sa_loop = zeros(1, n_Q5_sa);

% Run the simulated annealing algorithm a set number of times
for i = 1:n_Q5_sa
    [r_Q5_sa_loop(i,:), TTS_Q5_sa_loop(i), exitflag_Q5_sa_loop(i)] = simulannealbnd(@TTS, rand(1, N), lb, ub, options);
end

% Find the lowest value for TTS and its corresponding r(k) values
[TTS_Q5_sa, index] = min(TTS_Q5_sa_loop);
r_Q5_sa = r_Q5_sa_loop(index,:);
exitflag_Q5_sa = exitflag_Q5_sa_loop(index);

% Output answers
str = ['The optimal value of the TTS is ',num2str(TTS_Q5_sa, '%.2f'), ' when using a multistart simulated annealing algorithm with ', num2str(n_Q5_sa), ' uniformly distributed starting values', '\n'];
fprintf(str)
clear str

%% Q5: Run genetic algorithm
% Define options for the genetic algorithm
options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf, 'Display', 'off', 'MaxGenerations', 1000);

% Initialise vectors
r_Q5_ga_loop = zeros(n_Q5_ga, N);
TTS_Q5_ga_loop = zeros(1, n_Q5_ga);
exitflag_Q5_ga_loop = zeros(1, n_Q5_ga);

% Run the genetic algorithm a set number of times
for i = 1:n_Q5_ga
    [r_Q5_ga_loop(i,:), TTS_Q5_ga_loop(i), exitflag_Q5_ga_loop(i)] = ga(@TTS, N, [], [], [], [], lb, ub,[],[], options);
end

% Find the lowest value for TTS and its corresponding r(k) values
[TTS_Q5_ga, index] = min(TTS_Q5_ga_loop);
r_Q5_ga = r_Q5_ga_loop(index,:);
exitflag_Q5_ga = exitflag_Q5_ga_loop(index);

% Output answers
str = ['The optimal value of the TTS is ',num2str(TTS_Q5_ga, '%.2f'), ' when using a genetic algorithm with ', num2str(n_Q5_ga), ' runs', '\n\n'];
fprintf(str)
clear str

%% Q5: Obtain values
% Simulate the model with the optimal values of r(k) found by the
% multistart SQP algorithm
[~, x_Q5, V_Q5, ~, ~] = TTS(r_Q5_ms);

%% Q6: Run algorithms with extra constraint
% Define options for the multistart SQP algorithm
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1e6, 'StepTolerance', 1e-16, 'Display', 'off');
problem = createOptimProblem('fmincon','objective', @TTS, 'x0', rand(1,N),'lb',lb,'ub',ub,'options',options, 'nonlcon', @TTS_con);
ms = MultiStart('Display', 'off');

% Run the multistart SQP algorithm
[r_Q6_ms, TTS_Q6_ms, exitflag_Q6_ms] = run(ms, problem, RandomStartPointSet('NumStartPoints',n_Q6_ms));

% Output answers
str = ['For Question 6:\n', 'The optimal value of the TTS is ',num2str(TTS_Q6_ms, '%.2f'), ' when using a multistart SQP algorithm with ',num2str(n_Q6_ms),' uniformly distributed starting values', '\n'];
fprintf(str)
clear str


%% Q6: Genetic algorithm
if n_Q6_ga ~= 0
    % Define options for the genetic algorithm
    options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf, 'Display', 'off', 'MaxGenerations', 1000);
    
    % Initialise vectors
    r_Q6_ga_loop = zeros(n_Q6_ga, N);
    TTS_Q6_ga_loop = zeros(1, n_Q6_ga);
    exitflag_Q6_ga_loop = zeros(1, n_Q6_ga);
    
    % Run the genetic algorithm a set number of times
    for i = 1:n_Q6_ga
        [r_Q6_ga_loop(i,:), TTS_Q6_ga_loop(i), exitflag_Q6_ga_loop(i)] = ga(@TTS, N, [], [], [], [], lb, ub,@TTS_con,[], options);
    end
    
    % Find the lowest value for TTS and its corresponding r(k) values
    [TTS_Q6_ga, index] = min(TTS_Q6_ga_loop);
    r_Q6_ga = r_Q6_ga_loop(index,:);
    exitflag_Q6_ga = exitflag_Q6_ga_loop(index);
    
    % Output answers
    str = ['The optimal value of the TTS is ',num2str(TTS_Q6_ga, '%.2f'), ' when using a genetic algorithm with ', num2str(n_Q6_ga), ' run', '\n'];
    fprintf(str)
    clear str
end

%% Q6: Obtain values
% Simulate the model with the optimal values of r(k) found by the
% multistart SQP algorithm
[~, x_Q6, V_Q6, ~, ~] = TTS(r_Q6_ms);

%% Q7: Obtain values for the no control case
% Simulate the model with r(k) equal to one, i.e., no control on the
% on-ramp flow
[TTS_Q7, x_Q7, V_Q7, ~, ~] = TTS(ones(1,N));

% Output answers
str = ['\n', 'For Question 7:\n', 'The optimal value of the TTS is ',num2str(TTS_Q7, '%.2f'), ' when there is no control on the ramp metering rate', '\n\n'];
fprintf(str)
clear str

%% Plot figures
PlotTitles = 1;         % Set value to 1 for displaying figure titles

% Plot r(k)
figure; plot(r_Q4_1)
xlabel('Time interval'); ylabel('Value of r')
if PlotTitles == 1; title('Value of r for starting value r(k) = 0'); end
xlim([0 N]) ;ylim([-.1 1.1])

figure; plot(r_Q4_2)
xlabel('Time interval'); ylabel('Value of r')
if PlotTitles == 1; title('Value of r for starting value r(k) = 0.99'); end
xlim([0 N]); ylim([-.1 1.1])

figure; plot(r_Q6_ms)
xlabel('Time interval'); ylabel('Value of r')
if PlotTitles == 1; title('Value of r with constraint'); end
xlim([0 N]); ylim([-.1 1.1])

figure; plot(ones(1,N))
xlabel('Time interval'); ylabel('Value of r')
if PlotTitles == 1; title('Value of r without control'); end
xlim([0 N]); ylim([-.1 1.1])


% Plot densities
figure; hold on; plot(x_Q5(1,:)); plot(x_Q5(2,:)); plot(x_Q5(3,:)); plot(x_Q5(4,:));
xlabel('Time interval'); ylabel('Traffic density per segment'); legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4')
if PlotTitles == 1; title('Traffic density without constraint'); end
xlim([0 N]); ylim([0 115])

figure; hold on; plot(x_Q6(1,:)); plot(x_Q6(2,:)); plot(x_Q6(3,:)); plot(x_Q6(4,:));
xlabel('Time interval'); ylabel('Traffic density per segment'); legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4')
if PlotTitles == 1; title('Traffic density with constraint'); end
xlim([0 N]); ylim([0 115])

figure; hold on; plot(x_Q7(1,:)); plot(x_Q7(2,:)); plot(x_Q7(3,:)); plot(x_Q7(4,:));
xlabel('Time interval'); ylabel('Traffic density per segment'); legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4')
if PlotTitles == 1; title('Traffic density without control'); end
xlim([0 N]); ylim([0 115])


% Plot waiting queue
figure; plot(x_Q5(9,:));
xlabel('Time interval'); ylabel('Waiting queue')
if PlotTitles == 1; title('Waiting queue without constraint'); end
xlim([0 N]); ylim([0 max(x_Q5(9,:))*1.1])

figure; plot(x_Q6(9,:));
xlabel('Time interval'); ylabel('Waiting queue')
if PlotTitles == 1; title('Waiting queue with constraint'); end
xlim([0 N]); ylim([0 max(x_Q5(9,:))*1.1])

figure;plot(x_Q7(9,:));
xlabel('Time interval'); ylabel('Waiting queue')
if PlotTitles == 1; title('Waiting queue without control'); end
xlim([0 N]); ylim([0 max(x_Q5(9,:))*1.1])


% Plot mean velocities
figure; hold on; plot(x_Q5(5,:)); plot(x_Q5(6,:)); plot(x_Q5(7,:)); plot(x_Q5(8,:));
xlabel('Time interval'); ylabel('Mean velocity per segment [km/h]'); legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4')
if PlotTitles == 1; title('Mean velocity per segment without constraint'); end
xlim([0 N]); ylim([-5 max(x_Q6(5,:))*1.1])

figure; hold on; plot(x_Q6(5,:)); plot(x_Q6(6,:)); plot(x_Q6(7,:)); plot(x_Q6(8,:));
xlabel('Time interval'); ylabel('Mean velocity per segment [km/h]'); legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4')
if PlotTitles == 1; title('Mean velocity per segment with constraint'); end
xlim([0 N]); ylim([-5 max(x_Q6(5,:))*1.1])

figure; hold on; plot(x_Q7(5,:)); plot(x_Q7(6,:)); plot(x_Q7(7,:)); plot(x_Q7(8,:));
xlabel('Time interval'); ylabel('Mean velocity per segment [km/h]'); legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4')
if PlotTitles == 1; title('Mean velocity per segment without control'); end
xlim([0 N]); ylim([-5 max(x_Q6(5,:))*1.1])


%% Functions to calculate TTS and constraints
function [sum_y, x, V, c, ceq] = TTS(r)

StudentNumber1 = [4 3 6 7 1 4 6];               % Student number Bart de Jong
StudentNumber2 = [4 3 0 0 5 6 4];               % Student number Mike Pesselse
E1 = StudentNumber1(5) + StudentNumber2(5);     % Calculate variable E1
E2 = StudentNumber1(6) + StudentNumber2(6);     % Calculate variable E2
E3 = StudentNumber1(7) + StudentNumber2(7);     % Calculate variable E3

% METANET parameters
tau     = 18/3600;  % Model parameter [h]
mu      = 60;       % Model parameter [km^2 / h]
Cr      = 2000;     % Model parameter [veh / h]
rho_m   = 180;      % Model parameter [veh / km]
K       = 40;       % Model parameter [veh / km]
a       = 1.87;     % Model parameter []
v_f     = 100;      % Free-flow speed [km / h]
rho_c   = 33.5;     % Critical density [veh / km]

% Other parameters
N       = 120;      % number of intervals of length T
T       = 10/3600;  % Simulation time step [h]
L       = 1;        % Length of each segment [km]
lambda  = 4;        % Number of lanes [#]
Dr          = ones(1, N+1)*1500;    % Ramp demand [veh / h]
q0(1:30)    = 8000 + 100*E1;        % Flow entering mainline for k <  30 [veh / h]
q0(31:N+1)  = 4000 + 100*E2;        % Flow entering mainline for k >= 30 [veh / h]

% Initial conditions
rho1(1) = 30;       % Initial density first segment [veh / km]
rho2(1) = 30;       % Initial density second segment [veh / km]
rho3(1) = 30;       % Initial density third segment [veh / km]
rho4(1) = 30;       % Initial density fourth segment [veh / km]
v1(1)   = 80;       % Initial speed first segment [km / h]
v2(1)   = 80;       % Initial speed second segment [km / h]
v3(1)   = 80;       % Initial speed third segment [km / h]
v4(1)   = 80;       % Initial speed fourth segment [km / h]
w_r(1) 	= 0;        % Initial ramp queue [#]

x(:, 1) = [rho1(1) rho2(1) rho3(1) rho4(1) v1(1) v2(1) v3(1) v4(1) w_r(1)];

% Initialise matrices
qr = zeros(1, N);
V  = zeros(4, N);
y  = zeros(1, N);

% Model dynamics
for k = 1:120
    
    % Traffic flow that enters the fourth segment from the on-ramp
    qr(k) = min([r(k)*Cr, Dr(k) + x(9,k)/T, Cr*(rho_m - x(4,k))/(rho_m - rho_c)]);
    
    % Desired speed V per segment
    V(1, k) = v_f*exp(-1/a * (x(1, k)/rho_c)^a);
    V(2, k) = v_f*exp(-1/a * (x(2, k)/rho_c)^a);
    V(3, k) = v_f*exp(-1/a * (x(3, k)/rho_c)^a);
    V(4, k) = v_f*exp(-1/a * (x(4, k)/rho_c)^a);
    
    % Traffic density rho per segment
    x(1,k+1) = x(1,k) + T/(lambda*L)*(q0(k)                - lambda*x(1,k)*x(5,k)        );
    x(2,k+1) = x(2,k) + T/(lambda*L)*(lambda*x(1,k)*x(5,k) - lambda*x(2,k)*x(6,k)        );
    x(3,k+1) = x(3,k) + T/(lambda*L)*(lambda*x(2,k)*x(6,k) - lambda*x(3,k)*x(7,k)        );
    x(4,k+1) = x(4,k) + T/(lambda*L)*(lambda*x(3,k)*x(7,k) - lambda*x(4,k)*x(8,k) + qr(k));
    
    % Mean speed v per segment
    x(5,k+1) = x(5,k) + (T/tau)*(V(1, k) - x(5,k))                                  - (mu*T)/(tau*L)*(x(2,k) - x(1,k))/(x(1,k) + K);
    x(6,k+1) = x(6,k) + (T/tau)*(V(2, k) - x(6,k)) + (T/L)*x(6,k)*(x(5,k) - x(6,k)) - (mu*T)/(tau*L)*(x(3,k) - x(2,k))/(x(2,k) + K);
    x(7,k+1) = x(7,k) + (T/tau)*(V(3, k) - x(7,k)) + (T/L)*x(7,k)*(x(6,k) - x(7,k)) - (mu*T)/(tau*L)*(x(4,k) - x(3,k))/(x(3,k) + K);
    x(8,k+1) = x(8,k) + (T/tau)*(V(4, k) - x(8,k)) + (T/L)*x(8,k)*(x(7,k) - x(8,k))                                                ;
    
    % Set a minimum velocity of 4 km/h
    x(5, k+1) = max(4, x(5, k+1));
    x(6, k+1) = max(4, x(6, k+1));
    x(7, k+1) = max(4, x(7, k+1));
    x(8, k+1) = max(4, x(8, k+1));
    
    % Queue on the on-ramp w_r
    x(9,k+1) = x(9,k) + T*(Dr(k) - qr(k));
    
    
end

% calculate TTS for each time interval and sum them up
for i = 0:N-1
    y(i+1) = T*x(9, i+1) + T*L*lambda*x(1, i+1) + T*L*lambda*x(2, i+1) + T*L*lambda*x(3, i+1) + T*L*lambda*x(4, i+1);
end
sum_y = sum(y);

% Define the inequality constraint on the on-ramp queue
c = x(9,:)- 100 + E3;
ceq = [];

end

% Function outputting the constraints for implementation in fmincon
function [c, ceq] = TTS_con(r)
[~, ~, ~, c, ceq] = TTS(r);
end
