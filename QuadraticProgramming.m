% Second Assignment Quadratic Programming for the SC42055 Otimization in Systems and
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

RunExercise4 = 0;           % Set to 1 for running Exercise 4 (NOT RECOMMENDED!)
% This algorithm will run for at least 2.5 hours, without the possibility
% of stopping. Therefore, we recommend not to set RunExercise4 to 1.


%% Read data
filename    = 'qp_2019_measurements.csv';
data        = readtable(filename,'Delimiter',';', 'PreserveVariableNames', true);
N           = height(data);
dt          = 3600;         % Number of seconds between measurements
dT          = 273.15;       % Conversion for Celsius to Kelvin

% Store data in variables
T_c         = data(:,1);
q_dot_in    = data(:,2);
q_dot_solar = data(:,3);
q_dot_out   = data(:,4);
T_amb       = data(:,5);
I           = data(:,6);
phi_eur     = data(:,7);

% Convert to SI-units
T_c         = table2array(T_c)          +dT;
q_dot_in    = table2array(q_dot_in)     *1e6;
q_dot_solar = table2array(q_dot_solar)  *1e6;
q_dot_out   = table2array(q_dot_out)    *1e6;
T_amb       = table2array(T_amb)        +dT;
I           = table2array(I);
phi_eur     = table2array(phi_eur)      /3600;



%% Exercise 2
% Define model matrices
Y	= T_c(2:end);
phi	= [T_c(1:end-1) q_dot_solar(1:end-1) q_dot_out(1:end-1) q_dot_in(1:end-1) T_amb(1:end-1)];
H   = 2*(phi.'*phi);
c   = -2*phi.'*Y;
d  	= Y.'*Y;

% Define constraints
lb  = [-0.99 -Inf -Inf -Inf -Inf]';
ub  = [0.99 +Inf +Inf +Inf +Inf]';
Aeq = [0 0 1 1 0; 1 0 0 0 1];
beq = [0; 1];

% Run algorithm
options = optimoptions('quadprog','Algorithm','interior-point-convex', 'Display', 'off');
[x,~,exitflag_Q2,~]=quadprog(H,c,[],[],Aeq,beq,lb,ub,[],options);

% Extract variables a1, a2 and a3
a1 = x(4)/dt;
a2 = x(2)/dt;
a3 = x(5)/dt;

% Output answers
str = ['For Question 2:\na1 = ', num2str(a1), ' ', char(176), 'C/J', '\n', 'a2 = ', num2str(a2), ' ', char(176), 'C/J', '\na3 = ', num2str(a3),'  1/J', '\n\n'];
fprintf(str)
clear str


%% Exercise 3, first method
Tmin = 390 + dT;
T_1  = 331 + dT;

% Defining the minimization function:
c_fun = zeros(1, N);                        % Pre-define size
for i = 1:N
    c_fun(i) = phi_eur(i)*I(i)*dt;          % Calculate cost for each timestep
end

% Inequality constraint number one, defining the maximum collector area:
A3_ineq1 = zeros(1, N);                     % Make a zero-vector
A3_ineq1(1,1) = 1;                          % Let the first entry equal 1
b3_ineq1 = 4e6/max(I);                      % Define the b-vector

% Equality constraint number two, ensuring all entries in x3 are equal:
A3_eq = zeros(N-1, N);                      % Pre-define size
for i = 1:N-1
    for j = 1:N
        if i == j
            A3_eq(i,j) = 1;                 % Make diagonal entries equal to 1
        elseif i == j-1
            A3_eq(i, j) = -1;               % Make the entries to the right of the diagonal equal to -1
        end
    end
end
b3_eq = zeros(N-1, 1);                      % Define the b-vector

% Inequality constraint number two, ensuring that the temperature is not below 390 degrees celcius for certain timesteps:
a_Q3 = 1 - a3*dt;
b_Q3 = zeros(1, N);
for i = 1:N
    b_Q3(i) = a2*dt*I(i);
end
c_Q3 = zeros(1,N);
for i = 1:N
    c_Q3(i) = -a1*dt*q_dot_out(i) + a1*dt*q_dot_in(i) + a3*dt*T_amb(i);
end

A3_ineq2 = zeros(N-1, N);               % Construct A-matrix
for j=1:N-1
    k = j+1;
    for i = 1:j
        A3_ineq2(j) = A3_ineq2(j) + a_Q3^(i-1)*b_Q3(k-i);
    end
end

b3_ineq2 = zeros(N-1, 1);               % Construct b-matrix
sum = 0;
for j = 1:N-1
    k = j+1;
    if I(k) >= 500
        for i = 1:j
            sum = sum + a_Q3^(i-1)*c_Q3(k-i);
        end
        b3_ineq2(j) = Tmin - a_Q3^(j)*T_1 - sum;
        sum = 0;
    else
        b3_ineq2(j) = -1e10;
    end
end

% Combine constraint matrices
A3_ineq = [A3_ineq1; -A3_ineq2];
b3_ineq = [b3_ineq1; -b3_ineq2];

% Run optimization algorithm and extract optimal collector area
options = optimoptions('linprog','Algorithm', 'interior-point', 'Display', 'off');
lb      = zeros(1, N);
ub      = [];
[x3, ~, exitflag_Q3_1, ~] = linprog(c_fun, A3_ineq, b3_ineq, A3_eq, b3_eq, lb, ub, [], options);
Ac      = x3(1);


% Calculate temperature
T = zeros(N, 1);
T(1) = T_1;
for k=2:N
    T(k) = (1-a3*dt)*T(k-1) + a2*dt*I(k-1)*Ac - a1*dt*q_dot_out(k-1) + a1*dt*q_dot_in(k-1) + a3*dt*T_amb(k-1);
end

% Calculate cost
cost=zeros(1, N);
for k=2:N
    cost(k) = cost(k-1) + phi_eur(k)*I(k)*dt*Ac;
end


% Output answers
str = ['For Question 3, method 1:\n', 'The optimal collector area equals ', num2str(Ac), ' m^2', '\n', 'The cost is ', num2str(cost(end)/1e9, '%.2f'), ' billion euros', '\n\n'];
fprintf(str)
clear str


%% Exercise 3, second method
% Gather A and B matrices from question 2 and define constants
a = [a1 a2 a3];
A = (1 - a(3)*dt);
B = [a(2)*dt -a(1)*dt a(1)*dt a(3)*dt];

Qmax = 4*10^6;

% Build Equality constraint Aeq*X = Beq
Aeq1 = [-B(1)*I(1) 1];
Aeq1 = kron(eye(N),Aeq1);
Aeq2 = zeros(N,2*N);
Beq2 = zeros(N,1);

for n = 0:8758
    Aeq1(2+n,3+2*n) = -B(1)*I(n+2);
    Aeq1(2+n,2+2*n) = -A;
    Aeq2(1+n,1+2*n) = 1;
    Aeq2(1+n,3+2*n) = -1;
end

Beq1 = zeros(N, 1);
Beq1 = A*T_1+B(2)*q_dot_out(1)+B(3)*q_dot_in(1)+B(4)*T_amb(1);
for k = 2:8760
    Beq1 = [Beq1; B(2)*q_dot_out(k)+B(3)*q_dot_in(k)+B(4)*T_amb(k)];
end

Aeq = [Aeq1; Aeq2];
Beq = [Beq1; Beq2];

% Build first Inequality constraint Aineq*X = Bineq
Aineq = [1 0;-1 0;0 1];
Aineq = kron(-1*eye(8760), Aineq);
for n2 = 0:8759
    Aineq(1+3*n2,1+2*n2) = I(n2+1);
    Aineq(2+3*n2,1+2*n2) = -I(n2+1);
end
Aineq(end-2:end, :) = [];

Bineq = zeros(3*(N-1), 1);
Bineq = [Qmax; 0; 1e10];
for i = 2:8760-1
    if I(i+1) >=500
        Bineq = [Bineq; Qmax; 0; -Tmin];
    else
        Bineq = [Bineq; Qmax; 0; 1e10];
    end
end

% Build const function c^T*X
C = zeros(1, 2*N);
C = [phi_eur(1)*dt*I(1) 0];
for i2 = 2:N
    C = [C, phi_eur(i2)*dt*I(i2) 0];
end

%Upper and lowerbound
lb = zeros(size(C));
ub = inf(size(C));

options = optimoptions('linprog','Algorithm','dual-simplex', 'Display', 'off');
[x,fval,flag,output] = linprog(C,Aineq,Bineq,Aeq,Beq,lb,ub,[],options);

X = reshape(x,2,[]);

cost2=zeros(1, N);
for k=2:N
    cost2(k) = cost2(k-1) + phi_eur(k)*I(k)*dt*x(1);
end

% Output answers
str = ['For Question 3, method 2:\n', 'The optimal collector area equals ', num2str(x(1)), ' m^2', '\n', 'The cost is ', num2str(cost2(end)/1e9, '%.2f'), ' billion euros', '\n\n'];
fprintf(str)
clear str



%% Exercise 4
% This algorithm will run for at least 2.5 hours, without the possibility
% of stopping. Therefore, we recommend not to set RunExercise4 to 1.
if RunExercise4 == 1
    Tref = 475+dT;
    add_cost = (0.1+E2/10);
    Tmax = 525+dT;
    
    %compute new C^T
    C4 = [phi_eur(1)*dt*I(1) -2*add_cost*Tref];
    for i2 = 2:N
        C4 = [C4, phi_eur(i2)*dt*I(i2) -2*add_cost*Tref];
    end
    
    
    %compute new inequality constraints
    Aineq4 = [1 0;-1 0;0 -1; 0 1];
    Aineq4 = kron(eye(8760), Aineq4);
    for n2 = 0:8759
        Aineq4(1+4*n2,1+2*n2) = I(n2+1);
        Aineq4(2+4*n2,1+2*n2) = -I(n2+1);
    end
    Aineq4(end-3:end, :) = [];
    
    Bineq4 = [Qmax; 0; 1e10; Tmax];
    for i = 2:8760-1
        if I(i+1) >=500
            Bineq4 = [Bineq4; Qmax; 0; -Tmin; Tmax];
        else
            Bineq4 = [Bineq4; Qmax; 0; 1e10; Tmax];
        end
    end
    
    %compute Hessian
    H4 = zeros(2*N,2*N);
    for n = 0:8758
        H4(2+2*n,2+2*n) = 2*add_cost;
    end
    
    options = optimoptions('quadprog','Algorithm','interior-point-convex');
    [x4,~,exitflag_Q4,~]=quadprog(H4,C4,Aineq4,Bineq4,Aeq,Beq,lb,ub,[],options);
end

% Output answers
str = ['For Question 4:\n', 'The problem is unfeasible and yields, therefore, no results', '\n\n'];
fprintf(str)
clear str


%% Figure 1
subplot (3, 1, 1)
hold on
plot (T);
xlabel('Time in hours')
ylabel('Internal temperature T_c (K)')
axis ([0 N 0 max(T)])
plot ([0 , N] ,[Tmin ,Tmin] , '--')
title('Internal temperature T_c')
legend('Internal temperature', 'minimum temperature for I>500')

subplot (3, 1, 2)
plot (cost);
xlabel('Time in hours')
ylabel('Costs (euro/Wh)')
axis([0 N 0 2.0846e+09])
title('Cumulative costs')
legend('Cumulative costs')

subplot(3, 1, 3)
plot(I.*Ac)
hold on
plot ([0 , N] ,[4e6 ,4e6] , '--')
axis([0 N 0 4.5e6])
title('Solar heat input')
xlabel('Time in hours')
legend('Solar heat input', 'Maximum solar heat input')
ylabel('Solar heat input (W)')