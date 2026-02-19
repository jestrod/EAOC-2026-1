% .-------------------------------------------------.
% |                                                 |
% |     __  __      _ ___              __           |
% |    / / / /___  (_)   |  ____  ____/ /__  _____  |
% |   / / / / __ \/ / /| | / __ \/ __  / _ \/ ___/  |
% |  / /_/ / / / / / ___ |/ / / / /_/ /  __(__  )   |
% |  \____/_/ /_/_/_/  |_/_/ /_/\__,_/\___/____/    |
% |                                                 |
% '-------------------------------------------------'
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ~~~~~~~~~~~~~~~~~~ POINT 1 ~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    Author  : Esteban Rodríguez
%    Date    : 02.2026
%    Purpose : Fuzzy controller for a inverted pendullum
%    Notes   : EAOC 2026-1
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% --------------------------------- PARAM --------------------------------
clc; clear; close all;

%Integration step
step = 0.001;
% Controller period
T = 0.1;
% Control action counter
counter = round(T/step);
% Simulation time
tfinal = 5;
% Samples numbers
N = round(tfinal/step);

% Gains
% e
g1 = 1.5;
% de
g2 = 0.39;
% u
g0 = 9.8;  

% Setpoint angle
r= 0;

%% ------------------------------ INIT STATES -----------------------------

x = zeros(3,1);
% e(0)
x(1) = 0.4;   
% de(0)
x(2) = 0;
% u(0)
x(3) = 0;

%% --------------------------------- VARS ---------------------------------

% Discrete time vector
time = zeros(N,1);
% Pendulum angles vector
y = zeros(N,1);
% Actuator's control action vector
u = zeros(N,1);
% e(k-1) init
eold = 0;
%% ------------------------------ SIMULATION ------------------------------
for k = 1:N
    
    time(k) = (k-1)*step;
    y(k)    = x(1);
    
    % COMPUTE FUZZY CONTROLLER
    % Compute controller only if a control period is reached
    if mod(k,counter)==0
        
        % Inputs computation
        e  = r-x(1);
        
        de = (e - eold)/T;
        
        eold = e;
        
        % Input scaling
        e  = g1*e;
        de = g2*de;
        
        % Fuzzy controller computation for the actual state
        uf = fuzzy_controller(e,de);
        % Control action scaling
        u(k) = g0*uf;
        
    else
        u(k) = u(max(k-1,1));
    end
        
    % Compute the systems next state using RK4

    k1 = step*RK4_pendulum(x,u(k));
    k2 = step*RK4_pendulum(x+0.5*k1,u(k));
    k3 = step*RK4_pendulum(x+0.5*k2,u(k));
    k4 = step*RK4_pendulum(x+k3,u(k));
    
    x = x + (k1+2*k2+2*k3+k4)/6;
    
end

%% -------------------------------- GRAPHS --------------------------------       
    
    % y vs t
    figure
    plot(time,y,'LineWidth',1.5)
    xlabel('TIME [s]')
    ylabel('PENDULUM ANGLE y(t) [rad]')
    title('INVERTED PENDULUM RESPONSE')
    grid on
    
    % u vs y
    figure
    plot(time,u,'LineWidth',1.5)
    xlabel('TIME [s]')
    ylabel('u(t) [N]')
    title('ACTUATOR RESPONSE OVER THE PLANT')
    grid on
    
    % Non-lineal control surface

    % Evaluate different inputs of e and de on the controller
    e_vals  = linspace(-1,1,41);
    de_vals = linspace(-1,1,41);
    
    [E,DE] = meshgrid(e_vals,de_vals);
    
    % Compute the control surface values for the fuzzy controller
    U = zeros(size(E));
    
    for i = 1:length(e_vals)
        for j = 1:length(de_vals)
            U(j,i) = fuzzy_controller(E(j,i), DE(j,i));
        end
    end
    
    % Visualization
    figure
    surf(E,DE,U)
    xlabel('e')
    ylabel('de')
    zlabel('u')
    title('Control surface implemented by the fuzzy controller')
    colormap(white)
    grid on


%% --------------------------- LOCAL FUNCTIONS ----------------------------  
% ~~~~~~~~~~~~~~~~~~~~~ dx = RK4_pendulum(x,delta) ~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    Author  : Esteban Rodríguez
%    Date    : 02.2026
%    Purpose : RK4 implementation for inverted
%              pendulum plant.
%    Notes   : EAOC 2026-1
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function dx = RK4_pendulum(x,delta)

dx = zeros(3,1);

dx(1) = x(2);

dx(2) = ( 9.8*sin(x(1)) ...
         + cos(x(1)) * ((-x(3) - 0.25*x(2)^2*sin(x(1)))/1.5) ) ...
         / ( 0.5*((4/3)-(1/3)*(cos(x(1))^2)) );

dx(3) = -100*x(3) + 100*delta;

end

% ~~~~~~~~~~~~~~~~~~~~ uf = fuzzy_controller(e,de) ~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    Author  : Esteban Rodríguez
%    Date    : 02.2026
%    Purpose : Fuzzy controller implementation
%    Notes   : EAOC 2026-1
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function uf = fuzzy_controller(e,de)

% Fuzzification interface for e and de
% e and de certainly is obtained
mu_e  = fuzzify(e);
mu_de = fuzzify(de);

% Fuzzy controller rules definition
% Rows :    e[-1 -0.5 0 0.5 1].T
% Colums :  de[-1 -0.5 0 0.5 1]
% If e<0 -> The actuator must "push" the car
% If e>0 -> The actuator must "pull" the car
rules = [...
    1  1  0.5 0.5  0;
    1  0.5 0.5  0    -0.5;
    0.5 0.5  0    -0.5  -0.5;
    0.5  0    -0.5  -0.5  -1;
     0    -0.5  -0.5  -1    -1]; %e[-1 -0.5 0 0.5 1] vs de[-1 -0.5 0 0.5 1].T  

num = 0;
den = 0;

% Inference interface and defuzzification with CE method
for i=1:5
    for j=1:5
        
        % AND operator to compute implication certainty
        mu_premise = mu_e(i)*mu_de(j);

        %CE method
        mu_out = rules(i,j)* mu_premise;
        num = num + mu_out;
        den = den + mu_premise;
        
    end
end

uf = num/den;

end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~ mu = fuzzify(x) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    Author  : Esteban Rodríguez
%    Date    : 02.2026
%    Purpose : Variable fuzzification 
%    Notes   : EAOC 2026-1
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function mu = fuzzify(x)

% Membership functions centers
centers = [-1 -0.5 0 0.5 1];
% Membership width
width   = 0.5;

mu = zeros(1,5);

% Compute u for each fuzzy variable
% The MF are implemented as triangular functions
for i=1:5
    mu(i) = max(1 - abs(x-centers(i))/width , 0);
end

end

% ~~~~~~~~~~~~~~ mu_area = area_implied_fs(mu_premise,w) ~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    Author  : Esteban Rodríguez
%    Date    : 02.2026
%    Purpose : Compute the area of the truncated triangle with
%              height mu_premise
%    Notes   : EAOC 2026-1
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function mu_area = area_implied_fs(mu_premise,w) 
    mu_area = w * (mu_premise - (mu_premise^2)/(2)); 
end