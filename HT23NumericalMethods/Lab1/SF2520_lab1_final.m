%Lab 1


%% Part 1a)

a = (1/4)*[1 sqrt(11) 2]'; % Fixed vector
u0 = [0 0 1]'; %Initial condition
h_a = 0.1;

[result_a, time_a] = rk3_a(h_a, u0, 50);

subplot(1,2,1)
plot(time_a, result_a)
title("Plot of m relative the time")
xlabel("time")
ylabel("Magnetization")
legend("m(1)", "m(2)", "m(3)")
subplot(1,2,2)
plot3(result_a(1,:), result_a(2,:), result_a(3,:))
axis equal
hold on
quiver3(0,0,0,a(1), a(2), a(3))
title("Trajectory of m")
legend("Trajectory of m", "Plot of a", 'Location','north')
hold off




%% Part 1b)

N_b = [50, 100, 200, 400, 800, 1600];
T = 50;
h_b = T./N_b;
diff_m = zeros(1, length(N_b) - 1); % vector for saving error in each step
u0 = [0 0 1]'; %Initial condition
for i = 1:length(N_b)-1
    h1 = h_b(i);
    h2 = h_b(i + 1);
    [result_b1, time_b1] = rk3_a(h1, u0, T);
    [result_b2, time_b2] = rk3_a(h2, u0, T);
    diff_m(i) = norm(result_b1(:, end) - result_b2(:, end));
end
    
    loglog(h_b(1, 1:end-1), diff_m)
    hold on
    loglog(h_b(:,1:5),h_b(:,1:5).^3);
    hold off
    legend("error", "y=3x")
    grid on
    title("Error vs time step")
    xlabel("Stepsize, h")
    ylabel("Error difference")
    

%% Part 1c)

a1 = 1/4;
a2 = (1/4)*sqrt(11);
a3 = 2/4;

A1 = [0 -a3 a2] + 0.07*[-(a2^2+a3^2) a1*a2 a1*a3];
A2 = [a3 0 -a1] + 0.07*[a1*a2 -(a3^2+a1^2) a3*a2];
A3 = [-a2 a1 0] + 0.07*[a1*a3 a2*a3 -(a1^2 + a2^2)];

A = [A1 ; A2 ; A3];

eigA = eig(A);

h_stab = 0.1; % initial guess (a known stable step size)
s_vec = [-1,-1,-1];

% this loop looks for the smallest timestep h (with 0.01 
% steps) such that h*e_k (where e_k are the eigenvalues of A) 
% is a point outside the stability region

while (s_vec(1)<=0)&&(s_vec(2)<=0)&&(s_vec(3)<=0)
    for i=[1,2,3]
        s = stab_reg(h_stab * eigA(i,1));
        s_vec(i) = s;
    end
    h_stab = h_stab+0.01;
end



%% Part 1d)

u0 = [0 0 1]'; %Initial condition
[unstab_sol , unstab_tim] = rk3_a(2.07, u0, 500);
[stab_sol , stab_tim] = rk3_a(2.00, u0, 500);

plot(unstab_tim, unstab_sol(1,:));
title("Stability of the magnetization function")
xlabel("time")
ylabel("m(t)")
hold on 
plot(stab_tim ,stab_sol(1,:));
legend("h = 2.07", "h = 2.00")
hold off



%% Part 2a)

T_end = 10;
h = 0.0001; % Step length
r0 = [1.15, 0]'; % Initial relative position
dr0 = [0, -0.975]'; % Initial relative velocity
mu = 1/82.45; %  Ratio of moon/Earth mass
c0 = [-mu, 0].'; % Earth centre
c1 = [1 - mu, 0].'; %Moon center

% Letting Runge-Kutta method make the three initial positions
u1 = rk3_b(h, r0, dr0, h); 
u2 = rk3_b(h, u1(1:2), u1(3:4), h); 
u3 = rk3_b(h, u2(1:2), u2(3:4), h); 
u0 = [r0; dr0];
u_init = [u0, u1, u2, u3]; % first four positions and velocities


[sat_pos, sat_vel, tim] = Adam_Bash(u_init, h, T_end);

subplot(1,2,1)
plot(tim, sat_pos)
title("Position vs time")
legend("x", "y", "Location","northwest")

xlabel("time")
ylabel("Position")
subplot(1,2,2)
plot(tim, sat_vel)
title("Velocity vs time")
legend("x'", "y'", "Location","northwest")
%xlim([0 11])
xlabel("time")
ylabel("Velocity")
figure
plot(sat_pos(1,:), sat_pos(2,:))
title("The positions x(t) vs y(t)")
xlabel("x(t)")
ylabel("y(t)")
axis equal
hold on
plot(c1(1), c1(2), "ro",...
    'LineWidth',2,...
    'MarkerSize',7,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.7,0.5,0.7])
plot(c0(1), c0(2), "go",...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
legend("Satellite trajectory", "Moon centre", "Earth centre")
hold off

%% Part 2b)
r0 = [1.15, 0]'; % Initial relative position
dr0 = [0, -0.975]'; % Initial relative velocity
N_ab40 = 4400;
N_rk40 = 16500;
N_exp5 = 200000;
h_ab = 1/N_ab40;
u1 = rk3_b(h_ab, r0, dr0, h_ab) ;
u2 = rk3_b(h_ab, u1(1:2), u1(3:4), h_ab);
u3 = rk3_b(h_ab, u2(1:2), u2(3:4), h_ab);
u0 = [r0; dr0];
u_init = [u0, u1, u2, u3]; % first four positions and velocities


% N_rk3 = [705, 2100,16500];
% N_ab4 = [1035, 1246, 4400];
% N_ie = 2000000;

exact_vals = [0.4681, -0.2186, -1.4926; ...
              0.6355, -0.2136, -0.3339];

% below are three examples of how we calculated the correct value of N for
% each method. 



final_ab40 = Adam_Bash_b(u_init, 1/N_ab40, 40);
norm_ab40 = norm(final_ab40(1:2,:) - exact_vals(:,3));

final_rk40 = rk3_b(1/N_rk40, r0, dr0, 40);
norm_rk = norm(final_rk40(1:2) - exact_vals(:,3));

final_exp5 = expl_Eul(r0, dr0, 1/N_exp5, 5);
norm_exp = norm(final_exp5 - exact_vals(:,1));


%% Part 2c)
    r0 = [1.15, 0]'; % Initial relative position
    dr0 = [0, -0.975]'; % Initial relative velocity
    Y0 = [r0; dr0];
    mu = 1/82.45; %  Ratio of moon/Earth mass
    c0 = [-mu, 0].'; % Earth centre
    c1 = [1 - mu, 0].'; %Moon center
    B = [0, 1; -1, 0];
    %The 2nd derivative as defined in the Problem task
    dY = @(t, Y) [Y(3:4); -(1 - mu)*(Y(1:2) - c0)/(norm(Y(1:2) - c0))^3 - mu*(Y(1:2) - c1)/(norm(Y(1:2) - c1))^3 + 2*B*Y(3:4) + Y(1:2)];
    options = odeset('RelTol',1e-4);
    
    [TOUT20, YOUT20] = ode23(dY, [0 20],  Y0, options);
    
    %plot(TOUT, YOUT)
    %legend("x(t)", "y(t)", "x'(t)", "y'(t)")
    ex_val = [0.4681, -0.2186, -1.4926; ...
                0.6355, -0.2136, -0.3339];
    
    % Time T_end = 20
    
    r_out20 = YOUT20(end, 1:2);
    ex_val20 = norm(r_out20.' - ex_val(:,2));
    
    %largest step
    step_max = max(diff(TOUT20));
    step_min = min(diff(TOUT20));

    plot(TOUT20(1:end-1), diff(TOUT20))
    title("Step Size as a Function of Time" )
    xlabel("Time")
    ylabel("Step size")




%% Part 3a)

h=0.0001; %need timestep 0.0001
T=10; % end time
xa0=1; % initial values
xb0=0;
xc0=0;
x0 = [xa0;xb0;xc0];

[resk,timk, eigk] = rk3_c(h, x0,T);

subplot(1,2,1)
plot(timk,resk)
title("Plot of RK3 solution");
legend("xa","xb","xc");
xlabel("time");
ylabel("concentration");
hold on
subplot(1,2,2)
semilogy(timk,resk)
title("Semilogy of RK3 solution");
xlabel("time");
ylabel("concentration");
legend("xa","xb","xc");
hold off

%% part 3b)

% find jacobian and compute eigenvalues in each timestep

h=0.0007; %need approx. timestep 0.0007 for a stable solution
T=10; 

% initial values
xa0=1;
xb0=0;
xc0=0;
x0 = [xa0;xb0;xc0];
[res10,tim10, eig10] = rk3_c(h, x0, T); % this version of RK3 also keeps track of eigenvalues of the jacobian

plot(tim10, eig10, "linewidth", 2)

% plot(tim10, eig10(1,:))
% hold on
% plot(tim10, eig10(3,:))
% hold off
legend("eig1", "eig3");
title("Eigenvalues in every timestep")
xlabel("timesteps");
ylabel("eigenvalues");

eig_end=eig10(1,end); % this is the largest eigenvalue (in magnitude)
s = stab_reg(h*(eig_end)); % if s is negative, we are within the stability region for RK3

%% part 3c)

h=0.0001; %need approx. timestep 0.0007 for a stable solution
T=1000; 

% initial values
xa0=1;
xb0=0;
xc0=0;
x0 = [xa0;xb0;xc0];
tic; %measure time it takes to solve RK3 for T=1000
[res1000,tim1000, eig1000] = rk3_c(h, x0, T); 
tend = toc;

plot(tim1000, res1000)
legend("XA", "XB", "XC");
title("Solution using RK3 for T=1000")
xlabel("time");
ylabel("concentration");


%% part 3d)

[tout,yout] = impeuler(1000, x0, 1); % uses explicit euler to solve the Robertson problem
plot(tout,yout)
title("Implicit Euler solution to Robertson problem")
legend("xa(t)", "xb(t)", "xc(t)")
xlabel("time")
ylabel("concentration")




%% part 3e)

% for this calculation, the impeuler was temporarily modified so it didn't 
% save any values for plotting. it only returns the final values.
xa0=1;
xb0=0;
xc0=0;
x0 = [xa0;xb0;xc0];

xA_true = 0.293414227164;
xB_true = 0.000001716342048;
xC_true = 0.706584056494;
x_true = [xA_true; xB_true; xC_true];

tic;
u = impeuler(1000, x0, 1); 
imp_end = toc;

tic;
rk3_1000 = rk3_result(0.0001, x0, 1000);
rk3_end = toc;
rk3_norm = norm(rk3_1000 - x_true);
imp_norm = norm(u - x_true);



%% Function definitions part 1
function dm = dm(m) % Calculating the given Differential of m
    alpha = 0.07; % Damping parameter
    a1 = (1/4)*[1 sqrt(11) 2]'; % Fixed vector
    am = cross(a1, m);
    dm = am + alpha*cross(a1,am); %Defined differential of m
end

function [res, tim] = rk3_a(h, u0, T) %The defined Runge-Kutta 3 method for problem 1
    
    t0 = 0; %Init time
    t = t0;    
    l = length(u0);
    result = [u0 zeros(l, ceil(T/h))];
    time = [t0 zeros(1, ceil(T/h))];
    u = u0;  
    i = 2;    
    while t < T
        
        k1 = dm(u);
        k2 = dm(u + h*k1);
        k3 = dm(u + h*k1/4 + h*k2/4);
        t = t + h;
        u = u + h/6*(k1 + k2 + 4*k3);
        result(:, i) = u;
        time(i) = t;
        i = i + 1;
    end
    res = result;
    tim = time;
end



function s = stab_reg(z)
% function returns negative number if z is in the stability region for RK3,
% otherwise returns a positive number
    s = abs(1 + z + z^2 / 2 + z^3 / 6) - 1;
end



%% Function definitions part 2

% Explicit Euler function
function last_val = expl_Eul(u0, du0, step, T_end)
    h = step;
    mu = 1/82.45; %  Ratio of moon/Earth mass
    c0 = [-mu, 0].'; % Earth centre
    c1 = [1 - mu, 0].'; %Moon center
    B = [0, 1; -1, 0];
    dY = @(Y) [Y(3:4); -(1 - mu)*(Y(1:2) - c0)/(norm(Y(1:2) - c0))^3 - mu*(Y(1:2)...
        - c1)/(norm(Y(1:2) - c1))^3 + 2*B*Y(3:4) + Y(1:2)];
    u = [u0; du0];
    t = 0;
    while t < T_end
        u = u + h*dY(u);
        t = t + h;
    end
    last_val = u(1:2);
end

% The Adam-Bash function takes 3 arguments, a matrix of the first 4 position 
% and velocity vectors, the steplength and the final time
function [sat_pos, sat_vel, tim_vec] = Adam_Bash(u, steplength, T_end)
    mu = 1/82.45; %  Ratio of moon/Earth mass
    c0 = [-mu, 0].'; % Earth centre
    c1 = [1 - mu, 0].'; %Moon center
    B = [0, 1; -1, 0];
    dY = @(Y) [Y(3:4); -(1 - mu)*(Y(1:2) - c0)/(norm(Y(1:2) - c0))^3 - mu*(Y(1:2)...
        - c1)/(norm(Y(1:2) - c1))^3 + 2*B*Y(3:4) + Y(1:2)];
    
    h = steplength;
    % Satellite positions are saved in a matrix, we ceil the T_end/h since
    % it might not be an integer every time
    sat_pos = [u(1:2,:) zeros(2, ceil(T_end/h) - 4)];
    % Satellite velocities are saved in a matrix
    sat_vel = [u(3:4,:), zeros(2, ceil(T_end/h) - 4)];
    t0 = 0; % Initial time
    t = t0; % Setting initial time
    tim_vec = [t0, zeros(1, length(sat_pos)-1)];
    % Defining our step functions, which lets us take steps of the velocity
    % in the direction of the acceleration and the position in the
    % direction of the velocity
    f0 = dY(u(:,1));
    f1 = dY(u(:,2));
    f2 = dY(u(:,3));
    f3 = dY(u(:,4));
    i = 5;
    % u consists of 4 components, x, y positions and x', y' velocities in a
    % given time
    u = u(:,4); % set u as the fourth initial value
    while t < T_end
        % The Adam-Bashforth step
        u = u + (h/24)*(55*f3 - 59*f2 + 37*f1 - 9*f0);
        % Assign the new values to the step functions
        [f0, f1, f2, f3] = deal(f1, f2, f3, dY(u));
        t = t + h;
        % Don't want to save values after T_end, since four first values
        % already are saved
        if i <= length(sat_pos)
            sat_pos(:, i) = u(1:2);
            sat_vel(:, i) = u(3:4);
            tim_vec(i) = t;
        end
        i = i + 1;
    end
end


% The Adam-Bash function takes 3 arguments, a matrix of the first 4
% position and velocity vectors, the steplength and end time.
% Returns last value, without saving values in each step
function last_val= Adam_Bash_b(u, steplength, T_end)
    mu = 1/82.45; %  Ratio of moon/Earth mass
    c0 = [-mu, 0].'; % Earth centre
    c1 = [1 - mu, 0].'; %Moon center
    B = [0, 1; -1, 0];
    dY = @(Y) [Y(3:4); -(1 - mu)*(Y(1:2) - c0)/(norm(Y(1:2) - c0))^3 - mu*(Y(1:2)...
        - c1)/(norm(Y(1:2) - c1))^3 + 2*B*Y(3:4) + Y(1:2)];
    
    h = steplength;
    % Satellite positions are saved in a matrix, we ceil the T_end/h since
    % it might not be an integer every time
    t0 = 0; % Initial time
    t = t0; % Setting initial time
    % Defining our step functions, which lets us take steps of the velocity
    % in the direction of the acceleration and the position in the
    % direction of the velocity
    f0 = dY(u(:,1));
    f1 = dY(u(:,2));
    f2 = dY(u(:,3));
    f3 = dY(u(:,4));
    % u consists of 4 components, x, y positions and x', y' velocities in a
    % given time
    u = u(:,4); % set u as the fourth initial value
    while t < T_end
        % The Adam-Bashforth step
        u = u + (h/24)*(55*f3 - 59*f2 + 37*f1 - 9*f0);
        % Assign the new values to the step functions
        [f0, f1, f2, f3] = deal(f1, f2, f3, dY(u));
        t = t + h;
    end
    last_val = u;
end




function u = rk3_b(h, u0, du0, T) 
% This version of rk3 is created to solve questions in part 2 of the lab.
% Values are not saved in each time step. Returns only the final value.
    t0 = 0; %Init time
    t = t0;    
    u = [u0; du0];
    mu = 1/82.45; %  Ratio of moon/Earth mass
    c0 = [-mu, 0].'; % Earth centre
    c1 = [1 - mu, 0].'; %Moon center
    B = [0, 1; -1, 0];
    % Here we create an anonymous function of Y = [r(t); r'(t)] and spits
    % out the direction vectors at both points, i.e. Y' = [r'(t); r''(t)];
    dY = @(Y) [Y(3:4); -(1 - mu)*(Y(1:2) - c0)/(norm(Y(1:2) - c0))^3 ...
        - mu*(Y(1:2) - c1)/(norm(Y(1:2) - c1))^3 + 2*B*Y(3:4) + Y(1:2)];
      
    while t < T 
        k1 = dY(u);
        u1 = u + h*k1;
        k2 = dY(u1);
        u2 = u + h*k1/4 + h*k2/4;
        k3 = dY(u2);
        t = t + h;
        dir = (k1 + k2 + 4*k3)/6; 
        u = u + h*dir;
        
    end
end

%% Additional functions for part 3




function dx = der(x)
% takes a specific 3*1 vector as argument and returns its derivative
    r1 = 5*10^(-2); %rate constants
    r2 = 1.2*10^(4);
    r3=4*10^(7);
    xa=x(1);
    xb=x(2);
    xc=x(3);
    dxa = -r1*xa + r2*xb*xc;
    dxb = r1*xa -r2*xb*xc -r3*xb*xb;
    dxc = r3*xb*xb;
    dx = [dxa; dxb; dxc];
end


function j = dfdu(xb, xc)
% returns the jacobian of the right hand side in part 3 of the lab
    r1 = 5*10^(-2); % rate constants
    r2 = 1.2*10^(4);
    r3=4*10^(7);
    dja = [-r1, r2*xc, r2*xb];
    djb = [r1, -r2*xc-2*r3*xb, -r2*xb];
    djc = [0, 2*r3*xb, 0];
    j = [dja; djb; djc];
end



function [res, tim, eigs] = rk3_c(h, u0, T_end) %The defined Runge-Kutta 3 method for problem 3
% this version of the Runge-Kutta3 method returns the result and time
% vectors, as well as a vector of the eigenvalues in each timestep.
    T = T_end; %Time period
    t0 = 0; %Init time
    t = t0;    
    l = length(u0);
    result = [u0 zeros(l, round(T/h))];
    time = [t0 zeros(1, round(T/h))];
    jac0 = dfdu(u0(2), u0(3)); % jacobian in initial step
    eig0 = eig(jac0); % eigenvalues in initial step
    eigvals = [eig0, zeros(3,round(T/h))]; % matrix to save eigenvalues in each step
    u = u0;  
    i = 2;    
    while t < T
        k1 = der(u);
        k2 = der(u + h*k1);
        k3 = der(u + h*k1/4 + h*k2/4);
        t = t + h;
        u = u + h/6*(k1 + k2 + 4*k3);
        result(:, i) = u;
        time(i) = t;
        jacn = dfdu(u(2),u(3)); % the jacobian of f
        eign = eig(jacn);
        eigvals(:,i) = eign;
        i = i + 1;
    end
    res = result;
    tim = time;
    eigs = eigvals;
end



function rk3_final = rk3_result(h, u0, T_end) %The modified Runge-Kutta 3 method for problem 3.
% this version of the Runge-Kutta3 method returns only the final value of
% the vector at time t=T_end
    T = T_end; %Time period
    t0 = 0; %Init time
    t = t0;    
    u = u0;     
    while t < T
        k1 = der(u);
        k2 = der(u + h*k1);
        k3 = der(u + h*k1/4 + h*k2/4);
        t = t + h;
        u = u + h/6*(k1 + k2 + 4*k3);
    end
   rk3_final = u;
end