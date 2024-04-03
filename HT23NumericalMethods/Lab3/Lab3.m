%% Lab3 Part 1
d = 0.21;
T_end = 5;
%Stable solution: h_x/h_T^2 < 1/(2d)

% h_x/h_T^2 <= d/2 stability condition
c = 1/(2*d);
stab_cond = @(dt, dx) dt/dx^2;

%Plot stable and unstable solution
M = 10000;
dt = T_end/M;
i = 1;
figure
for N = [84, 85]
    dx = 1/N;
    A = assemble_A(N + 1);
    u_final = expEul(A, dt, T_end, dx);
    subplot(1,2,i)
    mesh(0:dt:2, 0:dx:1,  u_final)
    title("c = " + stab_cond(dt, dx))
    i = i + 1;
    hold on
end
hold off
%% Plot the distribution at t = 1.1

d = 0.21;
T_end = 5;
%Plot stable and unstable solution
NStab = 84;
MStab = 10000;
dt = T_end/MStab;
dx = 1/NStab;
AStab = assemble_A(NStab + 1);
u_finalStab = expEul(AStab, dt, T_end, dx);
subplot(2,2,1)
mesh(0:dt:2, 0:dx:1,  u_finalStab)

NUnstab = 83;
MUnstab = 9620;
dtUnstab = T_end/MUnstab;
dxUnstab = 1/NUnstab;
AUnstab = assemble_A(NUnstab + 1);
u_finalUnstab = expEul(AUnstab, dtUnstab, T_end, dxUnstab);
subplot(2,2,2)
mesh(0:dtUnstab:2, 0:dxUnstab:1,  u_finalUnstab)

% To get the h_t of the plot
t_val = 1.5/2*MStab;
subplot(2,2,3)
grid on
plot(0:dx:1, u_finalStab(:, t_val + 1))

% Plot endpoints of rod

% To get the h_t of the plot
x_init = 1;
x_final = 1*NStab;
subplot(2,2,4)
plot(0:dt:2, u_finalStab(x_init, :))
hold on
plot(0:dt:2,u_finalStab(x_final + 1, :))
legend("x = 0", "x = 1")
hold off


%% Lab3 Part1 b using ode23

% A = assemble_A(N+1);
% [t, u] = ode23(@(t, u) u_ode23(t, u, A, N), [0 2], zeros(101, 1));
% mesh(t, 0:1/100:1, u')
N = [100, 200, 400]';
timesteps = zeros(size(N));
timesteps_s = zeros(size(N));
timesteps_sj = zeros(size(N));
comp_time = zeros(size(N));
comp_time_s = zeros(size(N));
comp_time_sj = zeros(size(N));


i = 1;
while i <= length(N)
    A = assemble_A(N(i)+1);
    tic
    [t, u] = ode23(@(t, u) u_ode23(t, u, A, N(i)+1), [0 2], zeros(N(i)+1, 1));
    tid = toc;

    tic
    [t_s, u_s] = ode23s(@(t, u) u_ode23(t, u, A, N(i)+1), [0 2], zeros(N(i)+1, 1));
    tids = toc;

    tic
    options = odeset("Jacobian", A);
    [t_sj, u_sj] = ode23s(@(t, u) u_ode23(t, u, A, N(i)+1), [0 2], zeros(N(i)+1, 1), options);
    tidsj = toc;

    comp_time(i) = tid;
    comp_time_s(i) = tids;
    comp_time_sj(i) = tidsj;

    timesteps(i) = length(t);
    timesteps_s(i) = length(t_s);
    timesteps_sj(i) = length(t_sj);
    i = i + 1;
end

table(N, timesteps,timesteps_s, timesteps_sj, comp_time, comp_time_s, comp_time_sj)



%% Lab3 Part2

% Code from Lab 2

N = 60; % 60, 240 are also used
Lx = 12;
Ly = 5;
h = Lx/N;
M = round(Ly/h); %round to make sure stepsize is equal in x and y directions
Text = 25;

func_d = @(X,Y) 100*exp((-1/2)*(X-4).^2-4*(Y-1).^2);
A = create_A(N,M,h);
f = create_f(N, M, Lx, Ly, func_d);
dt = 1/8;
t_end = 40;
u0 = Text*ones(size(A(:,1)));
tic
u_mat = u_trap(A, f, dt, u0, t_end);
toc
u_end = u_mat(:, end);
u_end_mat = reshape(u_end, N + 1, M); % Reshaping T to a matrix to use for mesh plot
test_val1 = u_end_mat(6/h+1, 2/h);

% Part 2 b): Plot the solution for fixed tau
tau_vec = [0,1,4,12,22,40];
i = 1;
for tau = tau_vec
    subplot(2,3,i)
    u_tau = u_mat(:, tau/dt + 1);
    u_tau = reshape(u_tau, N + 1, M);
    u_tau = u_tau';
    mesh(0:h:Lx, h:h:Ly, u_tau)
    title("\tau = " + tau)
    i = i + 1;
end

%Point (6,2) for all time values
uSixTwo = u_mat((2/h - 1)* (N + 1) + 6/h + 1, :);
figure
plot(0:dt:t_end, uSixTwo, "LineWidth", 1.5)
title("The point (6,2) for all time values to t = 40")
xlabel("time")
ylabel("Temperature")
hold on 
yline(47.2201, "LineWidth", 1.5)
legend("Temp at (6,2) over time", "Temp at (6,2) in Lab 2", "Location","southeast")
hold off

% figure
% for t = 0:dt:t_end
%     u_tau = u_mat(:, tau/dt + 1);
%     u_tau = reshape(u_tau, N + 1, M);
%     u_tau = u_tau';
%     mesh(0:h:Lx, h:h:Ly, u_tau)
%     drawnow
%     pause(0.01)
% end


%% LAB 2 functions
% This creates the matrix A
function A = create_A(N,M,h)
    oneVec = ones((N+1)*(M-1),1);
    twoVec = [ones((N+1)*(M-2),1); 2*ones(N+1,1)];
    
    s_maindiag = -4*eye(N+1);
    s_subdiag = diag([ones(N-1,1); 2], -1);
    s_supdiag = diag([2; ones(N-1,1)], 1);
    S = s_maindiag + s_subdiag + s_supdiag; % S is a N+1 * N+1 matrix
    one_M = eye(M);
    kron_matrix = kron(one_M, S); %kronecker product 
    A = 1/h^2*(kron_matrix + diag(oneVec, N+1) + diag(twoVec, -N-1));    
end

%takes the function that determines ans creates answer vector f
function f=create_f(N, M, Lx, Ly, func) 
    Text = 25;
    h=Lx/N;
    [X,Y] = discretize(N,M,Lx,Ly);
    f=func(X,Y);
    f = f';
    f(:,1) = f(:,1) + 1/h^2*Text; % subtract boundary condition from first column of f (N+1 elements)
    f = reshape(f, [], 1); % reshape f to a vector
end

% This function is used to discretize a mesh from 0 to Lx and 0 to Ly
function [X,Y] = discretize(N,M,Lx,Ly) % [X,Y] is the mesh grid of omega discretized
    h=Lx/N;
    xIntv = linspace(0, Lx, N+1); % discretize x from 0 to Lx, length N+1
    yIntv = linspace(h, Ly, M); %discretize y from h to Ly, length M
    [X,Y] = meshgrid(xIntv, yIntv);
end

%% Lab3 Part1 functions

%Trapezoidal function
function u_mat = u_trap(A, f, dt, u0, t_end)
    u = u0;
    id = eye(length(A));
    S = (id - 1/2*dt*A);
    [L,U] = lu(S); % Lower upper factorisation
    K = (id + 1/2*dt*A);
    u_mat = [u, zeros(length(u), t_end/dt)];
    i = 2;
    t = dt;
    while t <= t_end
        u = U\(L\(K*u + dt*f));

        u_mat(:, i) = u;
        i = i + 1;
        t = t + dt;
    end
end

%Assemble the matrix A
function A = assemble_A(N)
    d = 0.21;
    dx = 1/(N-1);
    A_maindiag = -2*eye(N);
    A_supdiag = diag(ones(N-1,1),1);
    A_subdiag = diag([ones(N-2,1); 2],-1);
    A = A_maindiag + A_subdiag + A_supdiag;
    A = d/dx^2*A;
end
% Explicit Euler function
function u_mat = expEul(A, dt, T_end, dx)
    a = 1.4;
    d = 0.21;
    t = dt;
    u = zeros(size(A(1,:)));
    u = u';
    b = zeros(size(A(1,:)));
    b = b';
    u_mat = [u, zeros(length(u), T_end/dt)];
    i = 2;
    while t <= T_end
        if t <= a
            b(1) = d/dx^2*sin(pi*t/a);     
        else
            b(1) = 0;
        end

        u = u + dt*(A*u + b);
        %u = u;
        u_mat(:, i) = u;
        t = t + dt;
        i = i +1;
    end
end

function du = u_ode23(t, u, A, N)
    dx = 1/(N-1);
    a = 1.4;
    d = 0.21;
    b = zeros(size(A(1,:)));
    if t <= a
            b(1) = d/dx^2*sin(pi*t/a);     
    else
            b(1) = 0;
    end
    du = A*u + b';
end