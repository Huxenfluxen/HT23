%% SF2520 Lab 2 Part 1


%% Part 1 a) 
% This code plots the temperature distribution of problem 1
NPlot = [10,20,40,80];
NVal = [80, 160, 320];
v = 1;
i = 1;
mid_val_vec = zeros(3,1);
for N = NPlot
    [Temp_plot, zIntv, mid_val] = tempDist(v,N);
    plot(zIntv, Temp_plot)
    hold on
end
hold off
title("Temperature of fluid in the cylinder between 0 and 1")
legend("N = " + NPlot)
xlabel("Horizontal position in cylinder")
ylabel("Temperature")

for N = NVal
    [Temp_val, zIntv, mid_val] = tempDist(v,N);
    mid_val_vec(i) = mid_val;
    i = i + 1;
end

%% Part 1b)
vVal = [1, 5, 15, 100];
NVal = [40, 40, 40, 200];
vDict = dictionary(vVal, NVal);

for v = vVal
    [Temp_plot, zIntv, mid_val] = tempDist(v,vDict(v));
    plot(zIntv, Temp_plot)
    hold on
end
hold off

title("Temperature of fluid in the cylinder between 0 and 1")
legend("v = " + vVal + ", N =" + NVal)
xlabel("Horizontal position in cylinder")
ylabel("Temperature")



%% LAB 2 part 2

N_a = 120; % 60, 240 are also used
Lx = 12;
Ly = 5;
h_a = Lx/N_a;
M_a = round(Ly/h_a); %round to make sure stepsize is equal in x and y directions
Text = 25;
func_a = @(X,Y) 2* ones(size(X)); % X and Y are of the same size
func_d = @(X,Y) 100*exp((-1/2)*(X-4).^2-4*(Y-1).^2);
A_a = create_A(N_a,M_a,h_a);
f_a = create_f(N_a, M_a, Lx, Ly, func_d);
T_a = A_a\f_a;
T_a = [Text*ones(N_a+1,1); T_a]; % Add the initial Dirichlet conditions
T_a = reshape(T_a, N_a +1, M_a +1); % Reshaping T to a matrix to use for mesh plot
T_a = T_a'; % Set T as the transpose for getting the values in the right position
test_val = T_a(2/h_a+1, 6/h_a+1);
mesh(0:h_a:Lx, 0:h_a:Ly, T_a)
title("Temperature of the metal block when N = 120 as a function of coordinates (x,y)")
xlabel("x-coordinates")
ylabel("y-coordinates")
zlabel("Temperature")
% Following plots are for func_d
figure
imagesc(T_a)
title("Distribution of the temperature of the metal block when N = 120, M = 50")
xlabel("x-coordinates")
ylabel("y-coordinates")
figure
contour(T_a)
title("Contour plot of the temperature of the metal block when N = 120, M = 50")
xlabel("x-coordinates")
ylabel("y-coordinates")

%% Functions Part 1

%This function takes the fluid velocity and N, and ouputs the temperature
%distribution, the discretized interval as well as the median temp
function [Temp_out, zIntv, mid_val] = tempDist(v, N)
        
    a = 0.1;
    b = 0.4;
    alpha0 = 50; % Heat transfer coefficient
    Tout = 25;
    T0 = 100;

    %N = 20; % number of discretized points  of the interval [0, 1]
    z0 = 0;
    zN = 1;
    h = 1/N;
    alpha = @(v) sqrt(v^2/4 + alpha0^2) - v/2; % Heat transfer function?
    d_nmin1 = 1;
    d_n = -2*h*alpha(v);
    d_nplus1 = -Tout*d_n;
    
    i = 1;
    zIntv = linspace(z0, zN, N + 1); % Interval between z0 and zN with n-1 inner points
    %A = zeros(N, N);
    f = zeros(N,1);
    % For j in [1, N]
    aj = -1/(h^2) - v/(2*h); 
    bj = 2/(h^2);
    cj = -1/(h^2) + v/(2*h);
    
    aVec = aj*ones(N-1,1);
    bVec = bj*ones(N,1);
    cVec = cj*ones(N-1,1);
    A = diag(aVec, -1) + diag(bVec) + diag(cVec, 1);
    
    for j = 1:N
        zj = j*h;
        Qj = driv(zj, a, b);
        f(j) = Qj;
    end
    
    f(1) = f(1) - T0*aj; % Position 1 of f
    f(N) = f(N) - cj*d_nplus1; % Assigns the last value of the f vector
    A(N, N-1) = aj + d_nmin1*cj; % d_(n-1) = 1
    A(N, N) = bj + cj*d_n;
    
    Temp_out = [T0; A\f]; % Adding T0 to the first position

    mid_val = Temp_out(round(length(zIntv)/2));

end


%Driving function taking a, b as the staring and end point
function Q = driv(z, a, b)
    Q0 = 7000;

    if (z >= 0 && z < a) || (z > b && z <= 1)
        Q = 0;
    else
        Q = Q0*sin((z - a)*pi/(b-a));
    end
end

%% functions Part 2
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
    A = -1/h^2*(kron_matrix + diag(oneVec, N+1) + diag(twoVec, -N-1));    
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