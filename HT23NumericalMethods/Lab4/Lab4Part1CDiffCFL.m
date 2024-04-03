%% Lab4 with 
tEnd = 4;
tau = 2;
a = 2;
D = 10;
N = 100;
dx = D/N; %Fixed
xIntv = 0:dx:D;
% Different CFL numbers less than 1/2
lambda = [1/1.99, 1/1.999, 1/2.0001, 1/2.001, 1/2.01];
dt = lambda.*dx;

%Solutione for the g_sin and g_sq condition
uExactSin = @(x) sin(2*pi/tau*(tEnd-x./a));
uExactSquare = @(x) sign(uExactSin(x));
  
% Upwind method
M = tEnd./dt;

%Functions defined on boundary
xBoundarySin = @(t) sin(2*pi*t./tau);
xBoundarySquare = @(t) sign(xBoundarySin(t));

%Functions defined on boundary
xBoundarySin = @(t) sin(2*pi*t./tau);
xBoundarySquare = @(t) sign(xBoundarySin(t));

%Upwind with sin function
plottar = figure;
subplot(2,4,1)
plot(xIntv, uExactSin(xIntv), "LineWidth", 1.5)
hold on
for deltat = dt
    umat = upwind(dx, deltat, tEnd, D, a, xBoundarySin);
    umat = umat';

    plot(xIntv, umat(end, :), "LineWidth", 1.5)
end
    title("Upwind g_{sin}")
    %legend("CFL = " + dt./dx)
    axis([0,10,-2, 2])
    hold off

%Upwind with square function
subplot(2,4,5)
plot(xIntv, uExactSquare(xIntv), "LineWidth", 1.5)
hold on
for deltat = dt
    umat = upwind(dx, deltat, tEnd, D, a, xBoundarySquare);
    umat = umat';

    plot(xIntv, umat(end, :), "LineWidth", 1.5)
end
    title("Upwind g_{sq}")
    %legend("CFL = " + dt./dx)
    axis([0,10,-2, 2])
    hold off


%Lax Friedrich with sin function
subplot(2,4,2)
plot(xIntv, uExactSin(xIntv), "LineWidth", 1.5)
hold on
for deltat = dt
    umat = LF(dx, deltat, tEnd, D, a, xBoundarySin);
    umat = umat';

    plot(xIntv, umat(end, :), "LineWidth", 1.5)
end
    title("Lax Friedrich g_{sin}")
    %legend("CFL = " + dt./dx)
    axis([0,10,-2, 2])
    hold off

%Lax Friedrich with square function
subplot(2,4,6)
plot(xIntv, uExactSquare(xIntv), "LineWidth", 1.5)
hold on
for deltat = dt
    umat = LF(dx, deltat, tEnd, D, a, xBoundarySquare);
    umat = umat';

    plot(xIntv, umat(end, :), "LineWidth", 1.5)
end
    title("Lax Friedrich g_{sq}")
    %legend("CFL = " + dt./dx)
    axis([0,10,-2, 2])
    hold off


%Lax Wendroff with sine function
subplot(2,4,3)
plot(xIntv, uExactSin(xIntv), "LineWidth", 1.5)
hold on
for deltat = dt
    umat = LW(dx, deltat, tEnd, D, a, xBoundarySin);
    umat = umat';

    plot(xIntv, umat(end, :), "LineWidth", 1.5)
end
    title("Lax Wendroff g_{sin}")
    %legend("CFL = " + dt./dx)
    axis([0,10,-2, 2])
    hold off


%Lax Wendroff with square function
subplot(2,4,7)
plot(xIntv, uExactSquare(xIntv), "LineWidth", 1.5);
hold on
i = 2;
for deltat = dt
    umat = LW(dx, deltat, tEnd, D, a, xBoundarySquare);
    umat = umat';

    plot(xIntv, umat(end, :), "LineWidth", 1.5);
    i = i + 1;
end
    title("Lax Wendroff g_{sq}")
    %legend("CFL = " + dt./dx)
    axis([0,10,-2, 2])
    hold off

hSub = subplot(2,4,4);
plot(1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 'r');
set(hSub, 'Visible', 'off');
legend(hSub, "Exact solution", "CFL = " + dt./dx, 'Location', 'east');


%% Functions
% Upwind method
function umat = upwind(dx, dt, tEnd, D, a, boundaryFunc)
    lambda = dt/dx;    
    tIntv = 0:dt:tEnd-dt;
    xBCvec = boundaryFunc(tIntv); % Boundary condition vector for x = 0
    xIntv = 0:dx:D;
    tBCvec = 0 * xIntv;
    
    umat = zeros(length(tBCvec), length(xBCvec));
    umat(:, 1) = tBCvec';
    umat(1, :) = xBCvec;
    % disp(umat(:, 1))
    
    if a > 0
        for j = 2:length(tBCvec) 
            for n = 2:length(xBCvec)
                umat(j, n+1) = umat(j, n) - a*lambda*(umat(j, n) - umat(j-1, n));
            end
        end            
    end
    %umat = umat(2:end, :);
end

% Function for the Lax-Fridriechs methods

function umat = LF(dx, dt, tEnd, D, a, boundaryFunc)
    lambda = dt/dx;
    tIntv = 0:dt:tEnd-dt;
    xBCvec = boundaryFunc(tIntv); % Boundary condition vector for x = 0
    xIntv = 0:dx:D;
    tBCvec = 0 * xIntv;
    
    umat = zeros(length(tBCvec), length(xBCvec));
    umat(:, 1) = tBCvec';
    umat(1, :) = xBCvec;

    for n = 1:length(xBCvec) 
        for j = 2:length(tBCvec)-1
            umat(j, n+1) = 1/2*(umat(j+1, n) + umat(j-1, n)) - a*lambda/2*(umat(j+1, n) - umat(j-1, n));
        end
    end
    %Set the boundary for x = D
    umat(end, :) = 2*umat(end-1, :) - umat(end-2, :);

end

% Function for the Lax-Wendroff methods

function umat = LW(dx, dt, tEnd, D, a, boundaryFunc)
    lambda = dt/dx;
    tIntv = 0:dt:tEnd-dt;
    xBCvec = boundaryFunc(tIntv); % Boundary condition vector for x = 0
    xIntv = 0:dx:D;
    tBCvec = 0 * xIntv;

    umat = zeros(length(tBCvec), length(xBCvec));
    umat(:, 1) = tBCvec';
    umat(1, :) = xBCvec;

    for n = 1:length(xBCvec) 
        for j = 2:length(tBCvec)-1
            umat(j, n+1) = umat(j, n) - a*lambda/2*(umat(j+1, n) - umat(j-1, n)) + a^2*lambda^2/2*(umat(j+1, n) - 2*umat(j, n) + umat(j-1, n));
        end
    end
    %Set the boundary for x = D
    umat(end, :) = 2*umat(end-1, :) - umat(end-2, :);
end