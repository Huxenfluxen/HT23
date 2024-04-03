%% Lab 4

tEnd = 4;
tau = 2;
a = 2;
D = 10;
N = 100;
dx = D/N;
xIntv = 0:dx:D;
% We choose a CFL less than the stability limit 1/2
lambda = 1/4;
dt = lambda*dx;
tIntv = 0:dt:tEnd;


%Solution for the g_sin and g_sq condition, respectively
uExactSin = @(x) sin(2*pi/tau*(tEnd-x./a));
uExactSquare = @(x) sign(uExactSin(x));
  
% Upwind method
M = tEnd/dt;

%Functions defined on boundary
xBoundarySin = @(t) sin(2*pi*t./tau);
xBoundarySquare = @(t) sign(xBoundarySin(t));

% Upwind g_sin condition
upmat = upwind(dx, dt, tEnd, D, a, xBoundarySin);
upmat = upmat';
subplot(2,3,1)
plot(xIntv, uExactSin(xIntv), "LineWidth", 1.5)
hold on
plot(xIntv, upmat(end, :), "r", "LineWidth", 1.5)
title("Upwind g_{sin}")
axis([0,10,-2, 2])
hold off

% Upwind g_square condition
subplot(2,3,4)
plot(xIntv, uExactSquare(xIntv), "LineWidth", 1.5)
hold on
upmat2 = upwind(dx, dt, tEnd, D, a, xBoundarySquare);
upmat2 = upmat2';
plot(xIntv, upmat2(end, :), "r", "LineWidth", 1.5)
title("Upwind g_{sq}")
axis([0,10,-2, 2])
hold off


%Lax Friedrichs with sine function
LFmat = LF(dx, dt, tEnd, D, a, xBoundarySin);
LFmat = LFmat';

subplot(2,3,2)
plot(xIntv, uExactSin(xIntv), "LineWidth", 1.5)
hold on
plot(xIntv, LFmat(end, :), "r", "LineWidth", 1.5)
title("Lax Friedrich g_{sin}")
axis([0,10,-2, 2])
hold off

%Lax Friderich with square function
LFmat2 = LF(dx, dt, tEnd, D, a, xBoundarySquare);
LFmat2 = LFmat2';

subplot(2,3,5)
plot(xIntv, uExactSquare(xIntv), "LineWidth", 1.5)
hold on
plot(xIntv, LFmat2(end, :), "r", "LineWidth", 1.5)
title("Lax Friedrich g_{sq}")
axis([0,10,-2, 2])
hold off

%Lax-Wendroff with sine function
LWmat = LW(dx, dt, tEnd, D, a, xBoundarySin);
LWmat = LWmat';

subplot(2,3,3)
plot(xIntv, uExactSin(xIntv), "LineWidth", 1.5)
hold on
plot(xIntv, LWmat(end, :), "r", "LineWidth", 1.5)
title("Lax Wendroff g_{sin}")
axis([0,10,-2, 2])
hold off

%Lax Wendroff with square function
LWmat2 = LW(dx, dt, tEnd, D, a, xBoundarySquare);
LWmat2 = LWmat2';

subplot(2,3,6)
plot(xIntv, uExactSquare(xIntv), "LineWidth", 1.5)
hold on
plot(xIntv, LWmat2(end, :), "r", "LineWidth", 1.5)
title("Lax Wendroff g_{sq}")
axis([0,10,-2, 2])
hold off


%% Functions Part 1

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