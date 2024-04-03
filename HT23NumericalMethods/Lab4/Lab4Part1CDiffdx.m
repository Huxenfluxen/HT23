%% File for different N values

% produces multiple plots with different CFL numbers
Nvec = [10,50, 100, 500, 1000];
tEnd = 4;
t=tEnd;
tau = 2;
D = 10;
a = 2;
lambda = 1/2.2; %CFL number
% boundary conditions
xBoundarySin = @(t) sin(2*pi*t'./tau);
xBoundarySq = @(t) sign(sin(2*pi*t'./tau));
uExactSin = @(x) sin(2*pi/tau*(t-x./a));
uExactSinVec = uExactSin(0:0.1:D);
uExactSqVec = sign(uExactSinVec);
% upwind method, different CFL
%plot(xIntv,uExactSinVec); %exact solution
figure
subplot(2,4,1)
for N = Nvec
    dx = D/N; 
    dt = lambda*dx;
    xIntv = 0:dx:D;
    tIntv = 0:dt:tEnd;
    uUpwindSin = upwind(dx,dt,tEnd,D,a,xBoundarySin);
    uUpwindSin4 = uUpwindSin(:, end);
    plot(xIntv, uUpwindSin4, LineWidth=2 ,LineStyle="-"); %upwind solution
    hold on
end
plot(0:0.1:D, uExactSinVec);
hold off
title("Upwind, BC g_{sin}(t)")
subplot(2,4,5)
for N = Nvec
    dx = D/N; 
    dt = lambda*dx;
    xIntv = 0:dx:D;
    tIntv = 0:dt:tEnd;
    uUpwindSq = upwind(dx,dt,tEnd,D,a,xBoundarySq);
    uUpwindSq4 = uUpwindSq(:, end);
    plot(xIntv, uUpwindSq4, LineWidth=2 ,LineStyle="-"); %upwind solution
    hold on
end
plot(0:0.1:D, uExactSqVec);
hold off
title("Upwind, BC g_{sq}(t)")
% LF method, different CFL
subplot(2,4,2)
for N = Nvec
    dx = D/N; 
    dt = lambda*dx;
    xIntv = 0:dx:D;
    tIntv = 0:dt:tEnd;
    uLaxFriSin = LF(dx,dt,tEnd,D,a,xBoundarySin);
    uLaxFriSin4 = uLaxFriSin(:, end);
    plot(xIntv, uLaxFriSin4, LineWidth=2 ,LineStyle="-"); %LF solution
    hold on
end
plot(0:0.1:D, uExactSinVec);
hold off
title("LF, BC g_{sin}(t)")
subplot(2,4,6)
for N = Nvec
    dx = D/N; 
    dt = lambda*dx;
    xIntv = 0:dx:D;
    tIntv = 0:dt:tEnd;
    uLaxFriSq = LF(dx,dt,tEnd,D,a,xBoundarySq);
    uLaxFriSq4 = uLaxFriSq(:, end);
    plot(xIntv, uLaxFriSq4, LineWidth=2 ,LineStyle="-"); %LF solution
    hold on
end
plot(0:0.1:D, uExactSqVec);
hold off
title("LF, BC g_{sq}(t)")
% LW method, different CFL
subplot(2,4,3)
for N = Nvec
    dx = D/N; 
    dt = lambda*dx;
    xIntv = 0:dx:D;
    tIntv = 0:dt:tEnd;
    uLaxWenSin = LW(dx,dt,tEnd,D,a,xBoundarySin);
    uLaxWenSin4 = uLaxWenSin(:, end);
    plot(xIntv, uLaxWenSin4, LineWidth=2 ,LineStyle="-"); %LF solution
    hold on
end
plot(0:0.1:D, uExactSinVec);
hold off
title("LW, BC g_{sin}(t)")
subplot(2,4,7)
for N = Nvec
    dx = D/N; 
    dt = lambda*dx;
    xIntv = 0:dx:D;
    tIntv = 0:dt:tEnd;
    uLaxWenSq = LW(dx,dt,tEnd,D,a,xBoundarySq);
    uLaxWenSq4 = uLaxWenSq(:, end);
    plot(xIntv, uLaxWenSq4, LineWidth=2 ,LineStyle="-"); %LF solution
    hold on
end
plot(0:0.1:D, uExactSqVec);
hold off
title("LW, BC g_{sq}(t)")
hSub = subplot(2,4,8);
plot(1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 'r');
set(hSub, 'Visible', 'off');
legend(hSub, "N = " + Nvec, 'Location', 'east');
% Lax-Wendroff approximation
function umat = LW(dx,dt,tEnd, D, a, boundaryFunc)
    lambda = dt/dx;
    
    tInt = 0:dt:tEnd-dt;
    xBCvec =  boundaryFunc(tInt);  %BC at x=0
    
    xInt = 0:dx:D;
    tBCvec = 0*xInt; %BC at t=0
    umat = zeros(length(tBCvec), length(xBCvec));
    umat(:,1) = tBCvec;
    umat(1,:) = xBCvec;
    for n = 1:length(xBCvec)
        for j = 2:length(tBCvec)-1
            umat(j,n+1) = umat(j,n)- a*lambda/2*(umat(j+1,n)-umat(j-1,n)) + a^2*lambda^2/2*(umat(j+1,n)-2*umat(j,n)+umat(j-1,n));
        end
        umat(end, n) = 2*umat(end-1, n)-umat(end-2, n); % at boundary x = D
    end
end
% Lax-Friedrichs approximation
function umat = LF(dx,dt,tEnd, D, a, boundaryFunc)
    lambda = dt/dx;
    
    tInt = 0:dt:tEnd-dt;
    xBCvec =  boundaryFunc(tInt);  %BC at x=0
    
    xInt = 0:dx:D;
    tBCvec = 0*xInt; %BC at t=0
    umat = zeros(length(tBCvec), length(xBCvec));
    umat(:,1) = tBCvec;
    umat(1,:) = xBCvec;
    for n = 1:length(xBCvec)
        for j = 2:length(tBCvec)-1
            umat(j,n+1) = 1/2*(umat(j+1,n) + umat(j-1,n)) - a*lambda/2*(umat(j+1,n) - umat(j-1,n));
        end
        umat(end, n) = 2*umat(end-1, n)-umat(end-2, n); % at boundary x = D
    end
end
% upwind approximation
function umat = upwind(dx,dt,tEnd, D, a, boundaryFunc)
    lambda = dt/dx;
    
    tInt = 0:dt:tEnd-dt;
    xBCvec =  boundaryFunc(tInt);  %BC at x=0
    
    xInt = 0:dx:D;
    tBCvec = 0*xInt; %BC at t=0
    umat = zeros(length(tBCvec), length(xBCvec));
    umat(:,1) = tBCvec;
    umat(1,:) = xBCvec;
    for j = 2:length(tBCvec)
        for n = 2:length(xBCvec)
            umat(j, n+1) = umat(j,n) - a*lambda*(umat(j,n)-umat(j-1,n));
        end
    end
end