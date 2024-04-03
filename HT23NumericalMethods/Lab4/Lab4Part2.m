%% Lab 2 Part b (Minor changes needed in upwind and lax-wendroff methods!)

%
Fr = 0.35;
alpha = 1/Fr;

lambda = 1.0/(1*(alpha + 1));
dx = 0.002;
dt = (lambda)*dx;
tEnd = 0.15;
tIntv = 0:dt:tEnd;
L0 = -0.4;
L1 = 0.7;
A = [1, alpha; alpha, 1];


xIntv = L0:dx:L1;

%Defining sources and the system of equations
s = @(x) sin(20*pi.*x).*(abs(x) < 1/20);
r = @(t) (sin(40*pi.*t + pi/6) > 0.5);
f = @(x, t) s(x).*r(t);

% tTest = linspace(0, 0.15, 1000);
% xTest = linspace(-0.4, 0.7, 1000);
% [X, T] = meshgrid(xTest, tTest);
% figure
% mesh(xTest, tTest, f(X, T))

% Plotting perturbations and velocities of upwind and Lax Wendroff
% umat1 = upwind(dx, dt, tEnd, L0, L1, Ap, Am, f);
% umat2 = LW(dx, dt, tEnd, L0, L1, A, f);

[upwindMat, laxMat] = UpwindLax(dx, dt, tEnd, L0, L1, A, f);

figure
% Perturbations for both upwind and Lax-Wendroff
subplot(2,2,1)
plot(xIntv, upwindMat(1, :, end), "LineWidth", 1.5)
hold on
plot(xIntv, laxMat(1, :, end), "LineWidth", 1.5)
hold off
title("Perturbations")
legend("Upwind", "Lax-Wendroff")
xlabel("Distance")
ylabel("Amplitude")
% Velocitites for both upwind and Lax-Wendroff
subplot(2,2,2)
plot(xIntv, upwindMat(2, :, end), "LineWidth", 1.5)
hold on
plot(xIntv, laxMat(2, :, end), "LineWidth", 1.5)       
hold off
title("Velocities")
legend("Upwind", "Lax-Wendroff")
xlabel("Distance")
ylabel("Amplitude")

% 3D plots of the upwind and Lax-Wendroff solutions, respectively
subplot(2,2,3)
mesh(tIntv, xIntv, squeeze(upwindMat(1, :, :)))
title("Upwind, perturbation")
ylabel("Distance")
xlabel("Time")
zlabel("Amplitude")

subplot(2,2,4)
mesh(tIntv, xIntv, squeeze(laxMat(1, :, :)))
title("Lax-Wendroff, perturbation")
ylabel("Distance")
xlabel("Time")
zlabel("Amplitude")

%% Functions Part 2


% Function using Lax-Wendroff and Upwind method and outputs two matrices
% corresponding to Upwind and Lax-Wendroff, respectively.
% Need to change a few things in the methods!!! Wrong output!
function [upwindMat, laxMat] = UpwindLax(dx, dt, tEnd, L0, L1, A, sourceTerm)
    % Splitting up A for the Upwind method
    [V, Lambda] = eig(A);
    LambdaP = Lambda.*(Lambda > 0);
    LambdaM = Lambda.*(Lambda < 0);
    Ap = V*LambdaP*inv(V);
    Am = V*LambdaM*inv(V);
    lambda = dt/dx;
    tIntv = 0:dt:tEnd - dt;
    xBCvec = sourceTerm(0, tIntv); % Boundary condition vector for x = 0
    xBCvec = [0*xBCvec; xBCvec];
    xIntvPlus = 0:dx:L1;
    tBCvecPlus = [0; 0]*xIntvPlus;
    xIntvMinus = 0:-dx:L0;
    tBCvecMinus = [0; 0]*xIntvMinus;
    
    laxMatPlus = zeros(2, length(xIntvPlus), length(tIntv));
    laxMatMinus = zeros(2, length(xIntvMinus), length(tIntv));
    laxMatPlus(:, :, 1) = tBCvecPlus;
    laxMatMinus(:, :, 1) = tBCvecMinus;
    laxMatPlus(:, 1, :) = xBCvec;
    laxMatMinus(:, 1, :) = xBCvec;
    upwindMatPlus = laxMatPlus;
    upwindMatMinus = laxMatMinus;

    for j = 2:length(tBCvecPlus)-1
        for n = 2:length(xBCvec)  
            FLax = [0; (sourceTerm(dx*j, dt*n) + sourceTerm(dx*j, dt*(n+1)))/2] -...
                lambda/4*A*[0; (sourceTerm(dx*(j+1), dt*n) - sourceTerm(dx*(j-1), dt*n))];
            FUp = [0; sourceTerm(dx*j, dt*n)]; 
            
            laxMatPlus(:, j, n+1)  = laxMatPlus(:, j, n) - A*lambda/2*...
                (laxMatPlus(:, j+1, n) - laxMatPlus(:, j-1, n)) + A^2*lambda^2/2*...
                (laxMatPlus(:, j+1, n) - 2*laxMatPlus(:, j, n) + laxMatPlus(:, j-1, n)) + dt*FLax;
                      
            upwindMatPlus(:, j, n+1)  = upwindMatPlus(:, j, n) - Ap*lambda*...
                (upwindMatPlus(:, j, n) - upwindMatPlus(:, j-1, n)) -...
                Am*lambda*(upwindMatPlus(:, j+1, n) - upwindMatPlus(:, j, n)) + dt*FUp;
        end
    end
    % Loop for the negative x-interval
    for j = 2:length(tBCvecMinus)-1
        for n = 2:length(xBCvec)
            FLax = [0; (sourceTerm(-dx*j, dt*n) + sourceTerm(-dx*j, dt*(n+1)))/2] -...
                lambda/4*A*[0; (sourceTerm(-dx*(j+1), dt*n) - sourceTerm(-dx*(j-1), dt*n))];
            FUp = [0; sourceTerm(-dx*j, dt*n)];
           
            laxMatMinus(:, j, n+1)  = laxMatMinus(:, j, n) - A*lambda/2*...
                (laxMatMinus(:, j+1, n) - laxMatMinus(:, j-1, n)) + A^2*lambda^2/2*...
                (laxMatMinus(:, j+1, n) - 2*laxMatMinus(:, j, n) + laxMatMinus(:, j-1, n)) + dt*FLax;

            upwindMatMinus(:, j, n+1) = upwindMatMinus(:, j, n) - Ap*lambda*(upwindMatMinus(:, j, n) -...
                upwindMatMinus(:, j-1, n)) - Am*lambda*(upwindMatMinus(:, j+1, n) - upwindMatMinus(:, j, n)) + dt*FUp;
        end
    end   
% Reversing the order since we need the matrices from L0:0
    LaxMinus = zeros(size(laxMatMinus));
    upwindMinus = zeros(size(upwindMatMinus));
    for j = 2:length(tBCvecMinus)-1
        LaxMinus(:, j, :) = laxMatMinus(:, end-j+1, :);
        upwindMinus(:, j, :) = upwindMatMinus(:, end-j+1, :);
    end
% Extrapolation for the boundary 
    laxMat = cat(2, LaxMinus(:, 1:end-1, :), laxMatPlus);
    laxMat(:, 1, :) = 2*laxMat(:, 2, :) - laxMat(:, 3, :);
    laxMat(:, end, :) = 2*laxMat(:, end-1, :) - laxMat(:, end-2, :);

    upwindMat = cat(2, upwindMinus(:, 1:end-1, :), upwindMatPlus);
    upwindMat(:, 1, :) = 2*upwindMat(:, 2, :) - upwindMat(:, 3, :);
    upwindMat(:, end, :) = 2*upwindMat(:, end-1, :) - upwindMat(:, end-2, :);

end

