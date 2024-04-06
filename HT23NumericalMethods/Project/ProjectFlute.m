%% This section creates all the constants and vectors and matrices ought to be plotted. These have to be executed before the rest
% Initial constants
beta = 50;
delta = 0.01;
f0 = 523.25; %Fundamental frequency
c0 = 343; % Velocity of sound
L = c0/(2*f0); %Length of flute
w = 4*L/beta;
T = L/c0;
alpha = 3;
R = 3/100; %Radius of flute pipe

%Create a vector of positions of n number of holes of the same size
n =  7; % number of holes
S = 7*(1:n)/1000; %Diameter of fingerhole
s = L*alpha.*S./R^2; %Dimensionless finger hole
yHatVec = 0;
for i = 1:n
    yHatVec(i) = 1/2;%1/3 + i/3/n; % Location of the finger holes
end

%Time and space discretizations. CFL for this dimensionless problem
%is: dt/dx < 1
N = 500; %Spatial discretisation
M = 510; %Time discretisation
dt = 1/M; dy = 1/N; tEnd = 40;
yVec = linspace(0, 1, N);
% If using linspace for tVec don't forget to make sure dt < dx/2!
tVec = 0:dt:tEnd;

%Create D, A and u-matrix
D = createD(s, yVec, yHatVec, delta);
A = createA(length(yVec)-1, dy);
uMat = centralTime(A, D, yVec, tVec, beta);

%Find the dimensionless pressure over time on the point L/3, i.e. 1/3
strlk = size(uMat);
thirdIndx = round(strlk(1)/3) + 1; % +1 since index starts at 0

%% Section for mesh plot of uMat
% PLot the result
mesh(tVec, yVec, uMat)
title('Relative Pressure inside Recorder')
xlabel('Time (dimensionless)')
ylabel("Length of Flute (dimensionless)")
zlabel("Pressure (dimensionless)")


%% Section for plotting the pressure on u(tau, 1/3)
figure
plot(tVec, uMat(thirdIndx, :))
title("u(\tau, 1/3)")
xlabel('Time (dimensionless)')
ylabel("Pressure (dimensionless)")

%% Section for the fast Fourier transform of uMat
% Discrete (fast) Fourier Transform
fftRes = fft(uMat(thirdIndx, :));
tLen = length(tVec);
fVec = (0:tLen-1)/tLen/T/dt; % Some kind of Nyqvist division.......
fIndx = find(fVec <= 2000);
% [peak, peakLoc] = findpeaks(fVec(fIndx));
% fundFreq = peakLoc(1);

% Plot the amplitude of FFT vs frequency
figure
semilogy(fVec(fIndx), abs(fftRes(fIndx)));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of the Recorder Sound');

% hold on
% % Peaks in the frequency spectrum
% fs = length(uMat(1,:))/(tVec(end)*T);
% nr_peaks = [1 2 3];% Max nr of peaks
% for i = nr_peaks
%     resFreq = round(f0*i*length(fftRes)/fs + 1);
%     if fVec(resFreq) <= 2000% Dont want to plot above 2000 Hz
%         plot(fVec(resFreq), abs(fftRes(resFreq)),'*');
%     end
% end
% hold off

%% Section for computing eigenvalues and eigenvectors of -A + D
Afull = full(A); % Must take eigenvalues from non-sparse matrices
Dfull = full(D);
[eigVecs, ADeig] = eig(-Afull + Dfull);
tLen = length(tVec);
ADeig = sqrt(diag(ADeig))/2/pi; % Now it is a vector of fundamental frequencies and overtones
f0_2 = min(ADeig)/T; % The fundamental frequency f, first peak infrequency domain. Must divide by T since eigvals are normed
f0Ind = find(min(ADeig)); % Index of fundamental frequency
f0EigVec = eigVecs(:, f0Ind); % Not necessary to rescale as in fft since eigvecs are normed???
figure
plot(yVec(2:end-1), f0EigVec) % Plot w/o boundary points
% hold on
% plot(yVec, uMat(:,tVec(end)))
% hold off


%% Section of defined functions

%Central difference in time
function uMat = centralTime(A, D, yVec, tVec, beta)
    %Function g defined in project paper
    g = @(y) y.*exp(-beta.*y);
    
    uOld = zeros(length(yVec)-2, 1);
    u = uOld; % At t = 0
    dtau = diff(tVec(1:2));
    M = dtau^2*(A - D);
    gVec = dtau^2*g(yVec(2:end-1))';
    U = zeros(length(A), length(tVec));
    for n = 1:length(tVec)-1
        u = M*u + 2*u - uOld + gVec;
        U(:, n+1) = u;
        uOld = U(:, n);
    end
    %Create uMat with zeros 
    uMat = zeros(length(yVec), length(tVec)); 
    uMat(2:end-1, :) = U;
end

% Function for central difference in space
function A = createA(N, dx)
    Amain = -2*eye(N-1);
    Asub = diag(ones(N-2,1), -1);
    Asup = diag(ones(N-2,1), 1);
    Afinal = (1/dx^2)*(Amain + Asub + Asup);
    A = sparse(Afinal);
end

%Takes vector of dimensionless holes' sizes, vector of y values and vector of
% position of the holes. Outputs the diagonal matrix D of elements d(y_i).
function D = createD(sVec, yVec, yHatVec, delta)
    % The function h(y)
    h = @(y) (1/delta).*((1 + cos(2*pi.*y/delta)).*(abs(y)<=delta/2));
    d = 0;
    %Remove the points of Dirichlet condition
    yVec = yVec(2:end-1);
        for i = 1:length(sVec)
            s = sVec(i);
            yHat = yHatVec(i);
            d = d + s*h(yVec - yHat);
        end
    D = sparse(diag(d));
end   

