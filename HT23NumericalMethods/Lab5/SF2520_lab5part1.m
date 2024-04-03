%% LAB 5 PART 1

%% Part 1a
clear;

K = 5000; % max iterations

     %   [n,d]
ndVals = [45,2;...
          13,3;... 
          72,2;...
          200,2;
          34, 3;...
          ]; 


% check convergence using values below. when n=2n, 
% iterations k=2k for conjugate (O(n)), k=4k for jacobi (O(n^2)).

% ndVals = [10,2;...
%           20,2;... 
%           40,2;...
%           80,2;
%          ]; 


% select method as fist argument of function 

figure
ConvergenceComparison(@JacobiMethod, K, ndVals);
sgtitle("Jacobi Method");



figure
ConvergenceComparison(@ConGradMethod, K, ndVals);
sgtitle("Conjugate Gradient Method");




%% part 1b

% ConGrad is faster for smaller n...
clear;

ndVals = [10000,1;...
            100,2;...
             10,4];


ConGradTolTable = CGtab(ndVals);


% check the solutions are approximately the same, difference saved as vector d
A = lap(10,2);
b = rand(10^2,1);
c = ConGradTolMethod(A,b);
backslash = A\b;
d=c-backslash;




%% functions part  1


function ConGradTable = CGtab(ndVals)
% creates table comparing backslash to conjugate gradient methods in
% computational time
    nVals = ndVals(:,1);
    dVals = ndVals(:,2);
    NVals = nVals.^dVals;
    
    l = length(ndVals);
    
    ConGradTimes = zeros(l,1);
    BackSlashTimes = zeros(l,1);
    ConGradFaster = ones(l,1);

    for i =  1:l
        n= ndVals(i,1); d= ndVals(i,2);
        N = n^d;
        b1 = rand(N,1);
        A1 = lap(n,d);
        tic;
        xConGrad = ConGradTolMethod(A1,b1);
        ConGradTime = toc;
        tic;
        xBackslash = A1\b1;
        BackSlashTime = toc;
        BackSlashTimes(i,1) = BackSlashTime;
        ConGradTimes(i,1) = ConGradTime;
        ConGradFaster(i,1) = BackSlashTime > ConGradTime;
        
    end
    ConGradTable = table(NVals, nVals, dVals, ConGradTimes, BackSlashTimes, ConGradFaster);
end





function ConvergenceComparison(method, K, ndVals)
% function takes method (jacobi/ConGrad), max nr iterations (K) 
% and list of n and d values
    ln = length(ndVals);
    u = ceil(sqrt(ln));
    v = ceil(ln/u);

    for i = 1:ln
        
        n= ndVals(i,1); d= ndVals(i,2);
        

        N=n^d;
        A = lap(n,d);
        b = rand(N,1);
        
        x = method(A, b, K);

        rVec = zeros(K,1);
        for k = 1:K
            r = norm(A*x(:,k) - b)/ norm(b);
            rVec(k,1) = r;
        end
        
        
        
        semilogy(1:K, rVec, "LineWidth", 1.5);
        hold on
        
        xlabel("k");
        ylabel("r(k)")
        
    end
    n= ndVals(:,1); d= ndVals(:,2);
    N=n.^d;
    legend("n="+n + ", d="+d + ", N="+N, "location", "east");
    hold off
end




function x = ConGradTolMethod(A,b)
% Conjugate gradient method with tolerance 1e-10
% returns final x vector
    N=length(A);
    x = zeros(N,1); % initial guess
    
    beta0 = 0;
    pMin = x;
    rOld = b;
    p= rOld + beta0*pMin;
    
        
    tol = 1e-10;
    res = inf;
    i=1;
   
     while res > tol
          Ap = A*p;
          rrOld = rOld'*rOld;
          alpha= rrOld/(p'*Ap);
          x=x+alpha*p;
          rNew = rOld - alpha*Ap;
          beta= rNew'*rNew/(rrOld);
          p = rNew + beta*p;
          rOld=rNew; 
          res = norm(A*x - b)/ norm(b);
          i = i+1;
     end
    
     
end




function xMat = ConGradMethod(A,b,K)
% Conjugate gradient method

    N=length(A);
    x0 = zeros(N,1); % initial guess
    
    beta0 = 0;
    pMin = x0;
    rOld = b;
    p= rOld + beta0*pMin;
    
    xMat = zeros(N,K);
    x = x0;
        
     for k=1:K
          Ap = A*p;
          rrOld = rOld'*rOld;
          alpha= rrOld/(p'*Ap);
          x=x+alpha*p;
          rNew = rOld - alpha*Ap;
          beta= rNew'*rNew/(rrOld);
          p = rNew + beta*p;
          
          rOld=rNew; 
          xMat(:,k) = x;
     end
end






function xMat = JacobiMethod(A, b, K)
    
    D = diag(A); % elements on diagonal of A
    LpU = A - diag(D); % A without the diagonal elements (= L+U)
    
    N=length(A);
    nVec = ones(N,1); 
    
    xMat = zeros(N,K);

    x0 = 0 * nVec; % initial guess
    x = x0;
    
    for k = 1:K
        x = (b - LpU * x) ./ D;
        xMat(:,k) = x;
    end
end





function A = lap(n,d)

    % A = LAP(N,D) returns the system matrix corresponding to the
    % standard second order discretization of the Laplace operator
    % for D dimensions and N unknowns in each coordinate direction.
    % The size of the matrix is thus N^D times N^D.
    e = ones(n,1);
    A1 = -spdiags([e -2*e e], -1:1, n, n);
    I1 = speye(n,n);
    A=A1;
    I=I1;
    
    for k=2:d
        A = kron(A,I1)+kron(I,A1);
        I = kron(I,I1);
    end
    A = A*n^2;
end