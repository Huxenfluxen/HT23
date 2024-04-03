%% LAB 5 PART 2

%A2 is the matrix from cooling flange
%A3 is the matrix from convdiff

%% part 2a + 2b
clear;


% table of computation times for pcg method saved as pcgTable

cooling_flange = matfile('cooling_flange.mat');
A2 = cooling_flange.A;

pcgTable = pTable(A2);






% figure
% spy(A2);
% title("Sparsity of A");






%% part 2c
clear;

% table of computation times for gmres method saved as gmresTable

convdiff = matfile('convdiff.mat');
A3 = convdiff.A;

nonConvPlot(A3); % shows non-convergence of unsymmetric matrix

gmresTable = gTable(A3);




% figure
% spy(A3);
% title("Sparsity of A");



%% functions


function pcgTable = pTable(A2)
% plots convergence and returns pcgTable

    tol = 1e-4;
    maxit = 1000;
    
    N2 = length(A2);
    b2 = rand(N2,1);
    
    %   backslash
    tic;
    A2\b2;
    BackslashTime2a = toc;
    
    %   a) no precondition
    tic;
    [Xa,FLAGa,RELRESa,ITERa,RESVECa] = pcg(A2,b2, tol, maxit);
    pcgTime2a = toc;
    
    % part 2b
    
    M = diag(A2); M = diag(M);
    L = ichol(A2);
    
    %   bi) M = diag(A)
    tic;
    [Xbi,FLAGbi,RELRESbi,ITERbi,RESVECbi] = pcg(A2,b2, tol, maxit,M);
    pcgTime2bi = toc;
    
    %   bii) M = LL'
    tic;
    [Xbii,FLAGbii,RELRESbii,ITERbii,RESVECbii] = pcg(A2,b2, tol, maxit,L,L');
    pcgTime2bii = toc;
    
    figure
    plot(0:ITERa, RESVECa);
    hold on 
    plot(0:ITERbi, RESVECbi);
    plot(0:ITERbii, RESVECbii);
    hold off
    legend("No precondition", "M = diag(A)", "M = LL^T ");
    title("Convergence history of pcg solver")
    
    
    precondition = ["none"; "M=diag(A)"; "M=LL'"];
    iterations = [ITERa; ITERbi; ITERbii];
    relativeResidual = [RELRESa; RELRESbi; RELRESbii];
    PCGcompTimes = [pcgTime2a; pcgTime2bi; pcgTime2bii];
    backSlashTime = [BackslashTime2a; BackslashTime2a; BackslashTime2a ];
    
    pcgTable = table(precondition,iterations, relativeResidual, PCGcompTimes, backSlashTime);
end



function nonConvPlot(A3)
    tol = 1e-4;
    maxit = 2000;
    
    N3 = length(A3);
    b3 = rand(N3,1);

    [X2ci,FLAG2ci,RELRES2ci,ITER2ci,RESVEC2ci] = pcg(A3,b3, tol, maxit);
    figure
    plot(RESVEC2ci);
    title("(Non-)convergence history of non-symmetric matrix")
end


    

function gmresTable = gTable(A3)

    tol = 1e-4;
    maxit = 2000;
    
    N3 = length(A3);
    b3 = rand(N3,1);
    
    tic;
    A3\b3;
    BackslashTime2c = toc;
    
    
    
    tic;
    [X2cii,FLAG2cii,RELRES2cii,ITER2cii,RESVEC2cii] = gmres(A3,b3, [],tol, maxit);
    gmresTime2cii = toc;
    
    Mc = diag(A3); Mc = diag(Mc);
    tic;
    [X2ciii,FLAG2ciii,RELRES2ciii,ITER2ciii,RESVEC2ciii] = gmres(A3,b3, [],tol, maxit, Mc);
    gmresTime2ciii = toc;
    
    [L,U] = ilu(A3);
    tic;
    [X2civ,FLAG2civ,RELRES2civ,ITER2civ,RESVEC2civ] = gmres(A3,b3, [],tol, maxit, L, U);
    gmresTime2civ = toc;
    
    
    
    figure
    plot(RESVEC2cii);
    hold on
    plot(RESVEC2ciii);
    plot(RESVEC2civ);
    hold off
    legend("No precondition", "M = diag(A)", "M = LU")
    title("Convergence history of gmres solver")
    
    
    precondition = ["none"; "M=diag(A)"; "M=LU"];
    iterations = [ITER2cii(2); ITER2ciii(2); ITER2civ(2)];
    relativeResidual = [RELRES2cii; RELRES2ciii; RELRES2civ];
    GMREScompTimes = [gmresTime2cii; gmresTime2ciii; gmresTime2civ];
    backSlashTime = [BackslashTime2c; BackslashTime2c; BackslashTime2c];
    
    gmresTable = table(precondition,iterations, relativeResidual, GMREScompTimes, backSlashTime);

end





