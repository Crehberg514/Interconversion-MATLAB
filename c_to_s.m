function [S0, SMats, lambdas] = c_to_s(C0, CMats, rhos, coeffSize)
    % Converts prony series creep complance to relaxtion modulus
    % using Cholesky decomposition. Using algorithm 3 of "Interconversion of
    % linearly viscoelastic material functions expressed as Prony series"
    %
    %
    % Parameters
    % ----------
    % C0 : 2D matrix
    %     Equilibrium relaxation
    % CMats : 3D matrix
    %     Relaxtion modulus coefficient
    % rhos : Vector
    %     Relaxation time constants
    % coeffSize : int
    %     The size dim of the the matrix i.e. 6 if 6x6
    %
    % Returns
    % -------
    % S0 : 2D matrix
    %     Instantaneous creep modulus
    % SMats : 3D matrix
    %     Creep modulus coefficient
    % lambdas : Vector
    %     Creep time constants

    % Number of coefficent matrices
    [~, ~, numCoeff] = size(CMats);

    % Final amount of coefficient matrices returned
    % used to size different variables
    finalNumCoeff = coeffSize * numCoeff;

    L1 = zeros(6);
    
    L1 = L1 + C0;
    
    for i = 1:numCoeff
        L1 = L1 + CMats(:,:,i);
    end
    
    L2 = [];
    
    for i = 1:numCoeff
        tempMat = rhos(i) * CMats(:,:,i);
        tempMat = chol(tempMat);
        L2 = [L2, tempMat'];
    end
    
    L3 = eye(finalNumCoeff);
    L3Rhos = repelem(rhos, coeffSize);
    L3 = L3Rhos .* L3;
    
    B = eye(finalNumCoeff);
    
    A1 = inv(L1);
    
    A2 = L2' * A1;
    A2 = A2';
    
    A3 = L3 - L2' * A1 * L2;
    
    [P,~,~] = svd(A3);
    
    A3Star = P' * A3 * P;
    
    A2Star = P' * A2';
    A2Star = A2Star';
    
    S0 = A1;
    
    SMats = zeros(coeffSize, coeffSize, finalNumCoeff);
    lambdas = zeros(1, finalNumCoeff);
    
    for m = 1:finalNumCoeff
        for i = 1:coeffSize
            for j = 1:coeffSize
                SMats(i,j,m) = (A2Star(i,m) * A2Star(j,m)) / A3Star(m,m);
            end
        end
        lambdas(m) = A3Star(m,m) / B(m,m);
    end
    
    lambdas = fliplr(lambdas);
    
    end % of the function