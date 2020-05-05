function [C0, CMats, rhos] = s_to_c(S0, SMats, lambdas, coeffSize)
    % Converts prony series creep complance to relaxtion modulus
    % using Cholesky decomposition. Using algorithm 4 of "Interconversion of
    % linearly viscoelastic material functions expressed as Prony series"
    %
    %
    % Parameters
    % ----------
    % S0 : 2D matrix
    %     Instantaneous creep modulus
    % SMats : 3D matrix
    %     Creep modulus coefficient
    % lambdas : Vector
    %     Creep time constants
    % coeffSize : int
    %     The size dim of the the matrix i.e. 6 for 6x6
    %
    % Returns
    % -------
    % C0 : 2D matrix
    %     Equilibrium relaxation
    % CMats : 3D matrix
    %     Relaxtion modulus coefficient
    % rhos : Vector
    %     Relaxation time constants
    
    % Number of coefficent matrices
    [~, ~, numCoeff] = size(SMats);

    % Final amount of coefficient matrices returned
    % used to size different variables
    finalNumCoeff = coeffSize * numCoeff;

    A1 = S0;

    A2 = [];

    for i = 1:numCoeff
        tempMat = lambdas(i) * SMats(:,:,i);
        tempMat = chol(tempMat);
        A2 = [A2, tempMat'];

    end
    
    A3 = eye(finalNumCoeff);
    A3Lambdas = repelem(lambdas, coeffSize);
    A3 = A3Lambdas .* A3;
    
    B = eye(finalNumCoeff);
    
    L1 = inv(A1);
    
    L2 = A2' * L1;
    L2 = L2';
    
    L3 = A3 + A2' * L1 * A2;
    
    [P,~,~] = svd(L3);
    
    L3Star = P' * L3 * P;
    
    L2Star = P' * L2';
    L2Star = L2Star';
    
    CMats = zeros(coeffSize, coeffSize, finalNumCoeff);
    rhos = zeros(1, finalNumCoeff);
    
    for n = 1:finalNumCoeff
        for i = 1:coeffSize
            for j = 1:coeffSize
                
                CMats(i,j,n) = (L2Star(i, n) * L2Star(j, n)) / L3Star(n, n);
                
            end
        end
        rhos(n) = L3Star(n,n) / B(n,n);
    end

    C0 = L1;
    
    for i = 1:finalNumCoeff
        C0 = C0 - CMats(:,:,i);
    end
    
end % of the function