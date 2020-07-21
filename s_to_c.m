function [C0, CMats, rhos] = s_to_c(S0, SMats, lambdas)
    % Convert a prony series creep compliance to relaxation modulus.
    %
    % MATLAB numerical implementation of algorithm 4 from "Interconversion of
    % linearly viscoelastic material functions expressed as Prony series"[1].
    % Requiring the instantaneous creep modulus, the matrix creep coefficient
    % modulus, the inverted time series, and finally the size of the square
    % matrices.
    %
    % This version is of the algorithm is not limited to 6 by 6 matrices, but
    % has been minimally adjusted to accept any square matrices.
    %
    % Parameters
    % ----------
    % S0 : 2D matrix
    %     The instantaneous creep modulus as a 2D square array.
    % SMats : 3D matrix
    %     The creep modulus coefficient as a 3D array, with the third dimension
    %     being the coefficient matrices and the first and second being the row
    %     and columns respectfully.
    % lambdas : Vector
    %     The inverted creep time constants in a 1D array, in descending
    %     order.
    %
    % Returns
    % -------
    % C0 : 2D matrix
    %     The equilibrium relaxation modulus returned as a 2D square array.
    % CMats : 3D matrix
    %     The relaxation modulus coefficients returned as a 3D array, with the
    %     third dimension being the coefficient matrices and the first and
    %     second and being the row and columns respectfully.
    % rhos : Vector
    %     The inverted relaxation time constants in a 1D array.
    
    
    % Find the size of the square matrix and the number of coefficent matrices
    [~, coeffSize, numCoeff] = size(SMats);

    % Calculate amount of coefficient matrices created
    finalNumCoeff = coeffSize * numCoeff;

    % Follow the steps outlined in Algorithm 4 from [1]
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