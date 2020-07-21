function [S0, SMats, lambdas] = c_to_s(C0, CMats, rhos)
    % Convert a prony series relaxation modulus to creep compliance.
    %
    % MATLAB numerical implementation of algorithm 3 from "Interconversion of
    % linearly viscoelastic material functions expressed as Prony series"[1].
    % Requiring the equilibrium relaxation modulus, the relaxation coefficient
    % modulus matrices, the inverted time series, and finally the size of the
    % square matrices.
    %
    % This version is of the algorithm is not limited to 6 by 6 matrices, but
    % has been minimally adjusted to accept any square matrices.
    %
    % Parameters
    %----------
    % C0 : 2D matrix
    %     The equilibrium relaxation modulus as a 2D square array.
    % CMats : 3D matrix
    %     The relaxation modulus coefficients as a 3D array, with the third
    %     dimension being the coefficient matrices and the first and second
    %     being the row and columns respectfully.
    % rhos : Vector
    %     Relaxation time constants
    %
    % Returns
    % -------
    % S0 : 2D matrix
    %     The instantaneous creep modulus returned as a 2D square array.
    % SMats : 3D matrix
    %     The creep modulus coefficient returned as a 3D array, with the third
    %     dimension being the coefficient matrices and the first and second
    %     being the row and columns respectfully.
    % lambdas : Vector
    %     The inverted creep time constants in a 1D array.
    %
    % References
    % ----------
    % [1] Luk-Cyr, J., Crochon, T., Li, C. et al. Interconversion of linearly
    % viscoelastic material functions expressed as Prony series: a closure.
    % Mech Time-Depend Mater 17, 53–82 (2013). 
    % https://doi.org/10.1007/s11043-012-9176-y
    

    % Find the size of the square matrix and the number of coefficent matrices
    [~, coeffSize, numCoeff] = size(CMats);

    % Calculate amount of coefficient matrices created
    finalNumCoeff = coeffSize * numCoeff;

    % Follow the steps outlined in Algorithm 3 from [1]
    L1 = zeros(coeffSize);
    
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
      
    end % of the function