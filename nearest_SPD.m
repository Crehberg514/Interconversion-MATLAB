function SPD = nearest_SPD(A)
    % nearest_SPD - finds the nearest Symmetric Positive Definite matrix to
    % A
    %
    % From Higham's "Computing a nearest symmetric positive semidefinite
    % matrix"
    %
    % Parameters
    % ----------
    % A : 2D matrix
    %
    % Returns
    % -------
    % SPD : 2D nearest symmetric positive definite matrix
    
    % Symmetrize A into B
    B = (A + A')/2;
    
    % Compute the symmetric polar factor of B
    [~,S,V] = svd(B);
    H = V*S*V';
    
    % Calculate nearest SPD
    SPD = (B+H)/2;
    
    % Ensure symmetric
    SPD = (SPD + SPD')/2;
    
    % Test if SPD is pos def
    [~,p] = chol(SPD);
    k = 0;
    ident = eye(size(A));
    
    while p ~= 0
          
      % If not pos def, adjust  due to floating point errors
      mineig = max(eig(SPD));
      SPD = SPD + ident * (-mineig * k^2 + eps(mineig));
            
      [~,p] = chol(SPD);
      k = k + 1;
      
    end % of the function