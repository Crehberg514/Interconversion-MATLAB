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
    p = 1;
    k = 0;
    while p ~= 0
      [~,p] = chol(SPD);
      k = k + 1;
      
      % If not pos def, adjust  due to floating point errors
      if p ~= 0
        mineig = min(eig(SPD));
        SPD = SPD + (-mineig*k.^2 + eps(mineig))*eye(size(A));
      end
      
    end % of the function