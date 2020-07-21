function SPD = nearest_SPD(A)
    % Find the nearest positive-definite matrix to the input matrix.
    % 
    % Finds the nearest semi-positive definite matrix using the procedure
    % outlined in Higham's "Computing a nearest symmetric positive semidefinite
    % matrix." [1]
    %
    % Parameters
    % ----------
    % A : Matrix
    %     A 2D numpy array.
    %
    % Returns
    % -------
    % SPD : Matrix
    %     Semi-positive-definite 2D numpy array.
    %
    % References
    % ----------
    % [1] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    % matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
    
    % Symmetrize matrix A into matrix B
    B = (A + A')/2;
    
    % Compute the symmetric polar factor of matrix B
    [~,S,V] = svd(B);
    H = V*S*V';
    
    % Calculate the nearest semi-positive-definite matrix
    SPD = (B+H)/2;
    
    % Ensure the matrix is symmetric
    SPD = (SPD + SPD')/2;
    
    % Test if the matrix is positive-definite
    if is_pos_def(SPD)
      return
    end
    
    % Calculate minor adjustments to the matrix if it is not positive-definite
    % This is because of the how computers handle and calculate numbers
    k = 0;
    ident = eye(size(A));
    
    while ~is_pos_def(SPD)
      
      % If not pos def, adjust  due to floating point errors
      mineig = max(eig(SPD));
      SPD = SPD + ident * (-mineig * k^2 + eps(mineig));
      k = k + 1;
      
    end % of the function