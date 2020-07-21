function isSPD = is_pos_def(Mat)
    % Checks if the matrix is positive-definite.
    %
    % Checks via Cholesky decomposition if the given matrix is
    % positive-definite. If the matrix is positive-definite then the function
    % returns True, otherwise, catches the thrown exception and returns False.
    %
    % Parameters
    % ----------
    % B : Matrix
    %     Square matrix.
    %
    % Returns
    % -------
    % Bool : Bool
    %     Returns True if the matrix is positive-definite and false otherwise.
    
    [~, flag] = chol(Mat);
    
    if flag == 0
        isSPD = true;
    else
        isSPD = false;
    end
    
    end % of the function