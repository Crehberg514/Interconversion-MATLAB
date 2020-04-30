function isSPD = is_pos_def(Mat)
    % Comment describing the function
    
    [~, flag] = chol(Mat);
    
    if flag == 0
        isSPD = true;
    else
        isSPD = false;
    end
    
    end % of the function