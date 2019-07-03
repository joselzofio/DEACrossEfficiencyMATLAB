function [ out ] = crossEff( X, Y, crossEffModel)
% DEA Cross-Efficiency function, faciltating the solution of one of the following five cross-efficiency models:
% 1) Sexton et al. (1986) classic model
% 2) Sexton et al. (1986) ratio model
% 3) Cook and Zhu (2014) multiplicative model
% 4) Liang et al. (2008) game theory 
% 5) Cook and Zhu (2006) ordinal model
%
% For more detail on the different models and their relative benefits, please refer to Balk, B.M.,  de Koster, M.B.M.,  
% Kaps, C. and Zofio, J.L. (2017). "An Evaluation of Cross-Efficiency Methods, Applied to Measuring Warehouse Performance," 
% ERIM Report Series Research in Management ERS-2017-015-LIS, Erasmus Research.
% Institute of Management (ERIM). https://repub.eur.nl/pub/103185.
% 
% Copyright 2019 Christian Kaps, José L. Zofío
% https://github.com/joselzofio/DEACrossEfficiencyMATLAB
%
% Version: 1.0
% LAST UPDATE: 1 July, 2019

    % Check Size
    if size(X,1) ~= size(Y,1);
        error('Number of rows in X must be equal to number of rows in Y')
    end    
    
    % Get number of DMUs (n), inputs (m) and outputs (s)
    [n, m] = size(X);
    s = size(Y,2);
    
    % Check Correct specification of crossEffModel
    if ~any(strcmp(["classic", "ratio", "multiplicative", "gameTheory", "ordinal"],crossEffModel))
        error('Cross efficiency model must be specified. Choose from classic, ratio, multiplicative, gameTheory and ordinal.')
    end    
    
    % Calculate input oriented technical-efficiency DEA model
    DEA = dea(X,Y,'orient','io');
    
    % Select which model to run
    switch(crossEffModel)
        case {"classic"}
            CEff = crossEffClassic(X,Y,DEA.eff,n,m,s);
        case {"ratio"}
            CEff = crossEffRatio(X,Y,DEA.eff,n,m,s, DEA.dual.X, DEA.dual.Y);
        case {"multiplicative"}
            CEff = crossEffMultiplicativeCRS(X,Y,DEA.eff,n,m,s);
        case {"gameTheory"}
            CEff = crossEffGameTheory(X,Y,DEA.eff,n,m,s, DEA.dual.X, DEA.dual.Y);
        case {"ordinal"}
            CEff = crossEffOrdinal(X,Y,DEA.eff,n,m,s)
    end     
    
    %Dea Package output
    out.DEAFull = DEA;
    out.DEAEffScore = DEA.eff;
    
    %Assign different efficiency scores to output variables
    try 
        out.CEPeerAppraisals =  CEff.CEPeerAppraisals;
    catch
        out.CEPeerAppraisals = [];
    end
    
    try 
        out.CEEffArIn =  CEff.CEEffArIn;
        out.CEEffArInAvg = mean(out.CEEffArIn);
    catch
        out.CEEffArIn =  [];
        out.CEEffArInAvg = [];
    end
    
    try
        out.CEEffArEx =  CEff.CEEffArEx;
        out.CEEffArExAvg =  mean(out.CEEffArEx);
    catch
        out.CEEffArEx =  [];
        out.CEEffArExAvg = [];
    end
    
    try
        out.CEEffGeEx =  CEff.CEEffGeEx;
        out.CEEffGeExAvg = mean(out.CEEffGeEx);
    catch
        out.CEEffGeEx =  [];
        out.CEEffGeExAvg = [];
    end
    
    try
        out.CEEffGeIn =  CEff.CEEffGeIn;
        out.CEEffGeInAvg =  mean(out.CEEffGeIn);
    catch
        out.CEEffGeIn =  [];
        out.CEEffGeInAvg = [];
    end
   
    %print average results
    fprintf('Average %s cross-efficiency of %d DMUs, based on selected averaging method: \n', crossEffModel, n)
    fprintf("%6.4f - Arithmetic Mean including self appraisal. \n", out.CEEffArInAvg)
    fprintf("%6.4f - Arithmetic Mean excluding self appraisal. \n", out.CEEffArExAvg)
    fprintf("%6.4f - Geometric Mean including self appraisal. \n", out.CEEffGeInAvg)
    fprintf("%6.4f - Geometric Mean excluding self appraisal. \n",  out.CEEffGeExAvg)

end

