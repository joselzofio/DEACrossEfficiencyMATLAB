function out = crossEffOrdinal(X,Y,DEAEff,n,m,s)
% DEA Cross-Efficiency model using the Cook and Zhu (2006) ordinal model:
% 
% This function solves the standard cross-efficiency model as proposed by Cook, W. D., & Zhu, J. (2006) "Rank order data in DEA: 
% A general framework", in European Journal of Operational Research, 174(2), 1021-1038. For more detail on cross-efficiency calculation, please
% refer to Balk, B.M.,  de Koster, M.B.M.,  Kaps, C. and Zofio, J.L. (2017) "An Evaluation of Cross-Efficiency Methods, Applied 
% to Measuring Warehouse Performance," ERIM Report Series Research in Management ERS-2017-015-LIS, Erasmus Research.
% Institute of Management (ERIM). https://repub.eur.nl/pub/103185.
% Agressive weights are used for the optimization.
% 
% Copyright 2019 Christian Kaps, José L. Zofío
% https://github.com/joselzofio/DEACrossEfficiencyMATLAB
%
% Version: 1.0
% LAST UPDATE: 1 July, 2019

    minWeight =  1e-7; %Constraints value for DEA and XEff
    minWeight2 = 1e-7;

    %Create Limit Vector
    lb = zeros(1,s+m);
    A = [Y -X];
    for OutNIn = 1:(s+m)
        if sum(sum(A(:,OutNIn)==[0,1,-1])) ~= n %~= not equal
            lb(OutNIn) = minWeight;
        else
            lb(OutNIn) = minWeight2;
        end
    end

    DEAWeightsSame = zeros(n,s+m);

    for DMUvalue = 1:n
        %Concatenate all values of DMU row (output positive, input irrelevant, because forced to be one...)
        f = [Y(DMUvalue,:) zeros(1,m)];

        %Appraisal Scores have to be smaller or equal than 1
        A = [Y -X];
        b = zeros((n),1);

        %Inputs have to be 1
        Aeq = [zeros(1,s) X(DMUvalue,:)];
        beq = 1;

        %options to surpress the "different algorithm" warning
        options = optimoptions('linprog','Algorithm','dual-simplex', 'Display','off','OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-7, 'MaxIter',1e7);

        %optimize and store weights in output matrix DEA
        [z, ~, exitflag, ~, dual] = linprog(-f, A, b, Aeq, beq, lb, [],options);

        DEAWeightsSame(DMUvalue,:) = z';

    end

    DEAEffSameIntercept = sum(Y .* DEAWeightsSame(:,1:s),2);

    %Create Matrix of objective function weights (sum of all in/outputs excluding DMU)
    objFunValues = zeros(n,s+m);

    %Fill with output values
    for columnCount = 1:s
        allValuesThisColumn = sum(Y(:,columnCount));
        for DMUvalue = 1:n
            objFunValues(DMUvalue,columnCount)=allValuesThisColumn - Y(DMUvalue,columnCount);
        end
    end

    %Fill with input values (rightaftertheoutputvalues)
    for columnCount = (s+1):(s+m)
        allValuesThisColumn = sum(X(:,columnCount-s));
        for DMUvalue = 1:n
            objFunValues(DMUvalue,columnCount)=allValuesThisColumn - X(DMUvalue,columnCount-s);
        end
    end

    XeffWeights = zeros(n,s+m);
    for DMUvalue = 1:n
        %Concatenate all values of DMU row (output positive, input negative)
        %(spacebar just puts the values next to each other)
        f = [objFunValues(DMUvalue,1:s) -objFunValues(DMUvalue,(s+1):(s+m))];

        %Appraisal Scores have to be smaller or equal than 1
        A = [Y -X];
        b = zeros((n),1);

        %Inputs have to be 1
        Aeq = [[zeros(1,s) X(DMUvalue,:)];
            [Y(DMUvalue,:) zeros(1,m)]];
        beq = [1;
            DEAEffSameIntercept(DMUvalue)];

        %All Weights have lower bound of 0
        lb = (zeros(1,(s+m))+minWeight);

        %options to surpress the "different algorithm" warning
        options = optimoptions('linprog','Algorithm','dual-simplex', 'Display','off','OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-7, 'MaxIter',1e7);

        %optimize and store weights in output matrix XEFF
        [z, ~, exitflag, ~, dual] = linprog(-f, A, b, Aeq, beq, lb, [],options);

        XeffWeights(DMUvalue,:) = z';

    end

    %Calculate appraisal score for every DMU-DMU combination
    ApprScores = ones(n,n)* -1;
    for currentCalc = 1:n
        %Cycle through every DMU, except the calculated one
        for peerDMU = 1:n
            %Calculate Cross-Efficiency / Store in matrix
            %Every appraisal for DMU1 is stored in column1
            ApprScores(currentCalc,peerDMU) = dot(XeffWeights(peerDMU,1:s),Y(currentCalc,:))/dot(XeffWeights(peerDMU,(s+1):(s+m)),X(currentCalc,:));
        end
    end
       
    %Initialize variables to store results into later
    agressivepeerApprScores = zeros(n,n);
    agressivegeoEffIncl = zeros(n,1);
    agressivegeoEffExcl = zeros(n,1);
    agressiveariEffIncl = zeros(n,1);
    agressiveariEffExcl = zeros(n,1);
    
    for currentCalc = 1:n

        %Calculate Geo mean WITH own appraisal
        agressivegeoEffIncl(currentCalc,1) = geomean(ApprScores(currentCalc,:));
        %Calculate Geo mean WITHOUT own appraisal
        agressivegeoEffExcl(currentCalc,1) = geomean(ApprScores(currentCalc,[1:currentCalc-1,currentCalc+1:n]));
        %Calculate Arithmetic mean WITH own appraisal
        agressiveariEffIncl(currentCalc,1) = mean(ApprScores(currentCalc,:));
        %Calculate Arithmetic mean WITHOUT own appraisal
        agressiveariEffExcl(currentCalc,1) = mean(ApprScores(currentCalc,[1:currentCalc-1,currentCalc+1:n]));

    end

    out.CEPeerAppraisals =  ApprScores;
    out.CEEffArIn = agressiveariEffIncl;
    out.CEEffArEx = agressiveariEffExcl;
    out.CEEffGeIn = agressivegeoEffIncl;
    out.CEEffGeEx = agressivegeoEffExcl;

end

