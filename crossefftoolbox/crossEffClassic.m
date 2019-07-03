function out = crossEffClassic(X,Y,DEAEff,n,m,s) 
% DEA Cross-Efficiency model using the Sexton et al. (1986) classic model:
% 
% This function solves the standard cross-efficiency model as proposed by Sexton, T. R., Silkman, R. H., and Hogan, A. J. (1986) “Data 
% Envelopment Analysis: Critique and Extensions", in Measuring Efficiency: An Assessment of Data Envelopment Analysis, edited by R. H. 
% Silkman, New Directions for Program Evaluation 32 (Jossey-Bass, SanFrancisco/London). It corresponds to the linear surrogate as 
% presented in expression (9) in Balk, B.M.,  de Koster, M.B.M.,  Kaps, C. and Zofio, J.L., (2017) "An Evaluation of Cross-Efficiency
% Methods, Applied to Measuring Warehouse Performance," ERIM Report Series Research in Management ERS-2017-015-LIS, Erasmus Research.
% Institute of Management (ERIM). https://repub.eur.nl/pub/103185.
% Agressive weights are used for the optimization.
% 
% Copyright 2019 Christian Kaps, José L. Zofío
% https://github.com/joselzofio/DEACrossEfficiencyMATLAB
%
% Version: 1.0
% LAST UPDATE: 1 July, 2019

    %Create Matrix of objective function weights (sum of all in/outputs excluding DMU)
    objFunValues = zeros(n,m+s);

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

    %Set up the optimization for each DMU
    agressiveWeights = zeros(n,m+s);

    for DMUvalue = 1:n

        %Concatenate all values of DMU row (output positive, input negative)
        f = [objFunValues(DMUvalue,1:s) -objFunValues(DMUvalue,(s+1):(s+m))];

        %All calculated efficiencies with those weights have to be below 1 -
        A = [Y -X];
        b = zeros((n),1);

        %Force Efficiency to be the same as before
        Aeq = [zeros(1,s) X(DMUvalue,:);
            Y(DMUvalue,:) -DEAEff(DMUvalue,1)*X(DMUvalue,:)];
        beq = [1;0];

        %All Weights have lower bound of 0
        lb = zeros(1,(m+s));

        %Options to surpress the "different algorithm" warning
        options = getDEAoptions(n);
        linprogoptions = options.optimopts;

        %optimize and store weights in output matrix
        [z, ~, exitflag, ~, dual] = linprog(-f, A, b, Aeq, beq, lb, [],linprogoptions);

        agressiveWeights(DMUvalue,:) = z';

    end

    %Initialize variables to store results into later
    agressivepeerApprScores = zeros(n,n);
    agressivegeoEffIncl = zeros(n,1);
    agressivegeoEffExcl = zeros(n,1);
    agressiveariEffIncl = zeros(n,1);
    agressiveariEffExcl = zeros(n,1);

    %Calculate appraisal score for every DMU-DMU combination
    for currentCalc = 1:n
        %Cycle through every DMU, except the calculated one
        for peerDMU = 1:n

            %Calculate Cross-Efficiency / Store in matrix
            agressivepeerApprScores(currentCalc,peerDMU) = dot(agressiveWeights(peerDMU,1:s),Y(currentCalc,:))/dot(agressiveWeights(peerDMU,(s+1):(s+m)),X(currentCalc,:));
        end

    end

    for currentCalc = 1:n

        %Calculate Geo mean WITH own appraisal
        agressivegeoEffIncl(currentCalc,1) = geomean(agressivepeerApprScores(currentCalc,:));
        %Calculate Geo mean WITHOUT own appraisal
        agressivegeoEffExcl(currentCalc,1) = geomean(agressivepeerApprScores(currentCalc,[1:currentCalc-1,currentCalc+1:n]));
        %Calculate Arithmetic mean WITH own appraisal
        agressiveariEffIncl(currentCalc,1) = mean(agressivepeerApprScores(currentCalc,:));
        %Calculate Arithmetic mean WITHOUT own appraisal
        agressiveariEffExcl(currentCalc,1) = mean(agressivepeerApprScores(currentCalc,[1:currentCalc-1,currentCalc+1:n]));

    end

    out.CEPeerAppraisals =  agressivepeerApprScores;
    out.CEEffArIn = agressiveariEffIncl;
    out.CEEffArEx = agressiveariEffExcl;
    out.CEEffGeIn = agressivegeoEffIncl;
    out.CEEffGeEx = agressivegeoEffExcl;

end

