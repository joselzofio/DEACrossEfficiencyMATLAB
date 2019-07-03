function out = crossEffMultiplicativeCRS(X,Y,DEAEff,n,m,s)
% DEA Cross-Efficiency model using the Cook and Zhu (2014) multiplicative model:
%
% This function solves the multiplicative cross-efficiency model as proposed by Cook, W. D., & Zhu, J. (2014) "DEA Cobb–Douglas 
% frontier and cross-efficiency", in Journal of the Operational Research Society, 65(2), 265-268. It corresponds to the model 
% in expression (13) in Balk, B.M.,  de Koster, M.B.M.,  Kaps, C. and Zofio, J.L. (2017) "An Evaluation of Cross-Efficiency
% Methods, Applied to Measuring Warehouse Performance," ERIM Report Series Research in Management ERS-2017-015-LIS, Erasmus Research.
% Institute of Management (ERIM). https://repub.eur.nl/pub/103185.
% Agressive weights are used for the optimization.
% 
% Copyright 2019 Christian Kaps, José L. Zofío
% https://github.com/joselzofio/DEACrossEfficiencyMATLAB
%
% Version: 1.0
% LAST UPDATE: 1 July, 2019

    %Replace 0 data with e, to make the multiplicative model workable.
    X(X==0) = 0.0001;
    Y(Y==0) = 0.0001;

    % Take logarithms of all values to reduce skewness
    lnX = log(X);
    lnY = log(Y);

    %lower & upper bound bounds for variables
    DEAlb = ones(1,(s+m));
    DEAub = repelem(10,(s+m));

    %Inequality constraints
    DEAA = [lnY -lnX];
    DEAb = zeros(n,1);

    %options to surpress the "different algorithm" warning
    options = getDEAoptions(n);
    linprogoptions = options.optimopts;

    %Initialize weight matrix
    geometricDEAWeights = zeros(n,m+s);
    geometricDEAEfficiencies = zeros(n,1);

    %DEA Optimization function
    for DMUvalue = 1:n

        DEAf = [lnY(DMUvalue,:) -lnX(DMUvalue,:)];
        [z, ~, exitflag, ~, dual] = linprog(-DEAf, DEAA, DEAb,[],[], DEAlb, DEAub,linprogoptions);
        %Safe weights
        geometricDEAWeights(DMUvalue,:) = z';
        %Calculate efficiencies
        geometricDEAEfficiencies(DMUvalue,1) = prod(Y(DMUvalue,:).^geometricDEAWeights(DMUvalue,1:s))/prod(X(DMUvalue,:).^geometricDEAWeights(DMUvalue,1+s:s+m));

    end

    %log of Efficiencies
    lnEff = log(geometricDEAEfficiencies);

    %Create Cross-Efficiency Output Inequality Matrix
    outputEquationsTable = zeros(n*n,n*s);
    for outputCount = 1:s
        for kValue = 1:n
            for jValue = 1:n
                outputEquationsTable((((kValue-1)*n)+jValue),((outputCount-1)*n)+kValue) = lnY(jValue,outputCount);
            end
        end
    end

    %Create Cross-Efficiency Input Inequality Matrix
    inputEquationsTable = zeros(n*n,n*m);
    for inputCount = 1:m
        for kValue = 1:n
            for jValue = 1:n
                inputEquationsTable((((kValue-1)*n)+jValue),((inputCount-1)*n)+kValue) = -lnX(jValue,inputCount);
            end
        end
    end

    A = [outputEquationsTable, inputEquationsTable];
    b = zeros(size(A,1),1);

    %Create Cross-Efficiency Output Equality Matrix
    outputEquationsTable2 = zeros(n,n*s);
    for outputCount = 1:s
        for kValue = 1:n
            outputEquationsTable2(kValue,((outputCount-1)*n)+kValue) = lnY(kValue,outputCount);
        end
    end

    %Create Cross-Efficiency Input Equality Matrix
    inputEquationsTable2 = zeros(n,n*m);
    for inputCount = 1:m
        for kValue = 1:n
            inputEquationsTable2(kValue,((inputCount-1)*n)+kValue) = -lnX(kValue,inputCount);
        end
    end

    Aeq = [outputEquationsTable2, inputEquationsTable2];
    beq = lnEff;

    %Create lower and upper bounds
    lb = ones(1,n*(m+s));
    ub = repelem(20, n*(m+s));

    %Initialize results matrices
    geometricXeffWeights = zeros(n,size(A,2));
    standardAppraisals = zeros(n,n);
    maximumAppraisals = zeros(n,n);

    %Create Objective Function
    for DMUvalue = 1:n

        f = [repelem(lnY(DMUvalue,:),n) -repelem(lnX(DMUvalue,:),n)];
        [z, ~, exitflag, ~, dual] = linprog(-f, A, b, Aeq, beq, lb, ub,linprogoptions);
        %Safe weights
        geometricXeffWeights(DMUvalue,:) = z';

        %Calculate X-efficiencies
        for peerDMU = 1:n
            standardAppraisals(DMUvalue,peerDMU) = prod(Y(peerDMU,:).^geometricXeffWeights(DMUvalue,DMUvalue:n:(s*n)))/prod(X(peerDMU,:).^geometricXeffWeights(DMUvalue,((n*s)+DMUvalue):n:(s*n+m*n)));
            maximumAppraisals(peerDMU,DMUvalue) = prod(Y(DMUvalue,:).^geometricXeffWeights(DMUvalue,peerDMU:n:(s*n)))/prod(X(DMUvalue,:).^geometricXeffWeights(DMUvalue,((n*s)+peerDMU):n:(s*n+m*n)));
        end

    end
    
    %Initialize variables to store results into later
    agressivegeoEffIncl = zeros(n,1);
    agressivegeoEffExcl = zeros(n,1);
    agressiveariEffIncl = zeros(n,1);
    agressiveariEffExcl = zeros(n,1);
    
    
    %After termination, calculate efficiency scores
    for currentCalc = 1:n

        %Calculate Geo mean WITH own appraisal
        agressivegeoEffIncl(currentCalc,1) = geomean(standardAppraisals(currentCalc,:));
        %Calculate Geo mean WITHOUT own appraisal
        agressivegeoEffExcl(currentCalc,1) = geomean(standardAppraisals(currentCalc,[1:currentCalc-1,currentCalc+1:n]));
        %Calculate Arithmetic mean WITH own appraisal
        agressiveariEffIncl(currentCalc,1) = mean(standardAppraisals(currentCalc,:));
        %Calculate Arithmetic mean WITHOUT own appraisal
        agressiveariEffExcl(currentCalc,1) = mean(standardAppraisals(currentCalc,[1:currentCalc-1,currentCalc+1:n]));

    end
    
    out.CEPeerAppraisals =  standardAppraisals;
    out.CEEffArIn = agressiveariEffIncl;
    out.CEEffArEx = agressiveariEffExcl;
    out.CEEffGeIn = agressivegeoEffIncl;
    out.CEEffGeEx = agressivegeoEffExcl;

end

