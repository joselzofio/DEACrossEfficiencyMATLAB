function out = crossEffClassic(X, Y, DEAEff, n, m, s, dualX, dualY)
% DEA Cross-Efficiency model using the Sexton et al. (1986) ratio model:
% 
% This function solves the ratio cross-efficiency model as proposed by Sexton, T. R., Silkman, R. H., and Hogan, A. J. (1986) “Data 
% Envelopment Analysis: Critique and Extensions", in Measuring Efficiency: An Assessment of Data Envelopment Analysis, edited by R. H. 
% Silkman, New Directions for Program Evaluation 32 (Jossey-Bass, SanFrancisco/London). It corresponds to the exact ratio approach as 
% presented in expression (8) in Balk, B.M.,  de Koster, M.B.M.,  Kaps, C. and Zofio, J.L. (2017) "An Evaluation of Cross-Efficiency
% Methods, Applied to Measuring Warehouse Performance," ERIM Report Series Research in Management ERS-2017-015-LIS, Erasmus Research.
% Institute of Management (ERIM). https://repub.eur.nl/pub/103185.
% Agressive weights are used for the optimization.
% 
% Copyright 2019 Christian Kaps, José L. Zofío
% https://github.com/joselzofio/DEACrossEfficiencyMATLAB
%
% Version: 1.0
% LAST UPDATE: 1 July, 2019

    %Initialize weight matrix
    agressiveWeights = zeros(n,m + s);

    %Define global variables, so they can be used in the obj. function
    global funN funM funS funX funY funDMU funEff;
    funN = n;
    funM = m;
    funS = s;
    funX = X;
    funY = Y;
    funEff = DEAEff;
    
    %Switch off singular Matrix warnings that might pop up because of the
    %ratio nature of the set-up and the optimization procedure involved
    warning('off','MATLAB:nearlySingularMatrix')
    
    for DMUvalue = 1:n

        %Set global variable funDMU value to current DMU value
        funDMU = DMUvalue;

        %Decision Vector of 1
        objFunVector = ones(1,funN);

        %Set Function
        f= @funcRatioScore;

        %All calculated efficiencies with those weights have to be below 1 -
        A = [Y -X];
        b = zeros(n,1);

        %Force Efficiency to be the same as before
        Aeq = [zeros(1,s) X(DMUvalue,:);
            Y(DMUvalue,:) -DEAEff(DMUvalue,1)*X(DMUvalue,:)];
        beq = [1;0];

        %All Weights have lower bound of 0
        lb = zeros(1,(m + s));

        %Starting weights equal to DEA weights
        x0(1:m+s) = [dualY(DMUvalue,:) dualX(DMUvalue,:)];

        %Set options
        options = optimoptions(@fmincon,'MaxFunEvals',10^50);
        options = optimoptions(@fmincon,'Display','off');

        %Optimize and store weights in output matrix
        [z, fval, exitflag, output] = fmincon(f,x0,A, b, Aeq, beq, lb, [],[],options);
        agressiveWeights(DMUvalue,:) = z';

    end

    %Switch warning back on
    warning('on','MATLAB:nearlySingularMatrix')
    
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
            %Every appraisal for DMU1 is stored in column1
            agressivepeerApprScores(currentCalc,peerDMU) = dot(agressiveWeights(peerDMU,1:s),Y(currentCalc,:))/dot(agressiveWeights(peerDMU,(s+1):(s+m)),X(currentCalc,:));
        end

    end

    for currentCalc = 1:n

        %Calculate Geo mean WITH own appraisal
        agressivegeoEffIncl(currentCalc,1) = geomean(agressivepeerApprScores(currentCalc,:));
        %Calculate GEo mean WITHOUT own appraisal
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

