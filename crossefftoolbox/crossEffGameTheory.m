function out = crossEffGameTheory(X,Y,DEAEff,n,m,s, dualX, dualY)
% DEA Cross-Efficiency model using the Liang et al. (2008) game theory model:
% 
% This function solves the game theoretic, nash-equilibrium cross-efficiency model as proposed by Liang, L., Wu, J., Cook, W. D., & 
% Zhu, J. (2008) "The DEA game cross-efficiency model and its Nash equilibrium", in Operations Research, 56(5), 1278-1288. It corresponds 
% to the model presented in expression (15) in Balk, B.M.,  de Koster, M.B.M.,  Kaps, C. and Zofio, J.L. (2017) "An Evaluation 
% of Cross-Efficiency Methods, Applied to Measuring Warehouse Performance," ERIM Report Series Research in Management ERS-2017-015-LIS, 
% Erasmus Research. Institute of Management (ERIM). https://repub.eur.nl/pub/103185.
% Agressive weights are used for the optimization. 
%
% Copyright 2019 Christian Kaps, José L. Zofío
% https://github.com/joselzofio/DEACrossEfficiencyMATLAB
%
% Version: 1.0
% LAST UPDATE: 1 July, 2019

    %Run-time warning
    fprintf("Depending on size of input / output matrices, processor and speed of convergence, this method may very long. \n At the beginning of each optimization iteration, the current round number will be printed to the console \n")

    peerApprScores = zeros(n);

    %Calculate appraisal score for every DMU-DMU combination
    for currentCalc = 1:n
        %Cycle through every DMU, except the calculated one
        for peerDMU = 1:n
            %Calculate Cross-Efficiency / Store in matrix
            %Every appraisal for DMU1 is stored in column1
            peerApprScores(currentCalc,peerDMU) = dot(dualY(peerDMU,:),Y(currentCalc,:))/dot(dualX(peerDMU,:),X(currentCalc,:));
        end
    end

    %Calculate cross-efficiency scores
    for currentCalc = 1:n
        %Calculate Arithmetic mean WITH own appraisal
        alpha(currentCalc,1) = mean(peerApprScores(currentCalc,:));
    end

    %Initialize Appraisal Scores
    peerAppraisalScores = zeros(n,n);

    %Start iterations
    iteration = 1;
    condition = true;

    while condition
        fprintf("current iteration: %d \n", iteration)
        for DMUjValue = 1:n
            for DMUdValue = 1:n

                %Objective Function
                f = [zeros(1,m),-Y(DMUjValue,:)];

                %All results greater 0
                A = [-X Y;
                    alpha(DMUdValue,iteration)*X(DMUdValue,:) -Y(DMUdValue,:)];
                b = [zeros(n,1);
                    0];

                %Inputs = 1
                Aeq = [X(DMUjValue,:) zeros(1,s)];
                beq = 1;

                %All variables >=0
                lb = zeros(1,(m+s));
                ub = repelem(10,(m+s));

                %options to surpress the "different algorithm" warning
                options = getDEAoptions(n);
                linprogoptions = options.optimopts;

                %Optimization
                [z, ~, exitflag, ~, dual] = linprog(f, A, b, Aeq, beq, lb, ub,linprogoptions);
                if z == 0
                    z = [dualX(DMUjValue,:) dualY(DMUjValue,:)]';
                end
                if exitflag == 0
                    z = [dualX(DMUjValue,:) dualY(DMUjValue,:)]';
                end
                if isempty(z)
                    z = [dualX(DMUjValue,:) dualY(DMUjValue,:)]';
                end
                peerAppraisalScores(DMUjValue,DMUdValue) = dot(z(m+1:m+s,1)',Y(DMUjValue,:))/dot(z(1:m,1)',X(DMUjValue,:));
            end
            alpha(DMUjValue,iteration+1) = mean(peerAppraisalScores(DMUjValue,:));
        end

        %Check whether to end loop
        if  sum(abs((alpha(:,iteration+1)-alpha(:,iteration))) >= 0.0001) == 0
            condition = false;
        end
        iteration = iteration +1;
    end
    
    
    %Initialize variables to store results into later
    agressivegeoEffIncl = zeros(n,1);
    agressivegeoEffExcl = zeros(n,1);
    agressiveariEffIncl = zeros(n,1);
    agressiveariEffExcl = zeros(n,1);
    
    
    %After termination, calculate efficiency scores
    for currentCalc = 1:n

        %Calculate Geo mean WITH own appraisal
        agressivegeoEffIncl(currentCalc,1) = geomean(peerAppraisalScores(currentCalc,:));
        %Calculate Geo mean WITHOUT own appraisal
        agressivegeoEffExcl(currentCalc,1) = geomean(peerAppraisalScores(currentCalc,[1:currentCalc-1,currentCalc+1:n]));
        %Calculate Arithmetic mean WITH own appraisal
        agressiveariEffIncl(currentCalc,1) = mean(peerAppraisalScores(currentCalc,:));
        %Calculate Arithmetic mean WITHOUT own appraisal
        agressiveariEffExcl(currentCalc,1) = mean(peerAppraisalScores(currentCalc,[1:currentCalc-1,currentCalc+1:n]));

    end

    out.CEPeerAppraisals =  peerAppraisalScores;
    out.CEEffArIn = agressiveariEffIncl;
    out.CEEffArEx = agressiveariEffExcl;
    out.CEEffGeIn = agressivegeoEffIncl;
    out.CEEffGeEx = agressivegeoEffExcl;
   
end

