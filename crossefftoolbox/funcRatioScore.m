function f  = funcRatioScore(inputweights)
    
    %Call global variables
    global funN funM funS funX funY funDMU funRatioScore;
    
    %Initialize ratio vector
    ratioScore = zeros(funN,1);
    
    %Loop through ratio scores
    for DMUcounter =(1:funN)
               
        ratioScore(DMUcounter,1)= sum(inputweights(1:funS).*funY(DMUcounter,:)) / sum(inputweights((funS+1):(funS+funM)).*funX(DMUcounter,:));
            
    end
    
    f = -(sum(ratioScore)-ratioScore(funDMU));
    
    %Define funRatioScore, which will need to be able to simple efficiency
    funRatioScore = ratioScore(funDMU);
        
end