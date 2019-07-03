% Minimum Working Example of the “DEA Cross-efficiency toolbox for MATLAB”
%
% Copyright 2019 Christian Kaps, José L. Zofío
% https://github.com/joselzofio/DEACrossEfficiencyMATLAB
%
% Version: 1.0
% LAST UPDATE: 1 July, 2019

%Load sample data
load '../ExampleData/testData'
load '../ExampleData/testDataOrdinal'

%Call each model and store outputs
result1 = crossEff(Input12,Output12, "classic");
result2 = crossEff(Input12,Output12, "ratio");
result3 = crossEff(Input12,Output12, "multiplicative");
result4 = crossEff(Input12(1:20,:),Output12(1:20,:), "gameTheory"); %First 20 (out of 102) DMUs to speed up computation
result5 = crossEff(Input12Ordinal,Output12Ordinal, "ordinal");