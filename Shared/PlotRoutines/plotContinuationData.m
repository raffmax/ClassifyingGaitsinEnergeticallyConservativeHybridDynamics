function plotContinuationData(conData,continuationVector,PC_steps)
%plotContinuationData plots all available data in conData
if isfield(conData,'FloquetMultipliers')
    plotFloquetMultipliers(continuationVector,conData.FloquetMultipliers)
end

plotControlParameters(conData.xi,continuationVector)

plotTestFunctions(conData,continuationVector)

if nargin >2
    plotStepsInTangentSpace(continuationVector,PC_steps) 
end

if isfield(conData,'cost')
    if ~isempty(conData.cost)
        plotCostFunction(conData.cost,continuationVector)
    end
end

if ~isempty(conData.multiplier)
    plotLagrangeMultipliers(conData.multiplier,continuationVector)
end