function nTime   = getNumDomainTime(obj)
%getNumDomainTime counts domains of hybrid system

if strcmp(obj.integrationType,'TD')
    nTime = size(obj.sequence.FlowMapOrder,2);
else
    nTime = 0;
end

end

