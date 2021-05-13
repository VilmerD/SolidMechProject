classdef TOListener < handle
    
    properties
        statistics
    end
    
    methods
        function obj = TOListener()
            obj.statistics.g0 = [];
            obj.statistics.g1 = [];
            obj.statistics.designs = [];
            obj.statistics.factorizations = [];
            obj.statistics.lineqs = [];
        end
        
        function registerUpdate(obj, newStatistics)
            obj.statistics.g0(end + 1) = newStatistics.g0;
            obj.statistics.g1(end + 1) = newStatistics.g1;
            obj.statistics.designs(:, end+1) = newStatistics.design;
            obj.statistics.factorizations(end + 1) = newStatistics.factorizations;
            obj.statistics.lineqs(end + 1) = newStatistics.ncalls; 
        end
    end
end

