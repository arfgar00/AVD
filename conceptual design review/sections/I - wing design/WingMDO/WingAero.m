classdef WingAero
    properties
        circulation      % circulation distribution on each stripy
        dyArray          % width of each strip
        y                % spanwise location of strip center obtained from stripy

    end

    methods
        function init(obj,wingGeometry)
            obj.y = zeros(length(wingGeometry.stripy) + 1,1);
            obj.y(1) = wingGeometry.stripy(1) - wingGeometry.dyArray(1);
            obj.y(2:end) = wingGeometry.stripy + wingGeometry.dyArray;
        end

        function plotLiftDistribution()
	        
        end
    end
end