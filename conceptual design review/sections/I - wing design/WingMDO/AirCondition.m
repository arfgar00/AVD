classdef AirCondition
   properties
        M
        V
        h         %in meters
        rho
        T
        a
        P
        nu
        mu
        Re
   end
    
   methods
       function obj = init(obj,c)
            [obj.T, obj.a, obj.P, obj.rho, obj.nu, obj.mu] = atmosisa(obj.h);
            obj.V = obj.M * obj.a;
            obj.Re = obj.rho*obj.V*c/obj.mu;
       end

       function obj = calcM(obj)
           [obj.T, obj.a, obj.P, obj.rho, obj.nu, obj.mu] = atmosisa(obj.h);
            obj.M = obj.V/obj.a;
       end
   end
end