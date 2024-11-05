classdef KF_bank
    properties
        x = {}; 
        P = {};
        probs = [];
    end
    methods
        function obj = KF_bank(this_x,this_P,pr)
            obj.x = this_x;
            obj.P = this_P;
            obj.probs = pr;
        end
        function obj = set.x(obj,val)
            obj.x = val;
        end
        function obj = set.P(obj,val)
            obj.P = val;
        end
        function obj = set.probs(obj,val)
            obj.probs = val;
        end
        function val = get.x(obj)
            val = obj.x;
        end
        function val = get.P(obj)
            val = obj.P;
        end
        function val = get.probs(obj)
            val = obj.probs;
        end
    end
    
end