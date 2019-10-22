classdef LayerIO
    methods (Access = public)
        function obj = set_input(obj,input)
            obj.input = input;
        end
        function obj = set_output(obj,output)
            obj.output = output;
        end
    end
    properties
        input,output
    end
end