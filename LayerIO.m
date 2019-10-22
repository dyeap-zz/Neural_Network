classdef LayerIO
    methods (Access = public)
        function obj = set_input(obj,input)
            obj.input = input;
        end
        function obj = set_output(obj,output)
            obj.output = output;
        end
        function input = get_input(obj)
            input = obj.input;
        end
        function output = get_output(obj)
            output = obj.output;
        end
    end
    properties
        input,output
    end
end