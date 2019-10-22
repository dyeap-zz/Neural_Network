% if you're going to have multiple Neural Networks then you should use this
% it's a good idea to do this because each Neural Network will have
% different parameters

classdef ConvolutionalNeuralNetwork
    % methods accessible outside the class
    methods (Access = public)
        % Constructor
        function obj = ConvolutionalNeuralNetwork(num_epoch)
            obj.num_epoch = num_epoch;
        end
        % set methods
        function obj = append_layer_info(obj,layer)
            obj.layer_info(end+1) = {layer}; % append layer onto end  
        end
        % For debugging
        function obj = run_one_layer(obj)
            curr_layer = obj.layer_info{1,obj.curr_layer_num};
            
        end
        function obj = setup_begin_layer(obj,input)
           obj.layer_info. = obj.layer_info.setup_input(input)
        end
    end
    properties
        layer_info = {}; % Everything inside can only be of type layer
        num_epoch;
        curr_layer_num = 1;
    end   
end