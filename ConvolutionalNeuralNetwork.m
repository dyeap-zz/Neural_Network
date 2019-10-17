% if you're going to have multiple Neural Networks then you should use this
% it's a good idea to do this because each Neural Network will have
% different parameters

classdef ConvolutionalNeuralNetwork
    % methods accessible outside the class
    methods (Access = public)
        function this = append_layer_info(this,layer)
            this.layer_info(end+1) = {layer}; % append layer onto end
            
        end
    end
    properties
        layer_info = {}; % Everything inside can only be of type layer
        num_conv_layer;
    end   
end