classdef Error
    methods(Access = public)
        function [error,dw,b]= computeSoftMax(obj,prevOutput,currOutput)
            % calculate error on last layer
            if (bool_comp_error)
                    err1 = repmat(sum(error.*input, 1), [size(error,1) 1]);
                    tempOutput = -input.* (err1 -error);
            end
            cnn.layers{cnn.no_of_layers}.dW = cnn.layers{cnn.no_of_layers}.er{1} * ( cnn.layers{cnn.no_of_layers-1}.outputs)' / size(cnn.layers{cnn.no_of_layers}.er{1}, 2);
            cnn.layers{cnn.no_of_layers}.db = mean( cnn.layers{cnn.no_of_layers}.er{1}, 2);
        end
    end
    properties
    end
end