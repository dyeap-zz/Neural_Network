classdef Error
    methods(Access = public)
        function obj = computeSoftMax(obj,error,prevOutput,currOutput)
            % calculate error on last layer
            % applyactfunction
            err1 = repmat(sum(error.*currOutput, 1), [size(error,1) 1]);
            t_currLayerError = -currOutput.* (err1 - error); %2x5
            %
            t_dw = t_currLayerError * prevOutput' ./ size(t_currLayerError,2); % 2 x 15
            t_db = mean(t_currLayerError, 2);
            obj.currLayerError = t_currLayerError; % this is er{1}
            obj.dw = t_dw;
            obj.db = t_db;
            %cnn.layers{cnn.no_of_layers}.dW = cnn.layers{cnn.no_of_layers}.er{1} * ( cnn.layers{cnn.no_of_layers-1}.outputs)' / size(cnn.layers{cnn.no_of_layers}.er{1}, 2);
            %cnn.layers{cnn.no_of_layers}.db = mean( cnn.layers{cnn.no_of_layers}.er{1}, 2);
            % dw is 6x15 =  6x20 error after apply function * 15(15 from number of nodes in previous layer)x20 (prev_output)/ batch_size
            % cnn.layers{cnn.no_of_layers}.db = mean( cnn.layers{cnn.no_of_layers}.er{1}, 2);
            % db is 6x1
        end
        function obj = computeTanh(obj,error,prevOutput,currOutput)
            temp = 1-currOutput.*currOutput;
            t_currLayerError = temp.*error;
            
            t_dw = t_currLayerError * prevOutput' / size(t_currLayerError,2);
            t_db = mean(t_currLayerError, 2);
            
            % cnn.layers{i}.er{1} * ( cnn.layers{i-1}.outputs)' / size(cnn.layers{i}.er{1}, 2);
            obj.currLayerError = t_currLayerError; % this is er{1}
            obj.dw = t_dw;
            obj.db = t_db;
        end
        function currLayerError = get_currLayerError(obj)
            currLayerError = obj.currLayerError;            
        end
    end
    properties
        currLayerError; % This is er{1}
        dw;
        db;
    end
end