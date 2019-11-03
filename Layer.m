classdef Layer
    methods (Access = public)
        % Constructor
        function obj = Layer(name,num,num_filters,stride,filter_size,fc_numNodes,method,bool_finalFC)
            obj.name = name;
            obj.num = num;
            if(name == 'c') % convolutional layer
                obj.info = Convolution(stride,num_filters,filter_size);
            elseif(name == 'p') % pool layer
                obj.info = Pool(stride, filter_size);
            elseif(strcmp(name,'flat'))
                obj.info = Flatten();
            elseif(strcmp(name,'fc'))
                obj.info = FullyConnectedNeuralNetwork(fc_numNodes,method,bool_finalFC);
            elseif(name == 'r')
                % dont need to do anything
            else
                disp("Layer is not yet defined")
            end
        end
        function 
            % calculate error on last layer
            cnn.layers{cnn.no_of_layers}.dW = cnn.layers{cnn.no_of_layers}.er{1} * ( cnn.layers{cnn.no_of_layers-1}.outputs)' / size(cnn.layers{cnn.no_of_layers}.er{1}, 2);
            cnn.layers{cnn.no_of_layers}.db = mean( cnn.layers{cnn.no_of_layers}.er{1}, 2);
        end
        function obj = set_input(obj,input)
            obj.layer_io = LayerIO(input);
        end
    end
    properties
        % Needed for all layers
        name, num;
        % Define convolution layer
        info; % will either be a convolution, pool or relu
        io = LayerIO;
        error = Error;
    end
end