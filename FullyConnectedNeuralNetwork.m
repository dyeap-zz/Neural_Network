classdef FullyConnectedNeuralNetwork
    methods (Access = public)
        function obj = FullyConnectedNeuralNetwork(numNodes,method,bool_finalFC)
            obj.bool_finalFC = bool_finalFC;
            obj.numNodes = numNodes;
            obj.activation = ActivationFunction(method);
        end
        function output = computeFCLayer(obj,input)
            %cnn.layers{l}.W =0.5*rand([no_of_nodes cnn.layers{l}.no_of_inputs]) -0.25;
            %cnn.layers{l}.b = 0.5*rand([no_of_nodes 1]) - 0.25;
            % setup the weights and b for computation
            numNodesInput = size(input,1);
            
            obj.weight = randi([-99,99],obj.numNodes,numNodesInput)./100;
            obj.b = randi([-99,99],obj.numNodes,1)./100;
            batch_size = size(input,2);
            obj.b = repmat(obj.b,[1,batch_size]);
            Y = obj.weight*input + obj.b;
            tempOutput = obj.activation.compute(Y);
            output = tempOutput;
        end
        function bool = get_bool_finalFC(obj)
            bool = obj.bool_finalFC;
        end
        function obj = set_numNodes(obj,numNodes)
            obj.numNodes = numNodes;
        end
    end
    properties
        activation;
        bool_finalFC;
        numNodes;
        weight; %numNodes*batch_size
        b;
    end 
end