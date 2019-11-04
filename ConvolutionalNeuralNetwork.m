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
        function obj = append_layer(obj,layer)
            % initialize all iterations
            num_col = size(obj.layer,2) +1;
            for curr_batch = 1:obj.num_epoch
                obj.layer(curr_batch,num_col) = {layer}; % append layer onto end  
            end
        end
        % For debugging
        function obj = run_one_layer(obj)
            curr_batch = obj.curr_iteration;
            curr_layer_number = obj.curr_layer_num;
            curr_layer = obj.layer{curr_batch,curr_layer_number};
            input = obj.layer{curr_batch,curr_layer_number}.io.get_input();
            if (curr_layer.name == 'c')
                % used for second convolution layer
                %{
                if (iscell(input))
                    numFilters = size(input,2);
                    tempInput = [];
                    for currFilter=1:numFilters
                        tempInput = cat(3,tempInput,input{1,currFilter});
                    end
                    input = tempInput;
                end
                %}
                num_filters = curr_layer.info.get_num_filters();
                temp_output = cell([1,num_filters]);
                % compute multiple filters
                for i=1:num_filters
                    filter = curr_layer.info.get_filter(i);
                    if (iscell(input))
                        % compute first one only to init szie
                        [obj,temp] = obj.compute_forward_convolution(input{1,1},curr_layer,filter);
                        output = zeros(size(temp));
                        for prevFilterNum=1:size(input,2)
                            tempInput = input{1,prevFilterNum};
                            [obj,tempOutput] = obj.compute_forward_convolution(tempInput,curr_layer,filter);
                            output = output + tempOutput;
                        end
                    else
                        [obj,output] = obj.compute_forward_convolution(input,curr_layer,filter);
                    end
                    temp_output{1,i} = output;
                end
            elseif(curr_layer.name == 'p')
                % must iterate through all filters used
                num_filters = size(input,2);
                temp_output = cell([1,num_filters]);
                for i=1:num_filters
                    [obj,output] = obj.compute_forward_pool(input{1,i},curr_layer);
                    temp_output{1,i} = output;
                end
            elseif(curr_layer.name == 'r')
                temp_output = obj.apply_relu(input);
            elseif(strcmp(curr_layer.name,'flat'))
                temp_output = curr_layer.info.compute(input);
            elseif(strcmp(curr_layer.name,'fc'))
                %{
                if (curr_layer.info.get_bool_finalFC())
                    numNodes = unique(cellClassifications);
                    obj.layer{curr_batch,curr_layer_number}.info.set_numNodes(numNodes);
                end
                %}
                [temp_output,fc_obj] = curr_layer.info.computeFCLayer(input);
                obj.layer{curr_batch,curr_layer_number}.info = fc_obj;
            else
                disp("Error have not created such layer")
            end
            obj.layer{curr_batch,curr_layer_number}.io = obj.layer{curr_batch,curr_layer_number}.io.set_output(temp_output);
            if (obj.curr_layer_num < size(obj.layer,2))
                obj.curr_layer_num =  obj.curr_layer_num + 1;
                curr_layer_number = obj.curr_layer_num;
                obj.layer{curr_batch,curr_layer_number}.io = obj.layer{curr_batch,curr_layer_number}.io.set_input(temp_output);              
            end
            end
        function obj = setup_begin_layer(obj,input)
            obj.layer{obj.curr_iteration,1}.io = obj.layer{obj.curr_iteration,1}.io.set_input(input);
        end
        function [obj,output] = compute_forward_convolution(obj,input,curr_layer,curr_filter)
            %filter_rep = repmat(curr_filter,[1,1,size(input,3)]);
            stride = curr_layer.info.get_stride();
            % To perform correlation must rotate filter
            rotFilter = rot90(curr_filter,2);
            temp = convn(input,rotFilter,'valid');
            output = temp(1:stride:end,1:stride:end,:);
        end
        function [obj,output] = compute_forward_pool(obj,input,curr_layer)
            % remove number of samples and save 
            % X -- output activations of the previous layer, numpy array of shape (n_H_prev, n_W_prev) assuming input channels = 1
            % W -- Weights, numpy array of size (f, f) assuming number of filters = 1
            % Returns:
            % H -- conv output, numpy array of size (n_H, n_W)
            % cache -- cache of values needed for conv_backward() function
            
            stride = curr_layer.info.get_stride();
            % create a temp pool filter only used only to get size
            filter = cell(curr_layer.info.get_filter_size());    
            
            filter_height = size(filter,1);
            filter_width = size(filter,2);
            
            % left side curr iteration + right side of plus symbol computes the number of iteration
            width_iter = 1 + floor((size(input,2) - filter_width)/stride);
            height_iter = 1 + floor((size(input,1) - filter_height)/stride);
            % initialize new dot product variable
            temp_output = [];
            % for b = 1
            %temp_conv_layer = zeros(height_iter,width_iter);
            %row = 1;
            for sampleNum = 1:size(input,3)
                sampleInput = input(:,:,sampleNum);
                temp_pool_layer = zeros(height_iter,width_iter);
                row = 1;
                for i=1:stride:height_iter*stride
                    temp_intensity = sampleInput(i:(i-1)+filter_height,:);
                    col = 1;
                    for j=1:stride:width_iter*stride
                        sub_image = temp_intensity(:,j:(j-1)+filter_width);
                        % compute max pool
                        comp = obj.max_pool(sub_image);

                        % if batch_size = 1
                        %temp_conv_layer(row,col) = {comp};
                        temp_pool_layer(row,col) = comp;
                        col = col + 1;
                    end
                    row = row + 1;
                end
                temp_output = cat(3,temp_output,temp_pool_layer);
            end
            output = temp_output;
        end
        function double = dot_product(~,sub_image,filter)
            dot_product = dot(sub_image,filter,1);
            double = sum(dot_product);
            % should not be used because need to retain batch size
            %if(size(double,3)>1)
            %    double = sum(double,3);
            %end
        end      
        function double = max_pool(~,sub_image)
            double = max(max(sub_image));
        end
        function output = apply_relu(obj,input)
            num_filters = size(input,2);
            output = cell([1,num_filters]);
            for i=1:num_filters
                curr_image = input{1,i};
                neg_index = curr_image < 0;
                if (nnz(neg_index))
                    curr_image(neg_index) = 0;
                end                
                output{1,i}= curr_image;
            end
        end
        %%% Used to compute one layer back propogation
        function obj = run_one_layer_bp(obj,batch_training_labels)
            lastLayerNum = size(obj.layer,2);
            curr_batch = obj.curr_iteration;
            curr_layer_number = obj.curr_layer_num;
            if(curr_layer_number==size(obj.layer,2))
                lastLayerOutput = obj.layer{curr_batch,lastLayerNum}.io.get_output();
                obj.label_error = lastLayerOutput - batch_training_labels; % only used for last label error
                error_sq = obj.label_error.^2;
                obj.loss = sum(error_sq(:))/(2*size(error_sq,2)); %loss over all examples
            end
            
            curr_layer = obj.layer{curr_batch,curr_layer_number};
            prev_layer = obj.layer{curr_batch,curr_layer_number-1};
            curr_output = curr_layer.io.get_output();
            prev_output = prev_layer.io.get_output();
            
            if(strcmp(curr_layer.get_name(),'fc') && curr_layer_number==size(obj.layer,2))
                % compute dw and db
                obj.layer{curr_batch,curr_layer_number}.error = curr_layer.error.computeSoftMax(obj.label_error,prev_output,curr_output);             
            elseif(strcmp(curr_layer.get_name(),'fc') && curr_layer_number<size(obj.layer,2))
                % need to test out tanh layer
                front_layer = obj.layer{curr_batch,curr_layer_number+1};
                tempError = front_layer.info.get_weights()' * front_layer.error.get_currLayerError();% 15x5 = 2x15* front_error 2x5
                %cnn.layers{i}.er{1} = ( (cnn.layers{i+1}.W)' * cnn.layers{i+1}.er{1} ); %temp error because it gets overridden by the following line of code
                obj.layer{curr_batch,curr_layer_number}.error = curr_layer.error.computeTanh(tempError,prev_output,curr_output);
                %obj = curr_layer.error.computeFCLayer(obj.label_error,prev_output,curr_output);
            end
            obj.curr_layer_num = curr_layer_number - 1;
        end
    end
    properties
        loss;
        label_error; % only used for the last layer error
        layer = {}; % Everything inside can only be of type layer of fcnn
        num_epoch;
        curr_layer_num = 1;
        curr_iteration = 1;
        
        %layerio; % rows represent iteration and the col represents the layer
    end   
end