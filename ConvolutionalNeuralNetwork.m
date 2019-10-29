% if you're going to have multiple Neural Networks then you should use this
% it's a good idea to do this because each Neural Network will have
% different parameters

classdef ConvolutionalNeuralNetwork
    % methods accessible outside the class
    methods (Access = public)
        % Constructor
        function obj = ConvolutionalNeuralNetwork(num_epoch,num_conv_layers)
            obj.num_epoch = num_epoch;    
        end
        % set methods
        function obj = append_layer(obj,layer)
            % initialize all iterations
            num_col = size(obj.layer,2) +1;
            for curr_iter = 1:obj.num_epoch
                obj.layer(curr_iter,num_col) = {layer}; % append layer onto end  
            end
        end
        % For debugging
        function obj = run_one_layer(obj)
            curr_iter = obj.curr_iteration;
            curr_layer_number = obj.curr_layer_num;
            curr_layer = obj.layer{curr_iter,curr_layer_number};
            input = obj.layer{curr_iter,curr_layer_number}.io.get_input;
            if (curr_layer.name == 'c')
                num_filters = curr_layer.info.get_num_filters();
                temp_output = cell([1,num_filters]);
                % compute multiple filters
                for i=1:num_filters
                    filter = curr_layer.info.get_filter(i);
                    [obj,output] = obj.compute_forward_convolution_or_pool(input,curr_layer,filter);
                    temp_output{1,i} = output;
                end
            elseif(curr_layer.name == 'p')
                % must iterate through all filters used
                num_filters = size(input,2);
                temp_output = cell([1,num_filters]);
                for i=1:num_filters
                    [obj,output] = obj.compute_forward_convolution_or_pool(input{1,i},curr_layer);
                    temp_output{1,i} = output;
                end
            elseif(curr_layer.name == 'r')
                temp_output = obj.apply_relu(input);
            else
                disp("Error have not created such layer")
            end
            obj.layer{curr_iter,curr_layer_number}.io = obj.layer{curr_iter,curr_layer_number}.io.set_output(temp_output);
            obj.curr_layer_num =  obj.curr_layer_num + 1;
            curr_layer_number = obj.curr_layer_num;
            obj.layer{curr_iter,curr_layer_number}.io = obj.layer{curr_iter,curr_layer_number}.io.set_input(temp_output);
        end
        function obj = setup_begin_layer(obj,input)
            obj.layer{obj.curr_iteration,1}.io = obj.layer{obj.curr_iteration,1}.io.set_input(input);
        end
        function [obj,output] = compute_forward_convolution_or_pool(obj,input,curr_layer,curr_filter)
            % remove number of samples and save 
            % X -- output activations of the previous layer, numpy array of shape (n_H_prev, n_W_prev) assuming input channels = 1
            % W -- Weights, numpy array of size (f, f) assuming number of filters = 1
            % Returns:
            % H -- conv output, numpy array of size (n_H, n_W)
            % cache -- cache of values needed for conv_backward() function
            
            bool_conv = 0;
            if (nargin == 4)
                stride = curr_layer.info.get_stride();
                filter = curr_filter;
                bool_conv = 1;
            else
                stride = curr_layer.info.get_stride();
                % create a temp pool filter only used only to get size
                filter = cell(curr_layer.info.get_filter_size());
            end
            
            filter_height = size(filter,1);
            filter_width = size(filter,2);
            
            % left side curr iteration + right side of plus symbol computes the number of iteration
            width_iter = 1 + floor((size(input,2) - filter_width)/stride);
            height_iter = 1 + floor((size(input,1) - filter_height)/stride);
            % initialize new dot product variable
            temp_conv_layer = cell(height_iter,width_iter);
            % for b = 1
            %temp_conv_layer = zeros(height_iter,width_iter);
            row = 1;
            for i=1:stride:height_iter*stride
                temp_intensity = input(i:(i-1)+filter_height,:,:);
                col = 1;
                for j=1:stride:width_iter*stride
                    sub_image = temp_intensity(:,j:(j-1)+filter_width,:);
                    if (bool_conv == 1)
                        dot_filter = repmat(filter,[1,1,size(sub_image,3)]);
                        comp = obj.dot_product(sub_image,dot_filter);
                    else
                        comp = obj.max_pool(sub_image);
                    end
                    % if batch_size = 1
                    %temp_conv_layer(row,col) = {comp};
                    temp_conv_layer(row,col) = {comp};
                    col = col + 1;
                end
                row = row + 1;
            end
            output = temp_conv_layer;
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
            num_col = size(input,2);
            output = cell([1,num_col]);
            for i=1:num_col
                curr_image = input{1,i};
                neg_index = curr_image < 0;
                if (nnz(neg_index))
                    curr_image(neg_index) = 0;
                end                
                output{1,i}= curr_image;
            end
        end
    end
    properties
        layer = {}; % Everything inside can only be of type layer
        num_epoch;
        curr_layer_num = 1;
        curr_iteration = 1;
        %layerio; % rows represent iteration and the col represents the layer
    end   
end