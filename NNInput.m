% Naming conventions
% https://visualgit.readthedocs.io/en/latest/pages/naming_convention.html
% Make a class for the data or from PredictGCDMS send the data to this
% class. This class is specifically for Predict GCDMS. Can rename all the
% data here. Constructor can rename the classes. This will help reduce the
% number of variables in workspace and make it clear what variable is used
% for what. Possibly make another class that inherits from this class to.
% PredictGCDMS will create an object of this class. and will make it easier
% to read what is what.  

% Access set both setAccess and getAccess
% SetAccess relates only to changing the variable with . operator outside
% of class
% Get Access relates only to scope of getting the variables with . operator
% outside of class
classdef NNInput < AimsInput
    % Create instance variable and makes them private so that other
    % function cannot access and modify them. Only memeber functions can
    % access.
    % Later change to Access. SetAccess here so I can view what object has
    % in other methods
    properties
        convolution_filter;
        cf_stride;
        max_pool_filter_size;
        mp_stride;
        layer = {}
        %conv_layer = {};
        %pool = {};
        %relu = {};
        %layer_desc = struct("conv",{},"pool",{},"relu",{}...
        %                   );
        layer_transformation;
        curr_epoch;
    end
    methods
        function [obj,output] = compute_forward_convolution_or_pool(obj,input,bool_conv)
            % remove number of samples and save 
            % X -- output activations of the previous layer, numpy array of shape (n_H_prev, n_W_prev) assuming input channels = 1
            % W -- Weights, numpy array of size (f, f) assuming number of filters = 1
            % Returns:
            % H -- conv output, numpy array of size (n_H, n_W)
            % cache -- cache of values needed for conv_backward() function
            
            
            
            if (bool_conv == 1)
                stride = obj.cf_stride;
                filter = obj.convolution_filter;
            else
                stride = obj.mp_stride;
                % create a temp pool filter only used only to get size
                filter = cell(obj.max_pool_filter_size);
            end
            
            filter_height = size(filter,1);
            filter_width = size(filter,2);
            
            % left side curr iteration + right side of plus symbol computes the number of iteration
            width_iter = 1 + floor((size(input,2) - filter_width)/stride);
            height_iter = 1 + floor((size(input,1) - filter_height)/stride);
            % initialize new dot product variable
            temp_conv_layer = zeros(height_iter,width_iter);
            row = 1;
            for i=1:stride:height_iter*stride
                temp_intensity = input(i:(i-1)+filter_height,:);
                col = 1;
                for j=1:stride:width_iter*stride
                    sub_image = temp_intensity(:,j:(j-1)+filter_width);
                    if (bool_conv == 1)
                        comp = obj.dot_product(sub_image,filter);
                    else
                        comp = obj.max_pool(sub_image);
                    end
                    temp_conv_layer(row,col) = comp;
                    col = col + 1;
                end
                row = row + 1;
            end
            output = temp_conv_layer;
            obj.layer_transformation{obj.get_curr_epoch(),size(obj.layer_transformation,2)+1} = temp_conv_layer;
        end
        
        function double = dot_product(~,sub_image,filter)
            dot_product = dot(sub_image,filter,1);
            double = sum(dot_product);
        end
        
        function double = max_pool(~,sub_image)
            double = max(max(sub_image));
        end
        
        function [obj,flat_input] = flatten(obj,input)
            num_rows = size(input,1);
            num_cols = size(input,2);
            num_elements = num_rows*num_cols;
            flat_input = reshape(input,[num_elements,1]);
        end
        function obj = compute_forward_fully_connected(obj,input)
            
        end
        function obj = apply_relu(obj)
            num_samples = size(obj.get_sample_names(),1);
            last_layer = size(obj.layer_transformation,2);
            for curr_sample=1:num_samples
                curr_layer = obj.layer_transformation{curr_sample,last_layer};
                neg_index = curr_layer < 0;
                if (nnz(neg_index))
                    curr_layer(neg_index) = 0;
                    obj.layer_transformation{curr_sample,last_layer} = curr_layer;
                end                
            end
        end
        function obj = compute_forward_propagation(obj)
            
        end
        function obj = set_curr_epoch(obj,epoch)
            obj.curr_epoch = epoch;
        end
        function obj = set_convolution_filter(obj,num_filters,filter_size,num_conv_layers)
            for curr_layer=1:num_conv_layer
                for i = 1:num_filters
                    convolution_filter{i,1} = rand(filter_size);
                end
                obj.layer_desc.conv()
            end
        end
        function obj = set_cf_stride(obj,stride)
            obj.cf_stride = stride;
        end
        function obj = set_max_pool_filter_size(obj,filter)
            obj.max_pool_filter_size = filter;
        end
        function obj = set_mp_stride(obj,stride)
            obj.mp_stride = stride;
        end
        function mat_double = get_latest_layer(obj)
            mat_double = obj.layer_transformation{end,end};
        end
        function int = get_curr_epoch(obj)
            int = obj.curr_epoch;
        end
    end
end