classdef Layer
    methods (Access = public)
        % Constructor
        function obj = Layer(name,num,num_filters,stride,size_filter)
            obj.name = name;
            obj.num = num;
            if(name == 'c') % convolutional layer
                obj.stride  = stride;
                obj.filters = obj.generate_conv_filter(num_filters,size_filter);
                obj.filter_size = size_filter;
            elseif(name == 'p') % pool layer
                obj.stride  = stride;
                obj.num_filters = num_filters;
                obj.filter_size = size_filter;
            elseif(name == 'r')
                % dont need to do anything
            else
                disp("Layer is not yet defined")
            end
        end
        function conv_filters = generate_conv_filter(~,num_filters,filter_size)
            conv_filters = cell(1,num_filters);
            for i = 1:num_filters
                conv_filters{1,i} = rand(filter_size);
            end 
        end
    end
    properties
        % Needed for all layers
        name;
        num;
        % used for convolution layer
        num_filters; 
        stride;
        filters;
        filter_size;
        % may or may not need
        input;
        output;      
    end
end