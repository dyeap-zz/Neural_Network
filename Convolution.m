classdef Convolution
    methods (Access = public)
        function obj = Convolution(stride,num_filters,filter_size)
            obj.stride = stride;
            obj.filters = obj.generate_conv_filter(num_filters,filter_size);
            obj.filter_size = filter_size;       
        end
        function conv_filters = generate_conv_filter(~,num_filters,filter_size)
            conv_filters = cell(1,num_filters);
            for i = 1:num_filters
                conv_filters{1,i} = rand(filter_size);
            end 
        end
    end
    properties
        stride, filters, num_filters, filter_size;
    end
end