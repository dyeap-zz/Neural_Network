classdef Convolution
    methods (Access = public)
        function obj = Convolution(stride,num_filters,filter_size)
            obj.stride = stride;
            obj.filters = obj.generate_conv_filter(num_filters,filter_size);
            obj.filter_size = filter_size;    
            obj.num_filters = num_filters;
        end
        function conv_filters = generate_conv_filter(~,num_filters,filter_size)
            conv_filters = cell(1,num_filters);
            for i = 1:num_filters
                conv_filters{1,i} = rand(filter_size);
            end 
        end
        function stride = get_stride(obj)
            stride = obj.stride;
        end
        function filter = get_filter(obj,index)
            filter = obj.filters{1,index};
        end
        function num_filters = get_num_filters(obj)
            num_filters = obj.num_filters;
        end
    end
    properties
        stride, filters, num_filters, filter_size;
    end
end