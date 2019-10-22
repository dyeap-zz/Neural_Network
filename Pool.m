classdef Pool
    methods (Access = public)
        function obj = Pool(stride, filter_size)
            obj.stride = stride;
            obj.filter_size = filter_size;
        end
        function stride = get_stride(obj)
            stride = obj.stride;
        end
        function filter_size = get_filter_size(obj)
            filter_size = obj.filter_size;
        end
    end
    properties
        stride, filter_size
    end
end