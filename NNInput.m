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
        max_pool_filter;
        mp_stride;
        layer_transformation;
    end
    methods
        function obj = compute_convolution(obj)
            num_samples = size(obj.get_sample_names(),1);
            
            
        end
        function obj = set_convolution_filter(obj,filter)
            obj.convolution_filter = filter;
        end
        function obj = set_cf_stride(obj,stride)
            obj.cf_stride = stride;
        end
        function obj = set_max_pool_filter(obj,filter)
            obj.max_pool_filter = filter;
        end
        function obj = set_mp_stride(obj,stride)
            obj.mp_stride = stride;
        end
    end
end