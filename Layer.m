classdef Layer
    properties
        % Needed for all layers
        name = 'c';
        num = 1;
        % used for convolution layer
        filters = 2; % size of filter will give you number of filter
        stride = 3;
        % may or may not need
        input;
        output;      
    end
end