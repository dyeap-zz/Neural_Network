classdef Flatten
    methods (Access = public)
        function output = compute(~,input)
            tempOutput = [];
            for currFilter=1:size(input,2)
                currInput = input{1,currFilter};
                temp = reshape(currInput,size(currInput,1)*size(currInput,2),size(currInput,3));
                tempOutput = [tempOutput;temp];
            end
            output = tempOutput;
            %{
                for k=1:cnn.layers{i-1}.no_featuremaps
                   ss =size(cnn.layers{i-1}.featuremaps{k});
                   ss(3) =size(cnn.layers{i-1}.featuremaps{k},3);
                   if cnn.input_image_width == 1
                       ss(3) =ss(2);
                       ss(2)=1;
                   end
                   zz =[zz; reshape(cnn.layers{i-1}.featuremaps{k}, ss(1)*ss(2), ss(3))];
                   
                end
            %}    
        end
    end
    properties
    end
end