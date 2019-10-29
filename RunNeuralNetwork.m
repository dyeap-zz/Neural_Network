function RunNeuralNetwork
% Remove later
clear
clc
close all


% This is necessary to help get the data from the AnalyzeIMS workspace to
% the PredictGCDMS workspace. AnalyzeIMS will save .mat file in workspace
% and PredictGCDMS will load .mat file which requires all variables in the
% file to be predifened
global cellPlaylist cellData cellClassifications classifications cellCategories...
    cellCategoryInfo cellPreProcessing cellRawData cellSSAngleColorbar...
    numLV strBlank valCVMaxNeg valCVMaxPos valCVMinPos valCVMinNeg...
    valModelType valRTMaxNeg valRTMaxPos valRTMinPos valRTMinNeg...
    vecSSCurrAxes vecSSCurrShownIndices...

% Load the chemical data from AnalyzeIMS
sample_names_col = 2;
compensation_voltage_col = 1;
retention_time_col = 2;
intensity_col = 3;

load('but_hex_nn.mat');
nn_input = NNInput(cellPlaylist(:,sample_names_col),...
                  cellData(:,compensation_voltage_col),...
                  cellData(:,retention_time_col),...
                  cellData(:,intensity_col));

batch_size = 5;
num_epoch = 500;
num_conv_layers = 3;
num_samples = size(nn_input.get_sample_names(),1);
% Set up layer info
cnn = ConvolutionalNeuralNetwork(num_epoch,num_conv_layers);
% name,num,num_filters,stride,filters,size_filter
layer = Layer('c',1,4,2,[3,4]);
cnn = cnn.append_layer(layer);
layer = Layer('p',2,1,2,[2,2]);
cnn = cnn.append_layer(layer);
layer = Layer('r',3);
cnn = cnn.append_layer(layer);
% apply another layer
layer = Layer('c',1,4,2,[3,4]);
cnn = cnn.append_layer(layer);
layer = Layer('p',2,1,2,[2,2]);
cnn = cnn.append_layer(layer);
layer = Layer('r',3);
cnn = cnn.append_layer(layer);

% Try running one layer
for curr_epoch=1:num_epoch
    % go through all samples
    for curr_sample=1:num_samples
        % grab the number of samples for a batch
        % have for loop iterate 
        cnn = cnn.setup_begin_layer(nn_input.get_batch_intensity(curr_sample,curr_sample+batch_size-1));
        cnn = cnn.run_one_layer();
        cnn = cnn.run_one_layer();
        cnn = cnn.run_one_layer();
    end
end
cnn = cnn.append_layer(layer1);
cnn = cnn.append_layer(layer1);
% Initialize random weights for filter
num_filters = 3;
filter_size = [3,3];
num_conv_layers = 2;
nn = nn.set_convolution_filter(num_filters,filter_size,num_conv_layers);
stride = 1;
nn = nn.set_cf_stride(stride);

max_pool_width = 2;
max_pool_height = 2;
pool_stride = 2;
nn = nn.set_max_pool_filter_size([max_pool_width,max_pool_height]);
nn = nn.set_mp_stride(pool_stride);

for curr_epoch=1:num_epoch
    % go through all samples
    for curr_sample=1:num_samples
        nn = nn.set_curr_epoch(curr_epoch);
        curr_intensity = nn.get_intensity(curr_sample);
        bool_conv = 1;
        [nn,input_c] = nn.compute_forward_convolution_or_pool(curr_intensity,bool_conv);
        %nn = nn.apply_relu(); 
        bool_conv = 0;
        [nn,input_cp] = nn.compute_forward_convolution_or_pool(input_c,bool_conv);
        [nn,input_cpf] = nn.flatten(input_cp);
        
        
        
        
        
        
        
        
    end
end




%{
ax = figure
func_plot_graph(ax, nn_input.get_cv(1), nn_input.get_rt(1), nn_input.get_intensity(1),'bone');
ax2 = figure;
gray_image = mat2gray(nn_input.get_intensity(1));
func_plot_graph(ax2, nn_input.get_cv(1), nn_input.get_rt(1), gray_image,'bone');
%}
%z = cat(3,4,5)
%z = cat(3,double(nn_input.get_intensity(1)),double(nn_input.get_intensity(2)))
%z = double(nn_input.intensity)

%{
global cellPlaylist cellData cellClassifications classifications cellCategories...
    cellCategoryInfo cellPreProcessing cellRawData cellSSAngleColorbar...
    numLV strBlank valCVMaxNeg valCVMaxPos valCVMinPos valCVMinNeg...
    valModelType valRTMaxNeg valRTMaxPos valRTMinPos valRTMinNeg...
    vecSSCurrAxes vecSSCurrShownIndices...

% Load the chemical data from AnalyzeIMS
sample_names_col = 2;
compensation_voltage_col = 1;
retention_time_col = 2;
intensity_col = 3;

load('pent_but_hex_benz_mix.mat');

% Remove
row = [58,59,68,80,84,103,104,105];
cellData(row,:) = [];
cellPlaylist(row,:) = [];
cellClassifications(row,:) = [];

nn_input = NNInput(cellPlaylist(:,sample_names_col),...
                  cellData(:,compensation_voltage_col),...
                  cellData(:,retention_time_col),...
                  cellData(:,intensity_col));




nn_intensity = double(nn_input.get_intensity(1));
for i=1:size(nn_input.get_sample_names,1)-1
    nn_intensity = cat(3,nn_intensity,double(nn_input.get_intensity(i+1)));
end
nn_intensity = reshape(nn_intensity,[100,100,1,size(nn_input.get_sample_names,1)]);
nn_label = categorical(cellClassifications);

maxEpochs = 100;
miniBatchSize = 27;

layers = [ ...
    imageInputLayer([100 100 1])
    convolution2dLayer(12,25)
    reluLayer
    fullyConnectedLayer(size(unique(cellClassifications),1))
    softmaxLayer
    classificationLayer];

options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',false, ...
    'Plots','training-progress');

%net = trainNetwork(X,Y,layers,options);
net = trainNetwork(nn_intensity,nn_label,layers,options);
%}
%{
A = rand(2,3,4);
B = rand(2,3,5);
C = cat(3,A,B);
disp('a')
%}

end
