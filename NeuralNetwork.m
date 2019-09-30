function NeuralNetwork
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
nn = NNInput(cellPlaylist(:,sample_names_col),...
                  cellData(:,compensation_voltage_col),...
                  cellData(:,retention_time_col),...
                  cellData(:,intensity_col));
              
           
cvn_fltr = [1,0,1;...
             0,1,0;...
             1,0,1];              
nn.set_convolution_filter(cvn_fltr);
stride = 1;
nn.set_cf_stride(stride);
nn.compute_convolution();



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
nn_intensity = double(nn_input.get_intensity(1));
for i=1:size(nn_input.get_sample_names,1)-1
    nn_intensity = cat(3,nn_intensity,double(nn_input.get_intensity(i+1)));
end
nn_intensity = reshape(nn_intensity,[100,100,1,size(nn_input.get_sample_names,1)]);
nn_label = [1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0];

maxEpochs = 100;
miniBatchSize = 27;

layers = [ ...
    imageInputLayer([100 100 1])
    convolution2dLayer(12,25)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];

options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',false, ...
    'Plots','training-progress');

%net = trainNetwork(X,Y,layers,options);
net = trainNetwork(nn_intensity,nn_label,layers,options);

A = rand(2,3,4);
B = rand(2,3,5);
C = cat(3,A,B);
disp('a')


end
