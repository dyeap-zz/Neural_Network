function Prototpying
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

global cellPlaylist cellData cellClassifications classifications cellCategories...
    cellCategoryInfo cellPreProcessing cellRawData cellSSAngleColorbar...
    numLV strBlank valCVMaxNeg valCVMaxPos valCVMinPos valCVMinNeg...
    valModelType valRTMaxNeg valRTMaxPos valRTMinPos valRTMinNeg...
    vecSSCurrAxes vecSSCurrShownIndices...

%{
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
%}



nn_intensity = double(nn_input.get_intensity(1));
for i=1:size(nn_input.get_sample_names,1)-1
    nn_intensity = cat(3,nn_intensity,double(nn_input.get_intensity(i+1)));
end
nn_intensity = reshape(nn_intensity,[100,100,1,size(nn_input.get_sample_names,1)]);
nn_label = categorical(cellClassifications);

maxEpochs = 2;
miniBatchSize = 1;

layers = [ ...
    imageInputLayer([100 100 1])
    convolution2dLayer(3,18) % default stride 1 output is 100 - 5 + 1
    reluLayer
    fullyConnectedLayer(size(unique(cellClassifications),1))
    softmaxLayer
    classificationLayer];

 func = @(info)grabInfo(info,3);
options = trainingOptions('sgdm', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',true, ...
    'VerboseFrequency',1,...
    'Plots','none',...
    'OutputFcn', func);% training-progress

%net = trainNetwork(X,Y,layers,options);
net = trainNetwork(nn_intensity,nn_label,layers,options);
display(net);
%{
channels = 1:6 % pick what filters you want
%{
I = deepDreamImage(net,'conv',channels)
figure
for i = 1:max(channels)
    subplot(5,5,i)
    imshow(I(:,:,:,i))
end
%}
features = activations(net,nn_intensity,'conv_2'); % 96,96,6,18(number of images)

for filter_num=1:6
    figure
for i = 1:18
    subplot(5,5,i)% number of grid spots
    imshow(features(:,:,filter_num,i))
end
end
    %}
end
