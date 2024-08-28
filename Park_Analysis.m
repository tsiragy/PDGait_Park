%% Clean the workspace
clear all 
close all 
clc

%%  load V3D and C3D data
[filepath, name, ext] = fileparts(uigetfile({'*.mat'}));
load(name); %load the specified data into the workspace from trial of interest




[Rck_data] = Rck_Analysis(v3d_data);
[HL_data] = HL_Analysis(v3d_data);
[ML_data] = ML_Analysis(v3d_data);


data.HL_data = HL_data;
data.ML_data = ML_data;
data.Rck_data = Rck_data;

savedir = 'C:\Users\tsira\Documents\MATLAB\PD Project\Park Data\Park Data Gait Event Corrected\Processed';
save(fullfile(savedir, name), 'data');





