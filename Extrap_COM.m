function [ML_xCOM, AP_xCOM, w, l] = Extrap_COM(x, y, z) %enter COM position, COM velocity, and Lateral Heel Marker Trajectory Data

%UNTITLED Extrapolated Center of Mass Function: xCOM
%   Computed from formula given in McAndrew et al. 2012
%   Outputs a participant's instantaneous extrapolated Center of Mass 
%   Values must be in meters 
%   Written by: Tarique Siragy, University of Ottawa Human Kinetics 

pCOM = x; %COM position data  
vCOM = y; %COM velocity data
BOS_Marker = z; %Marker trajectories for the base of support from heel marker
pCOM_length = length(pCOM); 
vCOM_length = length(vCOM);
 
if pCOM_length == vCOM_length; %check if both signals have the same size
else
    disp('Both signals must be equal in length. Please recheck!');
    return;
end

%% Calculate Inverted Pendulum Length

difference = minus(pCOM,BOS_Marker); %difference between COM position and BOS position
squared = difference.^2; %square the differences 
total_sum = sum(squared,2); % sum the differences 
Sqrt_dis = sqrt(total_sum);
%average_sum = mean(total_sum);
l = mean(Sqrt_dis); %length of the inverted pendulum in meters based on pythagorem theorem 
w = sqrt(9.81/l); %Calculate eigenfrequency of the inverted pendulum 

%% Parse the position and velocity of the COM

pCOM_ML = pCOM(:,1); %Center of Mass ML direction
pCOM_AP = pCOM(:,2); %Center of Mass AP direction
pCOM_V = pCOM(:,3); %Center of Mass V direction

vCOM_ML = vCOM(:,1); %Center of Mass velocity in ML direction
vCOM_AP = vCOM(:,2); %Center of Mass velocity in AP direction 
vCOM_V = vCOM(:,3); %Center of Mass velocty in V direction 

%% Calculate the Extrapolated Center of Mass Variable in AP and ML

AP_xCOM = pCOM_AP + (vCOM_AP/w); %Calculate xCOM in the AP direction
ML_xCOM = pCOM_ML + (vCOM_ML/w); %Calculate xCOM in the ML direction

end
