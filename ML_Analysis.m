function [ML_data] = ML_Analysis(v3d_data)

%% Isolate Terrains 
Rock_ST = v3d_data.Rock_ST; %Rocky terrain start 
Rock_EN = v3d_data.Rock_EN; %Rocky terrain end
HL_ST = v3d_data.HL_ST; %Hilly terrain start 
HL_EN = v3d_data.HL_EN; %Hilly terrain end 
ML_ST = v3d_data.ML_ST; %Mediolateral start
ML_EN = v3d_data.ML_EN; %Mediolateral end

%% Hill Terrain Analysis %%
ML_t = minus(ML_EN, ML_ST); %Rocky terrain time duration in seconds
ML_Rind = round(v3d_data.ML_indR./10); %Right HS indices during rocky terrain
ML_Lind = round(v3d_data.ML_indL./10); %Left HS indices during rocky terrain 
R_StepT = (ML_Rind./100); %Right step times during rocky terrain
L_StepT = (ML_Lind./100); %Left step times during rocky terrain
Step_Times = [R_StepT; L_StepT]; %All right and left step times during rocky terrain 

%% Chronologically order the feet 
Rtemp_foot = v3d_data.RFootPos{:,:};
RFoot_loc = Rtemp_foot(ML_Rind,2);
Ltemp_foot = v3d_data.LFootPos{:,:};
LFoot_loc = Ltemp_foot(ML_Lind,2);
Step_loc = [RFoot_loc; LFoot_loc]; %Array of Right and Left Foot Locations 
Step_data = [Step_loc,Step_Times]; %Array of both Step Location (both feet) and Step Times (both Feet), first 125 in each column is right followed by left foot data
Ordered_Steps = sortrows(Step_data,2); %Chronologically order the right and left step times, right foot is odd numbers and left foot is even numbers 
[~,idx] = unique(Ordered_Steps(:,2)); %Find Duplicate steps 
Final_steps = Ordered_Steps(idx,:); 
%% Isolate the markers of interest at Right Heelstrike
RLHL = cell2mat(v3d_data.c3d_rlhl); %raw right lateral heel trajectories
RLHL_HS = RLHL(ML_Rind,:); %raw right lateral heel trajectories at heel strikes 
RHEEL = cell2mat(v3d_data.c3d_rheel); %raw right heel trajectories 
RHEEL_HS = RHEEL(ML_Rind,:); %raw right heel trajectories at heel strikes  
LLHL = cell2mat(v3d_data.c3d_llhl); %raw left lateral heel trajectories 
LLHL_HS = LLHL(ML_Lind,:); %raw left lateral heel trajectories at heel strikes
LHEEL = cell2mat(v3d_data.c3d_lheel); %raw left heel trajectories 
LHEEL_HS = LHEEL(ML_Lind,:); %raw left heel trajectories at heel strikes
HS_Foot_data_r = {RHEEL_HS, RLHL_HS, LHEEL_HS, LLHL_HS};


HS_Foot_data = cell(1,length(HS_Foot_data_r));
for i = 1:length(HS_Foot_data_r);
HS_Foot_data{:,i} = HS_Foot_data_r{i}(:,:)./1000; %combine all foot markers into a single array
end 
 

%% Create Butterworth Filter and Filter the data of interest 
Fe2 = 100; %sampling frequency
Fc2 = 12; %cut off frequency 
N2 = 2; %filter order
[B2,A2]=butter(N2,Fc2*2/Fe2,'low'); %filter parameter     

Filtered_MarkerData = cell(1,length(HS_Foot_data));

for i = 1:length(HS_Foot_data);
Filtered_MarkerData{:,i} = filtfilt(B2,A2,HS_Foot_data{i}(:,:)); %filter all foot trajectory data 
end


%% Parse the Filtered Data into their own variables for RHS 
Filtered_RHEEL = Filtered_MarkerData{1}; %XYZ trajectories of filtered right heel marker at right heel strike
Filtered_RLHL = Filtered_MarkerData{2}; %XYZ trajectories of filtered right lateral heel marker at right heel strike
Filtered_LHEEL = Filtered_MarkerData{3}; %XYZ trajectories of filtered left heel marker at right heel strike 
Filtered_LLHL = Filtered_MarkerData{4}; %XYZ trajectories of filtered left lateral heel marker at right heel strike 


RLHL_ML = Filtered_RLHL(:,1).*-1;
RLHL_AP = Filtered_RLHL(:,2).*-1;
RLHL_V = Filtered_RLHL(:,3);
RLHL_RHS = [RLHL_ML, RLHL_AP, RLHL_V];

RHEEL_ML = Filtered_RHEEL(:,1).*-1;
RHEEL_AP = Filtered_RHEEL(:,2).*-1;
RHEEL_V = Filtered_RHEEL(:,3);
RHEEL_RHS = [RHEEL_ML, RHEEL_AP, RHEEL_V];


LLHL_ML = Filtered_LLHL(:,1).*-1;
LLHL_AP = Filtered_LLHL(:,2).*-1;
LLHL_V = Filtered_LLHL(:,3);
LLHL_LHS = [LLHL_ML, LLHL_AP, LLHL_V];

LHEEL_ML = Filtered_LHEEL(:,1).*-1;
LHEEL_AP = Filtered_LHEEL(:,2).*-1;
LHEEL_V = Filtered_LHEEL(:,3);
LHEEL_LHS = [LHEEL_ML, LHEEL_AP, LHEEL_V];



%% Read in Center of Mass data (already filtered from V3D) 
COG_indeces = round(Final_steps(:,2).*100);
pCOG = v3d_data.COG{:,:}; %XYZ position of Center of Mass 
vCOG = v3d_data.COG_Velocity{:,:};%XYZ velocities of Center of Mass 
pCOG_RHS = pCOG(ML_Rind,:); %Center of Mass position at Right Heel Strikes 
pCOG_LHS = pCOG(ML_Lind,:); %Center of Mass position at Right Heel Strikes 
COG_Rind = ML_Rind - 1; %find the frames prior to right heel-strike to calculate instantenous acceleration 
COG_Lind = ML_Lind - 1; %find the frames prior to the left heel-strike to calculate instantaneous acceleration 
vCOG_prRHS = vCOG(COG_Rind,:);
vCOG_prLHS = vCOG(COG_Lind,:);
vCOG_RHS = vCOG(ML_Rind,:); %Center of Mass velocity at Right Heel Strikes 
vCOG_LHS = vCOG(ML_Lind,:); %Center of Mass velocity at Right Heel Strikes 


%% Calculate COG Acceleration
time = 1/100; %create time variable to calculate acceleration from velocity data 
vCOG_ML_Raw = vCOG(:,1);
vCOG_AP_Raw = vCOG(:,2);
vCOG_vert_Raw = vCOG(:,3); 
aCOG_Raw = diff(vCOG)./time;
aCOG_Indeces = COG_indeces - 1; %Indeces for COG acceleration heel-strikes, odd is right HS and even is left HS 
aCOG_ML = aCOG_Raw(aCOG_Indeces(1,1):aCOG_Indeces(end,1),1); %ML COG acceleration 
aCOG_AP = aCOG_Raw(aCOG_Indeces(1,1):aCOG_Indeces(end,1),2); %AP COG acceleration 
aCOG_V = aCOG_Raw(aCOG_Indeces(1,1):aCOG_Indeces(end,1),3); %Vert COG acceleration 
aCOG_RHS = minus(vCOG_RHS,vCOG_prRHS)./time; %calculate COG acclereation at right heel strikes 
aCOG_LHS = minus(vCOG_LHS,vCOG_prLHS)./time; %calculate COG acceleration at left heel strike


%% Angular Momentum 
Ang_Mom = v3d_data.ModelAngMmntm{:,:};
indext = COG_indeces(1:2:end);
idx =  cell(1,length(indext)-1);

for i = 1:(length(indext)-1);
    idx{:,i} = indext(i:i+1); 
end 


AngMomt =  cell(1,length(indext)-1);

for k = 1:length(AngMomt);
    AngMomt{:,k} = Ang_Mom(idx{k}(1):idx{k}(2),:);
end 



AngMmt_ML = cell(1,length(AngMomt)); %create blank cells to store data of each gait cycle in full (i.e. 2 steps or 1 full stride)
AngMmt_MLavgt = cell(1,length(AngMomt));
AngMmt_MLpkt = cell(1,length(AngMomt));
AngMmt_AP = cell(1,length(AngMomt));
AngMmt_APavgt = cell(1,length(AngMomt));
AngMmt_APpkt = cell(1,length(AngMomt));
AngMmt_V = cell(1,length(AngMomt)); 
AngMmt_Vavgt = cell(1,length(AngMomt));
AngMmt_Vpkt = cell(1,length(AngMomt));


for i = 1:length(AngMomt);
    AngMmt_ML{:,i} = AngMomt{i}(:,1);
    AngMmt_MLavgt{:,i} = mean(AngMmt_ML{i});
    AngMmt_MLpkt{:,i} = max(AngMmt_ML{i});
    AngMmt_AP{:,i} = AngMomt{i}(:,2);
    AngMmt_APavgt{:,i} = mean(AngMmt_AP{i});
    AngMmt_APpkt{:,i} = max(AngMmt_AP{i});
    AngMmt_V{:,i} = AngMomt{i}(:,3);
    AngMmt_Vavgt{:,i} = mean(AngMmt_V{i});
    AngMmt_Vpkt{:,i} = max(AngMmt_V{i});
end 

AngMmt_MLavg = abs(mean(cell2mat(AngMmt_MLavgt))); %average ML ang momentum throughout the gait trial 
AngMmt_MLstd = std(cell2mat(AngMmt_MLavgt)); %standard deviation of average ML ang momentum throughout the gait trial
AngMmt_MLpk = abs(mean(cell2mat(AngMmt_MLpkt))); %average ML ang momentum peak value 
AngMmt_MLpk_SD = std(cell2mat(AngMmt_MLpkt));  %standard deviation of ML ang momentum peak value
AngMmt_APavg = abs(mean(cell2mat(AngMmt_APavgt))); %average AP ang momentum throughout the gait trial
AngMmt_APstd = std(cell2mat(AngMmt_APavgt)); %standard deviation of average AP ang momentum throughout the gait trial
AngMmt_APpk = abs(mean(cell2mat(AngMmt_APpkt))); %average AP ang momentum peak value 
AngMmt_APpk_SD = std(cell2mat(AngMmt_APpkt)); %standard deviation of AP ang momentum peak value
AngMmt_Vavg = abs(mean(cell2mat(AngMmt_Vavgt))); %aveage Vrt ang momentum throughout the gait trial 
AngMmt_Vstd = std(cell2mat(AngMmt_Vavgt)); %standard deviation of average Vrt ang momentum throughout the gait trial 
AngMmt_Vpk = abs(mean(cell2mat(AngMmt_Vpkt))); %average Vertical ang momentum peak value 
AngMmt_Vpk_SD = std(cell2mat(AngMmt_Vpkt)); %standard deviation of Vertical ang momentum peak value

%% Trunk Linear Velocity
Trk_LV_R = v3d_data.TrunkLinVel{:,:}; %raw trunk linear velocity data

indext = COG_indeces(1:2:end);
idx =  cell(1,length(indext)-1);

for i = 1:(length(indext)-1);
    idx{:,i} = indext(i:i+1); 
end 


Trk_LV =  cell(1,length(indext)-1);

for k = 1:length(Trk_LV);
    Trk_LV{:,k} = Trk_LV_R(idx{k}(1):idx{k}(2),:);
end 





TrkLV_ML = cell(1,length(Trk_LV)); %create blank cells to store data of each gait cycle in full (i.e. 2 steps or 1 full stride)
TrkLV_MLavgt = cell(1,length(Trk_LV));
TrkLV_MLpkt = cell(1,length(Trk_LV));
TrkLV_AP = cell(1,length(Trk_LV));
TrkLV_APavgt = cell(1,length(Trk_LV));
TrkLV_APpkt = cell(1,length(Trk_LV));
TrkLV_V = cell(1,length(Trk_LV)); 
TrkLV_Vavgt = cell(1,length(Trk_LV));
TrkLV_Vpkt = cell(1,length(Trk_LV));


for i = 1:length(Trk_LV);
    TrkLV_ML{:,i} = Trk_LV{i}(:,1);
    TrkLV_MLavgt{:,i} = mean(TrkLV_ML{i});
    TrkLV_MLpkt{:,i} = max(TrkLV_ML{i});
    TrkLV_AP{:,i} = Trk_LV{i}(:,2);
    TrkLV_APavgt{:,i} = mean(TrkLV_AP{i});
    TrkLV_APpkt{:,i} = max(TrkLV_AP{i});
    TrkLV_V{:,i} = Trk_LV{i}(:,3);
    TrkLV_Vavgt{:,i} = mean(TrkLV_V{i});
    TrkLV_Vpkt{:,i} = max( TrkLV_V{i});
end 


TrkLV_MLavg = abs(mean(cell2mat(TrkLV_MLavgt))); %average ML trunk linear velocity throughout the gait trial 
TrkLV_MLstd = std(cell2mat(TrkLV_MLavgt)); %standard deviation of average ML trunk linear velocity throughout the gait trial
TrkLV_MLpk = abs(mean(cell2mat(TrkLV_MLpkt))); %average ML trunk linear velocity peak value 
TrkLV_MLpk_SD = std(cell2mat(TrkLV_MLpkt));  %standard deviation of ML trunk linear velocity peak value
TrkLV_APavg = abs(mean(cell2mat(TrkLV_APavgt))); %average AP trunk linear velocity throughout the gait trial
TrkLV_APstd = std(cell2mat(TrkLV_APavgt)); %standard deviation of average AP linear velocity throughout the gait trial
TrkLV_APpk = abs(mean(cell2mat(TrkLV_APpkt))); %average AP linear velocity peak value 
TrkLV_APpk_SD = std(cell2mat(TrkLV_APpkt)); %standard deviation of AP linear velocity peak value
TrkLV_Vavg = abs(mean(cell2mat(TrkLV_Vavgt))); %aveage Vrt linear velocity throughout the gait trial 
TrkLV_Vstd = std(cell2mat(TrkLV_Vavgt)); %standard deviation of average Vrt linear velocity throughout the gait trial 
TrkLV_Vpk = abs(mean(cell2mat(TrkLV_Vpkt))); %average Vertical linear velocity peak value 
TrkLV_Vpk_SD = std(cell2mat(TrkLV_Vpkt)); %standard deviation of Vertical linear velocity peak value



%% Trunk Angular Velocity
Trk_AV_R = v3d_data.TrunkAngVel{:,:}; %raw trunk angular velocity data

indext = COG_indeces(1:2:end);
idx =  cell(1,length(indext)-1);

for i = 1:(length(indext)-1);
    idx{:,i} = indext(i:i+1); 
end 


Trk_AV =  cell(1,length(indext)-1);

for k = 1:length(Trk_AV);
    Trk_AV{:,k} = Trk_AV_R(idx{k}(1):idx{k}(2),:);
end 

TrkAV_ML = cell(1,length(Trk_AV)); %create blank cells to store data of each gait cycle in full (i.e. 2 steps or 1 full stride)
TrkAV_MLavgt = cell(1,length(Trk_AV));
TrkAV_MLpkt = cell(1,length(Trk_AV));
TrkAV_AP = cell(1,length(Trk_AV));
TrkAV_APavgt = cell(1,length(Trk_AV));
TrkAV_APpkt = cell(1,length(Trk_AV));
TrkAV_V = cell(1,length(Trk_AV)); 
TrkAV_Vavgt = cell(1,length(Trk_AV));
TrkAV_Vpkt = cell(1,length(Trk_AV));


for i = 1:length(Trk_AV);
    TrkAV_ML{:,i} = Trk_AV{i}(:,1);
    TrkAV_MLavgt{:,i} = mean(TrkAV_ML{i});
    TrkAV_MLpkt{:,i} = max(TrkAV_ML{i});
    TrkAV_AP{:,i} = Trk_AV{i}(:,2);
    TrkAV_APavgt{:,i} = mean(TrkAV_AP{i});
    TrkAV_APpkt{:,i} = max(TrkAV_AP{i});
    TrkAV_V{:,i} = Trk_AV{i}(:,3);
    TrkAV_Vavgt{:,i} = mean(TrkAV_V{i});
    TrkAV_Vpkt{:,i} = max( TrkAV_V{i});
end 


TrkAV_MLavg = abs(mean(cell2mat(TrkAV_MLavgt))); %average ML trunk linear velocity throughout the gait trial 
TrkAV_MLstd = std(cell2mat(TrkAV_MLavgt)); %standard deviation of average ML trunk linear velocity throughout the gait trial
TrkAV_MLpk = abs(mean(cell2mat(TrkAV_MLpkt))); %average ML trunk linear velocity peak value 
TrkAV_MLpk_SD = std(cell2mat(TrkAV_MLpkt));  %standard deviation of ML trunk linear velocity peak value
TrkAV_APavg = abs(mean(cell2mat(TrkAV_APavgt))); %average AP trunk linear velocity throughout the gait trial
TrkAV_APstd = std(cell2mat(TrkAV_APavgt)); %standard deviation of average AP linear velocity throughout the gait trial
TrkAV_APpk = abs(mean(cell2mat(TrkAV_APpkt))); %average AP linear velocity peak value 
TrkAV_APpk_SD = std(cell2mat(TrkAV_APpkt)); %standard deviation of AP linear velocity peak value
TrkAV_Vavg = abs(mean(cell2mat(TrkAV_Vavgt))); %aveage Vrt linear velocity throughout the gait trial 
TrkAV_Vstd = std(cell2mat(TrkAV_Vavgt)); %standard deviation of average Vrt linear velocity throughout the gait trial 
TrkAV_Vpk = abs(mean(cell2mat(TrkAV_Vpkt))); %average Vertical linear velocity peak value 
TrkAV_Vpk_SD = std(cell2mat(TrkAV_Vpkt)); %standard deviation of Vertical linear velocity peak value

%% Left Arm Motion
Larm_raw = v3d_data.LShoulder{:,:}; %call trunk angular velocity data
indext = COG_indeces(1:2:end);
idx =  cell(1,length(indext)-1);

for i = 1:(length(indext)-1);
    idx{:,i} = indext(i:i+1); 
end 



Larm =  cell(1,length(indext)-1);

for k = 1:length(Larm);
    Larm{:,k} = Larm_raw(idx{k}(1):idx{k}(2),:);
end 

Larm_ML = cell(1,length(Larm)); %create blank cells to store data of each gait cycle in full (i.e. 2 steps or 1 full stride)
Larm_MLmaxt = cell(1,length(Larm));
Larm_MLmint = cell(1,length(Larm));
Larm_AP = cell(1,length(Larm));
Larm_APmaxt = cell(1,length(Larm));
Larm_APmint = cell(1,length(Larm));
Larm_V = cell(1,length(Larm)); 
Larm_Vmaxt = cell(1,length(Larm));
Larm_Vmint = cell(1,length(Larm)); 

for i = 1:length(Larm);
    Larm_AP{:,i} = Larm{i}(:,1);
    Larm_APmaxt{:,i} = max(Larm_AP{i});
    Larm_APmint{:,i} = min(Larm_AP{i});
    Larm_ML{:,i} = Larm{i}(:,2);
    Larm_MLmaxt{:,i} = max(Larm_ML{i});
    Larm_MLmint{:,i} = min(Larm_ML{i});
    Larm_V{:,i} = Larm{i}(:,3);
    Larm_Vmaxt{:,i} = max(Larm_V{i});
    Larm_Vmint{:,i} = min(Larm_V{i});
end 

Larm_MLmax = mean(cell2mat(Larm_MLmaxt));
Larm_MLmin = mean(cell2mat(Larm_MLmint));
Larm_APmax = mean(cell2mat(Larm_APmaxt));
Larm_APmin = mean(cell2mat(Larm_APmint));
Larm_Vmax = mean(cell2mat(Larm_Vmaxt));
Larm_Vmin = mean(cell2mat(Larm_Vmint));


Larm_ROM = Larm_APmax - Larm_APmin; 


%% Right Arm Motion

Rarm_raw = v3d_data.RShoulder{:,:}; %call trunk angular velocity data
indext = COG_indeces(1:2:end);
idx =  cell(1,length(indext)-1);

for i = 1:(length(indext)-1);
    idx{:,i} = indext(i:i+1); 
end



Rarm =  cell(1,length(indext)-1);

for k = 1:length(Rarm);
    Rarm{:,k} = Rarm_raw(idx{k}(1):idx{k}(2),:);
end

Rarm_ML = cell(1,length(Rarm)); %create blank cells to store data of each gait cycle in full (i.e. 2 steps or 1 full stride)
Rarm_MLmaxt = cell(1,length(Rarm));
Rarm_MLmint = cell(1,length(Rarm));
Rarm_AP = cell(1,length(Rarm));
Rarm_APmaxt = cell(1,length(Rarm));
Rarm_APmint = cell(1,length(Rarm));
Rarm_V = cell(1,length(Rarm)); 
Rarm_Vmaxt = cell(1,length(Rarm));
Rarm_Vmint = cell(1,length(Rarm)); 

for i = 1:length(Rarm);
    Rarm_AP{:,i} = Rarm{i}(:,1);
    Rarm_APmaxt{:,i} = max(Rarm_AP{i});
    Rarm_APmint{:,i} = min(Rarm_AP{i});
    Rarm_ML{:,i} = Rarm{i}(:,2);
    Rarm_MLmaxt{:,i} = max(Rarm_ML{i});
    Rarm_MLmint{:,i} = min(Rarm_ML{i});
    Rarm_V{:,i} = Rarm{i}(:,3);
    Rarm_Vmaxt{:,i} = max(Rarm_V{i});
    Rarm_Vmint{:,i} = min(Rarm_V{i});
end 

Rarm_MLmax = mean(cell2mat(Rarm_MLmaxt));
Rarm_MLmin = mean(cell2mat(Rarm_MLmint));
Rarm_APmax = mean(cell2mat(Rarm_APmaxt));
Rarm_APmin = mean(cell2mat(Rarm_APmint));
Rarm_Vmax = mean(cell2mat(Rarm_Vmaxt));
Rarm_Vmin = mean(cell2mat(Rarm_Vmint));


Rarm_ROM = Rarm_APmax - Rarm_APmin; 



%% Calculate the Margin of Stability in AP on Right Heel Strike(xCOM)
[ML_xCOM, AP_xCOM, w, l] = Extrap_COM(pCOG_RHS, vCOG_RHS,RLHL_RHS); %Call xCOM function 
temp_MOS_RAP = RHEEL_AP - AP_xCOM; %Calculate Anteroposterior Margin of Stability
MOS_avg_RAP = mean(temp_MOS_RAP); %temporary avg MOS AP 
MOS_std_RAP = std(temp_MOS_RAP); %temporary std MOS AP 
H_MOS_RAP = MOS_avg_RAP + 2*MOS_std_RAP;
L_MOS_RAP = MOS_avg_RAP - 2*MOS_std_RAP;
MOS_AP_Rindeces = find(temp_MOS_RAP>L_MOS_RAP & temp_MOS_RAP<H_MOS_RAP);
MOS_RAP = temp_MOS_RAP(MOS_AP_Rindeces);
MOS_avg_RAP = mean(MOS_RAP);
MOS_std_RAP = std(MOS_RAP);


%% Calculate the Margin of Stability in ML at Right Heel Strike (xCOM) 
temp_MOS_RML = RLHL_ML -ML_xCOM;%Calculate Mediolateral Margin of Stability
MOS_avg_RML = mean(temp_MOS_RML); %temporary avg MOS ML 
MOS_std_RML = std(temp_MOS_RML); %temporary std MOS ML 
H_MOS_RML = MOS_avg_RML + 2*MOS_std_RML; %data above 2std 
L_MOS_RML = MOS_avg_RML - 2*MOS_std_RML; %data below 2std 
MOS_ML_Rindeces = find(temp_MOS_RML>L_MOS_RML & temp_MOS_RML<H_MOS_RML); %find data after removing outliers 
MOS_RML = temp_MOS_RML(MOS_ML_Rindeces);
MOS_avg_RML = mean(MOS_RML);
MOS_std_RML = std(MOS_RML);
MOS_Rdata = [MOS_avg_RAP, MOS_avg_RML]; %store MOS_AP (first column) and MOS_ML (second column) in single matrix to save as a variable later

%% Calculate the Margin of Stability in AP on Left Heel Strike(xCOM)
[ML_xCOM, AP_xCOM, w, l] = Extrap_COM(pCOG_LHS, vCOG_LHS,LLHL_ML); %Call xCOM function 
temp_MOS_LAP = LHEEL_AP-AP_xCOM; %Calculate Anteroposterior Margin of Stability
MOS_avg_LAP = mean(temp_MOS_LAP); %temporary avg MOS AP 
MOS_std_LAP = std(temp_MOS_LAP); %temporary std MOS AP 
H_MOS_LAP = MOS_avg_LAP + 2*MOS_std_LAP;
L_MOS_LAP = MOS_avg_LAP - 2*MOS_std_LAP;
MOS_AP_Lindeces = find(temp_MOS_LAP>L_MOS_LAP & temp_MOS_LAP<H_MOS_LAP);
MOS_LAP = temp_MOS_LAP(MOS_AP_Lindeces);
MOS_avg_LAP = mean(MOS_LAP);
MOS_std_LAP = std(MOS_LAP);


%% Calculate the Margin of Stability in ML at Left Heel Strike (xCOM) 
temp_MOS_LML = LLHL_ML -ML_xCOM;%Calculate Mediolateral Margin of Stability
MOS_avg_LML = mean(temp_MOS_LML); %temporary avg MOS ML 
MOS_std_LML = std(temp_MOS_LML); %temporary std MOS ML 
H_MOS_LML = MOS_avg_LML + 2*MOS_std_LML; %data above 2std 
L_MOS_LML = MOS_avg_LML - 2*MOS_std_LML; %data below 2std 
MOS_ML_Lindeces = find(temp_MOS_LML>L_MOS_LML & temp_MOS_LML<H_MOS_LML); %find data after removing outliers 
MOS_LML = temp_MOS_LML(MOS_ML_Lindeces);
MOS_avg_LML = mean(MOS_LML);
MOS_std_LML = std(MOS_LML);
MOS_Ldata = [MOS_avg_LAP, MOS_avg_LML]; %store MOS_AP (first column) and MOS_ML (second column) in single matrix to save as a variable later



%% Calculate Right Step Interval CoV
StepT_diff = diff(Final_steps(:,2)); %Take the difference between time values to determine step time durations 
temp_RStep_Times = StepT_diff(1:2:end,:); %Vector of all right step times 
avg_RStepT_raw = mean(temp_RStep_Times); %Right Step Time average 
std_RStepT_raw = std(temp_RStep_Times); %Right Step Time standard deviation
CoV_RStepT_raw = (std_RStepT_raw/avg_RStepT_raw)*100; %before cropping step outliers 
H_Cutoff_RStepT = avg_RStepT_raw + 2*std_RStepT_raw; %calculate 2 standard deviations above the mean
L_Cutoff_RStepT = avg_RStepT_raw - 2*std_RStepT_raw; %calcualte 2 standard deviations below the mean
RStep_indeces = find(temp_RStep_Times>L_Cutoff_RStepT & temp_RStep_Times<H_Cutoff_RStepT); %Find indeces within 2 standard deviations from the mean
RStep_Times = temp_RStep_Times(RStep_indeces); % Right step times within 2 standard deviations from the mean
avg_RStepT = mean(RStep_Times); %calculate the new Right Step Time average
std_RStepT = std(RStep_Times); %calculate the new Right Step Time standard deviation
CoV_RStepT = (std_RStepT/avg_RStepT)*100; %Coefficient of Variation for Right Step Time after cropping outliers 

%% Calculate Left Step Interval CoV
temp_LStep_Times = StepT_diff(2:2:end,:); %Vector with all left Step Times
avg_LStepT_raw = mean(temp_LStep_Times); %Average of Left Step Times 
std_LStepT_raw = std(temp_LStep_Times); %Left Step Times Standard Deviation
CoV_LStepT_raw = (std_LStepT_raw/avg_LStepT_raw)*100; %before cropping step outliers
H_Cutoff_LStepT = avg_LStepT_raw + 2* std_LStepT_raw; %calculate 2 standard deviations above the mean
L_Cutoff_LStepT = avg_LStepT_raw - 2*std_LStepT_raw; %calcualte 2 standard deviations below the mean
LStepT_indeces = find(temp_LStep_Times>L_Cutoff_LStepT & temp_LStep_Times<H_Cutoff_LStepT); %Right step times within 2 standard deviations 
LStepTimes = temp_LStep_Times(LStepT_indeces);
avg_LStepT = mean(LStepTimes); %calculate the new Right Step Time average
std_LStepT = std(LStepTimes); %calculate the new Right Step Time standard deviation
CoV_LStepT = (std_LStepT/avg_LStepT)*100; %Coeffici



%% Calculate Right Step Length CoV 
Ltemp_foot = v3d_data.LFootPos{:,:}; % get left foot position data
LFoot_loc_AP = Ltemp_foot(ML_Rind,2); % Left foot AP position data during Right heel-strike
RFoot_loc_AP = Rtemp_foot(ML_Rind,2); % Right AP foot position data during Right heel-Strike 
temp_RSL = abs(LFoot_loc_AP - RFoot_loc_AP); %temporary difference in AP foot position
avg_RSL_raw = mean(temp_RSL); % temporary right step length average 
std_RSL_raw = std(temp_RSL); %temporary rihgt step length standard deviation 
CoV_RSL_raw = (std_RSL_raw/avg_RSL_raw)*100;
HCutoff_RSL = avg_RSL_raw + 2*std_RSL_raw; %calculate cutoff point above 2 standard deviations from average right step length
LCutoff_RSL = avg_RSL_raw - 2*std_RSL_raw; %calculate cutoff point below 2 standard deviations from average right step length
RSL_indeces = (temp_RSL>LCutoff_RSL & temp_RSL<HCutoff_RSL); %data within 2 standard deviations of the average right step length
RSL = temp_RSL(RSL_indeces);
avg_RSL_RHS = mean(RSL); %calculate new right step length average after removing outliers 
std_RSL_RHS = std(RSL); %calculate new rihgt step length standard deviation after removing outliers 
CoV_RSL = (std_RSL_RHS/avg_RSL_RHS)*100; %Coefficient of Variation for Right Step Length


%%  Calculate Left Step Length CoV 
Ltemp_foot = v3d_data.LFootPos{:,:}; % get left foot position data
LFoot_loc_AP_LHS = Ltemp_foot(ML_Lind,2); % AP left foot position data @ LHS
RFoot_loc_AP_LHS = Rtemp_foot(ML_Lind,2); % AP right foot position data @ LHS
temp_LSL = abs(LFoot_loc_AP_LHS - RFoot_loc_AP_LHS); %temporary difference in AP foot position, distance for step length at LHS 
avg_LSL_raw = mean(temp_LSL); % temporary left step length average 
std_LSL_raw = std(temp_LSL); %temporary left step length standard deviation 
CoV_LSL = (std_LSL_raw/avg_LSL_raw)*100;
HCutoff_LSL = avg_LSL_raw + 2*std_LSL_raw; %calculate cutoff point above 2 standard deviations from average step width
LCutoff_LSL = avg_LSL_raw - 2*std_LSL_raw; %calculate cutoff point below 2 standard deviations from average step width
LSL_indeces = (temp_LSL>LCutoff_LSL & temp_LSL<HCutoff_LSL); %data within 2 standard deviations of the average step width
LSL = temp_LSL(LSL_indeces);
avg_LSL_LHS = mean(LSL); %calculate new step width average after removing outliers 
std_LSL_LHS = std(LSL); %calculate new step width standard deviation after removing outliers 
CoV_LSL = (std_LSL_LHS/avg_LSL_LHS)*100; %Coeff



%% Calculate Step Width Cov at Right Heel Strike 
Ltemp_foot = v3d_data.LFootPos{:,:}; % get left foot position data
LFoot_loc_ML = Ltemp_foot(ML_Rind,1); % ML left foot position data
RFoot_loc_ML = Rtemp_foot(ML_Rind,1); % ML right foot position data
temp_SW = abs(LFoot_loc_ML - RFoot_loc_ML); %temporary difference in ML foot position, the step width
avg_SW_RHS_raw = mean(temp_SW); % temporary step width average 
std_SW_RHS_raw = std(temp_SW); %temporary step width standard deviation 
CoV_SW_RHS_raw = (std_SW_RHS_raw/avg_SW_RHS_raw)*100;
HCutoff_SW = avg_SW_RHS_raw + 2*std_SW_RHS_raw; %calculate cutoff point above 2 standard deviations from average step width
LCutoff_SW = avg_SW_RHS_raw - 2*std_SW_RHS_raw; %calculate cutoff point below 2 standard deviations from average step width
SW_indeces = (temp_SW>LCutoff_SW & temp_SW<HCutoff_SW); %data within 2 standard deviations of the average step width
SW = temp_SW(SW_indeces);
avg_SW_RHS = mean(SW); %calculate new step width average after removing outliers 
std_SW_RHS = std(SW); %calculate new step width standard deviation after removing outliers 
CoV_SW_RHS = (std_SW_RHS/avg_SW_RHS)*100; %step width Coefficient of Variation


%% Calculate COV Step Width at Left Heel-Strike
temp_Lfoot = v3d_data.LFootPos{:,:}; % get left foot position data
LFoot_loc_ML = temp_Lfoot(ML_Lind,1); % ML left foot position data at left HS
RFoot_loc_ML = Rtemp_foot(ML_Lind,1); % ML right foot position data at left HS 
temp_SW = abs(LFoot_loc_ML - RFoot_loc_ML); %temporary difference in ML foot position, the step width
avg_SW_LHS_raw = mean(temp_SW); % temporary step width average 
std_SW_LHS_raw = std(temp_SW); %temporary step width standard deviation 
CoV_SW_LHS_raw = (std_SW_LHS_raw/avg_SW_LHS_raw)*100; 
HCutoff_SW = avg_SW_LHS_raw + 2*std_SW_LHS_raw; %calculate cutoff point above 2 standard deviations from average step width
LCutoff_SW = avg_SW_LHS_raw - 2*std_SW_LHS_raw; %calculate cutoff point below 2 standard deviations from average step width
SW_indeces = (temp_SW>LCutoff_SW & temp_SW<HCutoff_SW); %data within 2 standard deviations of the average step width
SW_L = temp_SW(SW_indeces);
avg_SW_LHS = mean(SW_L); %calculate new step width average after removing outliers 
std_SW_LHS = std(SW_L); %calculate new step width standard deviation after removing outliers 
CoV_SW_LHS = (std_SW_LHS/avg_SW_LHS)*100; %step width Coefficient of Variation




%% Calculate ML Harmonic Ratio
[ML_pks, ML_frq] = Harmonic_RatioML1(aCOG_ML, 100); %call the harmonics function for ML COM acceleration
Har_ML_Od = sum(ML_pks(1:2:end)); %Sum of even Harmonic amplitudes for ML COM acceleration 
Har_ML_Ev= sum(ML_pks(2:2:end)); %Sum of odd Harmonic amplitudes for ML COM accelertation 
HR_ML = Har_ML_Od/Har_ML_Ev; %Calculate the Harmonic Ratio

%% AP Harmonic Ratio 
[AP_pks, AP_frq] = Harmonic_RatioAPV1(aCOG_AP,100);
Har_AP_Od = sum(AP_pks(1:2:end)); %Sum of even Harmonic amplitudes for AP COM acceleration 
Har_AP_Ev = sum(AP_pks(2:2:end)); %Sum of odd Harmonic amplitudes for AP COM accelertation 
HR_AP = Har_AP_Ev/Har_AP_Od; %Calculate the Harmonic Ratio
%% Vertical Harmonic Ratio 
[Vert_pks, Vert_frq] = Harmonic_RatioV(aCOG_V,100);
Har_Vert_Od = sum(Vert_pks(1:2:end)); %Sum of even Harmonic amplitudes for AP COM acceleration 
Har_Vert_Ev = sum(Vert_pks(2:2:end)); %Sum of odd Harmonic amplitudes for AP COM accelertation 
HR_Vert = Har_Vert_Ev/Har_Vert_Od; %Calculate the Harmonic Ratio

%% Save Variables 

ML_data.MOS_avg_RML = MOS_avg_RML;
ML_data.MOS_avg_RAP = MOS_avg_RAP;
ML_data.MOS_avg_LML = MOS_avg_LML;
ML_data.MOS_avg_LAP = MOS_avg_LAP;
ML_data.CoV_RStepT = CoV_RStepT;
ML_data.CoV_LStepT = CoV_LStepT;
ML_data.CoV_SW_RHS = CoV_SW_RHS;
ML_data.CoV_SW_LHS = CoV_SW_LHS;
ML_data.CoV_LSL = CoV_LSL;
ML_data.CoV_RSL = CoV_RSL;
ML_data.avg_SW_LHS = avg_SW_LHS;
ML_data.avg_SW_RHS = avg_SW_RHS;
ML_data.avg_LStepT = avg_LStepT;
ML_data.avg_RStepT = avg_RStepT;
ML_data.avg_LSL_LHS = avg_LSL_LHS;
ML_data.avg_RSL_RHS = avg_RSL_RHS;
ML_data.HR_AP = HR_AP;
ML_data.HR_ML = HR_ML;
ML_data.HR_Vert = HR_Vert;
ML_data.MOS_std_RML = MOS_std_RML;
ML_data.MOS_std_LML = MOS_std_LML;
ML_data.std_RStepT = std_RStepT;
ML_data.std_LStepT = std_LStepT;
ML_data.std_SW_RHS= std_SW_RHS;
ML_data.std_SW_LHS = std_SW_LHS;
ML_data.AngMmt_MLavg = AngMmt_MLavg; 
ML_data.AngMmt_MLstd = AngMmt_MLstd; 
ML_data.AngMmt_MLpk = AngMmt_MLpk;
ML_data.ANgMmt_MLpk_SD = AngMmt_MLpk_SD;
ML_data.AngMmt_APavg = AngMmt_APavg; 
ML_data.AngMmt_APstd = AngMmt_APstd;
ML_data.AngMmt_APpk = AngMmt_APpk;
ML_data.AngMmt_APpk_SD = AngMmt_APpk_SD;
ML_data.AngMmt_Vavg = AngMmt_Vavg;
ML_data.AngMmt_Vstd = AngMmt_Vstd;
ML_data.AngMmt_Vpk = AngMmt_Vpk;
ML_data.AngMmt_Vpk_SD = AngMmt_Vpk_SD;
ML_data.Trk_LV_MLavg = TrkLV_MLavg;
ML_data.Trk_LV_MLstd = TrkLV_MLstd;
ML_data.TrkLV_MLpk = TrkLV_MLpk;
ML_data.TrkLV_MLpk_SD = TrkLV_MLpk_SD;
ML_data.Trk_LV_APavg = TrkLV_APavg;
ML_data.Trk_LV_APstd = TrkLV_APstd;
ML_data.TrkLV_APpk = TrkLV_APpk;
ML_data.TrkLV_APpk_SD = TrkLV_APpk_SD;
ML_data.Trk_LV_Vavg = TrkLV_Vavg;
ML_data.Trk_LV_Vstd = TrkLV_Vstd;
ML_data.TrkLV_Vpk = TrkLV_Vpk;
ML_data.TrkLV_Vpk_SD = TrkLV_Vpk_SD;
ML_data.Trk_AV_MLavg = TrkAV_MLavg;
ML_data.Trk_AV_MLstd = TrkAV_MLstd;
ML_data.TrkAV_MLpk = TrkAV_MLpk;
ML_data.TrkAV_MLpk_SD = TrkAV_MLpk_SD;
ML_data.Trk_AV_APavg = TrkAV_APavg;
ML_data.Trk_AV_APstd = TrkAV_APstd;
ML_data.TrkAV_APpk = TrkAV_APpk;
ML_data.TrkAV_APpk_SD = TrkAV_APpk_SD;
ML_data.Trk_AV_Vavg = TrkAV_Vavg;
ML_data.Trk_AV_Vstd = TrkAV_Vstd;
ML_data.TrkAV_Vpk = TrkAV_Vpk;
ML_data.TrkAV_Vpk_SD = TrkAV_Vpk_SD;
ML_data.Hill_t = ML_t;
ML_data.Larm_ROM = Larm_ROM;
ML_data.Rarm_ROM = Rarm_ROM;

end 