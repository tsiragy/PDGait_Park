%% clean slate
clear all
close all
clc



%% save to excel

pathname = uigetdir;
allfiles = dir(fullfile(pathname,'*.mat'));
matfiles = allfiles([allfiles.isdir] == 0);   %probably not needed since it's unlikely that a xxx.mat is a directory
    for fileid = 1:numel(matfiles)
       [~, filename] = fileparts(matfiles(fileid).name);
       filetable(:,fileid) = struct2cell(load(fullfile(pathname, matfiles(fileid).name))); 
       %filetable.Properties.VariableNames{2} = filename;  %rename 2nd variable
       %fulltable = [fulltable, filetable(:, 2)];  %only copy second variable
    end

Left_COV_Time = filetable(1,:);
S1 = Left_COV_Time(1:6);
S2 = Left_COV_Time(7:12);
S3 = Left_COV_Time(13:18);
S4 = Left_COV_Time(19:24);
S5 = Left_COV_Time(25:30);
S6 = Left_COV_Time(31:36);
S7 = Left_COV_Time(37:42);
S8 = Left_COV_Time(43:48);
S9 = Left_COV_Time(49:54);
S10 = Left_COV_Time(55:60);
S11 = Left_COV_Time(61:66);
S12 = Left_COV_Time(67:72);
S13 = Left_COV_Time(73:78);
S14 = Left_COV_Time(79:84);
S15 = Left_COV_Time(85:90);

COV_LStepT = [S1; S2; S3; S4; S5; S6; S7; S8; S9; S10; S11; S12; S13; S14; S15];
xlswrite('COV_LeftStepT', COV_LStepT);


Left_COV_Length = filetable(2,:);
L_LS1 = Left_COV_Length(1:6);
L_LS2 = Left_COV_Length(7:12);
L_LS3 = Left_COV_Length(13:18);
L_LS4 = Left_COV_Length(19:24);
L_LS5 = Left_COV_Length(25:30);
L_LS6 = Left_COV_Length(31:36);
L_LS7 = Left_COV_Length(37:42);
L_LS8 = Left_COV_Length(43:48);
L_LS9 = Left_COV_Length(49:54);
L_LS10 = Left_COV_Length(55:60);
L_LS11 = Left_COV_Length(61:66);
L_LS12 = Left_COV_Length(67:72);
L_LS13 = Left_COV_Length(73:78);
L_LS14 = Left_COV_Length(79:84);
L_LS15 = Left_COV_Length(85:90);
Left_COV_SLength = [L_LS1; L_LS2; L_LS3; L_LS4; L_LS5; L_LS6; L_LS7; L_LS8; L_LS9; L_LS10; L_LS11; L_LS12; L_LS13; L_LS14; L_LS15];
xlswrite('COV_LStep_Length',Left_COV_SLength);


Right_COV_TIME = filetable(3,:);
R_S1 = Right_COV_TIME(1:6);
R_S2 = Right_COV_TIME(7:12);
R_S3 = Right_COV_TIME(13:18);
R_S4 = Right_COV_TIME(19:24);
R_S5 = Right_COV_TIME(25:30);
R_S6 = Right_COV_TIME(31:36);
R_S7 = Right_COV_TIME(37:42);
R_S8 = Right_COV_TIME(43:48);
R_S9 = Right_COV_TIME(49:54);
R_S10 = Right_COV_TIME(55:60);
R_S11 = Right_COV_TIME(61:66);
R_S12 = Right_COV_TIME(67:72);
R_S13 = Right_COV_TIME(73:78);
R_S14 = Right_COV_TIME(79:84);
R_S15 = Right_COV_TIME(85:90);
Right_COV_TIME = [R_S1; R_S2; R_S3; R_S4; R_S5; R_S6; R_S7; R_S8; R_S9; R_S10; R_S11; R_S12; R_S13; R_S14; R_S15];
xlswrite('COV_RStep_TIME',Right_COV_TIME);


Right_COV_SLength = filetable(4,:);
R_SL1 = Right_COV_SLength(1:6);
R_SL2 = Right_COV_SLength(7:12);
R_SL3 = Right_COV_SLength(13:18);
R_SL4 = Right_COV_SLength(19:24);
R_SL5 = Right_COV_SLength(25:30);
R_SL6 = Right_COV_SLength(31:36);
R_SL7 = Right_COV_SLength(37:42);
R_SL8 = Right_COV_SLength(43:48);
R_SL9 = Right_COV_SLength(49:54);
R_SL10 = Right_COV_SLength(55:60);
R_SL11 = Right_COV_SLength(61:66);
R_SL12 = Right_COV_SLength(67:72);
R_SL13 = Right_COV_SLength(73:78);
R_SL14 = Right_COV_SLength(79:84);
R_SL15 = Right_COV_SLength(85:90);
Right_COV_StepLength = [R_SL1; R_SL2; R_SL3; R_SL4; R_SL5; R_SL6; R_SL7; R_SL8; R_SL9; R_SL10; R_SL11; R_SL12; R_SL13; R_SL14; R_SL15];
xlswrite('COV_RStep_length',Right_COV_StepLength );


COV_SW = filetable(5,:);
SW1 = COV_SW(1:6);
SW2 = COV_SW(7:12);
SW3 = COV_SW(13:18);
SW4 = COV_SW(19:24);
SW5 = COV_SW(25:30);
SW6 = COV_SW(31:36);
SW7 = COV_SW(37:42);
SW8 = COV_SW(43:48);
SW9 = COV_SW(49:54);
SW10 = COV_SW(55:60);
SW11 = COV_SW(61:66);
SW12 = COV_SW(67:72);
SW13 = COV_SW(73:78);
SW14 = COV_SW(79:84);
SW15 = COV_SW(85:90);
COV_StepWidth = [SW1; SW2; SW3; SW4; SW5; SW6; SW7; SW8; SW9; SW10; SW11; SW12; SW13; SW14; SW15];
xlswrite('COV_StepWidth',COV_StepWidth);

avg_Left_StepLgth = filetable(6,:);
L_SL1 = avg_Left_StepLgth(1:6);
L_SL2 = avg_Left_StepLgth(7:12);
L_SL3 = avg_Left_StepLgth(13:18);
L_SL4 = avg_Left_StepLgth(19:24);
L_SL5 = avg_Left_StepLgth(25:30);
L_SL6 = avg_Left_StepLgth(31:36);
L_SL7 = avg_Left_StepLgth(37:42);
L_SL8 = avg_Left_StepLgth(43:48);
L_SL9 = avg_Left_StepLgth(49:54);
L_SL10 = avg_Left_StepLgth(55:60);
L_SL11 = avg_Left_StepLgth(61:66);
L_SL12 = avg_Left_StepLgth(67:72);
L_SL13 = avg_Left_StepLgth(73:78);
L_SL14 = avg_Left_StepLgth(79:84);
L_SL15 = avg_Left_StepLgth(85:90);
avg_Left_StepLgth = [L_SL1; L_SL2; L_SL3; L_SL4; L_SL5; L_SL6; L_SL7; L_SL8; L_SL9; L_SL10; L_SL11; L_SL12; L_SL13; L_SL14; L_SL15];
xlswrite('avg_Left_StepLgth', avg_Left_StepLgth);


avg_Left_StepTime = filetable(7,:);
L_ST1 = avg_Left_StepTime(1:6);
L_ST2 = avg_Left_StepTime(7:12);
L_ST3 = avg_Left_StepTime(13:18);
L_ST4 = avg_Left_StepTime(19:24);
L_ST5 = avg_Left_StepTime(25:30);
L_ST6 = avg_Left_StepTime(31:36);
L_ST7 = avg_Left_StepTime(37:42);
L_ST8 = avg_Left_StepTime(43:48);
L_ST9 = avg_Left_StepTime(49:54);
L_ST10 = avg_Left_StepTime(55:60);
L_ST11 = avg_Left_StepTime(61:66);
L_ST12 = avg_Left_StepTime(67:72);
L_ST13 = avg_Left_StepTime(73:78);
L_ST14 = avg_Left_StepTime(79:84);
L_ST15 = avg_Left_StepTime(85:90);
avg_Left_StepTime = [L_ST1; L_ST2; L_ST3; L_ST4; L_ST5; L_ST6; L_ST7; L_ST8; L_ST9; L_ST10; L_ST11; L_ST12; L_ST13; L_ST14; L_ST15];
xlswrite('avg_Left_StepTime', avg_Left_StepTime);


AP_MOS = filetable(8,:);
xcom1 = AP_MOS(1:6);
xcom2 = AP_MOS(7:12);
xcom3 = AP_MOS(13:18);
xcom4 = AP_MOS(19:24);
xcom5 = AP_MOS(25:30);
xcom6 = AP_MOS(31:36);
xcom7 = AP_MOS(37:42);
xcom8 = AP_MOS(43:48);
xcom9 = AP_MOS(49:54);
xcom10 = AP_MOS(55:60);
xcom11 = AP_MOS(61:66);
xcom12 = AP_MOS(67:72);
xcom13 = AP_MOS(73:78);
xcom14 = AP_MOS(79:84);
xcom15 = AP_MOS(85:90);
MOS_AP = [xcom1; xcom2; xcom3; xcom4; xcom5; xcom6; xcom7; xcom8; xcom9; xcom10; xcom11; xcom12; xcom13; xcom14; xcom15];
xlswrite('MOS_AP',MOS_AP);


ML_MOS = filetable(9,:);
MLxcom1 = ML_MOS(1:6);
MLxcom2 = ML_MOS(7:12);
MLxcom3 = ML_MOS(13:18);
MLxcom4 = ML_MOS(19:24);
MLxcom5 = ML_MOS(25:30);
MLxcom6 = ML_MOS(31:36);
MLxcom7 = ML_MOS(37:42);
MLxcom8 = ML_MOS(43:48);
MLxcom9 = ML_MOS(49:54);
MLxcom10 = ML_MOS(55:60);
MLxcom11 = ML_MOS(61:66);
MLxcom12 = ML_MOS(67:72);
MLxcom13 = ML_MOS(73:78);
MLxcom14 = ML_MOS(79:84);
MLxcom15 = ML_MOS(85:90);
MOS_ML = [MLxcom1; MLxcom2; MLxcom3; MLxcom4; MLxcom5; MLxcom6; MLxcom7; MLxcom8; MLxcom9; MLxcom10; MLxcom11; MLxcom12; MLxcom13; MLxcom14; MLxcom15];
xlswrite('MOS_ML',MOS_ML);



avg_Right_StepLgth = filetable(10,:);
R_SL1 = avg_Right_StepLgth(1:6);
R_SL2 = avg_Right_StepLgth(7:12);
R_SL3 = avg_Right_StepLgth(13:18);
R_SL4 = avg_Right_StepLgth(19:24);
R_SL5 = avg_Right_StepLgth(25:30);
R_SL6 = avg_Right_StepLgth(31:36);
R_SL7 = avg_Right_StepLgth(37:42);
R_SL8 = avg_Right_StepLgth(43:48);
R_SL9 = avg_Right_StepLgth(49:54);
R_SL10 = avg_Right_StepLgth(55:60);
R_SL11 = avg_Right_StepLgth(61:66);
R_SL12 = avg_Right_StepLgth(67:72);
R_SL13 = avg_Right_StepLgth(73:78);
R_SL14 = avg_Right_StepLgth(79:84);
R_SL15 = avg_Right_StepLgth(85:90);
avg_Right_StepLgth = [R_SL1; R_SL2; R_SL3; R_SL4; R_SL5; R_SL6; R_SL7; R_SL8; R_SL9; R_SL10; R_SL11; R_SL12; R_SL13; R_SL14; R_SL15];
xlswrite('avg_Right_StepLgth', avg_Right_StepLgth);


avg_Right_StepTime = filetable(11,:);
R_ST1 = avg_Right_StepTime(1:6);
R_ST2 = avg_Right_StepTime(7:12);
R_ST3 = avg_Right_StepTime(13:18);
R_ST4 = avg_Right_StepTime(19:24);
R_ST5 = avg_Right_StepTime(25:30);
R_ST6 = avg_Right_StepTime(31:36);
R_ST7 = avg_Right_StepTime(37:42);
R_ST8 = avg_Right_StepTime(43:48);
R_ST9 = avg_Right_StepTime(49:54);
R_ST10 = avg_Right_StepTime(55:60);
R_ST11 = avg_Right_StepTime(61:66);
R_ST12 = avg_Right_StepTime(67:72);
R_ST13 = avg_Right_StepTime(73:78);
R_ST14 = avg_Right_StepTime(79:84);
R_ST15 = avg_Right_StepTime(85:90);
avg_Right_StepTime = [R_ST1; R_ST2; R_ST3; R_ST4; R_ST5; R_ST6; R_ST7; R_ST8; R_ST9; R_ST10; R_ST11; R_ST12; R_ST13; R_ST14; R_ST15];
xlswrite('avg_Right_StepTime', avg_Right_StepTime);




avg_SW = filetable(12,:);
ASW1 = avg_SW(1:6);
ASW2 = avg_SW(7:12);
ASW3 = avg_SW(13:18);
ASW4 = avg_SW(19:24);
ASW5 = avg_SW(25:30);
ASW6 = avg_SW(31:36);
ASW7 = avg_SW(37:42);
ASW8 = avg_SW(43:48);
ASW9 = avg_SW(49:54);
ASW10 = avg_SW(55:60);
ASW11 = avg_SW(61:66);
ASW12 = avg_SW(67:72);
ASW13 = avg_SW(73:78);
ASW14 = avg_SW(79:84);
ASW15 = avg_SW(85:90);
avg_StepWidth = [ASW1; ASW2; ASW3; ASW4; ASW5; ASW6; ASW7; ASW8; ASW9; ASW10; ASW11; ASW12; ASW13; ASW14; ASW15];
xlswrite('avg_StepWidth', avg_StepWidth);
    
    










