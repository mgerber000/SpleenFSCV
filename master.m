%% Add function paths
p = pwd;
addpath(strcat(p,'\root'))
addpath(strcat(p,'\dependents'))
%% Create FSCV data storage object
rep = FCSVExperiment;

%% Trim FSCV Data and register with physiologic recordings
reg = PhysioFCSV_Reg(rep);

%% Quantify FSCV data and output as table
data = FCSVCompare(reg,[]);