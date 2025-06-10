clear
close all
clc

addpath('../../config')
addpath('../../data')
addpath('../../github_repo')
addpath('../../libs')
addpath('../../data/workspace_files')
addpath('../../apps/self_calibration')
addpath('../../libs/cdpr_model')
addpath('../../libs/export_utilities')
addpath('../../libs/numeric')
addpath('../../libs/orientation_geometry')
addpath('../../libs/under_actuated')
addpath('../../libs/over_actuated')
folder = '../../data';


%%% choose the desired robot
[cdpr_parameters, cdpr_variables, cdpr_ws_data ,cdpr_outputs,record,utilities] = ...
 LoadConfigAndInit("IRMA8_ws_analysis","IRMA8_ws_analysis");
cdpr_variables = UpdateIKZeroOrd([1.78352; 0.197381; -1.28125],...
  [0;0;0],cdpr_parameters,cdpr_variables);
record.SetFrame(cdpr_variables,cdpr_parameters);%
 ws_info = LoadWsInfo("8_cable_info");
 
%%% Final design update
% load("ipanema_grad_opt_ark.mat");
%%% Selection of force controlled cables
for cableCoupleIndex=8:8 % Update this till 28 (or 4 for planars)
    cableComb=nchoosek(1:cdpr_parameters.n_cables,cdpr_parameters.n_cables-cdpr_parameters.pose_dim);
    cablesForceControlled=cableComb(cableCoupleIndex,:);

    tic
    cdpr_outputs = CalcWorkspace(cdpr_parameters,cdpr_variables,...
        utilities,cdpr_outputs,folder,record,ws_info,cablesForceControlled);
    toc
    
    %%% UM conversions
    cdpr_outputs.rotCardouSens=180/pi*(cdpr_outputs.rotCardouSens)/1000;

    %%% Accuracy set limits
    for i = 1:cdpr_outputs.counter
        sigma_threshold=[8 1];
        if (cdpr_outputs.posCardouSens(i)<sigma_threshold(1) && cdpr_outputs.rotCardouSens(i)<sigma_threshold(2))
            cdpr_outputs.kaw(i)=1;
        else
            cdpr_outputs.kaw(i)=0;
        end
    end

    %%% Tension distribution sensitivity
    for i = 1:cdpr_outputs.counter
        if (cdpr_outputs.teiw(i)<3)
            cdpr_outputs.teiw(i)=1;
        else
            cdpr_outputs.teiw(i)=0;
        end
    end

    %%% pie chart contruction
    teiw_vol = nnz(cdpr_outputs.teiw)/cdpr_outputs.full_counter*cdpr_outputs.ws_volume;
    perc_not_wfw = 1-cdpr_outputs.counter/cdpr_outputs.full_counter;
    perc_err_sens = nnz(~cdpr_outputs.teiw)/cdpr_outputs.full_counter;
    perc_err_insens_inacc = nnz(cdpr_outputs.teiw-cdpr_outputs.kaw.*cdpr_outputs.teiw)/cdpr_outputs.full_counter;
    perc_ka_teiw = nnz(cdpr_outputs.kaw.*cdpr_outputs.teiw)/cdpr_outputs.full_counter;
    cdpr_outputs.X_pie = [perc_not_wfw perc_err_sens perc_err_insens_inacc perc_ka_teiw];

    %%% Compute wfw volume
    [k_points,wfw_vol]=boundary(cdpr_outputs.pose(1,:)',cdpr_outputs.pose(2,:)',cdpr_outputs.pose(3,:)',0.99);
    % trisurf(k_points,cdpr_outputs.pose(1,:),cdpr_outputs.pose(2,:),cdpr_outputs.pose(3,:));

    %%% save data
    avoidVariable = 'record';
    cableIdx=string(cablesForceControlled);
    if cdpr_parameters.pose_dim==3
        filename = strcat(folder,'/workspace_files/',record.figure_handle.FileName,cableIdx(1));
    else
        filename = strcat(folder,'/workspace_files/',record.figure_handle.FileName,cableIdx(1),cableIdx(end));
    end
    save(strcat(filename,'_WS','.mat'),'-regexp',['^(?!',avoidVariable,'$).'])
end
% DisplayAndSaveWorkspace(cdpr_parameters,cdpr_variables,cdpr_outputs,ws_info,folder,record);

%% Compile to show data
close all

% choose the correct title
sim_title="8_cable_HRPCable";
% sim_title="8_cable_IPAnema3_3";
% sim_title="Laser_engraver_planar";
config_name=sim_title; 
% cnt = 7; 
for cnt=7:7
    switch cnt
        case 1
            ws_info.display_criteria = 'TENSION_SENSITIVITYinputRatio';
            PlotWorkspaceIndex(cdpr_parameters,cdpr_variables,cdpr_outputs,ws_info,folder,sim_title,config_name,cablesForceControlled);
        case 2
            ws_info.display_criteria = 'TENSION_SENSITIVITY';
            PlotWorkspaceIndex(cdpr_parameters,cdpr_variables,cdpr_outputs,ws_info,folder,sim_title,config_name,cablesForceControlled);
        case 3
            ws_info.display_criteria = 'CABLE_SENSITIVITYrotCardou';
            PlotWorkspaceIndex(cdpr_parameters,cdpr_variables,cdpr_outputs,ws_info,folder,sim_title,config_name,cablesForceControlled);
        case 4
            ws_info.display_criteria = 'CABLE_SENSITIVITYposCardou';
            PlotWorkspaceIndex(cdpr_parameters,cdpr_variables,cdpr_outputs,ws_info,folder,sim_title,config_name,cablesForceControlled);
        case 5
            ws_info.display_criteria = 'TENSION_ERR_INSENSITIVE_WS';
            PlotWorkspaceIndex(cdpr_parameters,cdpr_variables,cdpr_outputs,ws_info,folder,sim_title,config_name,cablesForceControlled);
        case 6
            ws_info.display_criteria = 'KIN_ACCURATE_WS';
            PlotWorkspaceIndex(cdpr_parameters,cdpr_variables,cdpr_outputs,ws_info,folder,sim_title,config_name,cablesForceControlled);
        case 7
            ws_info.display_criteria = 'KA_TEI_WS';
            PlotWorkspaceIndex(cdpr_parameters,cdpr_variables,cdpr_outputs,ws_info,folder,sim_title,config_name,cablesForceControlled);
    end
end