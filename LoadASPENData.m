% Author: Victor Turpin Aguayo
% Course number: ASEN 3801
% File name: LoadASPENData
% Created: 1/20/26

function [t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = LoadASPENData(filename)
%
% Inputs: filename = name of the file that will be used
% 
% Outputs: t_vec =  1 x n time vector in seconds where n = the total number
%                   of frames from the dataset.
%          av_pos_inert = 3 x n matrix of position vectors in meters for 
%                         the aerospace vehicle in Frame E.
%          av_att = 3 x n matrix of attitude vectors listing the 3-2-1 
%                   Euler angles in radians for the aerospace vehicle 
%                   relative to Frame E.
%          tar_pos_inert = 3 x n matrix of position vectors in meters for 
%                          the target in Frame E.
%          tar_att = 3 x n matrix of attitude vectors listing the 3-2-1 
%                    Euler angles in radians for the target relative to 
%                    Frame E.
%
% Methodology: Convert the file into a matrix and clean the data. Assign
% columns to the input variables needed in the ConvertASPENData.m function.
% This function converts from N to E reference. Rename the outputs from
% ConvertASPENData.m to the necessary variable names.

M = readmatrix(filename);
M(1:3,:) = []; % Delete header rows
M(:,2) = []; % Delete subframe column

pos_av_aspen = M(:,11:13)'/1000; % Convert to m (from mm)
att_av_aspen = M(:,8:10)'; % Radians
pos_tar_aspen = M(:,5:7)'/1000; % Convert to m (from mm)
att_tar_aspen = M(:,2:4)';  % Radians

[pos_av_class, att_av_class, pos_tar_class, att_tar_class] = ConvertASPENData(pos_av_aspen, att_av_aspen,  pos_tar_aspen, att_tar_aspen);

av_pos_inert = pos_av_class;
av_att = att_av_class;
tar_pos_inert = pos_tar_class;
tar_att = att_tar_class;
t_vec = (M(:,1)-1)' * (1/100);

end