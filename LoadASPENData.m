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
t_vec = (M(:,1)-1)' * (1/100); % Converting frame numbers into time (1 Frame = 0.01 sec)

[pos_av_class, att_av_class, pos_tar_class, att_tar_class] = ConvertASPENData(pos_av_aspen, att_av_aspen,  pos_tar_aspen, att_tar_aspen);

% Parsing data to remove zero values
% Matrix concatenating each vector/matrix to make removal easier
parseMat = [t_vec;pos_av_class;att_av_class;pos_tar_class;att_tar_class];

% Searching through each index to find columns with exactly 0
parseMat(:,any(parseMat==0,1)) = [];

% Entering all ConvertASPENData outputs to wanted vectors
t_vec = parseMat(1,:); % Time
av_pos_inert = parseMat(2:4,:); % Position of vehicle
av_att = parseMat(5:7,:); % Attitude of vehicle
tar_pos_inert = parseMat(8:10,:); % Position of target
tar_att = parseMat(11:13,:); % Attitude of target

end