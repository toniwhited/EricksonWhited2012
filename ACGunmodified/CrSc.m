%PROGRAM:   CrSc.m
%Version:   4/10/2009
%Random variable assingnment (1 if Normal(0,1); 2 if lognormal(standardized); 3 if t-student - and you provide the degrees of freedom); 4 - Chi2 (standardized)
%INPUTS:    rep    - number of replications
%           obs    - number of observations in each cross section
%CALLS:     CSest.m
%%%%%%%%%%%%
echo off all
warning off all
clear; tic;
%%%%% INPUT  %%%%%
rep = 10000;
obs = 1000;
%%%%%%%%%%%%%%%%%%
Level = 1;

normal = 1; lognormal = 2; tstudent = 3; chi2 = 4; F = 5;

% % Normal (11.txt)
% dist_u  = lognormal;
% dist_xerror = lognormal;
% dist_chi = normal; t_chi = num2str(dist_chi);
% dist_XX = normal; t_XX = num2str(dist_XX);
% dist_fe = normal; 
% df = [1 1]; 
% 
% f_nameout = strcat(t_chi,t_XX,'CS.txt');
%    
% Dist = [dist_u dist_xerror dist_chi dist_XX dist_fe df];
% 
% status = CSest(f_nameout,Dist, rep, obs, Level)                
%

% lognormal (22.txt)
dist_u  = lognormal;
dist_xerror = lognormal;
dist_chi = lognormal; t_chi = num2str(dist_chi);
dist_XX = lognormal; t_XX = num2str(dist_XX);
dist_fe = normal; 
df = [1 1]; 

f_nameout = strcat(t_chi,t_XX,'CS.txt');

Dist = [dist_u dist_xerror dist_chi dist_XX dist_fe df];

status = CSest(f_nameout,Dist, rep, obs, Level)                
%

% % Chi2 (44_5.txt))
% dist_u  = lognormal;
% dist_xerror = lognormal;
% dist_chi = chi2; t_chi = num2str(dist_chi);
% dist_XX = chi2; t_XX = num2str(dist_XX);
% dist_fe = normal; 
% df = [5 1]; 
% 
% f_nameout = strcat(t_chi,t_XX,'CS.txt');
% 
% Dist = [dist_u dist_xerror dist_chi dist_XX dist_fe df];
% 
% status = CSest(f_nameout,Dist, rep, obs, Level)                
% %
%  
% % F (55_1040.txt)
% dist_u  = lognormal;
% dist_xerror = lognormal;
% dist_chi = F; t_chi = num2str(dist_chi);
% dist_XX = F; t_XX = num2str(dist_XX);
% dist_fe = normal; 
% df = [10 40]; 
% 
% f_nameout = strcat(t_chi,t_XX,'CS.txt');
%     
% Dist = [dist_u dist_xerror dist_chi dist_XX dist_fe df];
% 
% status = CSest(f_nameout,Dist, rep, obs, Level)                
% %

toc