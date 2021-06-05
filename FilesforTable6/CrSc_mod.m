%PROGRAM:   CrSc_mod.m
%Version:   6/15/2011
%Random variable assingnment (1 if Normal(0,1); 2 if lognormal(standardized); 3 if t-student - and you provide the degrees of freedom); 4 - Chi2 (standardized)
%INPUTS:    rep    - number of replications
%           obs    - number of observations in each cross section
%CALLS:     CSest.m
%Modified by Egor Matveyev and Toni Whited

%%%%%%%%%%%%
echo off all
warning off all
clear; clc; tic;
%%%%% INPUT  %%%%%
rep = 5000;
obs = 1000;
%%%%%%%%%%%%%%%%%%
Level = 1;

normal = 1; lognormal = 2; tstudent = 3; chi2 = 4; F = 5;


%%%%OUTPUT STUFF
outstr0  = 'number of trials, sample size, and seed %8.3f %8.3f %8.3f \n';
outstr1  = 'EW-GMM3 & Bias & %8.3f &  %8.3f &  %8.3f &  %8.3f   \\\\ \n';
outstr11 = 'EW-GMM4 & Bias & %8.3f &  %8.3f &  %8.3f &  %8.3f   \\\\ \n';
outstr12 = 'EW-GMM5 & Bias & %8.3f &  %8.3f &  %8.3f &  %8.3f   \\\\ \n';

outstr2  = ' & ACG     RMSE  & %8.3f &  %8.3f &  %8.3f &  %8.3f   \\\\ \n';
outstr3  = ' & Correct RMSE  & %8.3f &  %8.3f &  %8.3f &  %8.3f   \\\\ \n';
ti       = clock;
fname    = strcat(num2str(ti(1,1)),...
             num2str(ti(1,2)),...
             num2str(ti(1,3)),...
             num2str(ti(1,4)),...
             num2str(ti(1,5)),...
             num2str(ti(1,6),'%2.0f'));
fname1=strcat('output/',fname,'.out');
global fid;

jj=1;
while jj<=2;
    %SEEDS
    %rand('state',jj); randn('state',jj);
    %NOTE: If Matlab version is 7.7 and higher, replace the above line with
    %the following:
    seedj = RandStream('mt19937ar','Seed',jj); RandStream.setDefaultStream(seedj);
    %%%%%%


    %%%%%Almeida's Code
    dist_u  = lognormal;
    dist_xerror = lognormal;
    dist_chi = F; t_chi = num2str(dist_chi);
    dist_XX = F; t_XX = num2str(dist_XX);
    dist_fe = normal;
    df = [10 40];
    f_nameout = strcat(t_chi,t_XX,'CS.txt');
    Dist = [dist_u dist_xerror dist_chi dist_XX dist_fe df];
    [status biasgmm rmsegmm rmsegmm_pr] = CSest(f_nameout,Dist, rep, obs, Level);
    %%%%%
    fid = fopen(fname1, 'a');

    fprintf(fid, outstr0, [rep obs jj]);
    fprintf(fid, outstr1, biasgmm(1,2:5));
    fprintf(fid, outstr2, rmsegmm(1,2:5));
    fprintf(fid, outstr3, rmsegmm_pr(1,2:5));
    fprintf(fid, outstr11, biasgmm(2,2:5));
    fprintf(fid, outstr2, rmsegmm(2,2:5));
    fprintf(fid, outstr3, rmsegmm_pr(2,2:5));
    fprintf(fid, outstr12, biasgmm(3,2:5));
    fprintf(fid, outstr2, rmsegmm(3,2:5));
    fprintf(fid, outstr3, rmsegmm_pr(3,2:5));
    fprintf(fid,' \\hline \n');
    fclose('all');

    jj=jj+1;

end



toc
