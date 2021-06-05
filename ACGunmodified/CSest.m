function status = CSest(f_nameout, Dist_setup, rept, obs, Level)
%FUNCTION:  Monte Carlo simulations for different estimators (IV and Whitted)
%INPUTS:    f_nameout - name of the output file (depending on the errors and variables distribution
%           Dist_setup - vector indicating which distribution is used to generate errors and variables
%           rept - number of repetitions for the Monte Carlo procedure
%           obs - number of observations for each repetition of the Monte Carlo
%OUTPUTS:   status - indicator that no errors ocurred

%IMPORTANT: This function has its own 'inputs'. Here we control model specifications
%%%%%%%%%%%%%%%%%%

%Monte Carlo parameters
rep = rept;                         % # of replications 
n = obs;                            % # of observations in each estimation procedure
t = 1;

Xdim = 5;                                                      %Number of Variables 
mism = 1;                                                      %Number of mismeasured variables
Beta0 = [0 1 -1 1 -1];                                         %Number of parameters 

WD_order = 3;                                                  %EW estimation: GMM-order: 1->GMM3, 2->GG4,...

%%OUTPUTS
Wstat = zeros(rep,1); Wstat_Pval = zeros(rep,1);

BetahWL= zeros(rep,Xdim*(t+t*WD_order+WD_order)); BetaSEWL = zeros(rep,Xdim*(t+t*WD_order+WD_order));
J_L = zeros(rep,t*WD_order); J_PvalL = zeros(rep,t*WD_order);

BetahWW= zeros(rep,Xdim*(t+t*WD_order+WD_order)); BetaSEWW = zeros(rep,Xdim*(t+t*WD_order+WD_order));
J_W = zeros(rep,t*WD_order); J_PvalW = zeros(rep,t*WD_order);


%%START COMPUTATION
nt = n * t;
year = repmat((1:t)',n,1);

Constant = ones(nt,1);

[Xnames Zetadim] = genNames(Xdim, mism, WD_order);             %Creates names for the variables used by Whitted estimator

mu = 0;                                                        %--- Normal distrubution parameters ----
sig = 1;                                                       %using these parameters creates vector u and matrix Xmu = 0;


V = eye(Zetadim+mism) - 0.5*(eye(Zetadim+mism)-(ones(Zetadim+mism)));
% v_EW = chol(V);
% v_EW = 1;
[Qmx, vLmx] = eig(V);           
v_EW = Qmx*sqrt(vLmx)*Qmx';

for i=1:rep
    
    %Data Generation
    
    [u XX xerror IVerror fe] = genvars(nt, Zetadim, mism, mu, sig, Dist_setup); %Creates variables using specified distribution (XX without constant)
    
    XX = XX * v_EW;
    
    chi = XX(:,1);
    z1 = XX(:,2);
    z2 = XX(:,3);
    z3 = XX(:,4);     
    
    Y = Beta0(1) + Beta0(2) * chi + Beta0(3)* z1 + Beta0(4) * z2 + Beta0(5) * z3 + u;                             % IMPORTANT: # of parameters = Xdim 

    x = chi + xerror;
    
    %EW Estimation
    
    if Level

        D2 = [Y x Constant z1 z2 z3];
        
        % Pretest 
        
        [Wstat(i,:) Wstat_Pval(i,:)] = m_idtest1(D2, Xnames, t, n, year);
        
        % Estimation
                
        [BetahWL(i,:) BetaSEWL(i,:) J_L(i,:) J_PvalL(i,:)] = genBetaWD(n,t,Xdim,D2,WD_order,Xnames,Zetadim,year);
        
        clear u XX xerror IVerror fe Y BX D D2;
    else

        %Within transformation
        
        um = ones(t,1); U = um*um'; U = (1/t)*U;
        P = kron(eye(n), U);
        D = (eye(n*t) - P)*D;
        
        D = [D(:,1:1+mism) Constant D(:,mism+2:end)];             
        
        [BetahWW(i,:) BetaSEWW(i,:) J_W(i,:) J_PvalW(i,:)] = genBetaWD(n,t,Xdim,D,WD_order,Xnames,Zetadim,year);
        
        clear U P D;
    end
    %
end

clear U P D;

%%EW test
Wstat_m = mean(Wstat);
Wstat_Pm = mean(Wstat_Pval);
Wstat_Paux = Wstat_Pval < .05;
Wtest_size = mean(Wstat_Paux);

%%Averaging estimations results

Beta_mL = mean(BetahWL,1);
Beta_mSEL = mean(BetaSEWL,1);
Beta_mW = mean(BetahWW,1);
Beta_mSEW = mean(BetaSEWW,1);
clear BetaIV BetaIvar BetahWL BetaSEWL BetahWW BetaSEWW;
    
%OLS EW estimates - one estimation (natural order) for each t

Beta_mOLSL = Beta_mL(1:t*Xdim); Beta_mSEOLSL = Beta_mSEL(1:t*Xdim); 
Beta_mOLSW = Beta_mW(1:t*Xdim); Beta_mSEOLSW = Beta_mSEW(1:t*Xdim);
    
    %To order GMM EW estimates better - (t1: GMM3, GMM4, ... t2: GMM3, GMM4, ...); each GMMi in natural order
Beta_mGMML = zeros(WD_order, Xdim,t); Beta_mSEGMML = zeros(WD_order, Xdim,t);    
Beta_mGMMW = zeros(WD_order, Xdim,t); Beta_mSEGMMW = zeros(WD_order, Xdim,t);        
Beta_aux = Beta_mL( t*Xdim+1:t*Xdim*(1+WD_order) ); Beta_SEaux = Beta_mSEL( t*Xdim+1:t*Xdim*(1+WD_order) );
Beta_aux1 = Beta_mW( t*Xdim+1:t*Xdim*(1+WD_order) ); Beta_SEaux1 = Beta_mSEW( t*Xdim+1:t*Xdim*(1+WD_order) );
if WD_order > 1
    tdim = WD_order + WD_order + Zetadim*WD_order;
    for T=1:t
        Baux = Beta_aux((T-1)*tdim+1:T*tdim);
        BauxSE = Beta_SEaux((T-1)*tdim+1:T*tdim);
        Baux1 = Beta_aux1((T-1)*tdim+1:T*tdim);
        Baux1SE = Beta_SEaux1((T-1)*tdim+1:T*tdim);
        for i=1:WD_order
            Beta_mGMML(i,1,T) = Baux(i); Beta_mSEGMML(i,1,T) = BauxSE(i);
            Beta_mGMML(i,2,T) = Baux(WD_order + i); Beta_mSEGMML(i,2,T) = BauxSE(WD_order + i);
            Beta_mGMML(i,3:end,T) = Baux((2+(i-1)) * WD_order + 1:(2+(i-1)) * WD_order + Zetadim );
            Beta_mSEGMML(i,3:end,T) = BauxSE((2+(i-1)) * WD_order + 1:(2+(i-1)) * WD_order + Zetadim );
            
            Beta_mGMMW(i,1,T) = Baux1(i); Beta_mSEGMMW(i,1,T) = Baux1SE(i);
            Beta_mGMMW(i,2,T) = Baux1(WD_order + i); Beta_mSEGMMW(i,2,T) = Baux1SE(WD_order + i);
            Beta_mGMMW(i,3:end,T) = Baux1((2+(i-1)) * WD_order + 1:(2+(i-1)) * WD_order + Zetadim );
            Beta_mSEGMMW(i,3:end,T) = Baux1SE((2+(i-1)) * WD_order + 1:(2+(i-1)) * WD_order + Zetadim );            
        end
    end
end

%MD EW estimates

Beta_mMDL = Beta_mL( t*Xdim*(1+WD_order)+1:Xdim*(t+t*WD_order+WD_order) ); Beta_mSEMDL = Beta_mSEL( t*Xdim*(1+WD_order)+1:Xdim*(t+t*WD_order+WD_order) );  
Beta_mMDW = Beta_mW( t*Xdim*(1+WD_order)+1:Xdim*(t+t*WD_order+WD_order) ); Beta_mSEMDW = Beta_mSEW( t*Xdim*(1+WD_order)+1:Xdim*(t+t*WD_order+WD_order) );      

%J-statistics

JJL = mean(J_PvalL > .05,1);
JJW = mean(J_PvalW > .05,1);

%RMSE

Bias_OLSL = Beta_mOLSL - repmat(Beta0,1,t); 
RMSE_OLSL = sqrt(power(Beta_mSEOLSL,2) + power(Bias_OLSL,2));
Beta0Mat = repmat(Beta0,WD_order,1);
Bias_GMML = zeros(size(Beta_mGMML));
RMSE_GMML = zeros(size(Beta_mGMML));
for T=1:t
    Bias_GMML(:,:,T) = Beta_mGMML(:,:,T) - Beta0Mat;
    RMSE_GMML(:,:,T) = sqrt( power(Beta_mSEGMML(:,:,T),2) + power(Bias_GMML(:,:,T),2) );
end
Bias_MDL = Beta_mMDL - repmat(Beta0,1,WD_order);
RMSE_MDL = sqrt( power(Beta_mSEMDL,2) + power(Bias_MDL,2) );

Bias_OLSW = Beta_mOLSW - repmat(Beta0,1,t); 
RMSE_OLSW = sqrt(power(Beta_mSEOLSW,2) + power(Bias_OLSW,2));
Bias_GMMW = zeros(size(Beta_mGMMW));
RMSE_GMMW = zeros(size(Beta_mGMMW));
for T=1:t
    Bias_GMMW(:,:,T) = Beta_mGMMW(:,:,T) - Beta0Mat;
    RMSE_GMMW(:,:,T) = sqrt( power(Beta_mSEGMMW(:,:,T),2) + power(Bias_GMMW(:,:,T),2) );
end
Bias_MDW = Beta_mMDW - repmat(Beta0,1,WD_order);
RMSE_MDW = sqrt( power(Beta_mSEMDW,2) + power(Bias_MDW,2) );


%%File writing procedure

str = date;
fid = fopen(f_nameout, 'w');
fprintf(fid, '%s\t Rep:%5d\t Obs:%5d\n\n', str, rep, obs);
fprintf(fid, '\nRegression coefficients:\n\n');
fprintf(fid, 'Intercept, Mismeasured variable, and Zetas\n');
fprintf(fid, 'standard deviation\n\n');
fprintf(fid, 'RMSE\n\n---------------------------------------\n\n');

% Pre-test
fprintf(fid, 'EW idtest result:\n\nTest statistic\n(p-value)\n\nSize\n\n\n');
i = 1; printchar='';printchar_='';
while i <= t
   printchar = strcat(printchar,' %6.4f \t');
   printchar_ = strcat(printchar_,'(%6.4f)\t');
   i = i+1;
end
printchar = strcat(printchar,'\n'); printchar_ = strcat(printchar_,'\n\n');
fprintf(fid,printchar,Wstat_m);
fprintf( fid,printchar_, Wstat_Pm );
fprintf(fid,printchar, Wtest_size);

if Level
    
    % EW Results - Level
    fprintf(fid,'\n\n\n*EW Estimation - LEVEL\n');
    
    %OLS results
    fprintf(fid,'\n*OLS\n');
    i = 1; printchar='';printchar_='';
    while i <= Xdim
       printchar = strcat(printchar,' %6.4f \t');
       printchar_ = strcat(printchar_,'(%6.4f)\t');
       i = i+1;
    end
    printchar = strcat(printchar,'\n'); printchar_ = strcat(printchar_,'\n\n');
    for T=1:t
        fprintf(fid,'\nt=%d\n', T);
        fprintf(fid,printchar,Beta_mOLSL( (T-1)*Xdim+1:T*Xdim) );
        fprintf( fid,printchar_, Beta_mSEOLSL( (T-1)*Xdim+1:T*Xdim) );
        fprintf(fid,printchar,RMSE_OLSL( (T-1)*Xdim+1:T*Xdim ));
    end

    %GMM results
    fprintf(fid,'\n\n*GMM Estimation\n\n');
    for T=1:t
        fprintf(fid,'\nt=%d\n',T);   %t1: GMM3, GMM4, ...
        for i=1:WD_order
            fprintf(fid,'GMM%d\n', i+2);
            j = 1; printchar='';printchar_='';
            while j <= Xdim
               printchar = strcat(printchar,' %8.4f \t');
               printchar_ = strcat(printchar_,'(%8.4f)\t');
               j = j+1;
            end
            printchar = strcat(printchar,'\n'); printchar_ = strcat(printchar_,'\n\n');
            fprintf(fid,printchar, Beta_mGMML(i,:,T));
            fprintf( fid,printchar_, Beta_mSEGMML(i,:,T));
            printchar = strcat(printchar,'\n\n');
            fprintf(fid,printchar, RMSE_GMML(i,:,T));
        end
        fprintf(fid,'J-statistic\n(p-value)\n\nSize\n\n');
        i = 1; printchar='';printchar_='';
        while i <= WD_order
           printchar = strcat(printchar,' %6.4f \t');
           printchar_ = strcat(printchar_,'(%6.4f)\t');
           i = i+1;
        end
        printchar = strcat(printchar,'\n'); printchar_ = strcat(printchar_,'\n\n');
        fprintf(fid, printchar, JJL( (T-1)*WD_order+1:T*WD_order ));
    end

    %Minimum distance
    fprintf(fid,'\n\n*MD Estimation\n\n');
    for i=1:WD_order
        fprintf(fid,'MD%d\n', i+2);
        j = 1; printchar='';printchar_='';
        while j <= Xdim
           printchar = strcat(printchar,' %6.4f \t');
           printchar_ = strcat(printchar_,'(%6.4f)\t');
           j = j+1;
        end
        printchar = strcat(printchar,'\n'); printchar_ = strcat(printchar_,'\n\n');
        fprintf(fid,printchar, Beta_mMDL( (i-1)*Xdim+1:i*Xdim ) );
        fprintf( fid,printchar_, Beta_mSEMDL( (i-1)*Xdim+1:i*Xdim ) );
        printchar = strcat(printchar,'\n\n');
        fprintf(fid,printchar, RMSE_MDL( (i-1)*Xdim+1:i*Xdim ) );
    end
else

    % EW Results - Within
    fprintf(fid,'\n\n\n*EW Estimation - WHITHIN\n');
    
    %OLS results
    fprintf(fid,'\n*OLS\n');
    i = 1; printchar='';printchar_='';
    while i <= Xdim
       printchar = strcat(printchar,' %6.4f \t');
       printchar_ = strcat(printchar_,'(%6.4f)\t');
       i = i+1;
    end
    printchar = strcat(printchar,'\n'); printchar_ = strcat(printchar_,'\n\n');
    for T=1:t
        fprintf(fid,'\nt=%d\n', T);
        fprintf(fid,printchar,Beta_mOLSW( (T-1)*Xdim+1:T*Xdim) );
        fprintf( fid,printchar_, Beta_mSEOLSW( (T-1)*Xdim+1:T*Xdim) );
        fprintf(fid,printchar,RMSE_OLSW( (T-1)*Xdim+1:T*Xdim ));
    end

    %GMM results
    fprintf(fid,'\n\n*GMM Estimation\n\n');
    for T=1:t
        fprintf(fid,'\nt=%d\n',T);   %t1: GMM3, GMM4, ...
        for i=1:WD_order
            fprintf(fid,'GMM%d\n', i+2);
            j = 1; printchar='';printchar_='';
            while j <= Xdim
               printchar = strcat(printchar,' %8.4f \t');
               printchar_ = strcat(printchar_,'(%8.4f)\t');
               j = j+1;
            end
            printchar = strcat(printchar,'\n'); printchar_ = strcat(printchar_,'\n\n');
            fprintf( fid,printchar, Beta_mGMMW(i,:,T) );
            fprintf( fid,printchar_, Beta_mSEGMMW(i,:,T) );
            printchar = strcat(printchar,'\n\n');
            fprintf( fid,printchar, RMSE_GMMW(i,:,T) );
        end
        fprintf(fid,'J-statistic\n(p-value)\n\nSize\n\n');
        i = 1; printchar='';printchar_='';
        while i <= WD_order
           printchar = strcat(printchar,' %6.4f \t');
           printchar_ = strcat(printchar_,'(%6.4f)\t');
           i = i+1;
        end
        printchar = strcat(printchar,'\n'); printchar_ = strcat(printchar_,'\n\n');
        fprintf(fid, printchar, JJW( (T-1)*WD_order+1:T*WD_order ));
    end

    %Minimum distance results
    fprintf(fid,'\n\n*MD Estimation\n\n');
    for i=1:WD_order
        fprintf(fid,'MD%d\n', i+2);
        j = 1; printchar='';printchar_='';
        while j <= Xdim
           printchar = strcat(printchar,' %6.4f \t');
           printchar_ = strcat(printchar_,'(%6.4f)\t');
           j = j+1;
        end
        printchar = strcat(printchar,'\n'); printchar_ = strcat(printchar_,'\n\n');
        fprintf(fid,printchar, Beta_mMDW( (i-1)*Xdim+1:i*Xdim ) );
        fprintf( fid,printchar_, Beta_mSEMDW( (i-1)*Xdim+1:i*Xdim ) );
        printchar = strcat(printchar,'\n\n');
        fprintf(fid,printchar, RMSE_MDW( (i-1)*Xdim+1:i*Xdim ) );
    end
end
status = fclose(fid);
%%

