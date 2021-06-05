function [Betas, BetaVars, J_stat, J_Pval] = genBetaWD(n, t, Xdim, D, WD_order, Xnames, Zetadim, yrs)
%FUNCTION: genBetaWD
%Version: 9/15/2009
%Performs Whited estimation
%INPUTS:    n, t, Xdim, Zetadim  - cross-section and time dimensions, and # of mismeasured and perfectly measured variables, respectively  
%           D = [Y xi Constant Zetas];
%           WD_ORDER - WD_ORDER+2 is the max order moment used in Whited estimation (also WD_ORDER is the # of estimators)
%           Xnames - names of the variables
%           yr - elements of the time dimension (for each n); n times t
%IMPORTANT: WHITED PROGRAM HANDLES ONLY ONE MISMEASURED

if strcmp(Xnames{2}(1:4), 'xaux') 
    display('ERROR: Whited procedure can only handle one mismeasured variable')
    return;
end    
yaux = D(:,1); 
xaux1 = D(:,2);                  %xaux1 = xi+xerror
intercep = D(:,3);
zeta=zeros(size(D,1),Zetadim);  %Handles one mismeasured variable
for i=4:4+Zetadim-1
    zeta(:,i-3) = D(:,i);
end

yname  ='yaux';     %left-hand-side variable
xname = Xnames{1};  %mismeasured variable
zname = char( Xnames{2:end} );

%Creates output (Whitted estimates)
% olsmc= zeros(t*2,2,Xdim);
% gmmc = zeros(t*2,2,Xdim*WD_order);
% mdmc = zeros(WD_order*2,2,Xdim*WD_order*2);

year = yrs;


%step 5:------set the number of time periods.-------------------------------

nyr=t;

%step 6:------set the first time period.-----------------------------------

fyear=1;

%step 7:------set the number of observations in the largest cross-section.-

nco=n;

%step 8:------is your data a balanced panel?-------------------------------

panel=1;

%step 9:------set the number of estimators.--------------------------------

nestim=WD_order;

%step 10:-----set the starting value.--------------------------------------

cstart=0;

%step 11:-----set the output mode.-----------------------------------------

latex=0;

%step 12:-------set the maximum number of iterations.----------------------

maxiter=299;

%step 13:-----------------set the limit on the number of squeezes.---------

maxsqez=240;

%step 14:-----------------set the tolerance on the minimization routine.---

tol=1e-5;

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%                           End of user input.
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% 
% --------open data set.-------------------------------------------------

% fid = fopen(outname, 'w');

%---------number of perfectly measured regressors.----------------------

nz=size(zname,1)-1;

%-------set the flag for the standard error calculations.----------

bleh=2;

%bleh = 1 acts as if Ey and Ex known.  bleh = 2 optimally accounts for the


%-------------------parameters-----------------------

%ols

isave   = zeros(2,nyr);          % Intercept
bsave   = zeros(2,nyr);          % Coefficient on chi
if nz>0;
    zsave   = zeros(2*nz,nyr);       % Coefficients on perfectly measured regressors
    rzsave  = zeros(2*nz,nyr);       % Reverse perfectly measured regressors
end;
rbsave  = zeros(2,nyr);          % Reverse chi
rsave   = zeros(2,nyr);          % R2

%gmm

asave   = zeros(2,nyr*nestim);     % Intercept
csave   = zeros(2,nyr*nestim);     % Coefficient on chi
if nz>0;
    dsave   = zeros(2*nz,nyr*nestim);  % Coefficients on perfectly measured regressors
end;
tsave   = zeros(2,nyr*nestim);     % tau2
psave   = zeros(2,nyr*nestim);     % rho2
chisave = zeros(2,nyr*nestim);     % J-statistic

%-------------------minimum distance estimates-------

if panel==1;

%    gmm

    casave  = zeros(4,nestim);      % Intercept
    ccsave  = zeros(4,nestim);      % Coefficient on chi
    if nz>0;
        cdsave  = zeros(4*nz,nestim);   % Coefficients on perfectly measured regressors
    end;
    ctsave  = zeros(4,nestim);      % tau2
    cpsave  = zeros(4,nestim);      % rho2

%    ols

    bisave  = zeros(4,1);           % Intercept
    bbsave  = zeros(4,1);           % Coefficient on chi
    brsave  = zeros(4,1);           % R2
    if nz>0;
        bzsave  = zeros(4*nz,1);        % Coefficients on perfectly measured regressors
    end;

end;

%-------------------influence functions--------------

%gmm

iasave=zeros(nco,nyr*nestim);
icsave=zeros(nco,nyr*nestim);
if nz>0;
    idsave=zeros(nco*nz,nyr*nestim);
end;
ipsave=zeros(nco,nyr*nestim);
itsave=zeros(nco,nyr*nestim);

%ols


iisave=zeros(nco,nyr);
ibsave=zeros(nco,nyr);
irsave=zeros(nco,nyr);
if nz>0;
    izsave=zeros(nco*nz,nyr);
end;



y = eval(yname);
x = eval(xname);
xyz = [x,y];

if nz>0;

    z = zeros(size(y,1),nz);
    for j=2:(nz+1);
        z(:,j-1) = eval(deblank(zname(j,:)));
        
    end;
    xyz = [xyz,z];
end;

yrname='year';
yr = eval(yrname);

year=1;
while year<=nyr;

    %---------read in the data and define the variables.----------

    dta=xyz(yr==(year+fyear-1),:);

    n=size(dta,1);

    y=dta(:,2);
    x=dta(:,1);
    if nz>0;
        z=dta(:,3:(2+nz));
    end;


    icept=ones(n,1);
    if nz>0;
        za=z;
        iz=[icept za];
    else;
        iz=icept;
    end;

    mux=(iz'*iz)\(iz'*x);
    muy=(iz'*iz)\(iz'*y);
    y_d=y-iz*muy;
    x_d=x-iz*mux;


    y2_d=y_d.^2;
    yx_d=(x_d.*y_d);
    x2_d=x_d.^2;
    y2x_d=(x_d.*y2_d);
    yx2_d=(x2_d.*y_d);
    y3x_d=y2x_d.*y_d;
    y2x2_d=y2_d.*x2_d;
    yx3_d=yx2_d.*x_d;


    x3_d=x_d.*x2_d;
    y3_d=y_d.*y2_d;
    y4x_d=y3_d.*yx_d;
    y3x2_d=y3_d.*x2_d;
    y2x3_d=y2_d.*x3_d;
    yx4_d=yx_d.*x3_d;
    y4_d=y2_d.^2;
    x4_d=x2_d.^2;
    y5x_d=y_d.*y4x_d;
    y4x2_d=y_d.*y3x2_d;
    y3x3_d=y_d.*y2x3_d;
    y2x4_d=y_d.*yx4_d;
    yx5_d=yx4_d.*x_d;


    y5_d=y2_d.*y3_d;
    x5_d=x2_d.*x3_d;
    y6x_d=y5_d.*yx_d;
    y5x2_d=y5_d.*x2_d;
    y4x3_d=y4_d.*x3_d;
    y3x4_d=y3_d.*x4_d;
    y2x5_d=y2_d.*x5_d;
    yx6_d=yx_d.*x5_d;


    %    VERY IMPORTANT:  The following variables MUST be in the same order that
    %                        they appear in the subroutine deff.
    %


    mom={y2_d yx_d x2_d y2x_d yx2_d y3x_d y2x2_d yx3_d x3_d y3_d y4x_d ...
        y3x2_d y2x3_d yx4_d y4_d x4_d y5x_d y4x2_d y3x3_d y2x4_d yx5_d ...
        y5_d x5_d y6x_d y5x2_d y4x3_d y3x4_d y2x5_d yx6_d};


    ey2_d=mean(y2_d);
    eyx_d=mean(yx_d);
    ex2_d=mean(x2_d);
    ey2x_d=mean(y2x_d);
    eyx2_d=mean(yx2_d);
    ey3x_d=mean(y3x_d);
    ey2x2_d=mean(y2x2_d);
    eyx3_d=mean(yx3_d);


    ex3_d=mean(x3_d);
    ey3_d=mean(y3_d);
    ey4x_d=mean(y4x_d);
    ey3x2_d=mean(y3x2_d);
    ey2x3_d=mean(y2x3_d);
    eyx4_d=mean(yx4_d);
    ey4_d=mean(y4_d);
    ex4_d=mean(x4_d);
    ey5x_d=mean(y5x_d);
    ey4x2_d=mean(y4x2_d);
    ey3x3_d=mean(y3x3_d);
    ey2x4_d=mean(y2x4_d);
    eyx5_d=mean(yx5_d);
    ey5_d=mean(y5_d);
    ex5_d=mean(x5_d);
    ey6x_d=mean(y6x_d);
    ey5x2_d=mean(y5x2_d);
    ey4x3_d=mean(y4x3_d);
    ey3x4_d=mean(y3x4_d);
    ey2x5_d=mean(y2x5_d);
    eyx6_d=mean(yx6_d);


    emom={ey2_d eyx_d ex2_d ey2x_d eyx2_d ey3x_d ey2x2_d eyx3_d ex3_d ...
        ey3_d ey4x_d ey3x2_d ey2x3_d eyx4_d ey4_d ex4_d ey5x_d ey4x2_d ...
        ey3x3_d ey2x4_d eyx5_d ey5_d ex5_d ey6x_d ey5x2_d ey4x3_d ...
        ey3x4_d ey2x5_d eyx6_d};


    % Make the moments for the second stage of the partialling estimator.


    ezx=iz'*x/n;

    % Make the moments for the rest of the summary statistics.


    ey=mean(y);
    ex=mean(x);
    ey2=mean(y.^2);
    ex2=mean(x.^2);
    if nz>0;
        ez=mean(z);
        ez2=mean(z.^2);
    end;


    % Make the moments for the correct standard errors for beta.



    y6_d=y_d.^6;
    x6_d=x_d.^6;

    idx = ones(nz+1,1);
    ey2_z   = mean(y2_d(:,idx).*iz)';
    eyx_z   = mean(yx_d(:,idx).*iz)';
    ex2_z   = mean(x2_d(:,idx).*iz)';
    ey2x_z  = mean(y2x_d(:,idx).*iz)';
    eyx2_z  = mean(yx2_d(:,idx).*iz)';
    ey3x_z  = mean(y3x_d(:,idx).*iz)';
    ey2x2_z = mean(y2x2_d(:,idx).*iz)';
    eyx3_z  = mean(yx3_d(:,idx).*iz)';
    ey3_z   = mean(y3_d(:,idx).*iz)';
    ex3_z   = mean(x3_d(:,idx).*iz)';
    ey4x_z  = mean(y4x_d(:,idx).*iz)';
    ey3x2_z = mean(y3x2_d(:,idx).*iz)';
    ey2x3_z = mean(y2x3_d(:,idx).*iz)';
    eyx4_z  = mean(yx4_d(:,idx).*iz)';
    ey4_z   = mean(y4_d(:,idx).*iz)';
    ex4_z   = mean(x4_d(:,idx).*iz)';
    ey5x_z  = mean(y5x_d(:,idx).*iz)';
    ey4x2_z = mean(y4x2_d(:,idx).*iz)';
    ey3x3_z = mean(y3x3_d(:,idx).*iz)';
    ey2x4_z = mean(y2x4_d(:,idx).*iz)';
    eyx5_z  = mean(yx5_d(:,idx).*iz)';
    ey5_z   = mean(y5_d(:,idx).*iz)';
    ex5_z   = mean(x5_d(:,idx).*iz)';
    ey6x_z  = mean(y6x_d(:,idx).*iz)';
    ey5x2_z = mean(y5x2_d(:,idx).*iz)';
    ey4x3_z = mean(y4x3_d(:,idx).*iz)';
    ey3x4_z = mean(y3x4_d(:,idx).*iz)';
    ey2x5_z = mean(y2x5_d(:,idx).*iz)';
    eyx6_z  = mean(yx6_d(:,idx).*iz)';
    ey6_z   = mean(y6_d(:,idx).*iz)';
    ex6_z   = mean(x6_d(:,idx).*iz)';


    ezz = (iz'*iz)/n;
    inezz = inv(ezz);
    zy_d = iz.*y_d(:,idx);
    zx_d = iz.*x_d(:,idx);


    zmom = {inezz zy_d zx_d eyx_z ey2_z ex2_z ey3_z ey2x_z eyx2_z ex3_z ey4_z ...
    ey3x_z ey2x2_z eyx3_z ex4_z ey5_z ey4x_z ey3x2_z ey2x3_z eyx4_z ex5_z ...
    ey6_z ey5x_z ey4x2_z ey3x3_z ey2x4_z eyx5_z ex6_z};

    if nz>0; eza=mean(za); end;


    % The OLS estimator:


    if nz>0;
        des=[icept x za];
    else;
        des=[icept x];
    end;
    dd=inv((des'*des));
    allb=dd*des'*y;
    isave(1,year)=allb(1,1);
    bsave(1,year)=allb(2,1);
    if nz>0;
        for qq=1:nz; zsave((2*qq-1),year)=allb(2+qq,1); end;
    end;
    uhat=y-des*allb;
    wh=uhat(:,ones(size(des,2),1)).*des;
    seb=sqrt(diag(dd*(wh'*wh)*dd));
    naiveseb=sqrt(diag(dd*mean(uhat.^2)));
    isave(2,year)=seb(1,1);
    bsave(2,year)=seb(2,1);
    if nz>0;
        for qq=1:nz; zsave((2*qq),year)=seb(2+qq,1); end;
    end;
    rsave(1,year)=1-(uhat'*uhat)/(((y-ey)'*(y-ey)));


    olsinflnc=inv((des'*des)/n)*((des.*uhat(:,ones(size(des,2),1)))');
    iisave(1:n,year)=olsinflnc(1,:)';
    ibsave(1:n,year)=olsinflnc(2,:)';
    if nz>0;
        for qq=1:nz; izsave(((qq-1)*n+1):(qq*n),year)=olsinflnc(2+qq,:)';
        end;
    end;


    %----------first make the influence function for sigma_xz.--------------------
    if nz>0;


        xza=[x za];
        exza=mean(xza);
        nreg=nz+1;
        sigxz=((xza)-exza(ones(size(xza,1),1),:))'*((xza)-exza(ones(size(xza,1),1),:))/n;




        vecsigxz=zeros(nreg+nreg*(nreg-1)/2,1);
        vecsigxz(1:nreg,1)=diag(sigxz);


        counter=nreg+nreg*(nreg-1)/2;
        for qq=nreg:(-1):2;
            for ee=(qq-1):(-1):1;
                vecsigxz(counter,1)=sigxz(ee,qq);
                counter=counter-1;
            end;
        end;




        phixz=zeros(n,nreg+nreg*(nreg-1)/2);
        phixz(:,1:nreg)=((xza)-exza(ones(size(xza,1),1),:)).*((xza)-exza(ones(size(xza,1),1),:));


        counter=nreg+nreg*(nreg-1)/2;
        for qq=nreg:(-1):2;
            for ee=(qq-1):(-1):1;
                phixz(:,counter)=(xza(:,ee)-exza(ones(size(xza,1),1),ee)).*(xza(:,qq)-exza(ones(size(xza,1),1),qq));
                counter=counter-1;
            end;
        end;


        sigy=((y-ey)'*(y-ey))/n;
        phiy=(y-ey).^2-sigy;


        bigphi=[(olsinflnc(2:size(olsinflnc,1),:)); (phixz'); (phiy')];


        %--------------now make the derivative matrix.-------------------------------


        gee=zeros(nreg+(nreg+nreg*(nreg-1)/2)+1,1);




        %-------------first,the derivative wrt the ols coefficients.----------


        for qq=1:nreg;


            gee(qq,1)=(2*allb(2:nreg+1,1)'*sigxz(:,qq))/sigy;


        end;


        %------------derivatives wrt the first part of sigxz--


        for qq=1:nreg;


            gee(nreg+qq,1)=(allb(qq+1,1).^2)/sigy;


        end;


        %------------derivatives wrt the second part of sigxz--


        counter=nreg+(nreg+nreg*(nreg-1)/2);
        for qq=nreg:(-1):2;
            for ee=(qq-1):(-1):1;
                gee(counter,1)=2*allb(ee+1,1)*allb(qq+1,1)/sigy;
                counter=counter-1;
            end;
        end;


        %------------derivatives wrt sigy--


        counter=nreg+(nreg+nreg*(nreg-1)/2);


        gee(counter+1,1)=-rsave(1,year)/sigy;


        %-----------save the influence function.-------------


    else;


        %----------first make the influence function for sigma_xz.--------------------
        sigx=ex2-ex.^2;
        phix=(x-ex).^2-sigx;


        sigy=((y-ey)'*(y-ey))/n;
        phiy=(y-ey).^2-sigy;


        bigphi=[(olsinflnc(2:size(olsinflnc,1),:)); (phix'); (phiy')];


        %--------------now make the derivative matrix.--------------------------------


        gee=zeros(3,1);


        %-------------first,the derivative wrt the ols coefficients.----------


        gee(1,1)=(2*allb(2,1)*sigx)/sigy;


        %------------derivatives wrt the first part of sigxz--


        gee(2,1)=(allb(2,1).^2)/sigy;


        %------------derivatives wrt sigy--


        gee(3,1)=-rsave(1,year)/sigy;


        %-----------save the influence function.-------------


    end;


    irsave(1:n,year)=-bigphi'*gee;
    rsave(2,year)=sqrt(irsave(:,year)'*irsave(:,year)/(n.^2));




    if nz>0;
        des=[icept y za];
    else;
        des=[icept y];
    end;
    dd=inv((des'*des));
    allbr=dd*des'*x;
    rbsave(1,year)=1/allbr(2,1);
    if nz>0;
        for qq=1:nz; rzsave((2*qq-1),year)=-allbr(qq+2,1)/allbr(2,1);
        end;
    end;
    uhat=x-des*allbr;
    wh=uhat(:,ones(size(des,2),1)).*des;
    vc=(dd*(wh'*wh)*dd);
    gee=zeros(nz+2,nz+2);


    gee(1,1)=-1/allbr(2,1);
    gee(1,2)=allbr(1,1)/(allbr(2,1).^2);
    gee(2,2)=-allbr(2,1).^(-2);
    for qq=1:nz; gee(qq+2,qq+2)=-1/allbr(2,1); end;
    for qq=1:nz; gee(qq+2,2)=allbr(qq+2,1)/allbr(2,1); end;


    vcr2t=gee*vc*gee';
    serb=sqrt(diag(vcr2t));
    rbsave(2,year)=serb(2,1);
    if nz>0;
        for qq=1:nz; rzsave((2*qq),year)=serb(qq+1,1); end;
    end;

    % The GMM estimators:


    estim=1;
    while estim<=nestim;

        %--once exaclly identified and four iterative estimators.-----------------


        if estim==1;
            nc=5;
            neq=5;
        elseif estim==2;
            nc=6;
            neq=8;
        elseif estim==3;
            nc=9;
            neq=14;
        elseif estim==4;
            nc=12;
            neq=21;
        elseif estim==5;
            nc=15;
            neq=29;
        end;
        %----------here we input the starting values.----------------
        c=zeros(nc,1);


        c(1,1)=ey2x_d/eyx2_d;
        if cstart~=0 and estim>1; c(1,1)=cstart; end;
        c(2,1)=eyx_d/c(1,1);
        c(3,1)=ey2_d-(c(1,1).^2)*c(2,1);
        c(4,1)=ex2_d-c(2,1);
        c(5,1)=eyx2_d/c(1,1);
        if estim>1;
            c(6,1)=(eyx3_d-3*c(1,1)*c(2,1)*c(4,1))/c(1,1);
        end;
        if estim>2;
            c(7,1)=ex3_d-c(5,1);
            c(8,1)=ey3_d-(c(1,1).^3)*c(5,1);
            c(9,1)=(eyx4_d/c(1,1))-6*c(5,1)*c(4,1)-4*c(2,1)*c(7,1);
        end;
        if estim>3;
            c(10,1)=ey4_d-(c(1,1).^4)*c(6,1)-6*(c(1,1).^2)*c(2,1)*c(3,1);
            c(11,1)=ex4_d-c(6,1)-6*c(2,1)*c(4,1);
            c(12,1)=(eyx5_d/c(1,1))-10*c(6,1)*c(4,1)-10*c(5,1)*c(7,1)...
                -5*c(2,1)*c(11,1);
        end;
        if estim>4;
            c(13,1)=ey5_d-(c(1,1).^5)*c(9,1)-10*(c(1,1).^3)*c(5,1)*c(3,1)...
                -10*(c(1,1).^2)*c(2,1)*c(8,1);
            c(14,1)=ex5_d-c(9,1)-10*c(5,1)*c(4,1)-10*c(2,1)*c(7,1);
            c(15,1)=(eyx6_d/c(1,1))-15*c(9,1)*c(4,1)-20*c(6,1)*c(7,1)...
                -15*c(5,1)*c(11,1)-6*c(2,1)*c(14,1);
        end;


        if estim>1;


        %--------the program creates the weighting matrix.------------


            ff=optw(c,mom,emom,zmom,n,neq,estim,bleh);

            w=inv((ff'*ff)./n);


        %-------------start of iteration loop------------------------
            iter=1;
            dc=1;                      % Initialize the step length.
            while norm(dc,inf)>=tol;


                f=deff(c,emom,estim,neq,0);          % The program jumps to the subroutine that
                                                     %defines the f vector.


                if estim==2;
                    g=grad2(c,neq);                 % The program jumps to the subroutine that
                                                    %computes analytic derivatives.the matrix
                                                    %of partials of f with respect to the
                                                    %parameters is called g.
                elseif estim==3;
                    g=grad3(c,neq);
                elseif estim==4;
                    g=grad4(c,neq);
                elseif estim==5;
                    g=grad5(c,neq);
                end;


                obj = f'*w*f;               % This computes the value of the objective
                                            %function.


                gwg= g'*w*g;                % This uses the GAUSS-NEWTON method
                                            %to compute full step dc.
                gwf=g'*w*f;


                dc=gwg\gwf;


                if maxsqez > 0;               % This jumps to the subroutine that
                                              %adjusts the step length.
                    [c_new,sqz]=squeez(c,dc,emom,maxsqez,obj,estim,neq,w);
                else;
                    c_new=c-dc;
                end;


                dc=c_new-c;                      % Update variables for the next iteration.
                c=c_new;
                %--------------Print out the results for the current iteration------------
                %
                %"Number of Squeezes = "; ;  sqz;
                %"Objective Function = "; ;  obj*n;
                %
                %"            Value         Step";
                % format 14,6;
                % mm=1;
                % do until mm > nc;
                %  $"  "; ;  c[mm,1]; ;  "  "; ;  dc[mm,1];
                %  mm=mm+1;
                %  endo;
                % -----------------End of print out of current iteration-------------------
                %
                iter=iter+1;
                if iter >= maxiter;                      % Quit iterating if necessary.
                    break;
                end;
            end; %------end of iteration loop------------------------------------

        end; %---------end of exactly identified estimator condition-----------------


        %compute t-ratios,etc.


        ff=optw(c,mom,emom,zmom,n,neq,estim,bleh);
        f=deff(c,emom,estim,neq,0);


        w=inv((ff'*ff)./n);


        if estim==1;
            g=grad1(c,neq);
        elseif estim==2;
            g=grad2(c,neq);
        elseif estim==3;
            g=grad3(c,neq);
        elseif estim==4;
            g=grad4(c,neq);
        elseif estim==5;
            g=grad5(c,neq);
        end;


        gwg=g'*w*g;
        vc=inv(gwg)/n;
        stderr=sqrt(diag(vc));
        inflnc=-n*vc*g'*w*ff';


        csave(1,nestim*(year-1)+estim)=c(1,1);
        csave(2,nestim*(year-1)+estim)=stderr(1,1);
        icsave(1:n,nestim*(year-1)+estim)=inflnc(1,:)';


        if estim>1;
            chisq=obj*n;
            chisave(1,nestim*(year-1)+estim)=chisq;
            degf=neq-nc;
            chisave(2,nestim*(year-1)+estim)=1-chi2cdf(chisq,degf);
        end;


        % Make Delta


        c11=c(1,1);
        if nz>0;


            for qq=1:nz;
                dsave((qq*2-1),nestim*(year-1)+estim)=muy(qq+1,1)-c11*mux(qq+1,1);
            end;


            gee=zeros(((nz+1)*2+1),nz);


            gee(2:(nz+1),1:nz)=eye(nz);
            gee((nz+3):((nz+1)*2),1:nz)=-c11*eye(nz);


            for qq=1:nz; gee(((nz+1)*2+1),qq)=-mux(qq+1,1); end;


            %---correct standard error and influence function---


            bigphi=[([(inezz*zy_d'); (inezz*zx_d')]); (-inflnc(1,:))];
            avar=(bigphi*bigphi')./n.^2;
            dstd=sqrt(diag(gee'*avar*gee));
            omega=gee'*avar*gee;
            phidel=-bigphi'*gee;


            for qq=1:nz; dsave(qq*2,nestim*(year-1)+estim)=dstd(qq,1); end;


            for qq=1:nz;
                idsave(((qq-1)*n+1):(qq*n),nestim*(year-1)+estim)=phidel(:,qq);
            end;


        end;


        % Make the stupid intercept


        asave(1,nestim*(year-1)+estim)=muy(1,1)-c11*mux(1,1);


        gee=zeros(((nz+1)*2+1),1);


        gee(1,1)=1;
        gee((nz+2),1)=-c11;
        gee(((nz+1)*2+1),1)=-mux(1,1);


        %---correct standard error and influence function---


        bigphi=[([(inezz*zy_d'); (inezz*zx_d')]); (-inflnc(1,:))];
        avar=(bigphi*bigphi')./n.^2;
        astd=sqrt(diag(gee'*avar*gee));
        phidela=-bigphi'*gee;


        asave(2,nestim*(year-1)+estim)=astd(1,1);
        iasave(1:n,nestim*(year-1)+estim)=phidela(:,1);




        %--------------rho2 and tau2--------------------------------------------------




        if nz>0;


            %----------first make the influence function for sigma_z.---------------------


            sigz=(za-kron(ones(n,1),mean(za)))'*(za-kron(ones(n,1),mean(za)))/n;




            vecsigz=zeros(nz+nz*(nz-1)/2,1);
            vecsigz(1:nz,1)=diag(sigz);


            counter=nz+nz*(nz-1)/2;
            for qq=nz,(-1):2;
                for ee=(qq-1):(-1):1;
                    vecsigz(counter,1)=sigz(ee,qq);
                    counter=counter-1;
                end;
            end;




            phiz=zeros(n,nz+nz*(nz-1)/2);
            phiz(:,1:nz)=(za-kron(ones(n,1),eza)).*(za-kron(ones(n,1),eza));


            counter=nz+nz*(nz-1)/2;
            for qq=nz,(-1):2;
                for ee=(qq-1):(-1):1;
                    vecsigz(counter,1)=sigz(ee,qq);
                    phiz(:,counter)=(za(:,ee)-kron(ones(n,1),eza(1,ee))).*(za(:,qq)-kron(ones(n,1),eza(1,qq)));
                    counter=counter-1;
                end;
            end;




            phiz=phiz-kron(ones(n,1),vecsigz');


            phimuy=inezz*zy_d';
            phimux=inezz*zx_d';


            numer=(muy(2:nz+1,1)'*sigz*muy(2:nz+1,1)+c(1,1).^2*c(2,1));
            rho=numer/(numer+c(3,1));


            numer=(mux(2:nz+1,1)'*sigz*mux(2:nz+1,1)+c(2,1));
            tau=numer/(numer+c(4,1));


            %--make the influence functions for the standard errors for rho2 and tau2.


            bigphi=[phimux(2:nz+1,:); phimuy(2:nz+1,:); (phiz'); (-inflnc(1:4,:))];
            avar=(bigphi*bigphi')./n.^2;


            %--------------first column is for rho2 and the second is for tau2.


            gee = zeros(2*nz+(nz+nz*(nz-1)/2)+4,2);




            %-----------first,do rho2---------------------------




            numer=muy(2:nz+1,1)'*sigz*muy(2:nz+1,1)+c(1,1).^2*c(2,1);
            denom=numer+c(3,1); ;


            %------------derivatives wrt muy---------------------


            for qq=1:nz;


                gee(nz+qq,1)= (2*muy(2:nz+1,1)'*sigz(:,qq))/denom - numer*(2*muy(2:nz+1,1)'*sigz(:,qq))/(denom.^2);


            end;

            %------------derivatives wrt the first part of sigz--


            for qq=1:nz;


                gee(2*nz+qq,1)=(muy(qq+1,1).^2)/denom-numer/(denom.^2)*muy(qq+1,1).^2;


            end;


            %------------derivatives wrt the second part of sigz--


            counter=2*nz+(nz+nz*(nz-1)/2);
            for qq=nz,(-1):2;
                for ee=(qq-1):(-1):1;
                    gee(counter,1)= 2*muy(ee+1,1)*muy(qq+1,1)/denom-2*numer/(denom.^2)*muy(ee+1,1)*muy(qq+1,1);
                    counter=counter-1;
                end;
            end;


            %------------derivatives wrt c--


            counter=2*nz+(nz+nz*(nz-1)/2);


            gee(counter+1,1)=2*c(1,1)*c(2,1)/denom-2*numer/(denom.^2)*c(1,1)*c(2,1);


            gee(counter+2,1)=(c(1,1).^2)/denom-numer/(denom.^2)*c(1,1).^2;


            gee(counter+3,1)=-numer/(denom.^2);




            %-----------now for tau2------------------------




            numer=mux(2:nz+1,1)'*sigz*mux(2:nz+1,1)+c(2,1);
            denom=numer+c(4,1); ;




            %------------derivatives wrt mux---------------------


            for qq=1:nz;


                gee(qq,2)=...
                    (2*mux(2:nz+1)'*sigz(:,qq))/denom-numer*(2*mux(2:nz+1)'*sigz(:,qq))/(denom.^2);


            end;


            %------------derivatives wrt the first part of sigz--


            for qq=1:nz;


                gee(2*nz+qq,2)=(mux(qq+1,1).^2)/denom-numer/(denom.^2)*mux(qq+1,1).^2;


            end;


            %------------derivatives wrt the second part of sigz--


            counter=2*nz+(nz+nz*(nz-1)/2);
            for qq=nz,(-1):2;
                for ee=(qq-1):(-1):1;
                    gee(counter,2)=...
                        2*mux(ee+1,1)*mux(qq+1,1)/denom-2*numer/(denom.^2)*mux(ee+1,1)*mux(qq+1,1);
                    counter=counter-1;
                end;
            end;


            gee(2*nz+(nz+nz*(nz-1)/2)+2,2)=1/denom-numer/(denom.^2);
            gee(2*nz+(nz+nz*(nz-1)/2)+4,2)=-numer/(denom.^2);






        else;


            numer=(c(1,1).^2*c(2,1));
            rho=numer/(numer+c(3,1));


            numer=(c(2,1));
            tau=numer/(numer+c(4,1));


            %-----make the influence functions for the standard errors for rho2 and tau2.


            bigphi=(-inflnc(1:4,:));
            avar=(bigphi*bigphi')./n.^2;


            gee = zeros(4,2);  % First column is for rho2 and the second is for tau2.




            %-----------first,do rho2---------------------------




            numer=c(1,1).^2*c(2,1);
            denom=numer+c(3,1); ;


            %------------derivatives wrt c--




            gee(1,1)=2*c(1,1)*c(2,1)/denom-2*numer/(denom.^2)*c(1,1)*c(2,1);


            gee(2,1)=(c(1,1).^2)/denom-numer/(denom.^2)*c(1,1).^2;


            gee(3,1)=-numer/(denom.^2);




            %-----------now for tau2------------------------




            numer=c(2,1);
            denom=numer+c(4,1); ;




            gee(2,2)=1/denom-numer/(denom.^2);
            gee(4,2)=-numer/(denom.^2);




        end;


        vcrhotau=gee'*avar*gee;


        phirho=-bigphi'*gee(:,1);
        phitau=-bigphi'*gee(:,2);


        ipsave(1:n,nestim*(year-1)+estim)=phirho;
        itsave(1:n,nestim*(year-1)+estim)=phitau;




        % Now save both things.


        psave(1,nestim*(year-1)+estim)=rho;
        tsave(1,nestim*(year-1)+estim)=tau;


        psave(2,nestim*(year-1)+estim)=sqrt(vcrhotau(1,1));
        tsave(2,nestim*(year-1)+estim)=sqrt(vcrhotau(2,2));




        estim=estim+1;
    end;                        % The end of the estimator loop.


    year=year+1;


end;      %This ends the trial loop.


if panel==1;


%    This part does the classical minimum distance estimation.


    param = 1;  % 1 is beta, 2 is rho2, 3 is tau2,
                % 4 to nz+3 are the alphas, and the end is the intercept.
    while param<=nz+4;


        idx=(1:nestim:(nyr*nestim-1))-1;

        estim=1;
        while estim<=nestim;


            idx=idx+1;


            if param==1;
                saveit=csave;
                w=inv((icsave(:,idx)'*icsave(:,idx))/n);
            elseif param==2;
                saveit=psave;
                w=inv((ipsave(:,idx)'*ipsave(:,idx))/n);
            elseif param==3;
                saveit=tsave;
                w=inv((itsave(:,idx)'*itsave(:,idx))/n);
            elseif param==nz+4;
                saveit=asave;
                w=inv((iasave(:,idx)'*iasave(:,idx))/n);
            end;


            for qq=1:nz;
                if param==qq+3;
                    saveit=dsave(qq*2-1,:);
                    w=inv(idsave(((qq-1)*n+1):(qq*n),idx)'*idsave(((qq-1)*n+1):(qq*n),idx)/n);
                end;
            end;




            g=ones(nyr,1);
            theta=inv(g'*w*g)*g'*w*saveit(1,idx)';


            f=saveit(1,idx)'-theta;
            stderror=sqrt(inv(g'*w*g)/n);
            cmdtest=n*f'*w*f;


            if param==1; ccsave(1,estim)=theta; ccsave(2,estim)=cmdtest;
                ccsave(3,estim)=1-chi2cdf(cmdtest,3); ccsave(4,estim)=stderror;
            elseif param==2; cpsave(1,estim)=theta; cpsave(2,estim)=cmdtest;
                cpsave(3,estim)=1-chi2cdf(cmdtest,3); cpsave(4,estim)=stderror;
            elseif param==3; ctsave(1,estim)=theta; ctsave(2,estim)=cmdtest;
                ctsave(3,estim)=1-chi2cdf(cmdtest,3); ctsave(4,estim)=stderror;
            elseif param==nz+4; casave(1,estim)=theta; casave(2,estim)=cmdtest;
                casave(3,estim)=1-chi2cdf(cmdtest,3); casave(4,estim)=stderror;
            end;


            for qq=1:nz;
                if param==qq+3;
                    cdsave((qq-1)*4+1,estim)=theta;
                    cdsave((qq-1)*4+2,estim)=cmdtest;
                    cdsave((qq-1)*4+3,estim)=1-chi2cdf(cmdtest,3);
                    cdsave((qq-1)*4+4,estim)=stderror;
                end;
            end;


            estim=estim+1;
        end;


        param=param+1;
    end;






%    This part does the classical minimum distance estimation.
%    FOR OLS


    param = 1;  % 1 is beta, 2 through nz+1 are the alphas,
    %          nz+2 is r2, and the end is the intercept
    while param<=nz+3;


        if param==1;
            saveit=bsave(1,:);
            w=inv((ibsave'*ibsave)/n);
        elseif param==nz+2;
            saveit=rsave(1,:);
            w=inv((irsave'*irsave)/n);
        elseif param==nz+3;
            saveit=isave(1,:);
            w=inv((iisave'*iisave)/n);
        end;




        for qq=1:nz;
            if param==qq+1;
                saveit=zsave((2*qq-1),:);
                w=inv(izsave(((qq-1)*n+1):(qq*n),:)'*izsave(((qq-1)*n+1):(qq*n),:)/n);
            end;
        end;


        g=ones(nyr,1);
        theta=inv(g'*w*g)*g'*w*saveit';


        %compute t-ratios,etc.


        f=saveit'-theta;
        stderror=sqrt(inv(g'*w*g)/n);
        cmdtest=n*f'*w*f;


        if param==1;
            bbsave(1,1)=theta;
            bbsave(2,1)=cmdtest;
            bbsave(3,1)=1-chi2cdf(cmdtest,3);
            bbsave(4,1)=stderror;
        elseif param==nz+2;
            brsave(1,1)=theta;
            brsave(2,1)=cmdtest;
            brsave(3,1)=1-chi2cdf(cmdtest,3);
            brsave(4,1)=stderror;
        elseif param==nz+3;
            bisave(1,1)=theta;
            bisave(2,1)=cmdtest;
            bisave(3,1)=1-chi2cdf(cmdtest,3);
            bisave(4,1)=stderror;
        end;


        for qq=1:nz;
            if param==qq+1;
                bzsave((qq-1)*4+1,1)=theta;
                bzsave((qq-1)*4+2,1)=cmdtest;
                bzsave((qq-1)*4+3,1)=1-chi2cdf(cmdtest,3);
                bzsave((qq-1)*4+4,1)=stderror;
            end;
        end;


        param=param+1;
    end;


end;

%VERY IMPORTANT: OUTPUTS IN DIFFERENT ORDER FOR ZETAS
%OLS - natural order;intercept, xi, zetas for each time period
%GMM - intercepts (for GMM3, GMM4,...), xis(for GMM3, GMM4,...),BUT zetas for GMM3, THEN zetas for GMM4,...)
%MDMC - natural order; intercept, xi, zetas. ->MDMC3, MDMC4,...
for T=1:t
    seq = (T-1)*(WD_order)+1:T*(WD_order);
    olsmc( :,(T-1)*(Xdim)+1:T*(Xdim) ) = [isave(:,T),bsave(:,T),reshape(zsave(:,T),2,Zetadim)];        
    gmmc( :,(T-1)*(Xdim*WD_order)+1:T*(Xdim*WD_order) ) = [asave(:,seq),csave(:,seq),...
        reshape(dsave(:,seq),2,Zetadim*WD_order)];
end
% for i=1:WD_order
%     mdmc( :,(i-1)*(Xdim)+1:i*(Xdim) ) = [casave(1:3:end,i), ccsave(1:3:end,i),cdsave(1:3:4,i),cdsave(5:3:8,i),cdsave(9:3:12,i)];
% end
mdmc = zeros(2,WD_order + WD_order + Zetadim * WD_order);
for i=1:WD_order
    mdmc( 1,(i-1)*Xdim+1:i*Xdim ) = [casave(1,i), ccsave(1,i), cdsave(1,i),zeros(1,Xdim-3)];
    mdmc( 2,(i-1)*Xdim+1:i*Xdim ) = [casave(4,i), ccsave(4,i), cdsave(4,i),zeros(1,Xdim-3)];
    for j=1:Zetadim-1
        mdmc( 1,i*Xdim-(Zetadim-1)+j) = cdsave(j*4+1,i);
        mdmc( 2,i*Xdim-(Zetadim-1)+j) = cdsave(j*4+4,i);
    end
end

Betas = [olsmc(1,:) gmmc(1,:) mdmc(1,:)];
BetaVars = [olsmc(2,:) gmmc(2,:) mdmc(2,:)];
J_stat = chisave(1,:);
J_Pval = chisave(2,:);