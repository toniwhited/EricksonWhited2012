new ,60000;
ttt = time;
ddd = date;
call sysstate(14,1e-128);

fe = 0;
assets = 1;
nestim = 3;


nboot = 1000; @ number of bootstrap trials@

@ Step 1: ------Name of the Gauss data set.----------------------------------@

    dataset = "wdset";

@ Step 2: ------Name of the output file.-------------------------------------@

    outname = "ewboot";
    fmat="%lf";
    cnestim = nestim + 2;
    ce=ftos(cnestim,fmat,1,0);
    cf=ftos(fe,fmat,1,0);
    ca=ftos(assets,fmat,1,0);

    outname = outname $+ "_FE" $+ cf $+ "_ES" $+ ce $+ "_ASSET" $+ ca $+ ".out";
    output file = ^outname reset;


     if assets == 0;

@ Step 3: ------Name of the left-hand-side variable.-------------------------@

    let yname = ik;

@ Step 4: ------Name of the mismeasured variable.----------------------------@

    let xname = q;

@ Step 4: ------Names of the perfectly measured variables.-------------------@

    let zname = intercep cfk;

     elseif assets == 1;

@ Step 3: ------Name of the left-hand-side variable.-------------------------@

    let yname = ia;

@ Step 4: ------Name of the mismeasured variable.----------------------------@

    let xname = mb;

@ Step 4: ------Names of the perfectly measured variables.-------------------@

    let zname = intercep cfa;

    endif;


    let yrname = fyear;
@ Step 5: ------Set the number of time period.-------------------------------@

    nyr = 42;

@ Step 6: ------Set the first time period.-----------------------------------@

    fyear = 1967;


@ Step 6a: ---This shouldn't be changed, but I had to move the data reading
              because this is a bootstrap program.---------------------------@

    closeall f1;
    open f1=^dataset;
    aalldta=readr(f1,rowsf(f1));
    vnames = getname(dataset);
    aalldta = selif(aalldta,aalldta[.,indcv(yrname,vnames)].>=fyear .and aalldta[.,indcv(yrname,vnames)].<=(fyear+nyr-1));


@ Step 7: ------Set the number of observations in the largest cross-section.-@

    nco = 4100;

@ Step 8: ------Is your data a balanced panel?-------------------------------@

    panel = 0;

@ Step 9: ------Set the number of estimators.--------------------------------@

   @ I moved this up front so that it would be easier to rename the output file.@

@ Step 10: -----Set the starting value.--------------------------------------@

    cstart = 0;

@ Step 11: -----Set the output mode.-----------------------------------------@

    latex = 1;

@ Step 11: -----Set the output mode.-----------------------------------------@

    let idname = gvkey;
    isnumer    = 1;

@ Step 13: -------Set the maximum number of iterations.----------------------@

    maxiter = 299;

@ Step 14: -----------------Set the limit on the number of squeezes.---------@

    maxsqez=240;

@ Step 15: -----------------Set the tolerance on the minimization routine.---@

    tol=1e-9;


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*                           End of user input.                             */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


@--------Open GAUSS data set. -------------------------------------@

    outwidth  256;

@---------Number of perfectly measured regressors.----------------------@

    nz = rows(zname)-1;

@-------Set the flag for the standard error calculations.----------@

    bleh = 2;

/* bleh = 1 acts as if Ey and Ex known.  bleh = 2 optimally accounts for the
   fact that they are not known. */

@-------------------Parameters-----------------------@

   nsave   = zeros(1,nyr);
@OLS@

   isave   = zeros(2,nyr);         @ Intercept                             @
   bsave   = zeros(2,nyr);         @ Coefficient on chi                    @
   if nz>0;
   zsave   = zeros(2*nz,nyr);      @ Coefficients on
                                     perfectly measured regressors         @
   rzsave  = zeros(2*nz,nyr);      @ Reverse perfectly measured regressors @
   endif;
   rbsave  = zeros(2,nyr);         @ Reverse chi                           @
   rsave   = zeros(2,nyr);         @ R2                                    @

@GMM@

   asave   = zeros(2,nyr*nestim);    @ Intercept                           @
   csave   = zeros(2,nyr*nestim);    @ Coefficient on chi                  @
   if nz>0;
   dsave   = zeros(2*nz,nyr*nestim); @ Coefficients
                                       on perfectly measured regressors    @
   endif;
   tsave   = zeros(2,nyr*nestim);    @ tau2                                @
   psave   = zeros(2,nyr*nestim);    @ rho2                                @
   chisave = zeros(2,nyr*nestim);    @ J-statistic                         @

@-------------------Minimum Distance Estimates-------@


@GMM@

   casave  = zeros(4,nestim);     @ Intercept                              @
   ccsave  = zeros(4,nestim);     @ Coefficient on chi                     @
   if nz > 0;
   cdsave  = zeros(4*nz,nestim);  @ Coefficients on
                                    perfectly measured regressors          @
   endif;
   ctsave  = zeros(4,nestim);     @ tau2                                   @
   cpsave  = zeros(4,nestim);     @ rho2                                   @

@OLS@

   bisave  = zeros(4,1);          @ Intercept                              @
   bbsave  = zeros(4,1);          @ Coefficient on chi                     @
   brsave  = zeros(4,1);          @ R2                                     @
   if nz > 0;
   bzsave  = zeros(4*nz,1);       @ Coefficients on
                                    perfectly measured regressors          @
   endif;


@-------------------Influence Functions--------------@

@GMM@

   iasave  = zeros(nco,nyr*nestim);
   icsave  = zeros(nco,nyr*nestim);
   if nz > 0;
   idsave  = zeros(nco*nz,nyr*nestim);
   endif;
   ipsave  = zeros(nco,nyr*nestim);
   itsave  = zeros(nco,nyr*nestim);

@OLS@

   iisave  = zeros(nco,nyr);
   ibsave  = zeros(nco,nyr);
   irsave  = zeros(nco,nyr);
   if nz > 0;
   izsave  = zeros(nco*nz,nyr);
   endif;

@-----------------The saver for the "population" GMM objective function.-----------------------@
   effsave = zeros(29,nyr*nestim); @ 29 is the number of moment conditions for GMM7@


@-----------------The savers for the pivotal statistics.-----------------------@
   cpivot = zeros(nestim,nboot);
   apivot = zeros(nestim,nboot);
   dpivot = zeros(nz*nestim,nboot);
   ppivot = zeros(nestim,nboot);
   tpivot = zeros(nestim,nboot);

   bpivot = zeros(1,nboot);
   zpivot = zeros(nz,nboot);
   rpivot = zeros(1,nboot);
   ipivot = zeros(1,nboot);

cls;

/*   This part makes the index numbers for the different firms so that you can bootstrap by firm and still pick out the right observations. */
gvkeyall = aalldta[.,indcv(idname,vnames)];
if isnumer == 0;
  gvkeyall = stof(gvkeyall);
  isnumer = 1;
endif;
gvsall = unique(gvkeyall,isnumer);
nobsall = rows(gvsall);
pickco = miss(1,1)*zeros(nyr,nobsall);
nums = 1;
for iii(1,rows(gvsall),1);

dtai = selif(gvkeyall,gvkeyall.==gvsall[iii,1]);
pickco[1:rows(dtai),iii] = seqa(nums,1,rows(dtai));
nums = nums + rows(dtai);
endfor;

s1 = 1;
phat = ones(nobsall,1)/nobsall;
phatcdf=cumsumc(phat);
nums=seqa(1,1,rows(phatcdf));


bb = 0;
clear obj;
do while bb <= nboot;

output off;

"I am on bootstrap trial " bb;
@
"Bootstrap intermediate results";
print;
"Beta";
"     OLS:   ";
print;
"     GMM:   ";
print;
"Alpha_1";
"     OLS:   ";
print;
"     GMM:   ";
print;
print;
@
output on;
@----------------The saver for the identifiers.----@

bttdb:
   gvsave = -999*ones(nco,nyr);

   if bb == 0;
      alldta = aalldta;
   else;
      bootpick = rndus(nobsall,1,s1)';
      bootpick = maxc(ones(1,nobsall)|floor(nobsall*bootpick));

      rr=rndus(nobsall,1,s1);
      rr = subscat(rr,phatcdf,nums);

      pickout = pickco[.,rr];
      newgv = nums' + 0*pickout;
      pickout = packr(vec(pickout));
      newgv   = packr(vec(newgv));
      alldta = aalldta[pickout,.];
      alldta[.,1]=newgv;

   endif;

   gvkey = alldta[.,indcv(idname,vnames)];
   if isnumer == 0;
     gvkey = stof(gvkey);
     isnumer = 1;
   endif;
   gvs = unique(gvkey,isnumer);
   nobs = rows(gvs);




   if nz > 0;
   allvname = yname|xname|zname[2:rows(zname)];
   else;
   allvname = yname|xname;
   endif;
   indx = indcv(allvname,vnames);
   yindx = indx[1,1];
   xindx = indx[2,1];
   if nz > 0;
   zindx = indx[3:rows(indx),1];
   endif;
   yrindx = indcv(yrname,vnames);


   alldta = alldta[.,(1|yrindx|indx)];
   yrindx = 2;
   yindx = 3;
   xindx = 4;
   zindx = 5;
   let indx = 3 4 5;

   if fe==1;
   newdta = zeros(1,cols(alldta));
   gvs = unique(alldta[.,1],1);
   for ii(1,rows(gvs),1);
       dtai = selif(alldta,alldta[.,1].==gvs[ii,1]);
       if rows(dtai)>1;
         dtai[.,indx] = dtai[.,indx]-meanc(dtai[.,indx])';
         newdta = newdta|dtai;
       endif;
   endfor;
   alldta = newdta[2:rows(newdta),.];

   endif;



   startestim = nestim; @do not touch this@
   year = 1;
   do until year > nyr;

    @ --------- Read in the data and define the variables. ----------@

      dta=selif(alldta,(alldta[.,yrindx].==(fyear+year-1)));

      n=rows(dta);
      nsave[1,year] = n;

      gvsave[1:n,year] = dta[.,indcv(idname,vnames)];



    y = dta[.,yindx];
    x = dta[.,xindx];
    if nz > 0;
    z = dta[.,zindx];
    endif;


    icept=ones(n,1);
    if nz > 0;
    za = z;
    iz = icept~za;
    else;
    iz = icept;
    endif;

    mux = (invpd(moment(iz,0)))* (iz'x);
    muy = (invpd(moment(iz,0)))* (iz'y);     @   year+1972;;mux';;muy'; @
    y_d = y - iz*muy;
    x_d = x - iz*mux;

    y2_d = y_d^2;
    yx_d = (x_d .* y_d);
    x2_d = x_d^2;
    y2x_d = (x_d .* y2_d);
    yx2_d = (x2_d .* y_d);
    y3x_d = y2x_d .* y_d;
    y2x2_d = y2_d .* x2_d;
    yx3_d = yx2_d .* x_d;

    x3_d = x_d.*x2_d;
    y3_d = y_d.*y2_d;
    y4x_d = y3_d.*yx_d;
    y3x2_d = y3_d.*x2_d;
    y2x3_d = y2_d.*x3_d;
    yx4_d = yx_d.*x3_d;
    y4_d = y2_d^2;
    x4_d = x2_d^2;
    y5x_d = y_d.*y4x_d;
    y4x2_d = y_d.*y3x2_d;
    y3x3_d = y_d.*y2x3_d;
    y2x4_d = y_d.*yx4_d;
    yx5_d = yx4_d.*x_d;

    y5_d = y2_d.*y3_d;
    x5_d = x2_d.*x3_d;
    y6x_d = y5_d.*yx_d;
    y5x2_d = y5_d.*x2_d;
    y4x3_d = y4_d.*x3_d;
    y3x4_d = y3_d.*x4_d;
    y2x5_d = y2_d.*x5_d;
    yx6_d = yx_d.*x5_d;

/* VERY IMPORTANT:  The following variables MUST be in the same order that
                    they appear in the subroutine deff.
*/

    mom = { y2_d yx_d x2_d y2x_d yx2_d y3x_d y2x2_d yx3_d x3_d y3_d y4x_d
            y3x2_d y2x3_d yx4_d y4_d x4_d y5x_d y4x2_d y3x3_d y2x4_d yx5_d
            y5_d x5_d y6x_d y5x2_d y4x3_d y3x4_d y2x5_d yx6_d };

    Ey2_d = meanc(y2_d);
    Eyx_d = meanc(yx_d);
    Ex2_d = meanc(x2_d);
    Ey2x_d = meanc(y2x_d);
    Eyx2_d = meanc(yx2_d);
    Ey3x_d = meanc(y3x_d);
    Ey2x2_d = meanc(y2x2_d);
    Eyx3_d = meanc(yx3_d);

    Ex3_d     = meanc(x3_d);
    Ey3_d     = meanc(y3_d);
    Ey4x_d    = meanc(y4x_d);
    Ey3x2_d   = meanc(y3x2_d);
    Ey2x3_d   = meanc(y2x3_d);
    Eyx4_d    = meanc(yx4_d);
    Ey4_d     = meanc(y4_d);
    Ex4_d     = meanc(x4_d);
    Ey5x_d    = meanc(y5x_d);
    Ey4x2_d   = meanc(y4x2_d);
    Ey3x3_d   = meanc(y3x3_d);
    Ey2x4_d   = meanc(y2x4_d);
    Eyx5_d    = meanc(yx5_d);
    Ey5_d     = meanc(y5_d);
    Ex5_d     = meanc(x5_d);
    Ey6x_d    = meanc(y6x_d);
    Ey5x2_d   = meanc(y5x2_d);
    Ey4x3_d   = meanc(y4x3_d);
    Ey3x4_d   = meanc(y3x4_d);
    Ey2x5_d   = meanc(y2x5_d);
    Eyx6_d    = meanc(yx6_d);

    Emom = { Ey2_d Eyx_d Ex2_d Ey2x_d Eyx2_d Ey3x_d Ey2x2_d Eyx3_d Ex3_d
             Ey3_d Ey4x_d Ey3x2_d Ey2x3_d Eyx4_d Ey4_d Ex4_d Ey5x_d Ey4x2_d
             Ey3x3_d Ey2x4_d Eyx5_d Ey5_d Ex5_d Ey6x_d Ey5x2_d Ey4x3_d
             Ey3x4_d Ey2x5_d Eyx6_d };

/*  Make the moments for the second stage of the partialling estimator. */


    Ezx  = iz'x/n;

/* Make the moments for the rest of the summary statistics. */

   Ey  = meanc(y);
   Ex  = meanc(x);
   Ey2 = meanc(y^2);
   Ex2 = meanc(x^2);
   if nz > 0;
   Ez  = meanc(z);
   Ez2 = meanc(z^2);
   endif;

/* Make the moments for the correct standard errors for beta. */



    y6_d  = y_d^6;
    x6_d  = x_d^6;

    Ey2_z   =   meanc(y2_d   .* iz);
    Eyx_z   =   meanc(yx_d   .* iz);
    Ex2_z   =   meanc(x2_d   .* iz);
    Ey2x_z  =   meanc(y2x_d  .* iz);
    Eyx2_z  =   meanc(yx2_d  .* iz);
    Ey3x_z  =   meanc(y3x_d  .* iz);
    Ey2x2_z =   meanc(y2x2_d .* iz);
    Eyx3_z  =   meanc(yx3_d  .* iz);
    Ey3_z   =   meanc(y3_d   .* iz);
    Ex3_z   =   meanc(x3_d   .* iz);
    Ey4x_z    = meanc(y4x_d.*iz);
    Ey3x2_z   = meanc(y3x2_d.*iz);
    Ey2x3_z   = meanc(y2x3_d.*iz);
    Eyx4_z    = meanc(yx4_d.*iz);
    Ey4_z     = meanc(y4_d.*iz);
    Ex4_z     = meanc(x4_d.*iz);
    Ey5x_z    = meanc(y5x_d.*iz);
    Ey4x2_z   = meanc(y4x2_d.*iz);
    Ey3x3_z   = meanc(y3x3_d.*iz);
    Ey2x4_z   = meanc(y2x4_d.*iz);
    Eyx5_z    = meanc(yx5_d.*iz);
    Ey5_z     = meanc(y5_d.*iz);
    Ex5_z     = meanc(x5_d.*iz);
    Ey6x_z    = meanc(y6x_d.*iz);
    Ey5x2_z   = meanc(y5x2_d.*iz);
    Ey4x3_z   = meanc(y4x3_d.*iz);
    Ey3x4_z   = meanc(y3x4_d.*iz);
    Ey2x5_z   = meanc(y2x5_d.*iz);
    Eyx6_z    = meanc(yx6_d.*iz);
    Ey6_z     = meanc(y6_d.*iz);
    Ex6_z     = meanc(x6_d.*iz);

    Ezz     =   moment(iz,0)/n;
    inEzz   =   invpd(Ezz);
    zy_d    =   iz.*y_d;
    zx_d    =   iz.*x_d;

    if nz > 0; Eza   = meanc(za)'; endif;

/* The OLS estimator: */

    if nz > 0;
     des = icept~x~za;
    else;
     des = icept~x;
    endif;
    dd = invpd(moment(des,0));
    allb = dd*des'y;
    isave[1,year] = allb[1,1];
    bsave[1,year] = allb[2,1];
    if nz>0;
      for qq(1, nz, 1); zsave[(2*qq-1),year] = allb[2+qq,1]; endfor;
    endif;
    uhat = y - des*allb;
    wh = uhat.*des;
    seb = sqrt(diag(dd*(wh'wh)*dd));
    naiveseb = sqrt(diag(dd*meanc(uhat^2)));
    isave[2,year] = seb[1,1];
    bsave[2,year] = seb[2,1];
    if nz>0;
      for qq(1, nz, 1); zsave[(2*qq),year] = seb[2+qq,1]; endfor;
    endif;
    rsave[1,year] = 1 - moment(uhat,0)/(moment((y - Ey),0));

    olsinflnc = invpd(moment(des,0)/n)*((des.*uhat)');
    iisave[1:n,year] = olsinflnc[1,.]';
    ibsave[1:n,year] = olsinflnc[2,.]';
    if nz>0;
      for qq(1, nz, 1); izsave[((qq-1)*n+1):(qq*n),year] = olsinflnc[2+qq,.]';
      endfor;
    endif;

@----------First make the influence function for sigma_xz.--------------------@
if nz > 0;

xza = x~za;
Exza = meanc(xza);
nreg=nz+1;
sigxz = moment((xza)-Exza',0)/n;


vecsigxz = zeros(nreg+nreg*(nreg-1)/2,1);
vecsigxz[1:nreg,1] = diag(sigxz);

counter =  nreg+nreg*(nreg-1)/2;
for qq(nreg, 2, -1);
  for ee(qq-1, 1, -1);
    vecsigxz[counter,1] = sigxz[ee,qq];
    counter=counter-1;
  endfor;
endfor;


phixz = zeros(n,nreg+nreg*(nreg-1)/2);
phixz[.,1:nreg] = ((xza) - Exza').*((xza) - Exza');

counter =  nreg+nreg*(nreg-1)/2;
for qq(nreg, 2, -1);
  for ee(qq-1, 1, -1);
    phixz[.,counter] = (xza[.,ee] - Exza[ee,1]).*(xza[.,qq] - Exza[qq,1]);
    counter=counter-1;
  endfor;
endfor;

sigy = moment((y - Ey),0)/n;
phiy = (y - Ey)^2 - sigy;

bigphi = (olsinflnc[2:rows(olsinflnc),.])|(phixz')|(phiy');

@--------------Now make the derivative matrix.-------------------------------@

gee = zeros(nreg+(nreg+nreg*(nreg-1)/2)+1,1);


@-------------First, the derivative wrt the OLS coefficients.----------@

for qq(1, nreg, 1);

gee[qq,1] = (2*allb[2:nreg+1,1]'sigxz[.,qq])/sigy;

endfor;

@------------derivatives wrt the first part of sigxz--@

for qq(1, nreg, 1);

gee[nreg+qq,1] = (allb[qq+1,1]^2)/sigy;

endfor;

@------------derivatives wrt the second part of sigxz--@

counter = nreg+(nreg+nreg*(nreg-1)/2);
for qq(nreg, 2, -1);
  for ee(qq-1, 1, -1);
  gee[counter,1] = 2*allb[ee+1,1]*allb[qq+1,1]/sigy;
  counter=counter-1;
  endfor;
endfor;

@------------derivatives wrt sigy--@

counter = nreg+(nreg+nreg*(nreg-1)/2);

gee[counter+1,1] = -rsave[1,year]/sigy;

@-----------Save the influence function.-------------@

else;

@----------First make the influence function for sigma_xz.--------------------@
sigx = Ex2-Ex^2;
phix = (x - Ex)^2 - sigx;

sigy = moment((y - Ey),0)/n;
phiy = (y - Ey)^2 - sigy;

bigphi = (olsinflnc[2:rows(olsinflnc),.])|(phix')|(phiy');

@--------------Now make the derivative matrix.--------------------------------@

gee = zeros(3,1);

@-------------First, the derivative wrt the OLS coefficients.----------@

gee[1,1] = (2*allb[2,1]*sigx)/sigy;

@------------derivatives wrt the first part of sigxz--@

gee[2,1] = (allb[2,1]^2)/sigy;

@------------derivatives wrt sigy--@

gee[3,1] = -rsave[1,year]/sigy;

@-----------Save the influence function.-------------@

endif;

irsave[1:n,year] = -bigphi'gee;
 rsave[2,year] =  sqrt(moment(irsave[.,year],0)/n^2);


    if nz > 0;
      des = icept~y~za;
    else;
      des = icept~y;
    endif;
    dd = invpd(moment(des,0));
    allbr = dd*des'x;
    rbsave[1,year]  = 1/allbr[2,1];
    if nz > 0;
     for qq(1, nz, 1); rzsave[(2*qq-1),year] = -allbr[qq+2,1]/allbr[2,1];
     endfor;
    endif;
    uhat = x - des*allbr;
    wh = uhat.*des;
    vc=(dd*(wh'wh)*dd);
     gee = zeros(nz+2,nz+2);

     gee[1,1] = -1/allbr[2,1];
     gee[1,2] = allbr[1,1]/(allbr[2,1]^2);
     gee[2,2] = -allbr[2,1]^(-2);
     for qq(1, nz, 1); gee[qq+2,qq+2] = -1/allbr[2,1]; endfor;
     for qq(1, nz, 1); gee[qq+2,2] = allbr[qq+2,1]/allbr[2,1]; endfor;

     vcr2t = gee*vc*gee';
     serb = sqrt(diag(vcr2t));
     rbsave[2,year] = serb[2,1];
     if nz > 0;
       for qq(1, nz, 1); rzsave[(2*qq),year] = serb[qq+1,1]; endfor;
     endif;

/* The GMM estimators: */

    estim = startestim;
    do while estim <= nestim;

@--Once exaclly identified and four iterative estimators. -----------------@

if estim==1;
 nc = 5;
 neq = 5;
elseif estim==2;
 nc = 6;
 neq = 8;
elseif estim==3;
 nc = 9;
 neq = 14;
elseif estim==4;
 nc = 12;
 neq = 21;
elseif estim==5;
 nc = 15;
 neq = 29;
endif;
    @---------- Here we input the starting values.----------------@
    c = zeros(nc,1);

    c[1,1] = Ey2x_d/Eyx2_d;
    if cstart /= 0 and estim > 1; c[1,1] = cstart; endif;

    if bb == 0;

         if fe==0 and assets == 0;   @Note: these are the starting values that search.g spits out. This speeds things up substantially.@


             if (year+fyear-1)==1967;  if estim == 2; c[1,1] =     0.27575758   ; elseif estim == 3; c[1,1] =     -0.05151515     ;  endif;
         elseif (year+fyear-1)==1968;  if estim == 2; c[1,1] =    -0.05151515   ; elseif estim == 3; c[1,1] =      0.25757576     ;  endif;
         elseif (year+fyear-1)==1969;  if estim == 2; c[1,1] =    -0.06363636   ; elseif estim == 3; c[1,1] =     -0.04545455     ;  endif;
         elseif (year+fyear-1)==1970;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =     -0.05151515     ;  endif;
         elseif (year+fyear-1)==1971;  if estim == 2; c[1,1] =     0.25757576   ; elseif estim == 3; c[1,1] =      0.24545455     ;  endif;
         elseif (year+fyear-1)==1972;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.17272727     ;  endif;
         elseif (year+fyear-1)==1973;  if estim == 2; c[1,1] =     0.25757576   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1974;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1975;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1976;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.06363636     ;  endif;
         elseif (year+fyear-1)==1977;  if estim == 2; c[1,1] =     0.26363636   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1978;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.28787879     ;  endif;
         elseif (year+fyear-1)==1979;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1980;  if estim == 2; c[1,1] =     0.28181818   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1981;  if estim == 2; c[1,1] =     0.03939394   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1982;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25757576     ;  endif;
         elseif (year+fyear-1)==1983;  if estim == 2; c[1,1] =     0.27575758   ; elseif estim == 3; c[1,1] =      0.26969697     ;  endif;
         elseif (year+fyear-1)==1984;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1985;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1986;  if estim == 2; c[1,1] =     0.03939394   ; elseif estim == 3; c[1,1] =      0.14242424     ;  endif;
         elseif (year+fyear-1)==1987;  if estim == 2; c[1,1] =     0.06969697   ; elseif estim == 3; c[1,1] =      0.25757576     ;  endif;
         elseif (year+fyear-1)==1988;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.28181818     ;  endif;
         elseif (year+fyear-1)==1989;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.19090909     ;  endif;
         elseif (year+fyear-1)==1990;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1991;  if estim == 2; c[1,1] =     0.05151515   ; elseif estim == 3; c[1,1] =      0.03333333     ;  endif;
         elseif (year+fyear-1)==1992;  if estim == 2; c[1,1] =     0.23333333   ; elseif estim == 3; c[1,1] =      0.23939394     ;  endif;
         elseif (year+fyear-1)==1993;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.22727273     ;  endif;
         elseif (year+fyear-1)==1994;  if estim == 2; c[1,1] =     0.07575758   ; elseif estim == 3; c[1,1] =      0.10000000     ;  endif;
         elseif (year+fyear-1)==1995;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1996;  if estim == 2; c[1,1] =     0.19090909   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1997;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.21515152     ;  endif;
         elseif (year+fyear-1)==1998;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==1999;  if estim == 2; c[1,1] =     0.02727273   ; elseif estim == 3; c[1,1] =      0.26969697     ;  endif;
         elseif (year+fyear-1)==2000;  if estim == 2; c[1,1] =     0.32424242   ; elseif estim == 3; c[1,1] =      0.25757576     ;  endif;
         elseif (year+fyear-1)==2001;  if estim == 2; c[1,1] =     0.26363636   ; elseif estim == 3; c[1,1] =      0.24545455     ;  endif;
         elseif (year+fyear-1)==2002;  if estim == 2; c[1,1] =     0.37272727   ; elseif estim == 3; c[1,1] =      0.22121212     ;  endif;
         elseif (year+fyear-1)==2003;  if estim == 2; c[1,1] =     0.23939394   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==2004;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==2005;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.23939394     ;  endif;
         elseif (year+fyear-1)==2006;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.31818182     ;  endif;
         elseif (year+fyear-1)==2007;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         elseif (year+fyear-1)==2008;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =      0.25151515     ;  endif;
         endif;


         elseif fe==1 and assets == 0;


             if (year+fyear-1)==1967;  if estim == 2; c[1,1] =     -0.05757576     ; elseif estim == 3; c[1,1] =   -0.07575758     ;  endif;
         elseif (year+fyear-1)==1968;  if estim == 2; c[1,1] =      0.35454545     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1969;  if estim == 2; c[1,1] =     -0.05151515     ; elseif estim == 3; c[1,1] =    0.27575758     ;  endif;
         elseif (year+fyear-1)==1970;  if estim == 2; c[1,1] =      0.18484848     ; elseif estim == 3; c[1,1] =   -0.05757576     ;  endif;
         elseif (year+fyear-1)==1971;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.26969697     ;  endif;
         elseif (year+fyear-1)==1972;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =   -0.05757576     ;  endif;
         elseif (year+fyear-1)==1973;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.26363636     ;  endif;
         elseif (year+fyear-1)==1974;  if estim == 2; c[1,1] =      0.21515152     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1975;  if estim == 2; c[1,1] =      0.26363636     ; elseif estim == 3; c[1,1] =    0.26363636     ;  endif;
         elseif (year+fyear-1)==1976;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1977;  if estim == 2; c[1,1] =      0.26969697     ; elseif estim == 3; c[1,1] =   -0.05151515     ;  endif;
         elseif (year+fyear-1)==1978;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1979;  if estim == 2; c[1,1] =      0.27575758     ; elseif estim == 3; c[1,1] =   -0.05151515     ;  endif;
         elseif (year+fyear-1)==1980;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1981;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.26363636     ;  endif;
         elseif (year+fyear-1)==1982;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.28181818     ;  endif;
         elseif (year+fyear-1)==1983;  if estim == 2; c[1,1] =      0.31212121     ; elseif estim == 3; c[1,1] =    0.20303030     ;  endif;
         elseif (year+fyear-1)==1984;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.13636364     ;  endif;
         elseif (year+fyear-1)==1985;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1986;  if estim == 2; c[1,1] =      0.03333333     ; elseif estim == 3; c[1,1] =   -0.05151515     ;  endif;
         elseif (year+fyear-1)==1987;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1988;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1989;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1990;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1991;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1992;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.24545455     ;  endif;
         elseif (year+fyear-1)==1993;  if estim == 2; c[1,1] =      0.19696970     ; elseif estim == 3; c[1,1] =    0.22121212     ;  endif;
         elseif (year+fyear-1)==1994;  if estim == 2; c[1,1] =      0.27575758     ; elseif estim == 3; c[1,1] =    0.33636364     ;  endif;
         elseif (year+fyear-1)==1995;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25757576     ;  endif;
         elseif (year+fyear-1)==1996;  if estim == 2; c[1,1] =      0.30606061     ; elseif estim == 3; c[1,1] =    0.20909091     ;  endif;
         elseif (year+fyear-1)==1997;  if estim == 2; c[1,1] =      0.02727273     ; elseif estim == 3; c[1,1] =    0.26363636     ;  endif;
         elseif (year+fyear-1)==1998;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==1999;  if estim == 2; c[1,1] =      0.25757576     ; elseif estim == 3; c[1,1] =    0.25757576     ;  endif;
         elseif (year+fyear-1)==2000;  if estim == 2; c[1,1] =      0.43939394     ; elseif estim == 3; c[1,1] =    0.28181818     ;  endif;
         elseif (year+fyear-1)==2001;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.00303030     ;  endif;
         elseif (year+fyear-1)==2002;  if estim == 2; c[1,1] =      0.00909091     ; elseif estim == 3; c[1,1] =    0.00909091     ;  endif;
         elseif (year+fyear-1)==2003;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.25151515     ;  endif;
         elseif (year+fyear-1)==2004;  if estim == 2; c[1,1] =      0.02727273     ; elseif estim == 3; c[1,1] =    0.25757576     ;  endif;
         elseif (year+fyear-1)==2005;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =    0.21515152     ;  endif;
         elseif (year+fyear-1)==2006;  if estim == 2; c[1,1] =      0.00909091     ; elseif estim == 3; c[1,1] =    0.00303030     ;  endif;
         elseif (year+fyear-1)==2007;  if estim == 2; c[1,1] =      0.03333333     ; elseif estim == 3; c[1,1] =   -0.05151515     ;  endif;
         elseif (year+fyear-1)==2008;  if estim == 2; c[1,1] =      0.25151515     ; elseif estim == 3; c[1,1] =   -0.06969697     ;  endif;
         endif;



         elseif fe==0 and assets == 1;


             if (year+fyear-1)==1967;  if estim == 2; c[1,1] =    -0.05151515    ; elseif estim == 3; c[1,1] =       -0.02121212    ;  endif;
         elseif (year+fyear-1)==1968;  if estim == 2; c[1,1] =     0.37878788    ; elseif estim == 3; c[1,1] =       -0.03333333    ;  endif;
         elseif (year+fyear-1)==1969;  if estim == 2; c[1,1] =     0.39696970    ; elseif estim == 3; c[1,1] =       -0.05151515    ;  endif;
         elseif (year+fyear-1)==1970;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1971;  if estim == 2; c[1,1] =    -0.02727273    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1972;  if estim == 2; c[1,1] =    -0.02121212    ; elseif estim == 3; c[1,1] =        0.25757576    ;  endif;
         elseif (year+fyear-1)==1973;  if estim == 2; c[1,1] =    -0.05151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1974;  if estim == 2; c[1,1] =    -0.05151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1975;  if estim == 2; c[1,1] =    -0.05151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1976;  if estim == 2; c[1,1] =    -0.08181818    ; elseif estim == 3; c[1,1] =       -0.05757576    ;  endif;
         elseif (year+fyear-1)==1977;  if estim == 2; c[1,1] =     0.26969697    ; elseif estim == 3; c[1,1] =       -0.05151515    ;  endif;
         elseif (year+fyear-1)==1978;  if estim == 2; c[1,1] =    -0.06363636    ; elseif estim == 3; c[1,1] =       -0.06363636    ;  endif;
         elseif (year+fyear-1)==1979;  if estim == 2; c[1,1] =    -0.05151515    ; elseif estim == 3; c[1,1] =        0.25757576    ;  endif;
         elseif (year+fyear-1)==1980;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25757576    ;  endif;
         elseif (year+fyear-1)==1981;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1982;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1983;  if estim == 2; c[1,1] =     0.35454545    ; elseif estim == 3; c[1,1] =        0.26969697    ;  endif;
         elseif (year+fyear-1)==1984;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1985;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1986;  if estim == 2; c[1,1] =     0.17878788    ; elseif estim == 3; c[1,1] =        0.04545455    ;  endif;
         elseif (year+fyear-1)==1987;  if estim == 2; c[1,1] =     0.28787879    ; elseif estim == 3; c[1,1] =        0.26363636    ;  endif;
         elseif (year+fyear-1)==1988;  if estim == 2; c[1,1] =     0.24545455    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1989;  if estim == 2; c[1,1] =    -0.08787879    ; elseif estim == 3; c[1,1] =        0.28181818    ;  endif;
         elseif (year+fyear-1)==1990;  if estim == 2; c[1,1] =     0.26363636    ; elseif estim == 3; c[1,1] =        0.28787879    ;  endif;
         elseif (year+fyear-1)==1991;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =       -0.05757576    ;  endif;
         elseif (year+fyear-1)==1992;  if estim == 2; c[1,1] =     0.07575758    ; elseif estim == 3; c[1,1] =        0.24545455    ;  endif;
         elseif (year+fyear-1)==1993;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1994;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1995;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.30000000    ;  endif;
         elseif (year+fyear-1)==1996;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.17878788    ;  endif;
         elseif (year+fyear-1)==1997;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1998;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==1999;  if estim == 2; c[1,1] =     0.13636364    ; elseif estim == 3; c[1,1] =        0.22727273    ;  endif;
         elseif (year+fyear-1)==2000;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==2001;  if estim == 2; c[1,1] =     0.13636364    ; elseif estim == 3; c[1,1] =        0.31212121    ;  endif;
         elseif (year+fyear-1)==2002;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.25151515    ;  endif;
         elseif (year+fyear-1)==2003;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.14848485    ;  endif;
         elseif (year+fyear-1)==2004;  if estim == 2; c[1,1] =     0.23939394    ; elseif estim == 3; c[1,1] =        0.26969697    ;  endif;
         elseif (year+fyear-1)==2005;  if estim == 2; c[1,1] =     0.28181818    ; elseif estim == 3; c[1,1] =        0.25757576    ;  endif;
         elseif (year+fyear-1)==2006;  if estim == 2; c[1,1] =     0.25151515    ; elseif estim == 3; c[1,1] =        0.24545455    ;  endif;
         elseif (year+fyear-1)==2007;  if estim == 2; c[1,1] =     0.02727273    ; elseif estim == 3; c[1,1] =        0.26969697    ;  endif;
         elseif (year+fyear-1)==2008;  if estim == 2; c[1,1] =     0.19090909    ; elseif estim == 3; c[1,1] =        0.26363636    ;  endif;
         endif;



         elseif fe==1 and assets == 1;


             if (year+fyear-1)==1967;  if estim == 2; c[1,1] =    -0.05151515   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==1968;  if estim == 2; c[1,1] =     0.10000000   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==1969;  if estim == 2; c[1,1] =    -0.05757576   ; elseif estim == 3; c[1,1] =   0.33030303     ;  endif;
         elseif (year+fyear-1)==1970;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==1971;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==1972;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==1973;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==1974;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==1975;  if estim == 2; c[1,1] =    -0.05151515   ; elseif estim == 3; c[1,1] =  -0.09393939     ;  endif;
         elseif (year+fyear-1)==1976;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==1977;  if estim == 2; c[1,1] =     0.31212121   ; elseif estim == 3; c[1,1] =  -0.06969697     ;  endif;
         elseif (year+fyear-1)==1978;  if estim == 2; c[1,1] =    -0.05151515   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==1979;  if estim == 2; c[1,1] =    -0.01515152   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==1980;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.40303030     ;  endif;
         elseif (year+fyear-1)==1981;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =  -0.06363636     ;  endif;
         elseif (year+fyear-1)==1982;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.24545455     ;  endif;
         elseif (year+fyear-1)==1983;  if estim == 2; c[1,1] =     0.23939394   ; elseif estim == 3; c[1,1] =   0.27575758     ;  endif;
         elseif (year+fyear-1)==1984;  if estim == 2; c[1,1] =     0.28787879   ; elseif estim == 3; c[1,1] =   0.29393939     ;  endif;
         elseif (year+fyear-1)==1985;  if estim == 2; c[1,1] =     0.39090909   ; elseif estim == 3; c[1,1] =   0.34848485     ;  endif;
         elseif (year+fyear-1)==1986;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==1987;  if estim == 2; c[1,1] =     0.14848485   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==1988;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =  -0.06363636     ;  endif;
         elseif (year+fyear-1)==1989;  if estim == 2; c[1,1] =    -0.05151515   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==1990;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.17878788     ;  endif;
         elseif (year+fyear-1)==1991;  if estim == 2; c[1,1] =     0.33030303   ; elseif estim == 3; c[1,1] =   0.28787879     ;  endif;
         elseif (year+fyear-1)==1992;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==1993;  if estim == 2; c[1,1] =     0.18484848   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==1994;  if estim == 2; c[1,1] =     0.28787879   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==1995;  if estim == 2; c[1,1] =    -0.05757576   ; elseif estim == 3; c[1,1] =   0.25757576     ;  endif;
         elseif (year+fyear-1)==1996;  if estim == 2; c[1,1] =     0.26969697   ; elseif estim == 3; c[1,1] =   0.28181818     ;  endif;
         elseif (year+fyear-1)==1997;  if estim == 2; c[1,1] =     0.23333333   ; elseif estim == 3; c[1,1] =   0.27575758     ;  endif;
         elseif (year+fyear-1)==1998;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==1999;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==2000;  if estim == 2; c[1,1] =     0.26363636   ; elseif estim == 3; c[1,1] =   0.29393939     ;  endif;
         elseif (year+fyear-1)==2001;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.00909091     ;  endif;
         elseif (year+fyear-1)==2002;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.26969697     ;  endif;
         elseif (year+fyear-1)==2003;  if estim == 2; c[1,1] =    -0.05757576   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==2004;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==2005;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =  -0.05757576     ;  endif;
         elseif (year+fyear-1)==2006;  if estim == 2; c[1,1] =     0.25151515   ; elseif estim == 3; c[1,1] =   0.25151515     ;  endif;
         elseif (year+fyear-1)==2007;  if estim == 2; c[1,1] =    -0.05151515   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         elseif (year+fyear-1)==2008;  if estim == 2; c[1,1] =     0.05757576   ; elseif estim == 3; c[1,1] =  -0.05151515     ;  endif;
         endif;

         endif;

         c[2,1] = Eyx_d/c[1,1];
         c[3,1] = Ey2_d - (c[1,1]^2)*c[2,1];
         c[4,1] = Ex2_d - c[2,1];
         c[5,1] = Eyx2_d/c[1,1];
         if estim > 1;
             c[6,1] = (Eyx3_d - 3*c[1,1]*c[2,1]*c[4,1])/c[1,1];
         endif;
         if estim > 2;
             c[7,1] = Ex3_d - c[5,1];
             c[8,1] = Ey3_d - (c[1,1]^3)*c[5,1];
             c[9,1] = (Eyx4_d/c[1,1]) - 6*c[5,1]*c[4,1] - 4*c[2,1]*c[7,1];
         endif;
         if estim > 3;
             c[10,1] = Ey4_d - (c[1,1]^4)*c[6,1] - 6*(c[1,1]^2)*c[2,1]*c[3,1];
             c[11,1] = Ex4_d - c[6,1] - 6*c[2,1]*c[4,1];
             c[12,1] = (Eyx5_d/c[1,1]) - 10*c[6,1]*c[4,1] - 10*c[5,1]*c[7,1]
                       - 5*c[2,1]*c[11,1];
           endif;
           if estim > 4;
               c[13,1] = Ey5_d - (c[1,1]^5)*c[9,1] - 10*(c[1,1]^3)*c[5,1]*c[3,1]
                         - 10*(c[1,1]^2)*c[2,1]*c[8,1];
               c[14,1] = Ex5_d - c[9,1] - 10*c[5,1]*c[4,1] - 10*c[2,1]*c[7,1];
               c[15,1] = (Eyx6_d/c[1,1]) - 15*c[9,1]*c[4,1] - 20*c[6,1]*c[7,1]
                         - 15*c[5,1]*c[11,1] - 6*c[2,1]*c[14,1];
           endif;

          if estim > 1; { c, obj } = gmm(c); endif;  @Call the GMM routine.@

      else;

                         if estim == 1;

                             c[1,1] = Ey2x_d/Eyx2_d;
                             c[2,1] = Eyx_d/c[1,1];
                             c[3,1] = Ey2_d - (c[1,1]^2)*c[2,1];
                             c[4,1] = Ex2_d - c[2,1];
                             c[5,1] = Eyx2_d/c[1,1];


                         else;

                             ckeep = zeros(2,nc);
                             objkeep = zeros(2,1);

                             for mmm(1,2,1);

                                     if mmm==2;
                                        c[1,1] = Ey2x_d/Eyx2_d;
                                     else;
                                        c[1,1] = allb[2,1];
                                     endif;

                                     c[2,1] = Eyx_d/c[1,1];
                                     c[3,1] = Ey2_d - (c[1,1]^2)*c[2,1];
                                     c[4,1] = Ex2_d - c[2,1];
                                     c[5,1] = Eyx2_d/c[1,1];
                                     c[6,1] = (Eyx3_d - 3*c[1,1]*c[2,1]*c[4,1])/c[1,1];
                                     if estim > 2;
                                         c[7,1] = Ex3_d - c[5,1];
                                         c[8,1] = Ey3_d - (c[1,1]^3)*c[5,1];
                                         c[9,1] = (Eyx4_d/c[1,1]) - 6*c[5,1]*c[4,1] - 4*c[2,1]*c[7,1];
                                     endif;
                                     if estim > 3;
                                         c[10,1] = Ey4_d - (c[1,1]^4)*c[6,1] - 6*(c[1,1]^2)*c[2,1]*c[3,1];
                                         c[11,1] = Ex4_d - c[6,1] - 6*c[2,1]*c[4,1];
                                         c[12,1] = (Eyx5_d/c[1,1]) - 10*c[6,1]*c[4,1] - 10*c[5,1]*c[7,1]
                                                   - 5*c[2,1]*c[11,1];
                                     endif;
                                     if estim > 4;
                                         c[13,1] = Ey5_d - (c[1,1]^5)*c[9,1] - 10*(c[1,1]^3)*c[5,1]*c[3,1]
                                                   - 10*(c[1,1]^2)*c[2,1]*c[8,1];
                                         c[14,1] = Ex5_d - c[9,1] - 10*c[5,1]*c[4,1] - 10*c[2,1]*c[7,1];
                                         c[15,1] = (Eyx6_d/c[1,1]) - 15*c[9,1]*c[4,1] - 20*c[6,1]*c[7,1]
                                                   - 15*c[5,1]*c[11,1] - 6*c[2,1]*c[14,1];
                                     endif;

                                     { c, obj } = gmm(c);  @Call the GMM routine.@
                                     if obj == -999999 and bb > 0;  goto bttdb; endif;

                                     ckeep[mmm,.] = c';
                                     objkeep[mmm,1] = obj;
                               endfor;

                               cidx = selif((1|2),objkeep .==minc(objkeep));
                               if rows(cidx)>1; cidx = cidx[1,1];endif;
                               c = ckeep[cidx,.]';

                         endif;

     endif; @end of the bb = 0 conditional.@

    @ Compute t-ratios, etc.                                            @

         ff = optw(c,mom,Emom,bleh);
         f = deff(c,Emom,0,effsave);
         if bb == 0;                                                  @ The program jumps to the subroutine that defines the f vector @
           effsave[1:neq,(estim-1)*nestim+year] =  f;
         endif;

      w = invit(moment(ff,0)./n);
         if w== -999999 and bb > 0; goto bttdb; endif;

        if estim == 1;
        g=grad1(c);
        elseif estim == 2;
        g=grad2(c);
        elseif estim == 3;
        g=grad3(c);
        elseif estim == 4;
        g=grad4(c);
        elseif estim == 5;
        g=grad5(c);
        endif;

      gwg = g'w*g;
      vc=invit(gwg)/n;
        stderr=sqrt(diag(vc));
        inflnc=-n*vc*g'w*ff';

        csave[1,nestim*(year-1)+estim] = c[1,1];
        csave[2,nestim*(year-1)+estim] = stderr[1,1];
        icsave[1:n,nestim*(year-1)+estim] = inflnc[1,.]';

    if estim > 1;
    chisq = obj*n;
    chisave[1,nestim*(year-1)+estim] = chisq;
    degf=neq - nc;
    chisave[2,nestim*(year-1)+estim] = cdfchic(chisq,degf);
    endif;

    /* Make Delta */

    c11 = c[1,1];
    if nz > 0;

      for qq(1, nz, 1);
        dsave[(qq*2-1),nestim*(year-1)+estim] = muy[qq+1,1] - c11*mux[qq+1,1];
      endfor;

      gee = zeros(((nz+1)*2+1),nz);

      gee[2:(nz+1),1:nz] = eye(nz);
      gee[(nz+3):((nz+1)*2),1:nz] = -c11*eye(nz);

      for qq(1, nz, 1); gee[((nz+1)*2+1),qq] = -mux[qq+1,1]; endfor;

      @---Correct standard error and influence function---@

      bigphi = ( (inEzz*zy_d')|(inEzz*zx_d') )|(-inflnc[1,.]);
      avar = moment(bigphi',0)./n^2;
      dstd = sqrt(diag(gee'avar*gee));
      omega = gee'avar*gee;
      phidel = -bigphi'gee;

      for qq(1, nz, 1); dsave[qq*2,nestim*(year-1)+estim] = dstd[qq,1]; endfor;

      for qq(1, nz, 1);
        idsave[((qq-1)*n+1):(qq*n),nestim*(year-1)+estim] = phidel[.,qq];
      endfor;

    endif;

    /*  Make the stupid intercept */

    asave[1,nestim*(year-1)+estim] = muy[1,1] - c11*mux[1,1];

    gee = zeros(((nz+1)*2+1),1);

    gee[1,1] = 1;
    gee[(nz+2),1] = -c11;
    gee[((nz+1)*2+1),1] = -mux[1,1];

    @---Correct standard error and influence function---@

    bigphi = ( (inEzz*zy_d')|(inEzz*zx_d') )|(-inflnc[1,.]);
    avar = moment(bigphi',0)./n^2;
    aaastd = sqrt(diag(gee'avar*gee));
    phidela = -bigphi'gee;

    asave[2,nestim*(year-1)+estim] = aaastd[1,1];
    iasave[1:n,nestim*(year-1)+estim] = phidela[.,1];


    @--------------Rho2 and tau2--------------------------------------------------@


    if nz > 0;

    @----------First make the influence function for sigma_z.---------------------@

    sigz = moment(za-meanc(za)',0)/n;


    vecsigz = zeros(nz+nz*(nz-1)/2,1);
    vecsigz[1:nz,1] = diag(sigz);

    counter =  nz+nz*(nz-1)/2;
    for qq(nz, 2, -1);
      for ee(qq-1, 1, -1);
        vecsigz[counter,1] = sigz[ee,qq];
        counter=counter-1;
      endfor;
    endfor;


    phiz = zeros(n,nz+nz*(nz-1)/2);
    phiz[.,1:nz] = (za - Eza).*(za - Eza);

    counter =  nz+nz*(nz-1)/2;
    for qq(nz, 2, -1);
      for ee(qq-1, 1, -1);
        vecsigz[counter,1] = sigz[ee,qq];
        phiz[.,counter] = (za[.,ee] - Eza[1,ee]).*(za[.,qq] - Eza[1,qq]);
        counter=counter-1;
      endfor;
    endfor;


    phiz = phiz - vecsigz';

    phimuy = inEzz*zy_d';
    phimux = inEzz*zx_d';

    numer = (muy[2:nz+1,1]'sigz*muy[2:nz+1,1] + c[1,1]^2*c[2,1]);
    rho = numer/(numer + c[3,1]);

    numer = (mux[2:nz+1,1]'sigz*mux[2:nz+1,1] + c[2,1]);
    tau = numer/(numer + c[4,1]);

    @--Make the influence functions for the standard errors for rho2 and tau2.@

    bigphi = phimux[2:nz+1,.]|phimuy[2:nz+1,.]|(phiz')|(-inflnc[1:4,.]);
    avar = moment(bigphi',0)./n^2;

    @-------------- First column is for rho2 and the second is for tau2. @

    gee = zeros(2*nz+(nz+nz*(nz-1)/2)+4,2);  @fix this@


    @-----------First, do rho2---------------------------@


    numer = muy[2:nz+1,1]'sigz*muy[2:nz+1,1]+c[1,1]^2*c[2,1];
    denom =  numer+c[3,1]; ;

    @------------derivatives wrt muy---------------------@

    for qq(1, nz, 1);

    gee[nz+qq,1] =
    (2*muy[2:nz+1,1]'sigz[.,qq])/denom
    - numer*(2*muy[2:nz+1,1]'sigz[.,qq])/(denom^2);

    endfor;
    @------------derivatives wrt the first part of sigz--@

    for qq(1, nz, 1);

    gee[2*nz+qq,1] = (muy[qq+1,1]^2)/denom - numer/(denom^2)*muy[qq+1,1]^2;

    endfor;

    @------------derivatives wrt the second part of sigz--@

    counter = 2*nz+(nz+nz*(nz-1)/2);
    for qq(nz, 2, -1);
      for ee(qq-1, 1, -1);
      gee[counter,1] =
      2*muy[ee+1,1]*muy[qq+1,1]/denom - 2*numer/(denom^2)*muy[ee+1,1]*muy[qq+1,1];
      counter=counter-1;
      endfor;
    endfor;

    @------------derivatives wrt  c--@

    counter = 2*nz+(nz+nz*(nz-1)/2);

    gee[counter+1,1] = 2*c[1,1]*c[2,1]/denom - 2*numer/(denom^2)*c[1,1]*c[2,1];

    gee[counter+2,1] = (c[1,1]^2)/denom - numer/(denom^2)*c[1,1]^2;

    gee[counter+3,1] = -numer/(denom^2);


    @-----------Now for Tau2------------------------@


    numer = mux[2:nz+1,1]'sigz*mux[2:nz+1,1]+c[2,1];
    denom = numer+c[4,1]; ;


    @------------derivatives wrt mux---------------------@

    for qq(1, nz, 1);

    gee[qq,2] =
    (2*mux[2:nz+1]'sigz[.,qq])/denom - numer*(2*mux[2:nz+1]'sigz[.,qq])/(denom^2);

    endfor;

    @------------derivatives wrt the first part of sigz--@

    for qq(1, nz, 1);

    gee[2*nz+qq,2] = (mux[qq+1,1]^2)/denom - numer/(denom^2)*mux[qq+1,1]^2;

    endfor;

    @------------derivatives wrt the second part of sigz--@

    counter = 2*nz+(nz+nz*(nz-1)/2);
    for qq(nz, 2, -1);
      for ee(qq-1, 1, -1);
      gee[counter,2] =
      2*mux[ee+1,1]*mux[qq+1,1]/denom - 2*numer/(denom^2)*mux[ee+1,1]*mux[qq+1,1];
      counter=counter-1;
      endfor;
    endfor;

    gee[2*nz+(nz+nz*(nz-1)/2)+2,2] = 1/denom-numer/(denom^2);
    gee[2*nz+(nz+nz*(nz-1)/2)+4,2] = -numer/(denom^2);



    else;

    numer = (c[1,1]^2*c[2,1]);
    rho = numer/(numer + c[3,1]);

    numer = (c[2,1]);
    tau = numer/(numer + c[4,1]);

    @-----Make the influence functions for the standard errors for rho2 and tau2.@

    bigphi = (-inflnc[1:4,.]);
    avar = moment(bigphi',0)./n^2;

    gee = zeros(4,2); @ First column is for rho2 and the second is for tau2. @


    @-----------First, do rho2---------------------------@


    numer = c[1,1]^2*c[2,1];
    denom =  numer+c[3,1]; ;

    @------------derivatives wrt  c--@


    gee[1,1] = 2*c[1,1]*c[2,1]/denom - 2*numer/(denom^2)*c[1,1]*c[2,1];

    gee[2,1] = (c[1,1]^2)/denom - numer/(denom^2)*c[1,1]^2;

    gee[3,1] = -numer/(denom^2);


    @-----------Now for Tau2------------------------@


    numer = c[2,1];
    denom = numer+c[4,1]; ;


    gee[2,2] = 1/denom-numer/(denom^2);
    gee[4,2] = -numer/(denom^2);


    endif;

    vcrhotau = gee'avar*gee;

    phirho = -bigphi'gee[.,1];
    phitau = -bigphi'gee[.,2];

    ipsave[1:n,nestim*(year-1)+estim] = phirho;
    itsave[1:n,nestim*(year-1)+estim] = phitau;


    /*  Now save both things. */

        psave[1,nestim*(year-1)+estim]    = rho;
        tsave[1,nestim*(year-1)+estim]    = tau;

        psave[2,nestim*(year-1)+estim]    = sqrt(vcrhotau[1,1]);
        tsave[2,nestim*(year-1)+estim]    = sqrt(vcrhotau[2,2]);


        estim = estim + 1;
        endo;                       @ The end of the estimator loop. @

        year = year + 1;

    endo;     @This ends the trial loop.@

    /* This part does the classical minimum distance estimation. */
         sharefrac = zeros(nyr,nyr);
         for ff(1,nyr,1);
            for gg(1,nyr,1);
               isthere = selif(gvsave[.,ff],gvsave[.,ff]./=-999) .== selif(gvsave[.,gg],gvsave[.,gg]./=-999)';
               sharefrac[ff,gg] = sumc(sumc(isthere))/nsave[1,ff];
            endfor;
         endfor;
    param = 1; @ 1 is beta, 2 is rho2, 3 is tau2,
                 4 to nz+3 are the alphas, and the end is the intercept. @
    do while param <= nz+4;

      idxnot = seqa(1,nestim,nyr) - 1;

    estim = startestim;
    do while estim <= nestim;

      idx = idxnot + estim;
      allinfl = zeros(nobs,nyr);

      if panel == 1;

         if param == 1;
           saveit = csave;
           w = invit(moment(icsave[.,idx],0)/n^2);
         elseif param == 2;
           saveit = psave;
           w = invit(moment(ipsave[.,idx],0)/n^2);
         elseif param == 3;
           saveit = tsave;
           w = invit(moment(itsave[.,idx],0)/n^2);
         elseif param == nz+4;
           saveit = asave;
           w = invit(moment(iasave[.,idx],0)/n^2);
         endif;

         for qq(1, nz, 1);
         if param == qq+3;
           saveit = dsave[qq*2-1,.];
           w = invit(moment(idsave[((qq-1)*n+1):(qq*n),idx],0)/n^2);
         endif;
         endfor;


      else;

         for ffff(1,nyr,1);

           pickout = packr(indnv(gvsave[.,ffff],gvs));

           if param == 1;
             allinfl[pickout,ffff] = icsave[1:nsave[ffff],idx[ffff]];
             saveit = csave;
           elseif param == 2;
             saveit = psave;
             allinfl[pickout,ffff] = ipsave[1:nsave[ffff],idx[ffff]];
           elseif param == 3;
             saveit = tsave;
             allinfl[pickout,ffff] = itsave[1:nsave[ffff],idx[ffff]];
           elseif param == nz+4;
             saveit = asave;
             allinfl[pickout,ffff] = iasave[1:nsave[ffff],idx[ffff]];
           endif;


           for qq(1, nz, 1);
           if param == qq+3;
             saveit = dsave[qq*2-1,.]|dsave[qq*2,.];
                                     idd = idsave[((qq-1)*nco+1):(qq*nco),.];
             allinfl[pickout,ffff] = idd[1:nsave[ffff],idx[ffff]];
           endif;
           endfor;

         endfor;
         ninfl = substute(allinfl,allinfl./=0,1);
         divideby = moment(ninfl,0);

         w = invit(moment(allinfl,0)./(divideby^2).*(sharefrac^1).*((sharefrac')^1));   @ if param==4; diag(w)~((saveit[1,idx]./saveit[2,idx])');  endif;       @
         if w== -999999 and bb > 0; goto bttdb; endif;

      endif;


        g=ones(nyr,1);
        theta = invpd(g'w*g)*g'w*saveit[1,idx]';   @ if param == 4; diag(w)~saveit[1,idx]'; end; endif;   @

        f=saveit[1,idx]'-theta;
        stderror = sqrt(invpd(g'w*g));
        cmdtest = f'w*f;
        if param == 1;
                   if bb == 0;
                            ccsave[1,estim] = theta;
                            ccsave[2,estim] = cmdtest;
                            ccsave[3,estim] = cdfchic(cmdtest,nyr-1);
                            ccsave[4,estim] = stderror;
                   else;
                            cpivot[estim,bb] = (theta-ccsave[1,estim])/stderror;
                   endif;
        elseif param == 2;
                   if bb == 0;
                            cpsave[1,estim] = theta;
                            cpsave[2,estim] = cmdtest;
                            cpsave[3,estim] = cdfchic(cmdtest,nyr-1);
                            cpsave[4,estim] = stderror;
                   else;
                            ppivot[estim,bb] = (theta-cpsave[1,estim])/stderror;
                   endif;

        elseif param == 3;
                   if bb == 0;
                            ctsave[1,estim] = theta;
                            ctsave[2,estim] = cmdtest;
                            ctsave[3,estim] = cdfchic(cmdtest,nyr-1);
                            ctsave[4,estim] = stderror;
                   else;
                            tpivot[estim,bb] = (theta-ctsave[1,estim])/stderror;
                   endif;

        elseif param ==nz+4;
                   if bb == 0;
                            casave[1,estim] = theta;
                            casave[2,estim] = cmdtest;
                            casave[3,estim] = cdfchic(cmdtest,nyr-1);
                            casave[4,estim] = stderror;
                   else;
                            apivot[estim,bb] = (theta-casave[1,estim])/stderror;
                   endif;

        endif;

        for qq(1, nz, 1);
        if param == qq+3;

                   if bb == 0;
                            cdsave[(qq-1)*4+1,estim] = theta;
                            cdsave[(qq-1)*4+2,estim] = cmdtest;
                            cdsave[(qq-1)*4+3,estim] = cdfchic(cmdtest,nyr-1);
                            cdsave[(qq-1)*4+4,estim] = stderror;
                   else;
                            dpivot[(estim-1)*nz+qq,bb] = (theta-cdsave[(qq-1)*4+1,estim])/stderror;
                   endif;

        endif;
        endfor;


    estim = estim + 1;
    endo;

    param = param + 1;
    endo;



    /* This part does the classical minimum distance estimation. */
    /* FOR OLS */

    param = 1; @ 1 is beta, 2 through nz+1 are the alphas,
                 nz+2 is R2, and the end is the intercept @
    do while param <= nz+3;

       allinfl = zeros(nobs,nyr);
       if panel ==   1;
         if param == 1;
           saveit = bsave[1,.];
           w = invit(moment(ibsave,0)/n^2);
         elseif param == nz+2;
           saveit = rsave[1,.];
           w = invit(moment(irsave,0)/n^2);
         elseif param == nz+3;
           saveit = isave[1,.];
           w = invit(moment(iisave,0)/n^2);
         endif;


         for qq(1, nz, 1);
         if param == qq+1;
           saveit = zsave[(2*qq-1),.];
           w = invit(moment(izsave[((qq-1)*n+1):(qq*n),.],0)/n^2);
         endif;
         endfor;

       else;
         for ffff(1,nyr,1);

              pickout = packr(indnv(gvsave[.,ffff],gvs));

              if param == 1;
                saveit = bsave[1,.];
                allinfl[pickout,ffff] = ibsave[1:nsave[ffff],ffff];
              elseif param == nz+2;
                saveit = rsave[1,.];
                allinfl[pickout,ffff] = irsave[1:nsave[ffff],ffff];
              elseif param == nz+3;
                saveit = isave[1,.];
                allinfl[pickout,ffff] = iisave[1:nsave[ffff],ffff];
              endif;


              for qq(1, nz, 1);
              if param == qq+1;
                saveit = zsave[(2*qq-1),.];
                                        idd = izsave[((qq-1)*nco+1):(qq*nco),.];
                allinfl[pickout,ffff] = idd[1:nsave[ffff],ffff];
              endif;
              endfor;
         endfor;
         w = invit(moment(allinfl,0)./(divideby^2).*(sharefrac^1).*((sharefrac')^1));
         if w== -999999 and bb > 0; goto bttdb; endif;
       endif;




       g=ones(nyr,1);
         theta = invpd(g'w*g)*g'w*saveit';

    @ Compute t-ratios, etc.                                            @

        f=saveit'-theta;
        stderror = sqrt(invpd(g'w*g));
        cmdtest = f'w*f;

        if param == 1;
          if bb == 0;
             bbsave[1,1] = theta;
             bbsave[2,1] = cmdtest;
             bbsave[3,1] = cdfchic(cmdtest,nyr-1);
             bbsave[4,1] = stderror;
          else;
                   bpivot[1,bb] = (theta-bbsave[1,1])/stderror;
          endif;
        elseif param == nz+2;
          if bb == 0;
             brsave[1,1] = theta;
             brsave[2,1] = cmdtest;
             brsave[3,1] = cdfchic(cmdtest,nyr-1);
             brsave[4,1] = stderror;
          else;
                   rpivot[1,bb] = (theta-brsave[1,1])/stderror;
          endif;
        elseif param == nz+3;
          if bb == 0;
             bisave[1,1] = theta;
             bisave[2,1] = cmdtest;
             bisave[3,1] = cdfchic(cmdtest,nyr-1);
             bisave[4,1] = stderror;
           else;
                    ipivot[1,bb] = (theta-bisave[1,1])/stderror;
           endif;
        endif;

        for qq(1, nz, 1);
        if param == qq+1;
          if bb == 0;
             bzsave[(qq-1)*4+1,1] = theta;
             bzsave[(qq-1)*4+2,1] = cmdtest;
             bzsave[(qq-1)*4+3,1] = cdfchic(cmdtest,nyr-1);
             bzsave[(qq-1)*4+4,1] = stderror;
          else;
                   zpivot[qq,bb] = (theta-bzsave[(qq-1)*4+1,1])/stderror;
          endif;
        endif;
        endfor;

    param = param + 1;
    endo;


    if bb == 0;

    format /rd 8,3;
    if latex == 1;


    /*  This business makes a LaTeX table.   */


    print;
    "Intercept:";
    print;

    for jj(1, nyr, 1);
    "& " isave[1,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " asave[1,qq];;
    endfor; " \\" "\\";
    "&(" isave[2,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); ")&(" asave[2,qq];;
    endfor; ") \\" "\\";
    endfor;
    print;

    print;
    "Coefficient on the mismeasured regressor:";
    print;

    for jj(1, nyr, 1);
    jj;;"& " bsave[1,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " csave[1,qq];; endfor; " \\" "\\";
    jj;;"&(" bsave[2,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); ")&(" csave[2,qq];; endfor; ") \\" "\\";
    endfor;
    print;

    print;
    "Coefficients on the perfectly measured regressors:";
    print;

    for ee(1, nz, 1);

    print;
    for jj(1, nyr, 1);
    "& " zsave[(ee-1)*2+1,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " dsave[(ee-1)*2+1,qq];; endfor; " \\" "\\";
    "&(" zsave[ee*2,jj];;       for qq(((jj-1)*nestim+1), (nestim*jj), 1); ")&(" dsave[ee*2,qq];; endfor; ") \\" "\\";
    endfor;
    print;

    endfor;

    print;
    "R-squared: ";
    print;

    print;
    for jj(1, nyr, 1);
    "& " rsave[1,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " psave[1,qq];; endfor; " \\" "\\";
    "&(" rsave[2,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); ")&(" psave[2,qq];; endfor; ") \\" "\\";
    endfor;
    print; print;

    print;
    "tau-squared: ";
    print;

    print;
    for jj(1, nyr, 1);
         for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " tsave[1,qq];; endfor; " \\" "\\";
    " ";;for qq(((jj-1)*nestim+1), (nestim*jj), 1); "&(" tsave[2,qq];; ")";; endfor; " \\" "\\";
    endfor;
    print; print;


    print;
    "J-statistics: ";
    print;

    for jj(1, nyr, 1);
    for qq(((jj-1)*nestim+2), (nestim*jj), 1); " & " chisave[1,qq];; endfor; " \\" "\\";
    endfor;
    print;print;

    print;
    "P-values for the J-statistics: ";
    print;

    for jj(1, nyr, 1);
    for qq(((jj-1)*nestim+2), (nestim*jj), 1); " & " chisave[2,qq];; endfor; " \\" "\\";
    endfor;
    print;print;


    "----------------------------------------------------------------------------";
    "----------------------------------------------------------------------------";
    "----------------------------------------------------------------------------";
    "----------------------------------------------------------------------------";

    print;
    "OLS:";
    print;

    for yy(1,nyr,1);

    "& " isave[1,yy];; " & " bsave[1,yy];; for tt(1, nz, 1); " & " zsave[(tt*2-1),yy];; endfor; " & " Rsave[1,yy];; " \\" "\\";
    "&(" isave[2,yy];; ")&(" bsave[2,yy];; for tt(1, nz, 1); ")&(" zsave[(tt*2),yy];;   endfor; ")&(" Rsave[2,yy];; ") \\" "\\ [5pt]";

    endfor;
    print;print;

    for qq(1, nestim, 1);

    print;
    "GMM" qq+2;
    print;
    for yy(1,nyr,1);

    yy+fyear-1;;    "& " asave[1,nestim*(yy-1)+qq];; " & " csave[1,nestim*(yy-1)+qq];; for tt(1, nz, 1); " & " dsave[(tt*2-1),nestim*(yy-1)+qq];; endfor; " & " psave[1,nestim*(yy-1)+qq];; " & " tsave[1,nestim*(yy-1)+qq];; " & " chisave[1,nestim*(yy-1)+qq];;  " & "  nsave[1,yy];; " \\" "\\";
    @yy;; @"         &(" asave[2,nestim*(yy-1)+qq];; ")&(" csave[2,nestim*(yy-1)+qq];; for tt(1, nz, 1); ")&(" dsave[(tt*2),nestim*(yy-1)+qq];;endfor; ")&(" psave[2,nestim*(yy-1)+qq];; ")&(" tsave[2,nestim*(yy-1)+qq];; ")&(" chisave[2,nestim*(yy-1)+qq];;   ")&           \\" "\\ [5pt]";
    endfor;
    print;print;

    endfor;


    "---------------------------------------------------------------------------";
    "---------------------------------------------------------------------------";
    "---------------------------------------------------------------------------";
    "---------------------------------------------------------------------------";


       print;
       "CMD Intercept:";
       print;


       "& " bisave[1,1];; for qq(1, nestim, 1); " & " casave[1,qq];;
       endfor; " \\" "\\";
       "&(" bisave[4,1];; for qq(1, nestim, 1); ")&(" casave[4,qq];;
       endfor; ") \\" "\\";
       "& " bisave[2,1];; for qq(1, nestim, 1); " & " casave[2,qq];;
       endfor; " \\" "\\";
       "& " bisave[3,1];; for qq(1, nestim, 1); " & " casave[3,qq];;
       endfor; " \\" "\\";
       print; print;

       print;
       "CMD Coefficient on the mismeasured regressor:";
       print;

       "& " bbsave[1,1];; for qq(1, nestim, 1); " & " ccsave[1,qq];;
       endfor; " \\" "\\";
       "&(" bbsave[4,1];; for qq(1, nestim, 1); ")&(" ccsave[4,qq];;
       endfor; ") \\" "\\";
       "& " bbsave[2,1];; for qq(1, nestim, 1); " & " ccsave[2,qq];;
       endfor; " \\" "\\";
       "& " bbsave[3,1];; for qq(1, nestim, 1); " & " ccsave[3,qq];;
       endfor; " \\" "\\";
       print; print;

       print;
       "CMD Coefficients on the perfectly measured regressor:";
       print;

       for ee(1, nz, 1);
       "& " bzsave[(ee-1)*4+1,1];; for qq(1, nestim, 1); " & " cdsave[(ee-1)*4+1,qq];;
       endfor; " \\" "\\";
       "&(" bzsave[(ee-1)*4+4,1];; for qq(1, nestim, 1); ")&(" cdsave[(ee-1)*4+4,qq];;
       endfor; ") \\" "\\";
       "& " bzsave[(ee-1)*4+2,1];; for qq(1, nestim, 1); " & " cdsave[(ee-1)*4+2,qq];;
       endfor; " \\" "\\";
       "& " bzsave[(ee-1)*4+3,1];; for qq(1, nestim, 1); " & " cdsave[(ee-1)*4+3,qq];;
       endfor; " \\" "\\";
       print; print;
       endfor;

       print;
       "CMD R-squared:";
       print;

       "& " brsave[1,1];; for qq(1, nestim, 1); " & " cpsave[1,qq];;
       endfor; " \\" "\\";
       "&(" brsave[4,1];; for qq(1, nestim, 1); ")&(" cpsave[4,qq];;
       endfor; ") \\" "\\";
       "& " brsave[2,1];; for qq(1, nestim, 1); " & " cpsave[2,qq];;
       endfor; " \\" "\\";
       "& " brsave[3,1];; for qq(1, nestim, 1); " & " cpsave[3,qq];;
       endfor; " \\" "\\";
       print; print;

       print;
       "CMD tau-squared:";
       print;

       for qq(1, nestim, 1); " & " ctsave[1,qq];; endfor; " \\" "\\";
       " ";;for qq(1, nestim, 1); "&(" ctsave[4,qq];; ")";; endfor; " \\" "\\";
       for qq(1, nestim, 1); " & " ctsave[2,qq];; endfor; " \\" "\\";
       for qq(1, nestim, 1); " & " ctsave[3,qq];; endfor; " \\" "\\";

       pick = seqa(1,nestim,nyr);

       print;
       "J-Test Rejections:";
       print;

       for jj(1, nestim, 1);   " & " meanc(chisave[2,pick+jj-1]'.<0.05);;     endfor; " \\" "\\";


       print;
       "Fraction intercept significantly greater than zero:";
       print;

       "& ";; meanc((cdftc(isave[1,.]./isave[2,.],nyr).<0.025)');; for jj(1, nestim, 1);   " & " meanc((cdftc(asave[1,pick+jj-1]./asave[2,pick+jj-1],nyr).<0.025)');;     endfor; " \\" "\\";

       print;
       "Fraction mismeasured regressor coefficients significantly greater than zero:";
       print;

       "& ";; meanc((cdftc(bsave[1,.]./bsave[2,.],nyr).<0.025)');; for jj(1, nestim, 1);   " & " meanc((cdftc(csave[1,pick+jj-1]./csave[2,pick+jj-1],nyr).<0.025)');;     endfor; " \\" "\\";


       print;
       "Fraction perfectly measured regressor coefficients significantly greater than zero:";
       print;

       for ee(1, nz, 1);
       "& ";; meanc((cdftc(zsave[(ee-1)*2+1,.]./zsave[ee*2,.],nyr).<0.025)');; for jj(1, nestim, 1);   " & " meanc((cdftc(dsave[(ee-1)*2+1,pick+jj-1]./dsave[ee*2,pick+jj-1],nyr).<0.025)');;    endfor; " \\" "\\";
       endfor;

       print;
       "Fraction R2s significantly greater than zero:";
       print;

       "& ";; meanc((cdftc(rsave[1,.]./rsave[2,.],nyr).<0.025)');; for jj(1, nestim, 1);   " & " meanc((cdftc(psave[1,pick+jj-1]./psave[2,pick+jj-1],nyr).<0.025)');;     endfor; " \\" "\\";

       print;
       "Fraction tau-squareds significantly greater than zero:";
       print;

       for jj(1, nestim, 1);   " & " meanc((cdftc(tsave[1,pick+jj-1]./tsave[2,pick+jj-1],nyr).<0.025)');;     endfor; " \\" "\\";


    else;


    print;
    "Intercept:";
    print;

    for jj(1, nyr, 1);
    "  " isave[1,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); "   " asave[1,qq];;
    endfor; "   " "  ";
    " (" isave[2,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); ") (" asave[2,qq];;
    endfor; ")   " "  ";
    endfor;
    print;

    print;
    "Coefficient on the mismeasured regressor:";
    print;

    for jj(1, nyr, 1);
    "  " bsave[1,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); "   " csave[1,qq];;
    endfor; "   " "  ";
    " (" bsave[2,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); ") (" csave[2,qq];;
    endfor; ")   " "  ";
    endfor;
    print;

    print;
    "Coefficients on the perfectly measured regressors:";
    print;

    for ee(1, nz, 1);

    print;
    for jj(1, nyr, 1);
    "  " zsave[(ee-1)*2+1,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); "   " dsave[(ee-1)*2+1,qq];;
    endfor; "   " "  ";
    " (" zsave[ee*2,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); ") (" dsave[ee*2,qq];;
    endfor; ")   " "  ";
    endfor;
    print;

    endfor;

    print;
    "R-squared: ";
    print;

    print;
    for jj(1, nyr, 1);
    "  " rsave[1,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); "   " psave[1,qq];;
    endfor; "   " "  ";
    " (" rsave[2,jj];;
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); ") (" psave[2,qq];;
    endfor; ")   " "  ";
    endfor;
    print; print;

    print;
    "tau-squared: ";
    print;

    print;
    for jj(1, nyr, 1);
    for qq(((jj-1)*nestim+1), (nestim*jj), 1); "   " tsave[1,qq];;
    endfor; "   " "  ";
    " ";;for qq(((jj-1)*nestim+1), (nestim*jj), 1); " (" tsave[2,qq];; ")";;
    endfor; "   " "  ";
    endfor;
    print; print;


    print;
    "J-statistics: ";
    print;

    for jj(1, nyr, 1);
    for qq(((jj-1)*nestim+2), (nestim*jj), 1); "   " chisave[1,qq];;
    endfor; "   " "  ";
    endfor;
    print;print;

    print;
    "P-values for the J-statistics: ";
    print;

    for jj(1, nyr, 1);
    for qq(((jj-1)*nestim+2), (nestim*jj), 1); "   " chisave[2,qq];;
    endfor; "   " "  ";
    endfor;
    print;print;


    "---------------------------------------------------------------------------";
    "---------------------------------------------------------------------------";
    "---------------------------------------------------------------------------";
    "---------------------------------------------------------------------------";

    print;
    "OLS:";
    print;

    for yy(1,nyr,1);

    "  " isave[1,yy];; "   " bsave[1,yy];; for tt(1, nz, 1); "   " zsave[(tt*2-1),yy];;
    endfor; "   " Rsave[1,yy];; "   " "  ";
    " (" isave[2,yy];; ") (" bsave[2,yy];; for tt(1, nz, 1); ") (" zsave[(tt*2),yy];;
    endfor; ") (" Rsave[2,yy];; ")   " "   ";

    endfor;
    print;print;


    print;
    "GMM:";
    print;


    for qq(1, nestim, 1);

    for yy(1,nyr,1);

    "  " asave[1,nestim*(yy-1)+qq];; "   " csave[1,nestim*(yy-1)+qq];;
    for tt(1, nz, 1); "   " dsave[(tt*2-1),nestim*(yy-1)+qq];; endfor;
    "   " psave[1,nestim*(yy-1)+qq];; "   " tsave[1,nestim*(yy-1)+qq];;
    "   " chisave[1,nestim*(yy-1)+qq];;   "   " "  ";
    " (" asave[2,nestim*(yy-1)+qq];; ") (" csave[2,nestim*(yy-1)+qq];;
    for tt(1, nz, 1); ") (" dsave[(tt*2),nestim*(yy-1)+qq];;   endfor;
    ") (" psave[2,nestim*(yy-1)+qq];; ") (" tsave[2,nestim*(yy-1)+qq];;
    ") (" chisave[2,nestim*(yy-1)+qq];;   ")   " "   ";

    endfor;
    print;print;

    endfor;


    "---------------------------------------------------------------------------";
    "---------------------------------------------------------------------------";
    "---------------------------------------------------------------------------";
    "---------------------------------------------------------------------------";
    if panel == 1;

    print;
    "CMD Intercept:";
    print;


    "  " bisave[1,1];; for qq(1, nestim, 1); "   " casave[1,qq];; endfor; "   ";
    " (" bisave[4,1];; for qq(1, nestim, 1); ") (" casave[4,qq];; endfor; ")  ";
    "  " bisave[2,1];; for qq(1, nestim, 1); "   " casave[2,qq];; endfor; "   ";
    "  " bisave[3,1];; for qq(1, nestim, 1); "   " casave[3,qq];; endfor; "   ";
    print; print;

    print;
    "CMD Coefficient on the mismeasured regressor:";
    print;

    "  " bbsave[1,1];; for qq(1, nestim, 1); "   " ccsave[1,qq];; endfor; "   ";
    " (" bbsave[4,1];; for qq(1, nestim, 1); ") (" ccsave[4,qq];; endfor; ")  ";
    "  " bbsave[2,1];; for qq(1, nestim, 1); "   " ccsave[2,qq];; endfor; "   ";
    "  " bbsave[3,1];; for qq(1, nestim, 1); "   " ccsave[3,qq];; endfor; "   ";
    print; print;

    print;
    "CMD Coefficients on the perfectly measured regressor:";
    print;

    for ee(1, nz, 1);
    "  " bzsave[(ee-1)*4+1,1];; for qq(1, nestim, 1); "   " cdsave[(ee-1)*4+1,qq];;
    endfor; "   " "  ";
    " (" bzsave[(ee-1)*4+4,1];; for qq(1, nestim, 1); ") (" cdsave[(ee-1)*4+4,qq];;
    endfor; ")   " "  ";
    "  " bzsave[(ee-1)*4+2,1];; for qq(1, nestim, 1); "   " cdsave[(ee-1)*4+2,qq];;
    endfor; "   " "  ";
    "  " bzsave[(ee-1)*4+3,1];; for qq(1, nestim, 1); "   " cdsave[(ee-1)*4+3,qq];;
    endfor; "   " "  ";
    print; print;
    endfor;

    print;
    "CMD R-squared:";
    print;

    "  " brsave[1,1];; for qq(1, nestim, 1); "   " cpsave[1,qq];; endfor; "   ";
    " (" brsave[4,1];; for qq(1, nestim, 1); ") (" cpsave[4,qq];; endfor; ")  ";
    "  " brsave[2,1];; for qq(1, nestim, 1); "   " cpsave[2,qq];; endfor; "   ";
    "  " brsave[3,1];; for qq(1, nestim, 1); "   " cpsave[3,qq];; endfor; "   ";
    print; print;

    print;
    "CMD tau-squared:";
    print;

    for qq(1, nestim, 1); "   " ctsave[1,qq];; endfor; "   " "  ";
    " ";;for qq(1, nestim, 1); " (" ctsave[4,qq];; ")";;  endfor; "   " "  ";
    for qq(1, nestim, 1); "   " ctsave[2,qq];; endfor; "   " "  ";
    for qq(1, nestim, 1); "   " ctsave[3,qq];; endfor; "   " "  ";

    endif;

    endif;

    endif;

if bb == 0; cls; endif;
bb = bb + 1;
endo;






/* This business prints out the bootstrapped critical values and p-values for the CMD estimators. */


    pvalues    = zeros(nestim + 1,4+nz);
    locritvals = zeros(nestim + 1,4+nz);
    hicritvals = zeros(nestim + 1,4+nz);
    for estim(0,nestim,1);
    for ffff(1,4+nz,1);
       if ffff == 1;            if estim == 0;    pivotsort = sortc(ipivot[1,.]',1);                   distto = abs((bisave[1,estim]/bisave[4,estim])- pivotsort);
                                         else;    pivotsort = sortc(apivot[estim,.]',1);               distto = abs((casave[1,estim]/casave[4,estim])- pivotsort);
                                         endif;
       elseif ffff == 2;        if estim == 0;    pivotsort = sortc(bpivot[1,.]',1);                   distto = abs((bbsave[1,estim]/bisave[4,estim])- pivotsort);
                                         else;    pivotsort = sortc(cpivot[estim,.]',1);               distto = abs((ccsave[1,estim]/ccsave[4,estim])- pivotsort);
                                         endif;
       elseif ffff == 3;        if estim == 0;    pivotsort = sortc(rpivot[1,.]',1);                   distto = abs((brsave[1,estim]/bisave[4,estim])- pivotsort);
                                         else;    pivotsort = sortc(ppivot[estim,.]',1);               distto = abs((cpsave[1,estim]/cpsave[4,estim])- pivotsort);
                                         endif;
       elseif ffff == 4;        if estim == 0;    pivotsort = zeros(nboot,1);                          distto = zeros(nboot,1);
                                         else;    pivotsort = sortc(tpivot[estim,.]',1);               distto = abs((ctsave[1,estim]/ctsave[4,estim])- pivotsort);
                                         endif;
       elseif ffff >  4;        if estim == 0;    pivotsort = sortc(zpivot[ffff-4,.]',1);              distto = abs((bzsave[(ffff-5)*4+1,1]/bzsave[(ffff-5)*4+4,1])- pivotsort);
                                         else;    pivotsort = sortc(dpivot[(estim-1)*nz+ffff-4,.]',1); distto = abs((cdsave[(ffff-5)*4+1,estim]/cdsave[(ffff-5)*4+4,estim])- pivotsort);
                                         endif;
       endif;
       locritvals[estim+1,ffff] = pivotsort[maxc(1|floor(rows(pivotsort)*0.025)),1];
       hicritvals[estim+1,ffff] = pivotsort[ceil(rows(pivotsort)*0.975),1];
       pvalues[estim+1,ffff]    = 1-indnv(minc(distto),distto)/rows(distto);
    endfor;
    endfor;


"---------------------------------------------------------------------------";
"---------------------------------------------------------------------------";
"---------------------------------------------------------------------------";
"---------------------------------------------------------------------------";


   print;
   "CMD Intercept:";
   print;


   "& " bisave[1,1];;  for qq(1, nestim, 1); " & " casave[1,qq];;    endfor; "  \\" "\\";
   "&(" bisave[4,1];;  for qq(1, nestim, 1); ")&(" casave[4,qq];;    endfor; ") \\" "\\";
   "&[" pvalues[1,1];; for qq(1, nestim, 1); "]&[" pvalues[qq+1,1];; endfor; "] \\" "\\";
   "& " bisave[2,1];;  for qq(1, nestim, 1); " & " casave[2,qq];;    endfor; "  \\" "\\";
   "& " bisave[3,1];;  for qq(1, nestim, 1); " & " casave[3,qq];;    endfor; "  \\" "\\";
   "& " ((bisave[1,1]/bisave[4,1]) .> hicritvals[1,1] .or (bisave[1,1]/bisave[4,1]) .< locritvals[1,1]);;
   for qq(1, nestim, 1); " & " ((casave[1,qq]/casave[4,qq]) .> hicritvals[qq+1,1] .or (casave[1,qq]/casave[4,qq]) .< locritvals[qq+1,1]);; endfor; "  \\" "\\";
   print; print;

   print;
   "CMD Coefficient on the mismeasured regressor:";
   print;

   "& " bbsave[1,1];; for qq(1, nestim, 1); " & " ccsave[1,qq];; endfor; " \\" "\\";
   "&(" bbsave[4,1];; for qq(1, nestim, 1); ")&(" ccsave[4,qq];; endfor; ") \\" "\\";
   "&[" pvalues[1,2];; for qq(1, nestim, 1); "]&[" pvalues[qq+1,2];; endfor; "] \\" "\\";
   "& " bbsave[2,1];; for qq(1, nestim, 1); " & " ccsave[2,qq];; endfor; " \\" "\\";
   "& " bbsave[3,1];; for qq(1, nestim, 1); " & " ccsave[3,qq];; endfor; " \\" "\\";
   "& " ((bbsave[1,1]/bbsave[4,1]) .> hicritvals[1,2] .or (bbsave[1,1]/bbsave[4,1]) .< locritvals[1,2]);;
   for qq(1, nestim, 1); " & " ((ccsave[1,qq]/ccsave[4,qq]) .> hicritvals[qq+1,2] .or (ccsave[1,qq]/ccsave[4,qq]) .< locritvals[qq+1,2]);; endfor; " \\" "\\";
   print; print;

   print;
   "CMD Coefficients on the perfectly measured regressor:";
   print;

   for ee(1, nz, 1);
   "& " bzsave[(ee-1)*4+1,1];; for qq(1, nestim, 1); " & " cdsave[(ee-1)*4+1,qq];; endfor; " \\" "\\";
   "&(" bzsave[(ee-1)*4+4,1];; for qq(1, nestim, 1); ")&(" cdsave[(ee-1)*4+4,qq];; endfor; ") \\" "\\";
   "&[" pvalues[1,4+ee];; for qq(1, nestim, 1); "]&[" pvalues[qq+1,4+ee];; endfor; "] \\" "\\";
   "& " bzsave[(ee-1)*4+2,1];; for qq(1, nestim, 1); " & " cdsave[(ee-1)*4+2,qq];; endfor; " \\" "\\";
   "& " bzsave[(ee-1)*4+3,1];; for qq(1, nestim, 1); " & " cdsave[(ee-1)*4+3,qq];; endfor; " \\" "\\";
   "& " ((bzsave[(ee-1)*4+1,1]/bzsave[(ee-1)*4+4,1]) .> hicritvals[1,4+ee] .or (bzsave[(ee-1)*4+1,1]/bzsave[(ee-1)*4+4,1]) .< locritvals[1,4+ee]);;
   for qq(1, nestim, 1); " & " ((cdsave[(ee-1)*4+1,qq]/cdsave[(ee-1)*4+4,qq]) .> hicritvals[qq+1,4+ee] .or (cdsave[(ee-1)*4+1,qq]/cdsave[(ee-1)*4+4,qq]) .< locritvals[qq+1,4+ee]);; endfor; " \\" "\\";
   print; print;
   endfor;

   print;
   "CMD R-squared:";
   print;

   "& " brsave[1,1];; for qq(1, nestim, 1); " & " cpsave[1,qq];; endfor; " \\" "\\";
   "&(" brsave[4,1];; for qq(1, nestim, 1); ")&(" cpsave[4,qq];; endfor; ") \\" "\\";
   "&[" pvalues[1,3];; for qq(1, nestim, 1); "]&[" pvalues[qq+1,3];; endfor; "] \\" "\\";
   "& " brsave[2,1];; for qq(1, nestim, 1); " & " cpsave[2,qq];; endfor; " \\" "\\";
   "& " brsave[3,1];; for qq(1, nestim, 1); " & " cpsave[3,qq];; endfor; " \\" "\\";
   "& " ((brsave[1,1]/bisave[4,1]) .> hicritvals[1,3] .or (brsave[1,1]/bisave[4,1]) .< locritvals[1,3]);;
   for qq(1, nestim, 1); " & " ((cpsave[1,qq]/cpsave[4,qq]) .> hicritvals[qq+1,3] .or (cpsave[1,qq]/cpsave[4,qq]) .< locritvals[qq+1,3]);; endfor; " \\" "\\";
   print; print;

   print;
   "CMD tau-squared:";
   print;

   for qq(1, nestim, 1); " & " ctsave[1,qq];; endfor; " \\" "\\";
   " ";;for qq(1, nestim, 1); "&(" ctsave[4,qq];; ")";; endfor; " \\" "\\";
   " ";;for qq(1, nestim, 1); "&[" pvalues[qq+1,4];; "]";; endfor; " \\" "\\";
   for qq(1, nestim, 1); " & " ctsave[2,qq];; endfor; " \\" "\\";
   for qq(1, nestim, 1); " & " ctsave[3,qq];; endfor; " \\" "\\";
   for qq(1, nestim, 1); " & " ((ctsave[1,qq]/ctsave[4,qq]) .> hicritvals[qq+1,4] .or (ctsave[1,qq]/ctsave[4,qq]) .< locritvals[qq+1,4]);; endfor; " \\" "\\";











goto eop;

@ ----------------- Subroutines follow --------------- @

@ Subroutine to define the f vector for the optimal weighting matrix. @

proc optw(a,mmm,ooo,wf);
local f, ei;

f = zeros(n,neq);

f[.,1] = varget(mmm[1]) - varget(ooo[1]);
f[.,2] = varget(mmm[2]) - varget(ooo[2]);
f[.,3] = varget(mmm[3]) - varget(ooo[3]);
f[.,4] = varget(mmm[4]) - varget(ooo[4]);
f[.,5] = varget(mmm[5]) - varget(ooo[5]);
if estim > 1;
f[.,6] = varget(mmm[6]) - varget(ooo[6]);
f[.,7] = varget(mmm[7]) - varget(ooo[7]);
f[.,8] = varget(mmm[8]) - varget(ooo[8]);
endif;
if estim > 2;
f[.,9] = varget(mmm[9]) - varget(ooo[9]);
f[.,10] = varget(mmm[10]) - varget(ooo[10]);
f[.,11] = varget(mmm[11]) - varget(ooo[11]);
f[.,12] = varget(mmm[12]) - varget(ooo[12]);
f[.,13] = varget(mmm[13]) - varget(ooo[13]);
f[.,14] = varget(mmm[14]) - varget(ooo[14]);
endif;
if estim > 3;
f[.,15] = varget(mmm[15]) - varget(ooo[15]);
f[.,16] = varget(mmm[16]) - varget(ooo[16]);
f[.,17] = varget(mmm[17]) - varget(ooo[17]);
f[.,18] = varget(mmm[18]) - varget(ooo[18]);
f[.,19] = varget(mmm[19]) - varget(ooo[19]);
f[.,20] = varget(mmm[20]) - varget(ooo[20]);
f[.,21] = varget(mmm[21]) - varget(ooo[21]);
endif;
if estim > 4;
f[.,22] = varget(mmm[22]) - varget(ooo[22]);
f[.,23] = varget(mmm[23]) - varget(ooo[23]);
f[.,24] = varget(mmm[24]) - varget(ooo[24]);
f[.,25] = varget(mmm[25]) - varget(ooo[25]);
f[.,26] = varget(mmm[26]) - varget(ooo[26]);
f[.,27] = varget(mmm[27]) - varget(ooo[27]);
f[.,28] = varget(mmm[28]) - varget(ooo[28]);
f[.,29] = varget(mmm[29]) - varget(ooo[29]);
endif;
@-------------------------------------------------------@
@  This part makes the standard error adjustment.       @
@-------------------------------------------------------@


if wf == 2;

ei = zeros(n,neq);

ei[.,4] = (-2*Eyx_z'inEzz*zy_d' - Ey2_z'inEzz*zx_d')';
ei[.,5] = (-Ex2_z'inEzz*zy_d' - 2*Eyx_z'inEzz*zx_d')';

if estim > 1;
ei[.,6] = (-3*Ey2x_z'inEzz*zy_d' - Ey3_z'inEzz*zx_d')';
ei[.,7] = (-2*Eyx2_z'inEzz*zy_d' - 2*Ey2x_z'inEzz*zx_d')';
ei[.,8] = (-Ex3_z'inEzz*zy_d' - 3*Eyx2_z'inEzz*zx_d')';
endif;

if estim > 2;
ei[.,9]  = (-3*Ex2_z'inEzz*zx_d')';
ei[.,10] = (-3*Ey2_z'inEzz*zy_d')';
ei[.,11] = (-4*Ey3x_z'inEzz*zy_d' - Ey4_z'inEzz*zx_d')';
ei[.,12] = (-3*Ey2x2_z'inEzz*zy_d' - 2*Ey3x_z'inEzz*zx_d')';
ei[.,13] = (-2*Eyx3_z'inEzz*zy_d' - 3*Ey2x2_z'inEzz*zx_d')';
ei[.,14] = (-Ex4_z'inEzz*zy_d' - 4*Eyx3_z'inEzz*zx_d')';
endif;

if estim > 3;
ei[.,15] = (-4*Ey3_z'inEzz*zy_d')';
ei[.,16] = (-4*Ex3_z'inEzz*zx_d')';

ei[.,17] = (-5*Ey4x_z'inEzz*zy_d' -  Ey5_z'inEzz*zx_d')';
ei[.,18] = (-4*Ey3x2_z'inEzz*zy_d' - 2*Ey4x_z'inEzz*zx_d')';
ei[.,19] = (-3*Ey2x3_z'inEzz*zy_d' - 3*Ey3x2_z'inEzz*zx_d')';
ei[.,20] = (-2*Eyx4_z'inEzz*zy_d' - 4*Ey2x3_z'inEzz*zx_d')';
ei[.,21] = (-Ex5_z'inEzz*zy_d' -   5*Eyx4_z'inEzz*zx_d')';
endif;

if estim > 4;
ei[.,22] = (-5*Ey4_z'inEzz*zy_d')';
ei[.,23] = (-5*Ex4_z'inEzz*zx_d')';
ei[.,24] = (-6*Ey5x_z'inEzz*zy_d' -     Ey6_z'inEzz*zx_d')';
ei[.,25] = (-5*Ey4x2_z'inEzz*zy_d' -  2*Ey5x_z'inEzz*zx_d')';
ei[.,26] = (-4*Ey3x3_z'inEzz*zy_d' -  3*Ey4x2_z'inEzz*zx_d')';
ei[.,27] = (-3*Ey2x4_z'inEzz*zy_d' -  4*Ey3x3_z'inEzz*zx_d')';
ei[.,28] = (-2*Eyx5_z'inEzz*zy_d' -   5*Ey2x4_z'inEzz*zx_d')';
ei[.,29] = (-Ex6_z'inEzz*zy_d'   -    6*Eyx5_z'inEzz*zx_d')';
endif;

f = f + ei;
endif;

retp(f);
clear f;
endp;

@ Subroutine to define the f vector. @

proc deff(a,mmm,wf,effsave);
local f, ei;

f = zeros(1,neq);

f[.,1] = varget(mmm[1]) - (a[1,1]^2)*a[2,1] - a[3,1];
f[.,2] = varget(mmm[2]) - a[1,1]*a[2,1];
f[.,3] = varget(mmm[3]) - a[2,1] - a[4,1];
f[.,4] = varget(mmm[4]) - (a[1,1]^2)*a[5,1];
f[.,5] = varget(mmm[5]) - a[1,1]*a[5,1];

if estim > 1;

f[.,6] = varget(mmm[6]) - (a[1,1]^3)*a[6,1] - 3*a[1,1]*a[2,1]*a[3,1];
f[.,7] = varget(mmm[7]) - (a[1,1]^2)*( a[6,1] + a[2,1]*a[4,1] )
               - a[3,1]*( a[2,1] + a[4,1] );
f[.,8] = varget(mmm[8]) - a[1,1]*a[6,1] - 3*a[1,1]*a[2,1]*a[4,1];

endif;
if estim > 2;

f[.,9] = varget(mmm[9]) - a[5,1] - a[7,1];

f[.,10] = varget(mmm[10]) - (a[1,1]^3)*a[5,1] - a[8,1];

f[.,11] = varget(mmm[11]) - (a[1,1]^4)*a[9,1] - 6*(a[1,1]^2)*a[5,1]*a[3,1]
               - 4*a[1,1]*a[2,1]*a[8,1];

f[.,12] = varget(mmm[12]) - (a[1,1]^3)*( a[9,1] + a[5,1]*a[4,1] )
                - 3*a[1,1]*a[5,1]*a[3,1] - a[8,1]*( a[2,1] + a[4,1] );

f[.,13] = varget(mmm[13]) - (a[1,1]^2)*(a[9,1] +3*a[5,1]*a[4,1] + a[2,1]*a[7,1])
                - a[3,1]*( a[5,1] + a[7,1] );

f[.,14] = varget(mmm[14]) - a[1,1]*(a[9,1] + 6*a[5,1]*a[4,1] + 4*a[2,1]*a[7,1]);

endif;
if estim > 3;

f[.,15] = varget(mmm[15]) - (a[1,1]^4)*a[6,1] - 6*(a[1,1]^2)*a[2,1]*a[3,1]
                          - a[10,1];

f[.,16] = varget(mmm[16]) - a[6,1] - 6*a[2,1]*a[4,1] - a[11,1];

f[.,17] = varget(mmm[17]) - (a[1,1]^5)*a[12,1] - 10*(a[1,1]^3)*a[6,1]*a[3,1]
                     - 10*(a[1,1]^2)*a[5,1]*a[8,1] - 5*a[1,1]*a[2,1]*a[10,1];

f[.,18] = varget(mmm[18]) - (a[1,1]^4)*( a[12,1] + a[6,1]*a[4,1] )
                       - 6*(a[1,1]^2)*( a[6,1]*a[3,1] + a[2,1]*a[3,1]*a[4,1] )
                       - 4*a[1,1]*a[5,1]*a[8,1] - a[10,1]*( a[2,1] + a[4,1] );

f[.,19] = varget(mmm[19]) - (a[1,1]^3)*(a[12,1]+3*a[6,1]*a[4,1]+ a[5,1]*a[7,1])
                          - a[1,1]*( 3*a[6,1]*a[3,1] + 9*a[2,1]*a[3,1]*a[4,1] )
                          - a[8,1]*( a[5,1] + a[7,1] );

f[.,20] = varget(mmm[20]) - (a[1,1]^2)*( a[12,1] + 6*a[6,1]*a[4,1]
                             + 4*a[5,1]*a[7,1] + a[2,1]*a[11,1] )
                           - a[3,1]*( a[6,1] + 6*a[2,1]*a[4,1] + a[11,1] );

f[.,21] = varget(mmm[21]) - a[1,1]*( a[12,1] + 10*a[6,1]*a[4,1]
                                     + 10*a[5,1]*a[7,1] + 5*a[2,1]*a[11,1] );

endif;
if estim > 4;

f[.,22] = varget(mmm[22]) - (a[1,1]^5)*a[9,1] - 10*(a[1,1]^3)*a[5,1]*a[3,1]
              - 10*(a[1,1]^2)*a[2,1]*a[8,1] - a[13,1];


f[.,23] = varget(mmm[23])
         - a[9,1] - 10*a[5,1]*a[4,1] - 10*a[2,1]*a[7,1] - a[14,1];

f[.,24] = varget(mmm[24]) - (a[1,1]^6)*a[15,1]
               - 15*(a[1,1]^4)*a[9,1]*a[3,1]
               - 20*(a[1,1]^3)*a[6,1]*a[8,1]
               - 15*(a[1,1]^2)*a[5,1]*a[10,1]
               - 6*a[1,1]*a[2,1]*a[13,1];

f[.,25] = varget(mmm[25]) - (a[1,1]^5)*( a[15,1] + a[9,1]*a[4,1] )
                - 10*(a[1,1]^3)*( a[9,1]*a[3,1] + a[5,1]*a[3,1]*a[4,1] )
                - 10*(a[1,1]^2)*( a[6,1]*a[8,1] + a[2,1]*a[8,1]*a[4,1] )
                - 5*a[1,1]*a[5,1]*a[10,1]
                - a[13,1]*( a[2,1] + a[4,1] );


f[.,26] = varget(mmm[26])
- (a[1,1]^4)*( a[15,1] + 3*a[9,1]*a[4,1] + a[6,1]*a[7,1] )
- 6*(a[1,1]^2)*( a[9,1]*a[3,1] + 3*a[5,1]*a[3,1]*a[4,1] + a[2,1]*a[3,1]*a[7,1] )
- 4*a[1,1]*( a[6,1]*a[8,1] + 3*a[2,1]*a[8,1]*a[4,1] )
- a[10,1]*( a[5,1] + a[7,1] );


f[.,27] = varget(mmm[27])
- (a[1,1]^3)*( a[15,1] + 6*a[9,1]*a[4,1] + 4*a[6,1]*a[7,1] + a[5,1]*a[11,1] )
- a[1,1]*( 3*a[9,1]*a[3,1] + 18*a[5,1]*a[3,1]*a[4,1] + 12*a[2,1]*a[3,1]*a[7,1] )
-a[8,1]*( a[6,1] + 6*a[2,1]*a[4,1] + a[11,1] );

f[.,28] = varget(mmm[28])
- (a[1,1]^2)*( a[15,1] + 10*a[9,1]*a[4,1] + 10*a[6,1]*a[7,1] + 5*a[5,1]*a[11,1]
+ a[2,1]*a[14,1] )
- a[3,1]*( a[9,1] + 10*a[5,1]*a[4,1] + 10*a[2,1]*a[7,1] + a[14,1] );

f[.,29] = varget(mmm[29])
- a[1,1]*( a[15,1] + 15*a[9,1]*a[4,1] + 20*a[6,1]*a[7,1] +15*a[5,1]*a[11,1]
+ 6*a[2,1]*a[14,1] );

endif;

if wf == 0; f = f'; endif;

f = f -(bb>0)*effsave[1:neq,(estim-1)*nestim+year];

retp(f);
clear f;
endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad1(a);

local dvdc1,dvdc2,dvdc3,dvdc4,dvdc5;
dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];

dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;

dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;

dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;

dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];


retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5;
endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad2(a);

local dvdc1,dvdc2,dvdc3,dvdc4,dvdc5,dvdc6;
dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);
dvdc6  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];
dvdc1[6,1] = -3*(a[1,1]^2)*a[6,1] - 3*a[2,1]*a[3,1];
dvdc1[7,1] = -2*a[1,1]*( a[6,1] + a[2,1]*a[4,1] );
dvdc1[8,1] = -a[6,1] - 3*a[2,1]*a[4,1];

dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;
dvdc2[6,1] = -3*a[1,1]*a[3,1];
dvdc2[7,1] = -(a[1,1]^2)*a[4,1] - a[3,1];
dvdc2[8,1] = -3*a[1,1]*a[4,1];

dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;
dvdc3[6,1] = -3*a[1,1]*a[2,1];
dvdc3[7,1] = -( a[2,1] + a[4,1] );
dvdc3[8,1] = 0;

dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;
dvdc4[6,1] = 0;
dvdc4[7,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc4[8,1] = -3*a[1,1]*a[2,1];

dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];
dvdc5[6,1] = 0;
dvdc5[7,1] = 0;
dvdc5[8,1] = 0;

dvdc6[1,1] = 0;
dvdc6[2,1] = 0;
dvdc6[3,1] = 0;
dvdc6[4,1] = 0;
dvdc6[5,1] = 0;
dvdc6[6,1] = -(a[1,1]^3);
dvdc6[7,1] = -(a[1,1]^2);
dvdc6[8,1] = -a[1,1];


retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5~dvdc6);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5, dvdc6;
endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad3(a);

local dvdc,dvdc1,dvdc2,dvdc3,dvdc4,dvdc5,dvdc6,dvdc7,dvdc8,dvdc9;

dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);
dvdc6  = zeros(neq,1);
dvdc7  = zeros(neq,1);
dvdc8  = zeros(neq,1);
dvdc9  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];
dvdc1[6,1] = -3*(a[1,1]^2)*a[6,1] - 3*a[2,1]*a[3,1];
dvdc1[7,1] = -2*a[1,1]*( a[6,1] + a[2,1]*a[4,1] );
dvdc1[8,1] = -a[6,1] - 3*a[2,1]*a[4,1];
dvdc1[9,1] = 0;
dvdc1[10,1] = -3*(a[1,1]^2)*a[5,1];
dvdc1[11,1] = -4*(a[1,1]^3)*a[9,1] - 12*a[1,1]*a[5,1]*a[3,1]
              - 4*a[2,1]*a[8,1];
dvdc1[12,1] = -3*(a[1,1]^2)*( a[9,1] + a[5,1]*a[4,1] )
              - 3*a[5,1]*a[3,1];
dvdc1[13,1] = -2*a[1,1]*( a[9,1] + 3*a[5,1]*a[4,1] + a[2,1]*a[7,1] );
dvdc1[14,1] = -( a[9,1] + 6*a[5,1]*a[4,1] + 4*a[2,1]*a[7,1] );

dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;
dvdc2[6,1] = -3*a[1,1]*a[3,1];
dvdc2[7,1] = -(a[1,1]^2)*a[4,1] - a[3,1];
dvdc2[8,1] = -3*a[1,1]*a[4,1];
dvdc2[9,1] = 0;
dvdc2[10,1] = 0;
dvdc2[11,1] = -4*a[1,1]*a[8,1];
dvdc2[12,1] = -a[8,1];
dvdc2[13,1] = -(a[1,1]^2)*a[7,1];
dvdc2[14,1] = -4*a[1,1]*a[7,1];

dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;
dvdc3[6,1] = -3*a[1,1]*a[2,1];
dvdc3[7,1] = -( a[2,1] + a[4,1] );
dvdc3[8,1] = 0;
dvdc3[9,1] = 0;
dvdc3[10,1] = 0;
dvdc3[11,1] = -6*(a[1,1]^2)*a[5,1];
dvdc3[12,1] = -3*a[1,1]*a[5,1];
dvdc3[13,1] = -( a[5,1] + a[7,1] );
dvdc3[14,1] = 0;

dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;
dvdc4[6,1] = 0;
dvdc4[7,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc4[8,1] = -3*a[1,1]*a[2,1];
dvdc4[9,1] = 0;
dvdc4[10,1] = 0;
dvdc4[11,1] = 0;
dvdc4[12,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc4[13,1] = -3*(a[1,1]^2)*a[5,1];
dvdc4[14,1] = -6*a[1,1]*a[5,1];

dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];
dvdc5[6,1] = 0;
dvdc5[7,1] = 0;
dvdc5[8,1] = 0;
dvdc5[9,1] = -1;
dvdc5[10,1] = -(a[1,1]^3);
dvdc5[11,1] = -6*(a[1,1]^2)*a[3,1];
dvdc5[12,1] = -(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc5[13,1] = -3*(a[1,1]^2)*a[4,1] -a[3,1];
dvdc5[14,1] = -6*a[1,1]*a[4,1];

dvdc6[1,1] = 0;
dvdc6[2,1] = 0;
dvdc6[3,1] = 0;
dvdc6[4,1] = 0;
dvdc6[5,1] = 0;
dvdc6[6,1] = -(a[1,1]^3);
dvdc6[7,1] = -(a[1,1]^2);
dvdc6[8,1] = -a[1,1];
dvdc6[9,1] = 0;
dvdc6[10,1] = 0;
dvdc6[11,1] = 0;
dvdc6[12,1] = 0;
dvdc6[13,1] = 0;
dvdc6[14,1] = 0;

dvdc7[1,1] = 0;
dvdc7[2,1] = 0;
dvdc7[3,1] = 0;
dvdc7[4,1] = 0;
dvdc7[5,1] = 0;
dvdc7[6,1] = 0;
dvdc7[7,1] = 0;
dvdc7[8,1] = 0;
dvdc7[9,1] = -1;
dvdc7[10,1] = 0;
dvdc7[11,1] = 0;
dvdc7[12,1] = 0;
dvdc7[13,1] = -(a[1,1]^2)*a[2,1] -a[3,1];
dvdc7[14,1] = -4*a[1,1]*a[2,1];

dvdc8[1,1] = 0;
dvdc8[2,1] = 0;
dvdc8[3,1] = 0;
dvdc8[4,1] = 0;
dvdc8[5,1] = 0;
dvdc8[6,1] = 0;
dvdc8[7,1] = 0;
dvdc8[8,1] = 0;
dvdc8[9,1] = 0;
dvdc8[10,1] = -1;
dvdc8[11,1] = -4*a[1,1]*a[2,1];
dvdc8[12,1] = -( a[2,1] + a[4,1] );
dvdc8[13,1] = 0;
dvdc8[14,1] = 0;

dvdc9[1,1] = 0;
dvdc9[2,1] = 0;
dvdc9[3,1] = 0;
dvdc9[4,1] = 0;
dvdc9[5,1] = 0;
dvdc9[6,1] = 0;
dvdc9[7,1] = 0;
dvdc9[8,1] = 0;
dvdc9[9,1] = 0;
dvdc9[10,1] = 0;
dvdc9[11,1] = -(a[1,1]^4);
dvdc9[12,1] = -(a[1,1]^3);
dvdc9[13,1] = -(a[1,1]^2);
dvdc9[14,1] = -a[1,1];


retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5~dvdc6~dvdc7~dvdc8~dvdc9);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5, dvdc6, dvdc7, dvdc8, dvdc9;

endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad4(a);

local dvdc,dvdc1,dvdc2,dvdc3,dvdc4,dvdc5,dvdc6,dvdc7,dvdc8,dvdc9,
dvdc10,dvdc11,dvdc12;

dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);
dvdc6  = zeros(neq,1);
dvdc7  = zeros(neq,1);
dvdc8  = zeros(neq,1);
dvdc9  = zeros(neq,1);
dvdc10  = zeros(neq,1);
dvdc11  = zeros(neq,1);
dvdc12  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];
dvdc1[6,1] = -3*(a[1,1]^2)*a[6,1] - 3*a[2,1]*a[3,1];
dvdc1[7,1] = -2*a[1,1]*( a[6,1] + a[2,1]*a[4,1] );
dvdc1[8,1] = -a[6,1] - 3*a[2,1]*a[4,1];
dvdc1[9,1] = 0;
dvdc1[10,1] = -3*(a[1,1]^2)*a[5,1];
dvdc1[11,1] = -4*(a[1,1]^3)*a[9,1] - 12*a[1,1]*a[5,1]*a[3,1]
              - 4*a[2,1]*a[8,1];
dvdc1[12,1] = -3*(a[1,1]^2)*( a[9,1] + a[5,1]*a[4,1] )
              - 3*a[5,1]*a[3,1];
dvdc1[13,1] = -2*a[1,1]*( a[9,1] + 3*a[5,1]*a[4,1] + a[2,1]*a[7,1] );
dvdc1[14,1] = -( a[9,1] + 6*a[5,1]*a[4,1] + 4*a[2,1]*a[7,1] );

dvdc1[15,1] = -4*(a[1,1]^3)*a[6,1] -12*a[1,1]*a[2,1]*a[3,1];
dvdc1[16,1] = 0;
dvdc1[17,1] = -5*(a[1,1]^4)*a[12,1] - 30*(a[1,1]^2)*a[6,1]*a[3,1]
              -20*a[1,1]*a[5,1]*a[8,1] - 5*a[2,1]*a[10,1];
dvdc1[18,1] = -4*(a[1,1]^3)*( a[12,1] + a[6,1]*a[4,1] )
              -12*a[1,1]*( a[6,1]*a[3,1] + a[2,1]*a[3,1]*a[4,1] )
              -4*a[5,1]*a[8,1];
dvdc1[19,1] = -3*(a[1,1]^2)*( a[12,1] + 3*a[6,1]*a[4,1] + a[5,1]*a[7,1] )
              -( 3*a[6,1]*a[3,1] + 9*a[2,1]*a[3,1]*a[4,1] );
dvdc1[20,1] = -2*a[1,1]*( a[12,1] + 6*a[6,1]*a[4,1] + 4*a[5,1]*a[7,1]
                         + a[2,1]*a[11,1] );
dvdc1[21,1] = - ( a[12,1] + 10*a[6,1]*a[4,1] + 10*a[5,1]*a[7,1]
                  + 5*a[2,1]*a[11,1] );


dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;
dvdc2[6,1] = -3*a[1,1]*a[3,1];
dvdc2[7,1] = -(a[1,1]^2)*a[4,1] - a[3,1];
dvdc2[8,1] = -3*a[1,1]*a[4,1];
dvdc2[9,1] = 0;
dvdc2[10,1] = 0;
dvdc2[11,1] = -4*a[1,1]*a[8,1];
dvdc2[12,1] = -a[8,1];
dvdc2[13,1] = -(a[1,1]^2)*a[7,1];
dvdc2[14,1] = -4*a[1,1]*a[7,1];

dvdc2[15,1] = -6*(a[1,1]^2)*a[3,1];
dvdc2[16,1] = -6*a[4,1];
dvdc2[17,1] = -5*a[1,1]*a[10,1];
dvdc2[18,1] = -6*(a[1,1]^2)*a[3,1]*a[4,1] - a[10,1];
dvdc2[19,1] = -9*a[1,1]*a[3,1]*a[4,1];
dvdc2[20,1] = -(a[1,1]^2)*a[11,1] - 6*a[3,1]*a[4,1];
dvdc2[21,1] = -5*a[1,1]*a[11,1];



dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;
dvdc3[6,1] = -3*a[1,1]*a[2,1];
dvdc3[7,1] = -( a[2,1] + a[4,1] );
dvdc3[8,1] = 0;
dvdc3[9,1] = 0;
dvdc3[10,1] = 0;
dvdc3[11,1] = -6*(a[1,1]^2)*a[5,1];
dvdc3[12,1] = -3*a[1,1]*a[5,1];
dvdc3[13,1] = -( a[5,1] + a[7,1] );
dvdc3[14,1] = 0;

dvdc3[15,1] = -6*(a[1,1]^2)*a[2,1];
dvdc3[16,1] = 0;
dvdc3[17,1] = -10*(a[1,1]^3)*a[6,1];
dvdc3[18,1] = -6*(a[1,1]^2)*( a[6,1] + a[2,1] );
dvdc3[19,1] = -a[1,1]*( 3*a[6,1] + 9*a[2,1]*a[4,1] );
dvdc3[20,1] = -( a[6,1] + 6*a[2,1]*a[4,1] + a[11,1] );
dvdc3[21,1] = 0;



dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;
dvdc4[6,1] = 0;
dvdc4[7,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc4[8,1] = -3*a[1,1]*a[2,1];
dvdc4[9,1] = 0;
dvdc4[10,1] = 0;
dvdc4[11,1] = 0;
dvdc4[12,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc4[13,1] = -3*(a[1,1]^2)*a[5,1];
dvdc4[14,1] = -6*a[1,1]*a[5,1];

dvdc4[15,1] = 0;
dvdc4[16,1] = -6*a[2,1];
dvdc4[17,1] = 0;
dvdc4[18,1] = -(a[1,1]^4)*a[6,1] - 6*(a[1,1]^2)*a[2,1]*a[3,1] - a[10,1];
dvdc4[19,1] = -3*(a[1,1]^3)*a[6,1] - 9*a[1,1]*a[2,1]*a[3,1];
dvdc4[20,1] = -6*(a[1,1]^2)*a[6,1] - 6*a[2,1]*a[3,1];
dvdc4[21,1] = -10*a[1,1]*a[6,1];



dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];
dvdc5[6,1] = 0;
dvdc5[7,1] = 0;
dvdc5[8,1] = 0;
dvdc5[9,1] = -1;
dvdc5[10,1] = -(a[1,1]^3);
dvdc5[11,1] = -6*(a[1,1]^2)*a[3,1];
dvdc5[12,1] = -(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc5[13,1] = -3*(a[1,1]^2)*a[4,1] -a[3,1];
dvdc5[14,1] = -6*a[1,1]*a[4,1];

dvdc5[15,1] = 0;
dvdc5[16,1] = 0;
dvdc5[17,1] = -10*(a[1,1]^2)*a[8,1];
dvdc5[18,1] = -4*a[1,1]*a[8,1];
dvdc5[19,1] = -(a[1,1]^3)*a[7,1] - a[8,1];
dvdc5[20,1] = -4*(a[1,1]^2)*a[7,1];
dvdc5[21,1] = -10*a[1,1]*a[7,1];



dvdc6[1,1] = 0;
dvdc6[2,1] = 0;
dvdc6[3,1] = 0;
dvdc6[4,1] = 0;
dvdc6[5,1] = 0;
dvdc6[6,1] = -(a[1,1]^3);
dvdc6[7,1] = -(a[1,1]^2);
dvdc6[8,1] = -a[1,1];
dvdc6[9,1] = 0;
dvdc6[10,1] = 0;
dvdc6[11,1] = 0;
dvdc6[12,1] = 0;
dvdc6[13,1] = 0;
dvdc6[14,1] = 0;

dvdc6[15,1] = -(a[1,1]^4);
dvdc6[16,1] = -1;
dvdc6[17,1] = -10*(a[1,1]^3)*a[3,1];
dvdc6[18,1] = -(a[1,1]^4)*a[4,1] - 6*(a[1,1]^2)*a[3,1];
dvdc6[19,1] = -3*(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc6[20,1] = -6*(a[1,1]^2)*a[4,1] - a[3,1];
dvdc6[21,1] = -10*a[1,1]*a[4,1];



dvdc7[1,1] = 0;
dvdc7[2,1] = 0;
dvdc7[3,1] = 0;
dvdc7[4,1] = 0;
dvdc7[5,1] = 0;
dvdc7[6,1] = 0;
dvdc7[7,1] = 0;
dvdc7[8,1] = 0;
dvdc7[9,1] = -1;
dvdc7[10,1] = 0;
dvdc7[11,1] = 0;
dvdc7[12,1] = 0;
dvdc7[13,1] = -(a[1,1]^2)*a[2,1] -a[3,1];
dvdc7[14,1] = -4*a[1,1]*a[2,1];

dvdc7[15,1] = 0;
dvdc7[16,1] = 0;
dvdc7[17,1] = 0;
dvdc7[18,1] = 0;
dvdc7[19,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc7[20,1] = -4*(a[1,1]^2)*a[5,1];
dvdc7[21,1] = -10*a[1,1]*a[5,1];



dvdc8[1,1] = 0;
dvdc8[2,1] = 0;
dvdc8[3,1] = 0;
dvdc8[4,1] = 0;
dvdc8[5,1] = 0;
dvdc8[6,1] = 0;
dvdc8[7,1] = 0;
dvdc8[8,1] = 0;
dvdc8[9,1] = 0;
dvdc8[10,1] = -1;
dvdc8[11,1] = -4*a[1,1]*a[2,1];
dvdc8[12,1] = -( a[2,1] + a[4,1] );
dvdc8[13,1] = 0;
dvdc8[14,1] = 0;

dvdc8[15,1] = 0;
dvdc8[16,1] = 0;
dvdc8[17,1] = -10*(a[1,1]^2)*a[5,1];
dvdc8[18,1] = -4*a[1,1]*a[5,1];
dvdc8[19,1] = -( a[5,1] + a[7,1] );
dvdc8[20,1] = 0;
dvdc8[21,1] = 0;



dvdc9[1,1] = 0;
dvdc9[2,1] = 0;
dvdc9[3,1] = 0;
dvdc9[4,1] = 0;
dvdc9[5,1] = 0;
dvdc9[6,1] = 0;
dvdc9[7,1] = 0;
dvdc9[8,1] = 0;
dvdc9[9,1] = 0;
dvdc9[10,1] = 0;
dvdc9[11,1] = -(a[1,1]^4);
dvdc9[12,1] = -(a[1,1]^3);
dvdc9[13,1] = -(a[1,1]^2);
dvdc9[14,1] = -a[1,1];

dvdc9[15,1] = 0;
dvdc9[16,1] = 0;
dvdc9[17,1] = 0;
dvdc9[18,1] = 0;
dvdc9[19,1] = 0;
dvdc9[20,1] = 0;
dvdc9[21,1] = 0;


dvdc10[1,1] = 0;
dvdc10[2,1] = 0;
dvdc10[3,1] = 0;
dvdc10[4,1] = 0;
dvdc10[5,1] = 0;
dvdc10[6,1] = 0;
dvdc10[7,1] = 0;
dvdc10[8,1] = 0;
dvdc10[9,1] = 0;
dvdc10[10,1] = 0;
dvdc10[11,1] = 0;
dvdc10[12,1] = 0;
dvdc10[13,1] = 0;
dvdc10[14,1] = 0;

dvdc10[15,1] = -1;
dvdc10[16,1] = 0;
dvdc10[17,1] = -5*a[1,1]*a[2,1];
dvdc10[18,1] = -( a[2,1] + a[4,1] );
dvdc10[19,1] = 0;
dvdc10[20,1] = 0;
dvdc10[21,1] = 0;


dvdc11[1,1] = 0;
dvdc11[2,1] = 0;
dvdc11[3,1] = 0;
dvdc11[4,1] = 0;
dvdc11[5,1] = 0;
dvdc11[6,1] = 0;
dvdc11[7,1] = 0;
dvdc11[8,1] = 0;
dvdc11[9,1] = 0;
dvdc11[10,1] = 0;
dvdc11[11,1] = 0;
dvdc11[12,1] = 0;
dvdc11[13,1] = 0;
dvdc11[14,1] = 0;

dvdc11[15,1] = 0;
dvdc11[16,1] = -1;
dvdc11[17,1] = 0;
dvdc11[18,1] = 0;
dvdc11[19,1] = 0;
dvdc11[20,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc11[21,1] = -5*a[1,1]*a[2,1];


dvdc12[1,1] = 0;
dvdc12[2,1] = 0;
dvdc12[3,1] = 0;
dvdc12[4,1] = 0;
dvdc12[5,1] = 0;
dvdc12[6,1] = 0;
dvdc12[7,1] = 0;
dvdc12[8,1] = 0;
dvdc12[9,1] = 0;
dvdc12[10,1] = 0;
dvdc12[11,1] = 0;
dvdc12[12,1] = 0;
dvdc12[13,1] = 0;
dvdc12[14,1] = 0;

dvdc12[15,1] = 0;
dvdc12[16,1] = 0;
dvdc12[17,1] = -(a[1,1]^5);
dvdc12[18,1] = -(a[1,1]^4);
dvdc12[19,1] = -(a[1,1]^3);
dvdc12[20,1] = -(a[1,1]^2);
dvdc12[21,1] = -a[1,1];



retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5~dvdc6~dvdc7~dvdc8~dvdc9
          ~dvdc10~dvdc11~dvdc12);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5, dvdc6, dvdc7, dvdc8, dvdc9,
      dvdc10, dvdc11, dvdc12;

endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad5(a);

local dvdc,dvdc1,dvdc2,dvdc3,dvdc4,dvdc5,dvdc6,dvdc7,dvdc8,dvdc9,
dvdc10,dvdc11,dvdc12,dvdc13,dvdc14,dvdc15;

dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);
dvdc6  = zeros(neq,1);
dvdc7  = zeros(neq,1);
dvdc8  = zeros(neq,1);
dvdc9  = zeros(neq,1);
dvdc10  = zeros(neq,1);
dvdc11  = zeros(neq,1);
dvdc12  = zeros(neq,1);
dvdc13  = zeros(neq,1);
dvdc14  = zeros(neq,1);
dvdc15  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];
dvdc1[6,1] = -3*(a[1,1]^2)*a[6,1] - 3*a[2,1]*a[3,1];
dvdc1[7,1] = -2*a[1,1]*( a[6,1] + a[2,1]*a[4,1] );
dvdc1[8,1] = -a[6,1] - 3*a[2,1]*a[4,1];
dvdc1[9,1] = 0;
dvdc1[10,1] = -3*(a[1,1]^2)*a[5,1];
dvdc1[11,1] = -4*(a[1,1]^3)*a[9,1] - 12*a[1,1]*a[5,1]*a[3,1]
              - 4*a[2,1]*a[8,1];
dvdc1[12,1] = -3*(a[1,1]^2)*( a[9,1] + a[5,1]*a[4,1] )
              - 3*a[5,1]*a[3,1];
dvdc1[13,1] = -2*a[1,1]*( a[9,1] + 3*a[5,1]*a[4,1] + a[2,1]*a[7,1] );
dvdc1[14,1] = -( a[9,1] + 6*a[5,1]*a[4,1] + 4*a[2,1]*a[7,1] );
dvdc1[15,1] = -4*(a[1,1]^3)*a[6,1] -12*a[1,1]*a[2,1]*a[3,1];
dvdc1[16,1] = 0;
dvdc1[17,1] = -5*(a[1,1]^4)*a[12,1] - 30*(a[1,1]^2)*a[6,1]*a[3,1]
              -20*a[1,1]*a[5,1]*a[8,1] - 5*a[2,1]*a[10,1];
dvdc1[18,1] = -4*(a[1,1]^3)*( a[12,1] + a[6,1]*a[4,1] )
              -12*a[1,1]*( a[6,1]*a[3,1] + a[2,1]*a[3,1]*a[4,1] )
              -4*a[5,1]*a[8,1];
dvdc1[19,1] = -3*(a[1,1]^2)*( a[12,1] + 3*a[6,1]*a[4,1] + a[5,1]*a[7,1] )
              -( 3*a[6,1]*a[3,1] + 9*a[2,1]*a[3,1]*a[4,1] );
dvdc1[20,1] = -2*a[1,1]*( a[12,1] + 6*a[6,1]*a[4,1] + 4*a[5,1]*a[7,1]
                         + a[2,1]*a[11,1] );
dvdc1[21,1] = - ( a[12,1] + 10*a[6,1]*a[4,1] + 10*a[5,1]*a[7,1]
                  + 5*a[2,1]*a[11,1] );

dvdc1[22,1] = -5*(a[1,1]^4)*a[9,1] - 30*(a[1,1]^2)*a[5,1]*a[3,1]
              - 20*a[1,1]*a[2,1]*a[8,1];
dvdc1[23,1] = 0;
dvdc1[24,1] = -6*(a[1,1]^5)*a[15,1] - 60*(a[1,1]^3)*a[9,1]*a[3,1]
              -60*(a[1,1]^2)*a[6,1]*a[8,1] - 30*a[1,1]*a[5,1]*a[10,1]
              -6*a[2,1]*a[13,1];
dvdc1[25,1] = -5*(a[1,1]^4)*(a[15,1] + a[9,1]*a[4,1])
              -30*(a[1,1]^2)*(a[9,1]*a[3,1] + a[5,1]*a[3,1]*a[4,1])
              -20*a[1,1]*(a[6,1]*a[8,1] + a[2,1]*a[8,1]*a[4,1])
              -5*a[5,1]*a[10,1];

dvdc1[26,1] = -4*(a[1,1]^3)*(a[15,1] + 3*a[9,1]*a[4,1] + a[6,1]*a[7,1])
-12*a[1,1]*(a[9,1]*a[3,1] + 3*a[5,1]*a[3,1]*a[4,1] + a[2,1]*a[3,1]*a[7,1])
-4*(a[6,1]*a[8,1] + 3*a[2,1]*a[8,1]*a[4,1]);

dvdc1[27,1] =
-3*(a[1,1]^2)*(a[15,1] + 6*a[9,1]*a[4,1] + 4*a[6,1]*a[7,1] + a[5,1]*a[11,1])
-(3*a[9,1]*a[3,1] + 18*a[5,1]*a[3,1]*a[4,1] + 12*a[2,1]*a[3,1]*a[7,1]);

dvdc1[28,1] =
-2*a[1,1]*(a[15,1] + 10*a[9,1]*a[4,1] + 10*a[6,1]*a[7,1] + 5*a[5,1]*a[11,1]
+ a[2,1]*a[14,1]);

dvdc1[29,1] = -(a[15,1] + 15*a[9,1]*a[4,1] + 20*a[6,1]*a[7,1]
              + 15*a[5,1]*a[11,1] + 6*a[2,1]*a[14,1]);



dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;
dvdc2[6,1] = -3*a[1,1]*a[3,1];
dvdc2[7,1] = -(a[1,1]^2)*a[4,1] - a[3,1];
dvdc2[8,1] = -3*a[1,1]*a[4,1];
dvdc2[9,1] = 0;
dvdc2[10,1] = 0;
dvdc2[11,1] = -4*a[1,1]*a[8,1];
dvdc2[12,1] = -a[8,1];
dvdc2[13,1] = -(a[1,1]^2)*a[7,1];
dvdc2[14,1] = -4*a[1,1]*a[7,1];
dvdc2[15,1] = -6*(a[1,1]^2)*a[3,1];
dvdc2[16,1] = -6*a[4,1];
dvdc2[17,1] = -5*a[1,1]*a[10,1];
dvdc2[18,1] = -6*(a[1,1]^2)*a[3,1]*a[4,1] - a[10,1];
dvdc2[19,1] = -9*a[1,1]*a[3,1]*a[4,1];
dvdc2[20,1] = -(a[1,1]^2)*a[11,1] - 6*a[3,1]*a[4,1];
dvdc2[21,1] = -5*a[1,1]*a[11,1];

dvdc2[22,1] = -10*(a[1,1]^2)*a[8,1];
dvdc2[23,1] = -10*a[7,1];
dvdc2[24,1] = -6*a[1,1]*a[13,1];
dvdc2[25,1] = -10*(a[1,1]^2)*a[8,1]*a[4,1] - a[13,1];
dvdc2[26,1] = -6*(a[1,1]^2)*a[3,1]*a[7,1] - 12*a[1,1]*a[8,1]*a[4,1];
dvdc2[27,1] = -12*a[1,1]*a[3,1]*a[7,1] - 6*a[8,1]*a[4,1];
dvdc2[28,1] = -(a[1,1]^2)*a[14,1] - 10*a[3,1]*a[7,1];
dvdc2[29,1] = -6*a[1,1]*a[14,1];



dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;
dvdc3[6,1] = -3*a[1,1]*a[2,1];
dvdc3[7,1] = -( a[2,1] + a[4,1] );
dvdc3[8,1] = 0;
dvdc3[9,1] = 0;
dvdc3[10,1] = 0;
dvdc3[11,1] = -6*(a[1,1]^2)*a[5,1];
dvdc3[12,1] = -3*a[1,1]*a[5,1];
dvdc3[13,1] = -( a[5,1] + a[7,1] );
dvdc3[14,1] = 0;
dvdc3[15,1] = -6*(a[1,1]^2)*a[2,1];
dvdc3[16,1] = 0;
dvdc3[17,1] = -10*(a[1,1]^3)*a[6,1];
dvdc3[18,1] = -6*(a[1,1]^2)*( a[6,1] + a[2,1] );
dvdc3[19,1] = -a[1,1]*( 3*a[6,1] + 9*a[2,1]*a[4,1] );
dvdc3[20,1] = -( a[6,1] + 6*a[2,1]*a[4,1] + a[11,1] );
dvdc3[21,1] = 0;

dvdc3[22,1] = -10*(a[1,1]^3)*a[5,1];
dvdc3[23,1] = 0;
dvdc3[24,1] = -15*(a[1,1]^4)*a[9,1];
dvdc3[25,1] = -10*(a[1,1]^3)*(a[9,1] + a[5,1]*a[4,1]);
dvdc3[26,1] = -6*(a[1,1]^2)*(a[9,1] + 3*a[5,1]*a[4,1] + a[2,1]*a[7,1]);
dvdc3[27,1] = -a[1,1]*(3*a[9,1] + 18*a[5,1]*a[4,1] + 12*a[2,1]*a[7,1]);
dvdc3[28,1] = -(a[9,1] + 10*a[5,1]*a[4,1] + 10*a[2,1]*a[7,1] + a[14,1]);
dvdc3[29,1] = 0;



dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;
dvdc4[6,1] = 0;
dvdc4[7,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc4[8,1] = -3*a[1,1]*a[2,1];
dvdc4[9,1] = 0;
dvdc4[10,1] = 0;
dvdc4[11,1] = 0;
dvdc4[12,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc4[13,1] = -3*(a[1,1]^2)*a[5,1];
dvdc4[14,1] = -6*a[1,1]*a[5,1];
dvdc4[15,1] = 0;
dvdc4[16,1] = -6*a[2,1];
dvdc4[17,1] = 0;
dvdc4[18,1] = -(a[1,1]^4)*a[6,1] - 6*(a[1,1]^2)*a[2,1]*a[3,1] - a[10,1];
dvdc4[19,1] = -3*(a[1,1]^3)*a[6,1] - 9*a[1,1]*a[2,1]*a[3,1];
dvdc4[20,1] = -6*(a[1,1]^2)*a[6,1] - 6*a[2,1]*a[3,1];
dvdc4[21,1] = -10*a[1,1]*a[6,1];

dvdc4[22,1] = 0;
dvdc4[23,1] = -10*a[5,1];
dvdc4[24,1] = 0;
dvdc4[25,1] = -(a[1,1]^5)*a[9,1] - 10*(a[1,1]^3)*a[5,1]*a[3,1]
              -10*(a[1,1]^2)*a[2,1]*a[8,1] - a[13,1];
dvdc4[26,1] = -3*(a[1,1]^4)*a[9,1] - 18*(a[1,1]^2)*a[5,1]*a[3,1]
              -12*a[1,1]*a[2,1]*a[8,1];
dvdc4[27,1] = -6*(a[1,1]^3)*a[9,1] - 18*a[1,1]*a[5,1]*a[3,1]
              -6*a[8,1]*a[2,1];
dvdc4[28,1] = -10*(a[1,1]^2)*a[9,1] - 10*a[3,1]*a[5,1];
dvdc4[29,1] = -15*a[1,1]*a[9,1];



dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];
dvdc5[6,1] = 0;
dvdc5[7,1] = 0;
dvdc5[8,1] = 0;
dvdc5[9,1] = -1;
dvdc5[10,1] = -(a[1,1]^3);
dvdc5[11,1] = -6*(a[1,1]^2)*a[3,1];
dvdc5[12,1] = -(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc5[13,1] = -3*(a[1,1]^2)*a[4,1] -a[3,1];
dvdc5[14,1] = -6*a[1,1]*a[4,1];
dvdc5[15,1] = 0;
dvdc5[16,1] = 0;
dvdc5[17,1] = -10*(a[1,1]^2)*a[8,1];
dvdc5[18,1] = -4*a[1,1]*a[8,1];
dvdc5[19,1] = -(a[1,1]^3)*a[7,1] - a[8,1];
dvdc5[20,1] = -4*(a[1,1]^2)*a[7,1];
dvdc5[21,1] = -10*a[1,1]*a[7,1];

dvdc5[22,1] = -10*(a[1,1]^3)*a[3,1];
dvdc5[23,1] = -10*a[4,1];
dvdc5[24,1] = -15*(a[1,1]^2)*a[10,1];
dvdc5[25,1] = -10*(a[1,1]^3)*a[3,1]*a[4,1] - 5*a[1,1]*a[10,1];
dvdc5[26,1] = -18*(a[1,1]^2)*a[3,1]*a[4,1] - a[10,1];
dvdc5[27,1] = -(a[1,1]^3)*a[11,1] - 18*a[1,1]*a[3,1]*a[4,1];
dvdc5[28,1] = -5*(a[1,1]^2)*a[11,1] - 10*a[3,1]*a[4,1];
dvdc5[29,1] = -15*a[1,1]*a[11,1];



dvdc6[1,1] = 0;
dvdc6[2,1] = 0;
dvdc6[3,1] = 0;
dvdc6[4,1] = 0;
dvdc6[5,1] = 0;
dvdc6[6,1] = -(a[1,1]^3);
dvdc6[7,1] = -(a[1,1]^2);
dvdc6[8,1] = -a[1,1];
dvdc6[9,1] = 0;
dvdc6[10,1] = 0;
dvdc6[11,1] = 0;
dvdc6[12,1] = 0;
dvdc6[13,1] = 0;
dvdc6[14,1] = 0;
dvdc6[15,1] = -(a[1,1]^4);
dvdc6[16,1] = -1;
dvdc6[17,1] = -10*(a[1,1]^3)*a[3,1];
dvdc6[18,1] = -(a[1,1]^4)*a[4,1] - 6*(a[1,1]^2)*a[3,1];
dvdc6[19,1] = -3*(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc6[20,1] = -6*(a[1,1]^2)*a[4,1] - a[3,1];
dvdc6[21,1] = -10*a[1,1]*a[4,1];

dvdc6[22,1] = 0;
dvdc6[23,1] = 0;
dvdc6[24,1] = -20*(a[1,1]^3)*a[8,1];
dvdc6[25,1] = -10*(a[1,1]^2)*a[8,1];
dvdc6[26,1] = -(a[1,1]^4)*a[7,1] - 4*a[1,1]*a[8,1];
dvdc6[27,1] = -4*(a[1,1]^3)*a[7,1] - a[8,1];
dvdc6[28,1] = -10*(a[1,1]^2)*a[7,1];
dvdc6[29,1] = -20*a[1,1]*a[7,1];



dvdc7[1,1] = 0;
dvdc7[2,1] = 0;
dvdc7[3,1] = 0;
dvdc7[4,1] = 0;
dvdc7[5,1] = 0;
dvdc7[6,1] = 0;
dvdc7[7,1] = 0;
dvdc7[8,1] = 0;
dvdc7[9,1] = -1;
dvdc7[10,1] = 0;
dvdc7[11,1] = 0;
dvdc7[12,1] = 0;
dvdc7[13,1] = -(a[1,1]^2)*a[2,1] -a[3,1];
dvdc7[14,1] = -4*a[1,1]*a[2,1];
dvdc7[15,1] = 0;
dvdc7[16,1] = 0;
dvdc7[17,1] = 0;
dvdc7[18,1] = 0;
dvdc7[19,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc7[20,1] = -4*(a[1,1]^2)*a[5,1];
dvdc7[21,1] = -10*a[1,1]*a[5,1];

dvdc7[22,1] = 0;
dvdc7[23,1] = -10*a[2,1];
dvdc7[24,1] = 0;
dvdc7[25,1] = 0;
dvdc7[26,1] = -(a[1,1]^4)*a[6,1] - 6*(a[1,1]^2)*a[2,1]*a[3,1] - a[10,1];
dvdc7[27,1] = -4*(a[1,1]^3)*a[6,1] - 12*a[1,1]*a[2,1]*a[3,1];
dvdc7[28,1] = -10*(a[1,1]^2)*a[6,1] - 10*a[3,1]*a[2,1];
dvdc7[29,1] = -20*a[1,1]*a[6,1];



dvdc8[1,1] = 0;
dvdc8[2,1] = 0;
dvdc8[3,1] = 0;
dvdc8[4,1] = 0;
dvdc8[5,1] = 0;
dvdc8[6,1] = 0;
dvdc8[7,1] = 0;
dvdc8[8,1] = 0;
dvdc8[9,1] = 0;
dvdc8[10,1] = -1;
dvdc8[11,1] = -4*a[1,1]*a[2,1];
dvdc8[12,1] = -( a[2,1] + a[4,1] );
dvdc8[13,1] = 0;
dvdc8[14,1] = 0;
dvdc8[15,1] = 0;
dvdc8[16,1] = 0;
dvdc8[17,1] = -10*(a[1,1]^2)*a[5,1];
dvdc8[18,1] = -4*a[1,1]*a[5,1];
dvdc8[19,1] = -( a[5,1] + a[7,1] );
dvdc8[20,1] = 0;
dvdc8[21,1] = 0;

dvdc8[22,1] = -10*(a[1,1]^2)*a[2,1];
dvdc8[23,1] = 0;
dvdc8[24,1] = -20*(a[1,1]^3)*a[6,1];
dvdc8[25,1] = -10*(a[1,1]^2)*(a[6,1] + a[2,1]*a[4,1]);
dvdc8[26,1] = -4*a[1,1]*(a[6,1] + 3*a[2,1]*a[4,1]);
dvdc8[27,1] = -(a[6,1] + 6*a[2,1]*a[4,1] + a[11,1]);
dvdc8[28,1] = 0;
dvdc8[29,1] = 0;



dvdc9[1,1] = 0;
dvdc9[2,1] = 0;
dvdc9[3,1] = 0;
dvdc9[4,1] = 0;
dvdc9[5,1] = 0;
dvdc9[6,1] = 0;
dvdc9[7,1] = 0;
dvdc9[8,1] = 0;
dvdc9[9,1] = 0;
dvdc9[10,1] = 0;
dvdc9[11,1] = -(a[1,1]^4);
dvdc9[12,1] = -(a[1,1]^3);
dvdc9[13,1] = -(a[1,1]^2);
dvdc9[14,1] = -a[1,1];
dvdc9[15,1] = 0;
dvdc9[16,1] = 0;
dvdc9[17,1] = 0;
dvdc9[18,1] = 0;
dvdc9[19,1] = 0;
dvdc9[20,1] = 0;
dvdc9[21,1] = 0;

dvdc9[22,1] = -(a[1,1]^5);
dvdc9[23,1] = -1;
dvdc9[24,1] = -15*(a[1,1]^4)*a[3,1];
dvdc9[25,1] = -(a[1,1]^5)*a[4,1] - 10*(a[1,1]^3)*a[3,1];
dvdc9[26,1] = -3*(a[1,1]^4)*a[4,1] - 6*(a[1,1]^2)*a[3,1];
dvdc9[27,1] = -6*(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc9[28,1] = -10*(a[1,1]^2)*a[4,1] - a[3,1];
dvdc9[29,1] = -15*a[1,1]*a[4,1];



dvdc10[1,1] = 0;
dvdc10[2,1] = 0;
dvdc10[3,1] = 0;
dvdc10[4,1] = 0;
dvdc10[5,1] = 0;
dvdc10[6,1] = 0;
dvdc10[7,1] = 0;
dvdc10[8,1] = 0;
dvdc10[9,1] = 0;
dvdc10[10,1] = 0;
dvdc10[11,1] = 0;
dvdc10[12,1] = 0;
dvdc10[13,1] = 0;
dvdc10[14,1] = 0;
dvdc10[15,1] = -1;
dvdc10[16,1] = 0;
dvdc10[17,1] = -5*a[1,1]*a[2,1];
dvdc10[18,1] = -( a[2,1] + a[4,1] );
dvdc10[19,1] = 0;
dvdc10[20,1] = 0;
dvdc10[21,1] = 0;

dvdc10[22,1] = 0;
dvdc10[23,1] = 0;
dvdc10[24,1] = -15*(a[1,1]^2)*a[5,1];
dvdc10[25,1] = -5*a[1,1]*a[5,1];
dvdc10[26,1] = -(a[5,1] + a[7,1]);
dvdc10[27,1] = 0;
dvdc10[28,1] = 0;
dvdc10[29,1] = 0;



dvdc11[1,1] = 0;
dvdc11[2,1] = 0;
dvdc11[3,1] = 0;
dvdc11[4,1] = 0;
dvdc11[5,1] = 0;
dvdc11[6,1] = 0;
dvdc11[7,1] = 0;
dvdc11[8,1] = 0;
dvdc11[9,1] = 0;
dvdc11[10,1] = 0;
dvdc11[11,1] = 0;
dvdc11[12,1] = 0;
dvdc11[13,1] = 0;
dvdc11[14,1] = 0;
dvdc11[15,1] = 0;
dvdc11[16,1] = -1;
dvdc11[17,1] = 0;
dvdc11[18,1] = 0;
dvdc11[19,1] = 0;
dvdc11[20,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc11[21,1] = -5*a[1,1]*a[2,1];

dvdc11[22,1] = 0;
dvdc11[23,1] = 0;
dvdc11[24,1] = 0;
dvdc11[25,1] = 0;
dvdc11[26,1] = 0;
dvdc11[27,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc11[28,1] = -5*(a[1,1]^2)*a[5,1];
dvdc11[29,1] = -15*a[1,1]*a[5,1];



dvdc12[1,1] = 0;
dvdc12[2,1] = 0;
dvdc12[3,1] = 0;
dvdc12[4,1] = 0;
dvdc12[5,1] = 0;
dvdc12[6,1] = 0;
dvdc12[7,1] = 0;
dvdc12[8,1] = 0;
dvdc12[9,1] = 0;
dvdc12[10,1] = 0;
dvdc12[11,1] = 0;
dvdc12[12,1] = 0;
dvdc12[13,1] = 0;
dvdc12[14,1] = 0;
dvdc12[15,1] = 0;
dvdc12[16,1] = 0;
dvdc12[17,1] = -(a[1,1]^5);
dvdc12[18,1] = -(a[1,1]^4);
dvdc12[19,1] = -(a[1,1]^3);
dvdc12[20,1] = -(a[1,1]^2);
dvdc12[21,1] = -a[1,1];

dvdc12[22,1] = 0;
dvdc12[23,1] = 0;
dvdc12[24,1] = 0;
dvdc12[25,1] = 0;
dvdc12[26,1] = 0;
dvdc12[27,1] = 0;
dvdc12[28,1] = 0;
dvdc12[29,1] = 0;



dvdc13[1,1] = 0;
dvdc13[2,1] = 0;
dvdc13[3,1] = 0;
dvdc13[4,1] = 0;
dvdc13[5,1] = 0;
dvdc13[6,1] = 0;
dvdc13[7,1] = 0;
dvdc13[8,1] = 0;
dvdc13[9,1] = 0;
dvdc13[10,1] = 0;
dvdc13[11,1] = 0;
dvdc13[12,1] = 0;
dvdc13[13,1] = 0;
dvdc13[14,1] = 0;
dvdc13[15,1] = 0;
dvdc13[16,1] = 0;
dvdc13[17,1] = 0;
dvdc13[18,1] = 0;
dvdc13[19,1] = 0;
dvdc13[20,1] = 0;
dvdc13[21,1] = 0;
dvdc13[22,1] = -1;
dvdc13[23,1] = 0;
dvdc13[24,1] = -6*a[1,1]*a[2,1];
dvdc13[25,1] = -(a[2,1] + a[4,1]);
dvdc13[26,1] = 0;
dvdc13[27,1] = 0;
dvdc13[28,1] = 0;
dvdc13[29,1] = 0;



dvdc14[1,1] = 0;
dvdc14[2,1] = 0;
dvdc14[3,1] = 0;
dvdc14[4,1] = 0;
dvdc14[5,1] = 0;
dvdc14[6,1] = 0;
dvdc14[7,1] = 0;
dvdc14[8,1] = 0;
dvdc14[9,1] = 0;
dvdc14[10,1] = 0;
dvdc14[11,1] = 0;
dvdc14[12,1] = 0;
dvdc14[13,1] = 0;
dvdc14[14,1] = 0;
dvdc14[15,1] = 0;
dvdc14[16,1] = 0;
dvdc14[17,1] = 0;
dvdc14[18,1] = 0;
dvdc14[19,1] = 0;
dvdc14[20,1] = 0;
dvdc14[21,1] = 0;
dvdc14[22,1] = 0;
dvdc14[23,1] = -1;
dvdc14[24,1] = 0;
dvdc14[25,1] = 0;
dvdc14[26,1] = 0;
dvdc14[27,1] = 0;
dvdc14[28,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc14[29,1] = -6*a[1,1]*a[2,1];



dvdc15[1,1] = 0;
dvdc15[2,1] = 0;
dvdc15[3,1] = 0;
dvdc15[4,1] = 0;
dvdc15[5,1] = 0;
dvdc15[6,1] = 0;
dvdc15[7,1] = 0;
dvdc15[8,1] = 0;
dvdc15[9,1] = 0;
dvdc15[10,1] = 0;
dvdc15[11,1] = 0;
dvdc15[12,1] = 0;
dvdc15[13,1] = 0;
dvdc15[14,1] = 0;
dvdc15[15,1] = 0;
dvdc15[16,1] = 0;
dvdc15[17,1] = 0;
dvdc15[18,1] = 0;
dvdc15[19,1] = 0;
dvdc15[20,1] = 0;
dvdc15[21,1] = 0;
dvdc15[22,1] = 0;
dvdc15[23,1] = 0;
dvdc15[24,1] = -(a[1,1]^6);
dvdc15[25,1] = -(a[1,1]^5);
dvdc15[26,1] = -(a[1,1]^4);
dvdc15[27,1] = -(a[1,1]^3);
dvdc15[28,1] = -(a[1,1]^2);
dvdc15[29,1] = -a[1,1];

retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5~dvdc6~dvdc7~dvdc8~dvdc9
          ~dvdc10~dvdc11~dvdc12~dvdc13~dvdc14~dvdc15);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5, dvdc6, dvdc7, dvdc8, dvdc9,
      dvdc10, dvdc11, dvdc12, dvdc13, dvdc14, dvdc15;

endp;

@  ------------------------------------------------------- @
@ Subroutine to compute squeezes.

  This subroutine compares the values of the objective function at the
  points c+s*dc and c+0.5*s*dc, with dc = the proposed change in the vector
  of parameters, and with step length s initially at 1. s is halved until
  minus the objective function stops declining.
@
proc (2) = squeeze(s_c,s_dc,w,obj);
local s_c1,s_lm,s_itr,lc1,s_c2,lc2,s_f1,s_f2;


    s_c1=s_c - s_dc; s_lm=1/2; s_itr=1;
    s_f1 = deff(s_c1,Emom,0,effsave);
    lc1 = s_f1'w*s_f1;
    clear s_f1;
  do until s_itr > maxsqez;

    s_c2=s_c-s_lm*s_dc;
    s_f2 = deff(s_c2,Emom,0,effsave);
    lc2 = s_f2'w*s_f2;
    clear s_f2;

    if lc1 <= lc2 and lc1 <= obj;

       retp(s_c1,s_itr-1); goto eoproc;

    else;

       s_c1=s_c2; s_lm=s_lm/2; lc1=lc2;
       s_itr=s_itr+1;

    endif;

  endo;
retp(s_c2,s_itr-1);

eoproc:
endp;

@----------Subroutine to trap inversion errors-----------------@

proc invit(w);
local onemore,ww;
onemore = 1;
winv:
trap 1;
  ww = invpd(w);
trap 0;
  if scalerr(ww) > 0;
      if onemore == 1;
        call sysstate(14,1e-128);
        onemore = 0;
        goto winv;
      else;
        ww = -999999;
      endif;
  endif;

retp(ww);

endp;

@----------Subroutine to trap solving errors-----------------@

proc solvit(a,b);
local dc,onemore;
    onemore = 1;
    solveit:
    trap 1;
    dc = solpd(a,b);
    trap 0;
    if scalerr(dc) > 0;
      if onemore == 1;
        call sysstate(14,1e-128);
        onemore = 0;
        goto solveit;
      else;
        dc = -999999;
      endif;
    endif;

retp(dc);
endp;

proc (2) = gmm(c);

local ff, w, iter, dc, f, g, obj, gwg, gwf, c_new, sqz;


         @ --------The program creates the weighting matrix. ------------@

      ff = optw(c,mom,Emom,bleh);
      w = invit(moment(ff,0)./n);
      if w == -999999 and bb > 0; obj = -999999; goto escape; endif;
         @ ------------- Start of iteration loop ------------------------@
    iter = 1;
        dc=1;                     @ Initialize the step length. @
    do until abs(dc) < tol;

        f = deff(c,Emom,0,effsave);

        if estim == 2;
        g=grad2(c);                @ The program jumps to the subroutine that
                                    computes analytic derivatives. The matrix
                                    of partials of f with respect to the
                                    parameters is called g. @
        elseif estim == 3;
        g=grad3(c);
        elseif estim == 4;
        g=grad4(c);
        elseif estim == 5;
        g=grad5(c);
        endif;

        obj = f'w*f;              @ This computes the value of the objective
                                    function.@

        gwg= g'w*g;               @ This uses the GAUSS-NEWTON method
                                    to compute full step dc. @
        gwf = g'w*f;

        dc = solvit(gwf,gwg);
        if dc == -999999 and bb > 0; obj = -999999; goto escape; endif;

     if maxsqez > 0;              @ This jumps to the subroutine that
                                    adjusts the step length. @
       { c_new,sqz } = squeeze(c,dc,w,obj);
     else;
       c_new=c - dc;
     endif;


     dc=c_new-c;                     @ Update variables for the next iteration. @
     c=c_new;
    /*
    /*--------------Print out the results for the current iteration------------*/
    cls;
    "Number of Squeezes = ";; sqz;
    "Objective Function = ";; obj*n;

    "            Value         Step";
     format 14,6;
     mm=1;
     do until mm > nc;
      $"  ";; c[mm,1];; "  ";; dc[mm,1];
      mm=mm+1;
      endo;
    /*-----------------End of print out of current iteration--------------*/
    */
     iter=iter +1;
     if iter >= maxiter;  @ Quit iterating if necessary. @
       goto escape;
     endif;
    endo;
    @ ------ End of iteration loop ------------------------------------ @
    escape:

retp(c, obj);
endp;
@------------------------------------------------------------------------@

eop:

"Start time:  ";; timestr(ttt);
"End time:    ";; timestr(0);
"Start date:  ";; datestr(ddd);
"End date:    ";; datestr(0);

output off;
end;
