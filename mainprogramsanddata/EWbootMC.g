new ,150000;
ttt = time;
ddd = date;
/* Monte Carlo program to evaluate the 3rd--7th moment partialling
   estimators with the unrestricted weighting matrix, and the
   correct standard errors for the coefficient on the correctly measured
   regressor.

   Written by Tim Erickson and Toni Whited.

*/

load params[] = params.txt;

dofe       =   params[1,1]; @1 for fixed effects, 0 for levels@
nestim     =   params[2,1]; @1 for third moment estimstor, 2 for fourth moment estimator, 3 for fifth moment estimator, 4 for the sixth moment estimator, etc.@
samplesize =   params[3,1]; @1 for 1500, 2 for 4000@
kapital    =   params[4,1]; @1 for the capital stock DGP, 2 for the assets DGP @
trialno    =   params[5,1]; @we need something to reset the seed for each trial!!!!@



startestim = nestim;
format /rd 8,3;


outname = "bootmc";

    fmat="%lf";
    cnestim = nestim + 2;
    ce=ftos(cnestim,fmat,1,0);
    cf=ftos(dofe,fmat,1,0);
    ca=ftos(kapital,fmat,1,0);
    ct=ftos(trialno,fmat,1,0);
    outname = outname $+ "_FE" $+ cf $+ "_ES" $+ ce $+ "_KA" $+ ca $+ "_TR" $+ ct $+ ".txt";

    output file = ^outname reset;


outwidth 256;
rezat = 0;

      if samplesize == 1;       ncount = 1500;
      elseif samplesize == 2;   ncount = 4000;
      endif;
      @
      "----------------------------------------------------------------------------------------------------------------------------";
      "                                                sample size = ";; ncount;
      "----------------------------------------------------------------------------------------------------------------------------";
      @
      oid = 1;
      diffiv = 1;
      checkmoms = 0;
      mcount = 1;
      nboot  = 1000;
      betastart = 0;


if kapital == 1;


       zrho = 0.46;
       qrho = 0.72;
       erho = 0;
       urho = 0;

       /* bleh = 1 acts as if Ey and Ex known.  bleh = 2 optimally accounts for the
          fact that they are not known.
       */
       bleh = 2;
       nyr = 20;
       keepnyr = 10;

       ndum = 0;
       nz = ndum+1;
       beta = 0.02;
       acf = 0.05;
       delta = acf;
       coeffs = (beta|acf);

       Eofx = 2.279;
       Eofy = 0.129;
       EofZ = 0.176;

       let covmx[3,3] =
             0.0137       0.1573       0.0114
             0.1573      17.5878       0.4082
             0.0114       0.4082       0.1134;
       varY = covmx[1,1];
       varX = covmx[2,2];
       varZ = covmx[3,3];


       tau2 = 0.5;

else;

       zrho = 0.40;
       qrho = 0.67;
       erho = 0;
       urho = 0;

       /* bleh = 1 acts as if Ey and Ex known.  bleh = 2 optimally accounts for the
          fact that they are not known.
       */
       bleh = 2;
       nyr = 20;
       keepnyr = 10;

       ndum = 0;
       nz = ndum+1;
       beta = 0.03;
       acf = 0.08;
       delta = acf;
       coeffs = (beta|acf);

       Eofx = 1.492;
       Eofy = 0.068;
       EofZ = 0.087;

       let covmx[3,3] =
             0.0034       0.0089       0.0019
             0.0089       1.1191       0.0226
             0.0019       0.0226       0.0109;
       varY = covmx[1,1];
       varX = covmx[2,2];
       varZ = covmx[3,3];


       tau2 = 0.25;



endif;


Varchi = varX*tau2;
covtrumx = covmx[2:3,2:3];

covtrumx[1,1] = Varchi;

Vare   = varX - Varchi;
Varu   = VarY - coeffs'covtrumx*coeffs;
Sdchi  = sqrt(Varchi);
SDe    = sqrt(Vare);
SDu    = sqrt(Varu);


corrtrue = covtrumx./sqrt(diag(covtrumx))./sqrt(diag(covtrumx)');
{vLmx,Qmx} = eigrs2(corrtrue);
rtLmx = sqrt(diagrv(zeros(2,2),vLmx));
rtcov = Qmx*rtLmx*Qmx';


alpha0 = Eofy - beta'Eofx - delta'Eofz;
rho2 = 1 - (Sdu^2)/varY;
icept = ones(ncount,1);


      @--Print out parameter setting.--@
   @
      "beta   = " beta   ;
      "alphas = " acf'   ;
      "tau2   = " tau2   ;
      "rho2   = " rho2   ;
   @

      output off;
      cls;
      maxiter = 400;
      maxsqez = 40;

      SDtim = sqrt(96);
      rt3   = sqrt(3);
      rt5   = sqrt(5);
      rt17  = sqrt(17);
      rt50  = sqrt(50);
      rt2   = sqrt(2);
      rt5   = sqrt(5);
      rtp5  = sqrt(0.5);
      rt1p5 = sqrt(1.5);
      mln   = exp(0.5);
      sln   = sqrt(exp(1)*(exp(1)-1));

      @--Initialize the counters--@

      avAy = 0;
      avVy = 0;
      avSkewy = 0;
      avKurty = 0;
      avstd5y = 0;
      avstd6y = 0;
      avstd7y = 0;
      avAx = 0;
      avVx = 0;
      avSkewx = 0;
      avKurtx = 0;
      avstd5x = 0;
      avstd6x = 0;
      avstd7x = 0;
      avAchi = 0;
      avVchi = 0;
      avSkewchi = 0;
      avKurtchi = 0;
      avstd5chi = 0;
      avstd6chi = 0;
      avstd7chi = 0;
      avAz = 0;
      avVz = 0;
      avSkewz = 0;
      avKurtz = 0;
      avstd5z = 0;
      avstd6z = 0;
      avstd7z = 0;

      avAu = 0;
      avVu = 0;
      avSkewu = 0;
      avKurtu = 0;
      avstd5u = 0;
      avstd6u = 0;
      avstd7u = 0;

      meAy = 0;
      meVy = 0;
      meSkewy = 0;
      meKurty = 0;
      mestd5y = 0;
      mestd6y = 0;
      mestd7y = 0;
      meAx = 0;
      meVx = 0;
      meSkewx = 0;
      meKurtx = 0;
      mestd5x = 0;
      mestd6x = 0;
      mestd7x = 0;

      averhoy = 0;
      averhox = 0;
      averhoz = 0;


      aveEy2x = 0;
      aveEyx2 = 0;


      @-------------------Parameters-----------------------@


      @OLS@

      isave   = zeros(1,keepnyr);         @ Intercept                                     @
      bsave   = zeros(1,keepnyr);         @ Coefficient on chi                            @
      rsave   = zeros(1,keepnyr);         @ R2                                            @
      zsave   = zeros(nz,keepnyr);        @ Coefficients on perfectly measured regressors @
      asave   = zeros(1,keepnyr);         @ Stupid pointless thing from a long time ago.  @
      rtruesave = zeros(1,keepnyr);
      ttruesave = zeros(1,keepnyr);

      @GMM@

      csave   = zeros(1,keepnyr*nestim);    @ Coefficient on chi                            @
      dsave   = zeros(nz,keepnyr*nestim);   @ Coefficients on perfectly measured regressors @
      tsave   = zeros(1,keepnyr*nestim);    @ tau2                                          @
      psave   = zeros(1,keepnyr*nestim);    @ rho2                                          @

      @IVAB@

      bivsave = zeros(2,1);
      zivsave = zeros(nz*2,1);
      babsave = zeros(2,1);
      zabsave = zeros(nz*2,1);

      @-------------------Influence Functions--------------@

      @GMM@

      icsave  = zeros(ncount,keepnyr*nestim);
      idsave  = zeros(ncount*nz,keepnyr*nestim);
      ipsave  = zeros(ncount,keepnyr*nestim);
      itsave  = zeros(ncount,keepnyr*nestim);

      @OLS@

      iisave  = zeros(ncount,keepnyr);
      ibsave  = zeros(ncount,keepnyr);
      irsave  = zeros(ncount,keepnyr);
      izsave  = zeros(ncount*nz,keepnyr);
      iasave  = zeros(ncount,keepnyr);


      @------------------CMD Savers--------------------------@

      @OLS@

      bcmd = zeros(3,1); @ "1" is the coefficient. "2" is the standard error. "3" is the oid test.@
      rcmd = zeros(3,1);
      acmd = zeros(3,1);
      zcmd = zeros(3*nz,1);

      @GMM@

      ccmd = zeros(3,nestim);
      dcmd = zeros(3*nz,nestim);
      scmd = zeros(3,nestim);
      pcmd = zeros(3,nestim);
      tcmd = zeros(3,nestim);

      @------------------Monte Carlo Counters----------------@

      @OLS@

      sumb      = 0;
      sumerr2b  = 0;
      sumerrb   = 0;
      sumseb    = 0;
      pcountb   = zeros(1,1);
      bkeep     = zeros(mcount,1);

      sumz     =  zeros(nz,1);
      sumerr2z =  zeros(nz,1);
      sumerrz  =  zeros(nz,1);
      sumsez   =  zeros(nz,1);
      pcountz  =  zeros(1*nz,1);
      zkeep     = zeros(mcount,nz);

      suma      = 0;
      sumerr2a  = 0;
      sumerra   = 0;
      sumsea    = 0;
      pcounta   = zeros(1,1);

      sumARR    = 0;
      sumerr2R  = 0;
      sumerrR   = 0;
      sumseR    = 0;
      pcountR   = zeros(1,1);
      Rkeep     = zeros(mcount,1);


      sumtrueR = 0;
      sumtrueT = 0;

      @IV@

      biv         = zeros(10,1);
      sumbiv      = 0;
      sumerr2biv  = 0;
      sumerrbiv   = 0;
      sumsebiv    = 0;
      pcountbiv   = zeros(1,1);
      bivkeep     = zeros(mcount,1);

      sumziv     =  zeros(nz,1);
      sumerr2ziv =  zeros(nz,1);
      sumerrziv  =  zeros(nz,1);
      sumseziv   =  zeros(nz,1);
      pcountziv  =  zeros(1*nz,1);
      zivkeep    = zeros(mcount,1);

      sumARRiv    = 0;
      sumerr2Riv  = 0;
      sumerrRiv   = 0;
      sumseRiv    = 0;
      pcountRiv   = zeros(1,1);
      Rivkeep     = zeros(mcount,1);

      @AB@

      bnew        = zeros(10,1);
      sumbab      = 0;
      sumerr2bab  = 0;
      sumerrbab   = 0;
      sumsebab    = 0;
      pcountbab   = zeros(1,1);
      babkeep     = zeros(mcount,1);

      sumzab     =  zeros(nz,1);
      sumerr2zab =  zeros(nz,1);
      sumerrzab  =  zeros(nz,1);
      sumsezab   =  zeros(nz,1);
      pcountzab  =  zeros(1*nz,1);
      zabkeep     = zeros(mcount,nz);

      sumARRab    = 0;
      sumerr2Rab  = 0;
      sumerrRab   = 0;
      sumseRab    = 0;
      pcountRab   = zeros(1,1);
      Rabkeep     = zeros(mcount,1);



      @GMM@

      sumg = zeros(1,nestim);
      sumerr2g = zeros(1,nestim);
      sumerrg = zeros(1,nestim);
      sumseg = zeros(1,nestim);
      pcountg = zeros(1,nestim);
      gkeep   = zeros(mcount,nestim);

      sumd =      zeros(nz,nestim);
      sumerr2d =  zeros(nz,nestim);
      sumerrd =   zeros(nz,nestim);
      sumsed  =   zeros(nz,nestim);
      pcountd =   zeros(nz,nestim);
      dkeep   = zeros(mcount,nestim*nz);

      sump =  zeros(1,nestim);
      sumerr2p = zeros(1,nestim);
      sumerrp = zeros(1,nestim);
      sumsep = zeros(1,nestim);
      pcountp = zeros(1,nestim);
      pkeep   = zeros(mcount,nestim);

      sumt = zeros(1,nestim);
      sumerr2t = zeros(1,nestim);
      sumerrt = zeros(1,nestim);
      sumset = zeros(1,nestim);
      pcountt = zeros(1,nestim);
      tkeep   = zeros(mcount,nestim);


      @----------------ALL of the overidentification test counters.--------------------@
      hansenj  = zeros(1,nestim);
      EWjcount = zeros(1,nestim);
      IVjcount = 0;
      ABjcount = 0;


      @-----------------End of counters.----------------------------------------------------@


      numwh = 0;

      acritval = 1.645;
      twocrit  = 1.960;

      let oidcrit =
       0.000000
       5.987470
      11.065330
      16.912750
      23.677540;

      teedelta = zeros(nestim,2);
      teebeta  = zeros(nestim,2);
      teetau   = zeros(nestim,2);
      teerho   = zeros(nestim,2);




      bootc    = zeros(nestim,1);
      bootd    = zeros(nestim,nz);
      bootp    = zeros(nestim,1);
      boott    = zeros(nestim,1);
      bootab   = zeros(2,1);
      bootiv   = zeros(2,1);



      /* Make the seeds for the dgp. */
      sb  = 123456;      @NEVER CHANGE THIS SEED!@

      let allseed = 0 1 2 3 4 5 6 7 8 9 10 11 12;
      s1  = allseed[1,1]  + trialno*100;
      s2  = allseed[2,1]  + trialno*100;
      s3  = allseed[3,1]  + trialno*100;
      s4  = allseed[4,1]  + trialno*100;
      s5  = allseed[5,1]  + trialno*100;
      s10 = allseed[6,1]  + trialno*100;
      s11 = allseed[7,1]  + trialno*100;
      s12 = allseed[8,1]  + trialno*100;
      s23 = allseed[9,1]  + trialno*100;
      t21 = allseed[10,1] + trialno*100;
      t22 = allseed[11,1] + trialno*100;
      t23 = allseed[12,1] + trialno*100;
      t24 = allseed[13,1] + trialno*100;


      if dofe==1;
      icount = 0;
      else;
      icount = 1;
      endif;
      enditnow = 0;
      nrun = 0;
      do while icount <= mcount and enditnow == 0;

          if key==116; enditnow = 1; endif;

          bttdb:

          @------------This part makes the "years" of data.------------------@
          @ original@


          if kapital == 1;

          agamE = 0.023;
          bgamE = 1;

          agamR = 0.027;
          bgamR = 1;

          agamU = 0.25;
          bgamU = 1;

          agamZ = 0.8;
          bgamZ = 1;

          else;

          agamE = 0.2;
          bgamE = 1;

          agamR = 0.025;
          bgamR = 1;

          agamU = 0.45;
          bgamU = 1;

          agamZ = 0.45;
          bgamZ = 1;


          endif;



          if icount == 0;
             nncount = 10000000;
          else;
             nncount = ncount;
          endif;


          creg    = zeros(nncount*nyr,1);
          zreg    = zeros(nncount*nyr,1);
          ek      = zeros(nncount*nyr,1);
          uk      = zeros(nncount*nyr,1);


          { einnov, s4 } = rndLCgam(nncount*nyr, 1, agamE, s4);      einnov = bgamE*einnov;              einnov =  (einnov   - agamE*bgamE)/sqrt(agamE*bgamE^2);
          { uinnov, s5 } = rndLCgam(nncount*nyr, 1, agamU, s5);      uinnov = bgamU*uinnov;              uinnov =  (uinnov   - agamU*bgamU)/sqrt(agamU*bgamU^2);
          { rinnov, s1 } = rndLCgam(nncount*nyr, 1, agamR, s1);      rinnov = bgamR*rinnov;              rinnov =  (rinnov -   agamR*bgamR)/sqrt(agamR*bgamR^2);
          { zinnov, s3 } = rndLCgam(nncount*nyr, 1, agamZ, s3);      zinnov = bgamZ*zinnov;              zinnov = -(zinnov   - agamZ*bgamZ)/sqrt(agamZ*bgamZ^2);

          idx = seqa(1,nyr,nncount);                                @        This just indexes the different firms in the matrix. They are stacked on top of each other.    @
          { ek[idx,1],   s12   } = rndLCgam(nncount, 1, agamE,s12);    ek[idx,1]   = bgamE*ek[idx,1];      ek[idx,1]   =    (ek[idx,1]   - agamE*bgamE)/sqrt(agamE*bgamE^2);
          { uk[idx,1],   s11   } = rndLCgam(nncount, 1, agamU,s11);    uk[idx,1]   = bgamU*uk[idx,1];      uk[idx,1]   =    (uk[idx,1]   - agamU*bgamU)/sqrt(agamU*bgamU^2);
          { creg[idx,1], s10 }   = rndLCgam(nncount, 1, agamR,s10);    creg[idx,1] = bgamR*creg[idx,1];    creg[idx,1] =  (creg[idx,1]   - agamR*bgamR)/sqrt(agamR*bgamR^2);
          { zreg[idx,1], s2  }   = rndLCgam(nncount, 1, agamZ,s2 );    zreg[idx,1] = bgamZ*zreg[idx,1];    zreg[idx,1] = -(zreg[idx,1]   - agamZ*bgamZ)/sqrt(agamZ*bgamZ^2);
          yr = zeros(nncount*nyr,1);
          yr[idx,1] = ones(nncount,1);


           for ttt(0,nyr-2,1);

              yr[idx+ttt+1,1] = yr[idx+ttt,1] + 1;
              zreg[idx+ttt+1,1]     = zreg[idx+ttt,1]*zrho + zinnov[idx+ttt+1,1];
              creg[idx+ttt+1,1]     = creg[idx+ttt,1]*qrho + rinnov[idx+ttt+1,1];
              uk[idx+ttt+1,1] =         uk[idx+ttt,1]*urho + uinnov[idx+ttt+1,1];
              ek[idx+ttt+1,1] =         ek[idx+ttt,1]*erho + einnov[idx+ttt+1,1];

           endfor;

           /* This part truncates */

          if keepnyr>1 or nyr > 1;
              creg = reshape(creg,nncount,nyr); creg = creg[.,nyr-keepnyr+1:nyr]; creg = vec(creg');
              zreg = reshape(zreg,nncount,nyr); zreg = zreg[.,nyr-keepnyr+1:nyr]; zreg = vec(zreg');
              uk = reshape(uk,nncount,nyr); uk = uk[.,nyr-keepnyr+1:nyr];         uk = vec(uk');
              ek = reshape(ek,nncount,nyr); ek = ek[.,nyr-keepnyr+1:nyr];         ek = vec(ek');
              yr = reshape(yr,nncount,nyr); yr = yr[.,nyr-keepnyr+1:nyr];         yr = vec(yr');  yr = yr - nyr + keepnyr;
          endif;
          creg = creg*sqrt(1-qrho^2);
          zreg   = zreg*sqrt(1-zrho^2);
          czreg = (creg~zreg)*rtcov;
          creg = czreg[.,1];
          zreg = czreg[.,2];


          uk     = uk*sqrt(1-urho^2);
          ek     = ek*sqrt(1-erho^2);
          @ icount;; stdc(uk)^2;;   @


          allc     = creg*Sdchi + EofX;
          u        = uk*Sdu;
          epsilon  = ek*Sde;
          allz     = EofZ + sqrt(varz)*zreg;

          ally = alpha0 + beta*allc + delta*allz + u;
          allx = allc + epsilon;


          if icount == 0;
             @---------------------------This part makes true rho2 and tau2 for the within stuff------------------------------------@

               if dofe == 1;

                    allx = reshape(allx,nncount,keepnyr); allx = allx - meanc(allx'); allx = vec(allx');
                    ally = reshape(ally,nncount,keepnyr); ally = ally - meanc(ally'); ally = vec(ally');
                    allz = reshape(allz,nncount,keepnyr); allz = allz - meanc(allz'); allz = vec(allz');

                    idx = seqa(1,keepnyr,nncount);
                    year = 1;
                    do while year <= keepnyr;

                    y = ally[idx+year-1,1];
                    x = allx[idx+year-1,1];
                    z = allz[idx+year-1,1];
                    chi = allc[idx+year-1,1];

                    za = z;

                    /*----------------------------------OLS Regressions using chi-------------------------*/


                    des = @ones(rows(chi),1)~@chi~za;
                    dd = invpd(moment(des,0));
                    bbbb = dd*des'y;
                    uhat = y - des*bbbb;
                    rtruesave[1,year] = 1 - moment(uhat,0)/(moment((y - meanc(y)),0));

                    des = @ones(rows(chi),1)~@chi;
                    dd = invpd(moment(des,0));
                    bbbb = dd*des'x;
                    uhat = x - des*bbbb;
                    ttruesave[1,year] = 1 - moment(uhat,0)/(moment((x - meanc(x)),0));



                    year = year + 1;
                    endo;

                    rho2 = meanc(rtruesave');   "Within rho-squared = ";; rho2;
                    tau2 = meanc(ttruesave');   "Within tau-squared = ";; tau2;
                endif;

            else;



          /*===========================================This is where the bootstrap loop goes!========================================*/
          @-----------------The saver for the "population" GMM objective function.-----------------------@
          effsave = zeros(29,nyr*nestim); @ 29 is the number of moment conditions for GMM7@

          nzab = 4;
          abfsave = zeros(nzab*(keepnyr-2),1);
          nziv = 6;
          ivfsave = zeros(nziv,1);

          @-----------------The savers for the pivotal statistics.-----------------------@
          cpivot = zeros(nestim,nboot+(nboot==0));
          apivot = zeros(nestim,nboot+(nboot==0));
          dpivot = zeros(nz*nestim,nboot+(nboot==0));
          ppivot = zeros(nestim,nboot+(nboot==0));
          tpivot = zeros(nestim,nboot+(nboot==0));

          bpivot = zeros(1,nboot+(nboot==0));
          zpivot = zeros(nz,nboot+(nboot==0));
          rpivot = zeros(1,nboot+(nboot==0));
          ipivot = zeros(1,nboot+(nboot==0));

          abpivot = zeros(2,nboot+(nboot==0));
          ivpivot = zeros(2,nboot+(nboot==0));

          saveallx = reshape(allx,nncount,keepnyr);
          saveally = reshape(ally,nncount,keepnyr);
          saveallz = reshape(allz,nncount,keepnyr);

          rndseed 12345;


         pickidx = seqa(1,1,keepnyr);
         bb=0;
          do while bb<=nboot;

             bttfdb:
             if bb==0;
                allx = saveallx;
                ally = saveally;
                allz = saveallz;

             else;
                pickboot = ceil(rndus(nncount,1,sb)*nncount);
                pickboot = maxc(ones(1,rows(pickboot))|(pickboot'));
                allx = saveallx[pickboot,.];
                ally = saveally[pickboot,.];
                allz = saveallz[pickboot,.];

             endif;

             if dofe == 1;

                  allx = allx - meanc(allx');
                  ally = ally - meanc(ally');
                  allz = allz - meanc(allz');

             endif;
             allx = vec(allx');
             ally = vec(ally');
             allz = vec(allz');


             clear allb, rb, rd;

             if checkmoms == 0;
                 idx = seqa(1,keepnyr,nncount);                                           @This just indexes the different firms in the matrix. They are stacked on top of each other. @
                 year = 1;
                 do while year <= keepnyr;

                     y = ally[idx+year-1,1];
                     x = allx[idx+year-1,1];
                     z = allz[idx+year-1,1];
                     chi = allc[idx+year-1,1];

                     za = z;
                     iz = icept~z;

                     mux = (invpd(moment(iz,0)))* (iz'x);
                     muy = (invpd(moment(iz,0)))* (iz'y);


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


                     Ezx  = iz'x/nncount;

                     /* Make the moments for the rest of the summary statistics. */

                     Ey  = meanc(y);
                     Ex  = meanc(x);
                     Ez  = meanc(z);
                     Ey2 = meanc(y^2);
                     Ex2 = meanc(x^2);
                     Ez2 = meanc(z^2);


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

                     Ezz     =   moment(iz,0)/nncount;
                     inEzz   =   invpd(Ezz);
                     zy_d    =   iz.*y_d;
                     zx_d    =   iz.*x_d;

                     Eza   = meanc(za)';



                     /* The OLS estimator: */

                     des = icept~x~za;
                     dd = invpd(moment(des,0));
                     allb = dd*des'y;
                     uhat = y - des*allb;
                     wh = uhat.*des;
                     seb = sqrt(diag(dd*(wh'wh)*dd));
                     R2 = 1 - moment(uhat,0)/(moment((y - Ey),0));
                     olsinflnc = invpd(moment(des,0)/nncount)*((des.*uhat)');

                     bsave[1,year] = allb[2,1];
                     asave[1,year] = sumc(allb[3:3+ndum/2,1]);
                     for qq(1, nz, 1); zsave[qq,year] = allb[2+qq,1]; endfor;

                     ibsave[.,year] = olsinflnc[2,.]';
                     iasave[.,year] = sumc(olsinflnc[3:3+ndum/2,.]);
                     for qq(1, nz, 1); izsave[((qq-1)*nncount+1):(qq*nncount),year] = olsinflnc[2+qq,.]'; endfor;

                     rsave[1,year] = 1 - moment(uhat,0)/(moment((y - Ey),0));

                     des = icept~y~za;
                     dd = invpd(moment(des,0));
                     allbr = dd*des'x;
                     rb  = 1/allbr[2,1];
                     rd  = -allbr[3,1]/allbr[2,1];
                     @----------First make the influence function for sigma_xz.----------------------@
                     xza = x~za;
                     Exza = meanc(xza);
                     nreg=nz+1;
                     sigxz = moment((xza)-Exza',0)/nncount;


                     vecsigxz = zeros(nreg+nreg*(nreg-1)/2,1);
                     vecsigxz[1:nreg,1] = diag(sigxz);

                     counter =  nreg+nreg*(nreg-1)/2;
                     for qq(nreg, 2, -1);
                       for ee(qq-1, 1, -1);
                         vecsigxz[counter,1] = sigxz[ee,qq];
                         counter=counter-1;
                       endfor;
                     endfor;


                     phixz = zeros(nncount,nreg+nreg*(nreg-1)/2);
                     phixz[.,1:nreg] = ((xza) - Exza').*((xza) - Exza');

                     counter =  nreg+nreg*(nreg-1)/2;
                     for qq(nreg, 2, -1);
                       for ee(qq-1, 1, -1);
                         phixz[.,counter] = (xza[.,ee] - Exza[ee,1]).*(xza[.,qq] - Exza[qq,1]);
                         counter=counter-1;
                       endfor;
                     endfor;

                     sigy = moment((y - Ey),0)/nncount;
                     phiy = (y - Ey)^2 - sigy;

                     bigphi = (olsinflnc[2:rows(olsinflnc),.])|(phixz')|(phiy');

                     @--------------Now make the derivative matrix.--------------------------------------@

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

                     irsave[.,year] = -bigphi'gee;

                     /*----------------------------------OLS Regressions using chi-------------------------*/

                     if year == 1;

                     des = icept~chi~za;
                     dd = invpd(moment(des,0));
                     bbbb = dd*des'y;
                     uhat = y - des*bbbb;
                     rtruesave[1,year] = 1 - moment(uhat,0)/(moment((y - Ey),0));

                     des = icept~chi;
                     dd = invpd(moment(des,0));
                     bbbb = dd*des'x;
                     uhat = x - des*bbbb;
                     ttruesave[1,year] = 1 - moment(uhat,0)/(moment((x - Ex),0));

                     endif;

                     clear uhat, wh, des, dd, olsinflnc;


                     /* The GMM estimators: */

                     /* bleh = 1 acts as if Ey and Ex known.  bleh = 2 optimally accounts for the
                        fact that they are not known.
                     */

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


                                     @ --------The program creates the weighting matrix. ------------@
                                     clear w;
                                     ff = optw(c,mom,Emom,bleh);
                                     w = invit(moment(ff,0)./nncount);
                                     if w == -999999; if bb==0; goto bttdb; else;  goto bttfdb; endif; endif;

                                     @ ------------- Start of iteration loop ------------------------@
                                     iter = 1;
                                     tol=1e-9;                 @ Set the convergence criterion. @
                                     dc=1;                     @ Initialize the step length. @
                                     do until abs(dc) < tol;
                                         clear g, f, obj, sqz;

                                         f = deff(c,Emom,0,effsave);


                                         /*   g=gradtoni(c);*/
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
                                         if dc == -999999; if bb==0; goto bttdb; else; goto bttfdb; endif; endif;

                                         if maxsqez > 0;              @ This jumps to the subroutine that
                                                                        adjusts the step length. @
                                           { c_new,sqz } = squeeze(c,dc);

                                         else;
                                           c_new=c - dc;
                                         endif;


                                         dc=c_new-c;                     @ Update variables for the next iteration. @
                                         c=c_new;

                                        /*--------------Print out the results for the current iteration------------*/
                                        /*
                                        if bb>0;
                                        cls;
                                        "Number of Squeezes = ";; sqz;; bb;
                                        "Objective Function = ";; obj*nncount;

                                        "            Value         Step";
                                         format 14,6;
                                         mm=1;
                                         do until mm > nc;
                                          $"  ";; c[mm,1];; "  ";; dc[mm,1];
                                          mm=mm+1;
                                          endo;   pause(0.2);
                                        endif;
                                        */
                                        /*-----------------End of print out of current iteration--------------*/

                                         iter=iter +1;
                                         if iter >= maxiter;                     @ Quit iterating if necessary. @
                                         @  maxc(abs(dc));;
                                           "crap";; bb;;year;;estim; @
                                           goto escape;
                                       @   else;
                                           if abs(dc)<tol;
                                           "good";; bb;;year;;estim;
                                           endif; @
                                         endif;
                                     endo;
                                     @ ------ End of iteration loop ------------------------------------ @
                                     escape:
                                     ckeep[mmm,.] = c';
                                     objkeep[mmm,1] = obj;
                               endfor;

                               cidx = selif((1|2),objkeep .==minc(objkeep));
                               if rows(cidx)>1; cidx = cidx[1,1];endif;
                               c = ckeep[cidx,.]';

                         endif;
                         @ Compute t-ratios, etc.                                            @
                         f = deff(c,Emom,0,effsave);
                         ff = optw(c,mom,Emom,bleh);
                         if bb == 0;
                            effsave[1:neq,(estim-1)*nestim+year] =  f;
                         endif;

                         w = invit(moment(ff,0)./nncount);
                         if w == -999999; if bb==0; goto bttdb; else; goto bttfdb; endif; endif;
                         if estim > 1;
                            hansenj[1,estim] = nncount*f'w*f-oidcrit[estim,1];
                         endif;

                         if estim == 1;         g=grad1(c);
                         elseif estim == 2;     g=grad2(c);
                         elseif estim == 3;     g=grad3(c);
                         elseif estim == 4;     g=grad4(c);
                         elseif estim == 5;     g=grad5(c);
                         endif;

                         gwg = g'w*g;
                         vc=invit(gwg)/nncount;
                         if vc == -999999; if bb==0; goto bttdb; else; goto bttfdb; endif; endif;
                            if (not(diag(vc) > 0));
                               if bb==0; goto bttdb; else; goto bttfdb; endif;
                               endif;

                         stderr=sqrt(diag(vc));
                         inflnc=-nncount*vc*g'w*ff';

                         clear g, f, sqz, gwg, gwf, dc;

                         /* GMM estimator ends here. */

                         csave[1,nestim*(year-1)+estim] = c[1,1];
                         icsave[.,nestim*(year-1)+estim] = inflnc[1,.]';

                         @----This part substitutes in the trimmed estimates if necessary.-----------@
                         if rezat == 1;
                         b = allb[2,1];
                         if b > 0;
                           if csave[1,estim] < b;
                             csave[1,estim] = b;
                           elseif csave[1,estim] > rb;
                             csave[1,estim] = rb;
                           endif;
                         elseif b <= 0;
                           if csave[1,estim] > b;
                             csave[1,estim] = b;
                           elseif csave[1,estim] < rb;
                             csave[1,estim] = rb;
                           endif;
                         endif;
                         endif;
                         /* Make the coefficients on the perfectly measured regressors. */

                         /* Make Delta */

                         c11 = c[1,1];

                         gee = zeros(((nz+1)*2+1),nz);

                         gee[2:(nz+1),1:nz] = eye(nz);
                         gee[(nz+3):((nz+1)*2),1:nz] = -c11*eye(nz);

                         for qq(1, nz, 1); gee[((nz+1)*2+1),qq] = -mux[qq+1,1]; endfor;

                         @---Correct standard error and influence function---@

                         bigphi = ( (inEzz*zy_d')|(inEzz*zx_d') )|(-inflnc[1,.]);
                         avar = moment(bigphi',0)./nncount^2;
                         dstd = sqrt(diag(gee'avar*gee));
                         omega = gee'avar*gee;
                         phidel = -bigphi'gee;

                         for qq(1, nz, 1); dsave[qq,nestim*(year-1)+estim] = muy[qq+1,1] - c11*mux[qq+1,1]; endfor;
                         for qq(1, nz, 1); idsave[((qq-1)*nncount+1):(qq*nncount),nestim*(year-1)+estim] = phidel[.,qq]; endfor;



                         @----------First make the influence function for sigma_z.----------------------@

                         sigz = moment(za-meanc(za)',0)/nncount;

                         vecsigz = zeros(nz+nz*(nz-1)/2,1);
                         vecsigz[1:nz,1] = diag(sigz);

                         counter =  nz+nz*(nz-1)/2;
                         for qq(nz, 2, -1);
                           for ee(qq-1, 1, -1);
                             vecsigz[counter,1] = sigz[ee,qq];
                             counter=counter-1;
                           endfor;
                         endfor;


                         phiz = zeros(nncount,nz+nz*(nz-1)/2);
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

                         @---------Make the influence functions for the standard errors for rho2 and tau2.@

                         bigphi = phimux[2:nz+1,.]|phimuy[2:nz+1,.]|(phiz')|(-inflnc[1:4,.]);

                         gee = zeros(2*nz+(nz+nz*(nz-1)/2)+4,2); @ First column is for rho2 and the second is for tau2. @

                         @-----------First, do rho2---------------------------@


                         numer = muy[2:nz+1,1]'sigz*muy[2:nz+1,1]+c[1,1]^2*c[2,1];
                         denom =  numer+c[3,1]; ;

                         @------------derivatives wrt muy---------------------@

                         for qq(1, nz, 1);

                         gee[nz+qq,1] = (2*muy[2:nz+1,1]'sigz[.,qq])/denom - numer*(2*muy[2:nz+1,1]'sigz[.,qq])/(denom^2);

                         endfor;
                         @------------derivatives wrt the first part of sigz--@

                         for qq(1, nz, 1);

                         gee[2*nz+qq,1] = (muy[qq+1,1]^2)/denom - numer/(denom^2)*muy[qq+1,1]^2;

                         endfor;

                         @------------derivatives wrt the second part of sigz--@

                         counter = 2*nz+(nz+nz*(nz-1)/2);
                         for qq(nz, 2, -1);
                           for ee(qq-1, 1, -1);
                           gee[counter,1] = 2*muy[ee+1,1]*muy[qq+1,1]/denom - 2*numer/(denom^2)*muy[ee+1,1]*muy[qq+1,1];
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

                         gee[qq,2] = (2*mux[2:nz+1]'sigz[.,qq])/denom - numer*(2*mux[2:nz+1]'sigz[.,qq])/(denom^2);

                         endfor;

                         @------------derivatives wrt the first part of sigz--@

                         for qq(1, nz, 1);

                         gee[2*nz+qq,2] = (mux[qq+1,1]^2)/denom - numer/(denom^2)*mux[qq+1,1]^2;

                         endfor;

                         @------------derivatives wrt the second part of sigz--@

                         counter = 2*nz+(nz+nz*(nz-1)/2);
                         for qq(nz, 2, -1);
                           for ee(qq-1, 1, -1);
                           gee[counter,2] = 2*mux[ee+1,1]*mux[qq+1,1]/denom - 2*numer/(denom^2)*mux[ee+1,1]*mux[qq+1,1];
                           counter=counter-1;
                           endfor;
                         endfor;

                         gee[2*nz+(nz+nz*(nz-1)/2)+2,2] = 1/denom-numer/(denom^2);
                         gee[2*nz+(nz+nz*(nz-1)/2)+4,2] = -numer/(denom^2);


                         ipsave[.,nestim*(year-1)+estim] = -bigphi'gee[.,1];
                         itsave[.,nestim*(year-1)+estim] = -bigphi'gee[.,2];

                             psave[1,nestim*(year-1)+estim]    = rho;
                             tsave[1,nestim*(year-1)+estim]    = tau;


                         /* This trims rho2 and tau2 if necessary. */
                         if rezat ==1;
                         if psave[1,nestim*(year-1)+estim] > 1;
                         psave[1,nestim*(year-1)+estim] = 1;  endif;
                         if psave[1,nestim*(year-1)+estim] < 0;
                         psave[1,nestim*(year-1)+estim] = 0;  endif;

                         if tsave[1,nestim*(year-1)+estim] > 1;
                         tsave[1,nestim*(year-1)+estim] = 1;  endif;
                         if tsave[1,nestim*(year-1)+estim] < 0;
                         tsave[1,nestim*(year-1)+estim] = 0;  endif;
                         endif;


                     estim = estim + 1;
                     endo;                       @ The end of the estimator loop. @




                year = year + 1;
                endo;
             endif; /*end of the block out the calculations conditional */

             if checkmoms == 0;

                   /* This part does the classical minimum distance estimation of beta and alpha_1. */
                   /* For GMM. Optimal Weighting Matrix */

                   param = 1; @ 1 is beta, 2 is rho2, 3 is tau2, 4 is the sum, and 5 to nz+4 are the alphas. @
                   do while param <= nz+3;

                     idxnot = seqa(1,nestim,keepnyr) - 1;

                     estim = startestim;
                     do while estim <= nestim;

                         idx = idxnot + estim;

                         if param == 1;
                           saveit = csave;
                           w = invit(moment(icsave[.,idx],0)/nncount);
                         elseif param == 2;
                           saveit = psave;
                           w = invit(moment(ipsave[.,idx],0)/nncount);
                         elseif param == 3;
                           saveit = tsave;
                           w = invit(moment(itsave[.,idx],0)/nncount);
                         endif;

                         for qq(1, nz, 1);
                         if param == qq+3;
                           saveit = dsave[qq,.];
                           w = invit(moment(idsave[((qq-1)*nncount+1):(qq*nncount),idx],0)/nncount);

                         endif;
                         endfor;

                         if w == -999999; if bb==0; goto bttdb; else;  goto bttfdb; endif; endif;

                         g=ones(keepnyr,1);
                         theta = invpd(g'w*g)*g'w*saveit[1,idx]';

                         f=saveit[1,idx]'-theta;
                         stderror = sqrt(invpd(g'w*g)/nncount);
                         cmdtest = nncount*f'w*f;

                         if param == 1;
                                    if bb == 0;
                                             ccmd[1,estim] = theta;
                                             ccmd[3,estim] = cmdtest;
                                             ccmd[2,estim] = stderror;
                                    else;
                                             cpivot[estim,bb] = (theta-ccmd[1,estim])/stderror;
                                    endif;
                         elseif param == 2;
                                    if bb == 0;
                                             pcmd[1,estim] = theta;
                                             pcmd[3,estim] = cmdtest;
                                             pcmd[2,estim] = stderror;
                                    else;
                                             ppivot[estim,bb] = (theta-pcmd[1,estim])/stderror;
                                    endif;

                         elseif param == 3;
                                    if bb == 0;
                                             tcmd[1,estim] = theta;
                                             tcmd[3,estim] = cmdtest;
                                             tcmd[2,estim] = stderror;
                                    else;
                                             tpivot[estim,bb] = (theta-tcmd[1,estim])/stderror;
                                    endif;


                         endif;

                         for qq(1, nz, 1);
                         if param == qq+3;

                                    if bb == 0;
                                             dcmd[(qq-1)*3+1,estim] = theta;
                                             dcmd[(qq-1)*3+3,estim] = cmdtest;
                                             dcmd[(qq-1)*3+2,estim] = stderror;
                                    else;
                                             dpivot[(estim-1)*nz+qq,bb] = (theta-dcmd[(qq-1)*3+1,estim])/stderror;
                                    endif;

                         endif;
                         endfor;



                     estim = estim + 1;
                     endo;

                   param = param + 1;
                   endo;



                   /* This part does the classical minimum distance estimation of beta and alpha_1. */
                   /* FOR OLS. Optimal Weighting Matrix.*/

                   param = 1; @ 1 is beta, 2 is R2, 3 is the sum, and 4 through nz+1 are the alphas.@
                   do while param <= nz+3;

                     if param == 1;
                       saveit = bsave[1,.];
                       w = invit(moment(ibsave,0)/nncount);
                     elseif param == 2;
                       saveit = Rsave[1,.];
                       w = invit(moment(irsave,0)/nncount);
                     elseif param == 3;
                       saveit = asave[1,.];
                       w = invit(moment(iasave,0)/nncount);
                     endif;

                       for qq(1, nz, 1);
                       if param == qq+3;
                         saveit = zsave[qq,.];
                         w = invit(moment(izsave[((qq-1)*nncount+1):(qq*nncount),.],0)/nncount);
                       endif;
                       endfor;

                       g=ones(keepnyr,1);
                       theta = invpd(g'w*g)*g'w*saveit';

                   @ Compute t-ratios, etc.                                            @

                       f=saveit'-theta;
                       stderror = sqrt(invpd(g'w*g)/nncount);
                       cmdtest = nncount*f'w*f;

                       if param == 1;
                         if bb == 0;
                            bcmd[1,1] = theta;
                            bcmd[3,1] = cmdtest;
                            bcmd[2,1] = stderror;
                         else;
                                  bpivot[1,bb] = (theta-bcmd[1,1])/stderror;
                         endif;
                       elseif param == nz+2;
                         if bb == 0;
                            rcmd[1,1] = theta;
                            rcmd[3,1] = cmdtest;
                            rcmd[2,1] = stderror;
                         else;
                                  rpivot[1,bb] = (theta-rcmd[1,1])/stderror;
                         endif;
                       endif;

                       for qq(1, nz, 1);
                       if param == qq+1;
                         if bb == 0;
                            zcmd[(qq-1)*4+1,1] = theta;
                            zcmd[(qq-1)*4+3,1] = cmdtest;
                            zcmd[(qq-1)*4+2,1] = stderror;
                         else;
                                  zpivot[qq,bb] = (theta-zcmd[(qq-1)*4+1,1])/stderror;
                         endif;
                       endif;
                       endfor;

                   param = param + 1;
                   endo;


             endif;  /* End of the check moments conditional */




          /* This does the IV stuff. */

          if checkmoms == 2 or checkmoms == 0;
              biv = zeros(3,1);
              bnew = zeros(3,1);
              if keepnyr>2;
                   yt0 = selif(ally,yr.>3);
                   xt0 = selif(allx,yr.>3);
                   zt0 = selif(allz,yr.>3);

                   yt1 = selif(ally,yr.<(keepnyr-0) .and yr.>2);
                   xt1 = selif(allx,yr.<(keepnyr-0) .and yr.>2);
                   zt1 = selif(allz,yr.<(keepnyr-0) .and yr.>2);

                   yt2 = selif(ally,yr.<(keepnyr-1) .and yr.>1);
                   xt2 = selif(allx,yr.<(keepnyr-1) .and yr.>1);
                   zt2 = selif(allz,yr.<(keepnyr-1) .and yr.>1);

                   yt3 = selif(ally,yr.<(keepnyr-2));
                   xt3 = selif(allx,yr.<(keepnyr-2));
                   zt3 = selif(allz,yr.<(keepnyr-2));


                   des = ones(rows(yt1),1)~(xt0-xt1)~(zt0-zt1);
                      ivs = ones(rows(yt1),1)~(zt0-zt1)~xt2~zt2~(xt3)~(zt3);
                      zzzz = moment(ivs,0);
                      inezz = invpd(zzzz);
                      desivs = des'ivs;
                      biv = invpd(desivs*inezz*desivs')*desivs*inezz*(ivs'(yt0-yt1) - ivfsave);
                      ehat = (yt0-yt1) - des*biv;

                      wmx = clust((ehat.*ivs - ivfsave'),nncount,rows(yt0)/nncount);
                      inxzx = desivs*inezz*desivs';
                      inxzx = invpd(inxzx);


                      sbiv  = inxzx*desivs*inezz*wmx*inezz*desivs'inxzx;
                      sbiv  = sqrt(diag(sbiv));

                      jtest = ehat'ivs*inezz*ivs'ehat/meanc(ehat^2);


                   if bb == 0;
                        IVJcount = IVJcount + (cdfchic(jtest,cols(ivs)-cols(des)).<0.05);
                        bivsave[1,1] =  biv[2,1];
                        bivsave[2,1] = sbiv[2,1];
                        zivsave[1,1] =  biv[3,1];
                        zivsave[2,1] = sbiv[3,1];
                        ivfsave = ivs'ehat;
                   else;
                        ivpivot[1,bb] = (biv[2,1]-bivsave[1,1])/sbiv[2,1];
                        ivpivot[2,bb] = (biv[3,1]-zivsave[1,1])/sbiv[3,1];
                   endif;
                 @AB@
                 yrshape = reshape(yr,nncount,keepnyr);
                 yshape = reshape(ally,nncount,keepnyr);
                 xshape = reshape(allx,nncount,keepnyr);
                 zshape = reshape(allz,nncount,keepnyr);  @ yr~ally;
                                                           vec(yrshape)~vec(yshape);
                                                          @
                 touty = vec(yshape);

                 touty = touty[nncount*2+1:rows(touty),1] - touty[nncount+1:rows(touty)-nncount,1];


                 toutx = vec(xshape);
                 toutz = vec(zshape);
                 toutx = (toutx[nncount*2+1:rows(toutx),1] - toutx[nncount+1:rows(toutx)-nncount,1])~(toutz[nncount*2+1:rows(toutz),1] - toutz[nncount+1:rows(toutz)-nncount,1]);
                 toutx = ones(rows(toutx),1)~toutx;
                 nzab = 4;
                 toutz = zeros(nncount*(keepnyr-2),(keepnyr-2)*(nzab));

                 xz = vec(xshape)~vec(zshape)~vec(yshape);

                 for ii(1,keepnyr-2,1);
                  toutz[(ii-1)*nncount+1:ii*nncount,(ii-1)*(nzab)+1:ii*(nzab)]  = ones(nncount,1)~xz[(ii-1)*nncount+1:ii*nncount,1:2]~(xz[(ii+1)*nncount+1:(ii+2)*nncount,2] - xz[(ii-0)*nncount+1:(ii+1)*nncount,2] );
                 endfor;

                 inezz = invpd(moment(toutz,0));
                 bhat = invpd(toutx'toutz*inezz*toutz'toutx)*toutx'toutz*inezz*(toutz'touty - abfsave);

                 vhat = touty - toutx*bhat;
                 omega = zeros(cols(toutz),cols(toutz));

                 @
                 for ii(1,keepnyr-2,1);
                    fu = (toutz[(ii-1)*nncount+1:ii*nncount,(ii-1)*(nzab)+1:ii*(nzab)].*vhat[(ii-1)*nncount+1:ii*nncount,1]) - abfsave[(ii-1)*(nzab)+1:ii*(nzab),1]';
                    omega[(ii-1)*(nzab)+1:ii*(nzab),(ii-1)*(nzab)+1:ii*(nzab)]  =fu'fu;
                 endfor;
                 @

                 fu = toutz.*vhat;
                 omega = fu'fu;
                 iomega = invpd(omega);

                 bnew = invpd(toutx'toutz*iomega*toutz'toutx)*toutx'toutz*iomega*(toutz'touty - abfsave);
                 vhat = touty - toutx*bnew;

                 sbnew = invpd(toutx'toutz*iomega*toutz'toutx);
                 sbnew = sqrt(diag(sbnew));

                 jtest = vhat'toutz*iomega*toutz'vhat;

              endif;

              if bb == 0;
                   ABJcount = ABJcount + (cdfchic(jtest,cols(toutz)-cols(toutx)).<0.05);
                   babsave[1,1] = bnew[2,1];
                   babsave[2,1] = sbnew[2,1];
                   zabsave[1,1] = bnew[3,1];
                   zabsave[2,1] = sbnew[3,1];
                   abfsave = (vhat'toutz)';
              else;
                   abpivot[1,bb] = (bnew[2,1]-babsave[1,1])/sbnew[2,1];
                   abpivot[2,bb] = (bnew[3,1]-zabsave[1,1])/sbnew[3,1];
              endif;

           endif; /*End of the check the moments conditional */






            if ceil(bb/100)==(bb/100) and bb/=0@ and ceil(nrun/10)==(nrun/10)@;

                output off;
                "Baseline  Trial Number " nrun+1;; "               ";;
                "Bootstrap Trial Number " bb;
                @
                "CMD Estimator Intermediate Results";
                "Unconstrained Mean Absolute Deviations and biases";
                print;
                "Beta";
                "     OLS:   ";; (sumerrb/(nrun+1));; (sumb/(nrun+1));
                print;
                "     GMM:   ";; (sumerrg/(nrun+1));; (sumg/(nrun+1));
                print;
                "Alpha_1";
                "     OLS:   ";; (sumerrz/(nrun+1))';; (sumz/(nrun+1))';
                print;
                "     GMM:   ";; (sumerrd/(nrun+1));;  (sumd/(nrun+1));
                print;
                print; @
                output on;
           endif;
          bb = bb + 1;
          endo;                  /* end of the bootstrap loop */

          if checkmoms == 0;
             estim = startestim;
             do while estim<=nestim;

             pivotsort = sortc(cpivot[estim,.]',1);
             locritval = pivotsort[maxc(1|floor(rows(pivotsort)*0.025)),1];
             hicritval = pivotsort[ceil(rows(pivotsort)*0.975),1];
             tee = (((ccmd[1,estim]-beta)')./(ccmd[2,estim]'));
             bootc[estim,1] = bootc[estim,1] + (tee .> hicritval .or tee .< locritval);


             pivotsort = sortc(ppivot[1,.]',1);
             locritval = pivotsort[maxc(1|floor(rows(pivotsort)*0.025)),1];
             hicritval = pivotsort[ceil(rows(pivotsort)*0.975),1];
             tee = (((pcmd[1,estim]-rho2)')./(pcmd[2,estim]'));
             bootp[estim,1] = bootp[estim,1] + (tee .> hicritval .or tee .< locritval);


             pivotsort = sortc(tpivot[1,.]',1);
             locritval = pivotsort[maxc(1|floor(rows(pivotsort)*0.025)),1];
             hicritval = pivotsort[ceil(rows(pivotsort)*0.975),1];
             tee = (((tcmd[1,estim]-tau2)')./(tcmd[2,estim]'));
             boott[estim,1] = boott[estim,1] + (tee .> hicritval .or tee .< locritval);

             for qq(1, nz, 1);
             pivotsort = sortc(dpivot[(estim-1)*nz+qq,.]',1);
             locritval = pivotsort[maxc(1|floor(rows(pivotsort)*0.025)),1];
             hicritval = pivotsort[ceil(rows(pivotsort)*0.975),1];
             tee = (((dcmd[(qq-1)*3+1,estim]-delta[qq,1])')./(dcmd[(qq-1)*3+2,estim]'));
             bootd[estim,qq] = bootd[estim,qq] + (tee .> hicritval .or tee .< locritval);
             endfor;

             estim = estim + 1;
             endo;
          endif;

          if checkmoms == 2 or checkmoms == 0;

             pivotsort = sortc(abpivot[1,.]',1);
             locritval = pivotsort[maxc(1|floor(rows(pivotsort)*0.025)),1];
             hicritval = pivotsort[ceil(rows(pivotsort)*0.975),1];
             tee = (((babsave[1,1]-beta))./(babsave[2,1]));
             bootab[1,1] = bootab[1,1] + (tee .> hicritval .or tee .< locritval);


             pivotsort = sortc(abpivot[2,.]',1);
             locritval = pivotsort[maxc(1|floor(rows(pivotsort)*0.025)),1];
             hicritval = pivotsort[ceil(rows(pivotsort)*0.975),1];
             tee = (((zabsave[1,1]-delta))./(zabsave[2,1]));
             bootab[2,1] = bootab[2,1] + (tee .> hicritval .or tee .< locritval);


             pivotsort = sortc(ivpivot[1,.]',1);
             locritval = pivotsort[maxc(1|floor(rows(pivotsort)*0.025)),1];
             hicritval = pivotsort[ceil(rows(pivotsort)*0.975),1];
             tee = (((bivsave[1,1]-beta))./(bivsave[2,1]));
             bootiv[1,1] = bootiv[1,1] + (tee .> hicritval .or tee .< locritval);


             pivotsort = sortc(ivpivot[2,.]',1);
             locritval = pivotsort[maxc(1|floor(rows(pivotsort)*0.025)),1];
             hicritval = pivotsort[ceil(rows(pivotsort)*0.975),1];
             tee = (((zivsave[1,1]-delta))./(zivsave[2,1]));
             bootiv[2,1] = bootiv[2,1] + (tee .> hicritval .or tee .< locritval);
          endif;



          /* Update the counters. */

          EWjcount = EWjcount + (hansenj.>0);
          sumg     = sumg     + ccmd[1,.];
          sumseg   = sumseg   + ccmd[2,.];
          sumerr2g = sumerr2g + (ccmd[1,.] - beta)^2;
          sumerrg  = sumerrg  + abs(ccmd[1,.] - beta);
          gkeep[icount,.] = ccmd[1,.];

          idx=seqa(1,3,nz);
          sumd     = sumd     + dcmd[idx,.];
          sumsed   = sumsed   + dcmd[idx+1,.];
          sumerr2d = sumerr2d + (dcmd[idx,.] - acf)^2;
          sumerrd  = sumerrd  + abs(dcmd[idx,.] - acf);
          dkeep[icount,.] = dcmd[idx,.];

          sumt     = sumt     + tcmd[1,.];
          sumset   = sumset   + tcmd[2,.];
          sumerr2t = sumerr2t + (tcmd[1,.] - tau2)^2;
          sumerrt  = sumerrt  + abs(tcmd[1,.] - tau2);
          tkeep[icount,.] = tcmd[1,.];

          sump     = sump     + pcmd[1,.];
          sumsep   = sumsep   + pcmd[2,.];
          sumerr2p = sumerr2p + (pcmd[1,.] - rho2)^2;
          sumerrp  = sumerrp  + abs(pcmd[1,.] - rho2);
          pkeep[icount,.] = pcmd[1,.];

          testg = abs(ccmd[1,.] - beta);
          testd = abs(dcmd[idx,.] - acf);
          testp = abs(pcmd[1,.] - rho2);
          testt = abs(tcmd[1,.] - tau2);

          pcountg[1,.] =  pcountg[1,.]       + (testg .<= (0.2*beta));
          pcountd[1:nz,.] =  pcountd[1:nz,.] + (testd .<= (0.2*acf));
          pcountp[1,.] =  pcountp[1,.]       + (testp .<= (0.2*rho2));
          pcountt[1,.] =  pcountt[1,.]       + (testt .<= (0.2*tau2));

          sumtrueR = sumtrueR + rtruesave[1,1];
          sumtrueT = sumtrueT + ttruesave[1,1];

          sumb = sumb + bcmd[1,1];
          sumerr2b = sumerr2b + (bcmd[1,1] - beta)^2;
          sumerrb = sumerrb + abs(bcmd[1,1] - beta);
          sumseb = sumseb + bcmd[2,1];
          bkeep[icount,.] = bcmd[1,1];

          sumARR = sumARR + rcmd[1,1];
          sumerr2R = sumerr2R + (rcmd[1,1] - beta)^2;
          sumerrR = sumerrR + abs(rcmd[1,1] - beta);
          sumseR = sumseR + rcmd[2,1];
          rkeep[icount,.] = rcmd[1,.];

          idx=seqa(1,3,nz);
          sumz     = sumz     + zcmd[idx,.];
          sumsez   = sumsez   + zcmd[idx+1,.];
          sumerr2z = sumerr2z + (zcmd[idx,.] - acf)^2;
          sumerrz  = sumerrz  + abs(zcmd[idx,.] - acf);
          zkeep[icount,.] = zcmd[idx,.];

          testb = abs(bcmd[1,1] - beta);
          testa = abs(acmd[1,1] - sumc(acf[1:(1+ndum/2),1]));
          testz = abs(zcmd[idx,1] - acf);
          testR = abs(zcmd[idx,1] - rho2);

          pcountb[1,1] =  pcountb[1,1]       + (testb <= (0.2*beta));
          pcounta[1,1] =  pcounta[1,1]       + (testa <= (0.2));
          pcountR[1,1] =  pcountR[1,1]       + (testR <= (0.2*rho2));
          pcountz[1:nz,1] =  pcountz[1:nz,1] + (testz <= (0.2*acf));


          idx=seqa(1,2,nz);
          sumbiv = sumbiv + bivsave[1,1];
          sumerr2biv = sumerr2biv + (bivsave[1,1] - beta)^2;
          sumerrbiv = sumerrbiv + abs(bivsave[1,1] - beta);
          sumsebiv = sumsebiv + bivsave[2,1];
          bivkeep[icount,.] = bivsave[1,1];

          sumziv     = sumziv     + zivsave[idx,1];
          sumseziv   = sumseziv   + zivsave[idx+1,1];
          sumerr2ziv = sumerr2ziv + (zivsave[idx,1] - acf)^2;
          sumerrziv  = sumerrziv  + abs(zivsave[idx,1] - acf);
          zivkeep[icount,.] = zivsave[idx,1];

          testbiv = abs(bivsave[1,1] - beta);
          testziv = abs(zivsave[idx,1] - acf);

          pcountbiv[1,1] =  pcountbiv[1,1] + (testbiv <= (0.2*beta));
          pcountziv[1:nz,1] =  pcountziv[idx,1] + (testziv <= (0.2*acf));

          sumbab = sumbab + babsave[1,1];
          sumerr2bab = sumerr2bab + (babsave[1,1] - beta)^2;
          sumerrbab = sumerrbab + abs(babsave[1,1] - beta);
          sumsebab = sumsebab + babsave[2,1];
          babkeep[icount,.] = babsave[1,1];

          sumzab     = sumzab     + zabsave[idx,1];
          sumsezab   = sumsezab   + zabsave[idx+1,.];
          sumerr2zab = sumerr2zab + (zabsave[idx,1] - acf)^2;
          sumerrzab  = sumerrzab  + abs(zabsave[idx,1] - acf);
          zabkeep[icount,.] = zabsave[idx,1];

          testbab = abs(babsave[1,1] - beta);
          testzab = abs(zabsave[idx,1] - acf);

          pcountbab[1,1] =  pcountbab[1,1] + (testbab <= (0.2*beta));
          pcountzab[1:nz,1] =  pcountzab[1,1] + (testzab <= (0.2*acf));

          /*  Do the test that alpha_1 is zero. */

          teedelta[.,1] = teedelta[.,1] + ((((dcmd[1,.]-delta)')./(dcmd[2,.]')) .> acritval);
          teedelta[.,2] = teedelta[.,2] + (abs(((dcmd[1,.]-delta)')./(dcmd[2,.]')) .> twocrit);

endif; @ This is the end of the icount==0 conditional from way up there @
          nrun = icount;
          icount = icount + 1;


      endo;     @This ends the trial loop.@

      /* ------------------------------------------------------------*/
      /*
      cls;
      format 14,6;
      output on;
      print "  Mcmd Monte Carlo Program";
      print;
      print "sample size =" nncount;
      print "number of samples =" nrun;
      print;


      format /rd 8,3;



      print " Probability that the hypothesis (alpha_1 - alpha_1(true)) > 0 is rejected.     ";
      print " Columns are estimators. Rows are significance levels of   ";
      print " 0.05 ";

      for qq(1, nestim, 1); "&";; teedelta[qq,1]/nrun;; endfor; "\\" "\\";


      print " Probability that the hypothesis (alpha_1 - alpha_1(true)) /= 0 is rejected.     ";
      print " Columns are estimators. Rows are significance levels of   ";
      print " 0.05 ";

      for qq(1, nestim, 1); "&";; teedelta[qq,2]/nrun;; endfor; "\\" "\\";

      print " Probability that the overidentifying restrictions are rejected ";
      print " Columns are estimators. ";
      print;
      EWjcount/nrun;
      print;
      print "OLS-IV OID test ";
      IVJcount/nrun;
      print;
      print "AB OID test ";
      ABJcount/nrun;
      print;
      print;
      print "OLS estimates of rho2 and tau2 ";
      sumtrueR/nrun;
      sumtrueT/nrun;

      "===============================================================================================";
      "===============================================================================================";
      "===============================================================================================";
      "===============================================================================================";
      */

                     for qq(startestim, nestim, 1); bootc[qq,1]/nrun;; endfor;
      for ee(1,nz,1);for qq(startestim, nestim, 1); bootd[qq,ee]/nrun;; endfor;  endfor;
                       bootiv[1,1]/nrun;;
                       bootiv[2,1]/nrun;;
                       bootab[1,1]/nrun;;
                       bootab[2,1]/nrun;


      @ OLD OUTPUT
     @
     /*

      /*  This business makes a LaTeX table.   */

      print;
      format /rd 8,3;


      "------------------------------CMD in percent-----------------------------------";

      "&"  sumb/nrun/beta-1;;               "&"  sumbiv/nrun/beta-1;;               "&"  sumbab/nrun/beta-1;;                  for qq (1, nestim, 1); "&" sumg[1,qq]/nrun/beta-1;;          endfor; "\\" "\\";; " [2pt] ";
      "&"  median(bkeep[1:nrun,.])/beta-1;; "&"  median(bivkeep[1:nrun,.])/beta-1;; "&"  median(babkeep[1:nrun,.])/beta-1;;    for qq (1, nestim, 1); "&" median(gkeep[1:nrun,qq])/beta-1;; endfor; "\\" "\\";; " [2pt] ";
      "&"  sumerrb/nrun;;                   "&"  sumerrbiv/nrun;;                   "&"  sumerrbab/nrun;;                      for qq (1, nestim, 1); "&" sumerrg[1,qq]/nrun;;              endfor; "\\" "\\";; " [2pt] ";
      "&"  pcountb[1,1]/nrun;;              "&"  pcountbiv[1,1]/nrun;;              "&"  pcountbab[1,1]/nrun;;                 for qq (1, nestim, 1); "&" pcountg[1,qq]/nrun;;              endfor; "\\" "\\";; " [11pt] ";


      for ee(1, nz, 1);
      "&" sumz[ee,1]/nrun/delta[ee,1]-1;;          "&" sumziv/nrun/delta[ee,1]-1;;                "&" sumzab/nrun/delta[ee,1]-1;;               for qq (1, nestim, 1); "&" sumd[ee,qq]/nrun/delta[ee,1]-1;;                       endfor; "\\" "\\";; " [2pt] ";
      "&" median(zkeep[1:nrun,ee])/delta[ee,1]-1;; "&" median(zivkeep[1:nrun,.])/delta[ee,1]-1;;  "&" median(zabkeep[1:nrun,.])/delta[ee,1]-1;; for qq (1, nestim, 1); "&" median(dkeep[1:nrun,(ee-1)*nestim+qq])/delta[ee,1]-1;; endfor; "\\" "\\";; " [2pt] ";
      "&" sumerrz[ee,1]/nrun;;                     "&" sumerrziv/nrun;;                           "&" sumerrzab/nrun;;                          for qq (1, nestim, 1); "&" sumerrd[ee,qq]/nrun;;                                  endfor; "\\" "\\";; " [2pt] ";
      "&" pcountz[ee,1]/nrun;;                     "&" pcountziv[1,1]/nrun;;                      "&" pcountzab[1,1]/nrun;;                     for qq (1, nestim, 1); "&" pcountd[ee,qq]/nrun;;                                  endfor; "\\" "\\";; " [11pt] ";
      endfor;

      "&" sumARR/nrun/rho2-1;;              "&" "  -----  ";;     "&" "  -----  ";;     for qq (1, nestim, 1); "&" sump[1,qq]/nrun/rho2-1;;          endfor; "\\" "\\";; " [2pt] ";
      "&"  median(rkeep[1:nrun,.])/rho2-1;; "&" "  -----  ";;     "&" "  -----  ";;     for qq (1, nestim, 1); "&" median(pkeep[1:nrun,qq])/rho2-1;; endfor; "\\" "\\";; " [2pt] ";
      "&" sumerrr/nrun;;                    "&" "  -----  ";;     "&" "  -----  ";;     for qq (1, nestim, 1); "&" sumerrp[1,qq]/nrun;;              endfor; "\\" "\\";; " [2pt] ";
      "&" pcountr/nrun;;                    "&" "  -----  ";;     "&" "  -----  ";;     for qq (1, nestim, 1); "&" pcountp[1,qq]/nrun;;              endfor; "\\" "\\";; " [11pt] ";

      "&" "  -----  ";;      "&" "  -----  ";;     "&" "  -----  ";;         for qq (1, nestim, 1); "&" sumt[1,qq]/nrun/tau2-1;;          endfor; "\\" "\\";; " [2pt] ";
      "&" "  -----  ";;      "&" "  -----  ";;     "&" "  -----  ";;         for qq (1, nestim, 1); "&" median(tkeep[1:nrun,qq])/tau2-1;; endfor; "\\" "\\";; " [2pt] ";
      "&" "  -----  ";;      "&" "  -----  ";;     "&" "  -----  ";;         for qq (1, nestim, 1); "&" sumerrt[1,qq]/nrun;;              endfor; "\\" "\\";; " [2pt] ";
      "&" "  -----  ";;      "&" "  -----  ";;     "&" "  -----  ";;         for qq (1, nestim, 1); "&" pcountt[1,qq]/nrun;;              endfor; "\\" "\\";; " [11pt] ";

      "------------------------------CMD----------------------------------------------";

      "&"  sumb/nrun;;               "&"  sumbiv/nrun;;               "&"  sumbab/nrun;;                  for qq (1, nestim, 1); "&" sumg[1,qq]/nrun;;          endfor; "\\" "\\";; " [2pt] ";
      "&"  median(bkeep[1:nrun,.]);; "&"  median(bivkeep[1:nrun,.]);; "&"  median(babkeep[1:nrun,.]);;    for qq (1, nestim, 1); "&" median(gkeep[1:nrun,qq]);; endfor; "\\" "\\";; " [2pt] ";
      "&"  sumerrb/nrun;;            "&"  sumerrbiv/nrun;;            "&"  sumerrbab/nrun;;               for qq (1, nestim, 1); "&" sumerrg[1,qq]/nrun;;       endfor; "\\" "\\";; " [2pt] ";
      "&"  pcountb[1,1]/nrun;;       "&"  pcountbiv[1,1]/nrun;;       "&"  pcountbab[1,1]/nrun;;          for qq (1, nestim, 1); "&" pcountg[1,qq]/nrun;;       endfor; "\\" "\\";; " [11pt] ";


      for ee(1, nz, 1);
      "&" sumz[ee,1]/nrun;;          "&" sumziv/nrun;;                "&" sumzab/nrun;;               for qq (1, nestim, 1); "&" sumd[ee,qq]/nrun;;                       endfor; "\\" "\\";; " [2pt] ";
      "&" median(zkeep[1:nrun,ee]);; "&" median(zivkeep[1:nrun,.]);;  "&" median(zabkeep[1:nrun,.]);; for qq (1, nestim, 1); "&" median(dkeep[1:nrun,(ee-1)*nestim+qq]);; endfor; "\\" "\\";; " [2pt] ";
      "&" sumerrz[ee,1]/nrun;;       "&" sumerrziv/nrun;;             "&" sumerrzab/nrun;;            for qq (1, nestim, 1); "&" sumerrd[ee,qq]/nrun;;                    endfor; "\\" "\\";; " [2pt] ";
      "&" pcountz[ee,1]/nrun;;       "&" pcountziv[1,1]/nrun;;        "&" pcountzab[1,1]/nrun;;       for qq (1, nestim, 1); "&" pcountd[ee,qq]/nrun;;                    endfor; "\\" "\\";; " [11pt] ";
      endfor;

      "&" sumARR/nrun;;                     "&" "  -----  ";;     "&" "  -----  ";;     for qq (1, nestim, 1); "&" sump[1,qq]/nrun;;          endfor; "\\" "\\";; " [2pt] ";
      "&"  median(rkeep[1:nrun,.]);;        "&" "  -----  ";;     "&" "  -----  ";;     for qq (1, nestim, 1); "&" median(pkeep[1:nrun,qq]);; endfor; "\\" "\\";; " [2pt] ";
      "&" sumerrr/nrun;;                    "&" "  -----  ";;     "&" "  -----  ";;     for qq (1, nestim, 1); "&" sumerrp[1,qq]/nrun;;       endfor; "\\" "\\";; " [2pt] ";
      "&" pcountr/nrun;;                    "&" "  -----  ";;     "&" "  -----  ";;     for qq (1, nestim, 1); "&" pcountp[1,qq]/nrun;;       endfor; "\\" "\\";; " [11pt] ";

      "&" "  -----  ";;      "&" "  -----  ";;     "&" "  -----  ";;         for qq (1, nestim, 1); "&" sumt[1,qq]/nrun;;          endfor; "\\" "\\";; " [2pt] ";
      "&" "  -----  ";;      "&" "  -----  ";;     "&" "  -----  ";;         for qq (1, nestim, 1); "&" median(tkeep[1:nrun,qq]);; endfor; "\\" "\\";; " [2pt] ";
      "&" "  -----  ";;      "&" "  -----  ";;     "&" "  -----  ";;         for qq (1, nestim, 1); "&" sumerrt[1,qq]/nrun;;       endfor; "\\" "\\";; " [2pt] ";
      "&" "  -----  ";;      "&" "  -----  ";;     "&" "  -----  ";;         for qq (1, nestim, 1); "&" pcountt[1,qq]/nrun;;       endfor; "\\" "\\";; " [11pt] ";

      "------------------------------------------------IV and AB----------------------------------------------";


      "&"  sumbiv/nrun;;               "&"  sumbab/nrun;
      "&"  median(bivkeep[1:nrun,.]);; "&"  median(babkeep[1:nrun,.]);
      "&"  sumerrbiv/nrun;;            "&"  sumerrbab/nrun;
      "&"  pcountbiv[1,1]/nrun;;       "&"  pcountbab[1,1]/nrun;

      for ee(1, nz, 1);
      "&" sumziv/nrun;;                "&" sumzab/nrun;
      "&" median(zivkeep[1:nrun,.]);;  "&" median(zabkeep[1:nrun,.]);
      "&" sumerrziv/nrun;;             "&" sumerrzab/nrun;
      "&" pcountziv[1,1]/nrun;;        "&" pcountzab[1,1]/nrun;
      endfor;

      print; print;


      format /rd 9,3;
      "Summary Statistics ";

      "&" avAy/(keepnyr*nrun);; "&" avVy/(keepnyr*nrun);; "&" avskewy/(keepnyr*nrun);; "&" avkurty/(keepnyr*nrun);;
      "&" avstd5y/(keepnyr*nrun);; "&" avstd6y/(keepnyr*nrun);; "&" avstd7y/(keepnyr*nrun);; "\\" "\\";
      "&" avAx/(keepnyr*nrun);; "&" avVx/(keepnyr*nrun);; "&" avskewx/(keepnyr*nrun);; "&" avkurtx/(keepnyr*nrun);;
      "&" avstd5x/(keepnyr*nrun);; "&" avstd6x/(keepnyr*nrun);; "&" avstd7x/(keepnyr*nrun);; "\\" "\\";
      "&" avAz/(keepnyr*nrun);; "&" avVz/(keepnyr*nrun);; "&" avskewz/(keepnyr*nrun);; "&" avkurtz/(keepnyr*nrun);;
      "&" avstd5z/(keepnyr*nrun);; "&" avstd6z/(keepnyr*nrun);; "&" avstd7z/(keepnyr*nrun);; "\\" "\\";

      print; "CHI"; print;
      "&" avAchi/(keepnyr*nrun);; "&" avVchi/(keepnyr*nrun);; "&" avskewchi/(keepnyr*nrun);; "&" avkurtchi/(keepnyr*nrun);;
      "&" avstd5chi/(keepnyr*nrun);; "&" avstd6chi/(keepnyr*nrun);; "&" avstd7chi/(keepnyr*nrun);; "\\" "\\";
      print;

      print; "OLS residual"; print;
      "&" avAu/(keepnyr*nrun);; "&" avVu/(keepnyr*nrun);; "&" avskewu/(keepnyr*nrun);; "&" avkurtu/(keepnyr*nrun);;
      "&" avstd5u/(keepnyr*nrun);; "&" avstd6u/(keepnyr*nrun);; "&" avstd7u/(keepnyr*nrun);; "\\" "\\";
      print;

      "Summary Statistics for the Residuals";
      print;
      "&" meAy/(keepnyr*nrun);; "&" meVy/(keepnyr*nrun);; "&" meskewy/(keepnyr*nrun);; "&" mekurty/(keepnyr*nrun);;
      "&" mestd5y/(keepnyr*nrun);; "&" mestd6y/(keepnyr*nrun);; "&" mestd7y/(keepnyr*nrun);; "\\" "\\";
      "&" meAx/(keepnyr*nrun);; "&" meVx/(keepnyr*nrun);; "&" meskewx/(keepnyr*nrun);; "&" mekurtx/(keepnyr*nrun);;
      "&" mestd5x/(keepnyr*nrun);; "&" mestd6x/(keepnyr*nrun);; "&" mestd7x/(keepnyr*nrun);; "\\" "\\";


      "Summary Statistics for the Serial Correlations";
      "&"  averhoy/nrun;
      "&"  averhox/nrun;
      "&"  averhoz/nrun;

      "Summary Statistics for Ey2x and Eyx2";

      "&" aveEy2x/(keepnyr*nrun);
      "&" aveEyx2/(keepnyr*nrun);


      "End Time ";; timestr(0);
      */


@ ----------------- Subroutines follow --------------- @

@ Subroutine to define the f vector for the optimal weighting matrix. @

proc optw(a,mmm,ooo,wf);
local f, ei;

f = zeros(nncount,neq);

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

ei = zeros(nncount,neq);

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
f = f - (bb>0)*effsave[1:neq,(estim-1)*nestim+year];


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
proc (2) = squeeze(s_c,s_dc);
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
      @  sumbw = sumbw + b;@
        numwh = numwh + 1;
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
      @  sumbw = sumbw + b;@
        numwh = numwh + 1;
        dc = -999999;
      endif;
    endif;

retp(dc);
endp;

proc clust(h,capN,capT);
local ff,qq,hi;

    ff = 0;
    for qq(1,capN,1);
        hi = sumc(h[(qq-1)*(capT)+1:qq*(capT),.]);
        ff = ff + hi*hi';
    endfor;

    ff = ff/capN;

ff = invit(ff);
if scalerr(ff)>0;
  retp(-999999);
else;
  retp(ff);
endif;

endp;

eop:
@
"Start time:  ";; timestr(ttt);
"End time:    ";; timestr(0);
"Start date:  ";; datestr(ddd);
"End date:    ";; datestr(0);
@
output off;
end;
