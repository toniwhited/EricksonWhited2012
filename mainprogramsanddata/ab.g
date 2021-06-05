new ,60000;
ttt = time;
ddd = date;
call sysstate(14,1e-128);
assets = 1;


@ Step 1: ------Name of the Gauss data set.----------------------------------@

    dataset = "abdset";

@ Step 2: ------Name of the output file.-------------------------------------@

    output file = ab.out reset;


    open f1=^dataset;
    vnames = getname(dataset);
    dta = readr(f1,rowsf(f1));
    year = dta[.,indcv("fyear",vnames)];
    ydum = seqa(minc(year),1,(maxc(year)-minc(year)+1));
    ydum = year.==(ydum');


    doyd = 1;
    if doyd == 1;
      ifo = ones(rows(year),1)~ydum[.,2:cols(ydum)];
    else;
      ifo = ones(rows(year),1);
    endif;

    if assets == 0;

    let yname = Dik;                                        yidx = indcv(yname,vnames);
    let wname = Dq Dcfk;                                    widx = indcv(wname,vnames);
    let zname = LLq Dcfk@ LLLq@  LLcfk@ LLLcfk@;            zidx = indcv(zname,vnames);      nz = rows(zidx)+1;

    else;

    let yname = Dia;                                         yidx = indcv(yname,vnames);
    let wname = Dmb Dcfa;                                    widx = indcv(wname,vnames);
    let zname = LLmb Dcfa@ LLLq@  LLcfa@ LLLcfk@;            zidx = indcv(zname,vnames);      nz = rows(zidx)+1;

    endif;

    nyr  = maxc(year)-minc(year)+1;
    nco = rows(dta)/nyr;
    newdta = zeros(1,cols(dta));
    idx = seqa(1,nyr,nco);
    for ii(1,nyr,1);
       newdta  = newdta|dta[idx+ii-1,.];
    endfor;
    dta = newdta[2:rows(newdta),.];

    format /rd 8,3;

    n = rows(dta);
    icept=ones(nco,1);

    y = packr(dta[.,yidx]);
    w = packr(ifo~dta[.,widx]);
    ally = dta[.,yidx];
    allw = ifo~dta[.,widx];
    allz = dta[.,zidx];
    zomega = zeros(nco*nyr,nyr*nz);
    for ii(1,nyr,1);
    zomega[(ii-1)*nco+1:ii*nco,(ii-1)*nz+1:ii*nz]  = icept~allz[(ii-1)*nco+1:ii*nco,.];
    endfor;
    z = packr(zomega);
    zomega = missrv(zomega,0);

    inezz = invpd(moment(z,0));
    bhat = invpd(w'z*inezz*z'w)*w'z*inezz*z'y;

    vhat = missrv(ally,0) - missrv(allw,0)*bhat;  @    bhat;   @

    omega = zeros(cols(zomega),cols(zomega));
    for ii(1,nyr,1);
    for jj(1,nyr,1);@  ii;;jj;  @
      omega[(ii-1)*nz+1:ii*nz,(jj-1)*nz+1:jj*nz] = (vhat[(ii-1)*nco+1:ii*nco,1].*zomega[(ii-1)*nco+1:ii*nco,(ii-1)*nz+1:ii*nz])'
                                                   (vhat[(jj-1)*nco+1:jj*nco,1].*zomega[(jj-1)*nco+1:jj*nco,(jj-1)*nz+1:jj*nz]);

    endfor;
    endfor;

    iomega = invpd(omega);@iomega;@

    bnew = invpd(w'z*iomega*z'w)*w'z*iomega*z'y;
    avar = invpd(w'z*iomega*z'w);
    stderr = sqrt(diag(avar));

    "Unconstrained parameters and standard errors";
    "beta  ";; bnew[rows(bnew)-1];; stderr[rows(stderr)-1];
    "alpha1";; bnew[rows(bnew)];;   stderr[rows(stderr)];
 @   bnew~stderr;   print; @

    vhat = missrv(ally,0) - missrv(allw,0)*bnew;  @    bhat;   @

    for ii(1,nyr,1);
    for jj(1,nyr,1);@  ii;;jj;  @
      omega[(ii-1)*nz+1:ii*nz,(jj-1)*nz+1:jj*nz] = (vhat[(ii-1)*nco+1:ii*nco,1].*zomega[(ii-1)*nco+1:ii*nco,(ii-1)*nz+1:ii*nz])'
                                                   (vhat[(jj-1)*nco+1:jj*nco,1].*zomega[(jj-1)*nco+1:jj*nco,(jj-1)*nz+1:jj*nz]);

    endfor;
    endfor;

    iomega = invpd(omega);@iomega;@

    jstat = vhat'zomega*iomega*zomega'vhat;
    pval = cdfchic(jstat,rows(iomega)-cols(w));
    "J-test and p-value";
    jstat;; pval;


"Start time:  ";; timestr(ttt);
"End time:    ";; timestr(0);
"Start date:  ";; datestr(ddd);
"End date:    ";; datestr(0);

output off;
end;
