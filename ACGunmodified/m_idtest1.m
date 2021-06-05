function [Stat Stat_P] = m_idtest1(D, Xnames, t, n, year)


Y = D(:,1);
xaux1 = D(:,2);
if ~(sum( D(:,3)))
    display ('Error: No intercept. m_idtest1')
    return;
end
Intercep = D(:,3);
zeta = zeros(size(D,1), size(D,2)-3); 
for i=4:size(D,2)
    zeta(:,i-3) = D(:,i);
end


%step 1:------name of the data set.----------------------------------------
%step 2:------name of the output file.-------------------------------------
%step 3:------name of the left-hand-side variable.-------------------------
yname = 'Y';
%step 4:------name of the mismeasured variable.----------------------------
xname = Xnames{1};
%step 4:------names of the perfectly measured variables.-------------------
zname = char( Xnames{2:end} );
%step 5:------set the number of time periods.-------------------------------
nyr=t;
%step 6:------set the first time period.-----------------------------------
fyear=1;
%step 7:------set the number of observations in the largest cross-section.-
nco=n;
%step 8:-----set the output mode.-----------------------------------------
latex=0;

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%                           End of user input.
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%--------open data set.-------------------------------------------------
%---------number of perfectly measured regressors.----------------------
nz=size(zname,1)-1;
%---------------initialize the counters.---------------------------
csave=zeros(nyr,2);
%-------set the number of equations & parameters.----------------
nc=2;
neq=2;
%-------now we enter the main body of the program.-----------




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


    %---------construct the moments for subsequent estimation.----
    y2=y_d.^2;
    x2=x_d.^2;
    yx=x_d.*y_d;
    y2x=(y_d.^2).*x_d;
    yx2=(x_d.^2).*y_d;


    ey2x=mean(y2x);
    eyx2=mean(yx2);


    ezz = (iz'*iz)/n;
    inezz = inv(ezz);

    idx = ones(nz+1,1);

    ey2_z=mean(y2(:,idx).*iz);
    eyx_z=mean(yx(:,idx).*iz);
    ex2_z=mean(x2(:,idx).*iz);
    zy_d=iz.*y_d(:,idx);
    zx_d=iz.*x_d(:,idx);


    %-------initialize the parameters.---------------------------




    c=zeros(nc,1);


    c(1,1)=ey2x;
    c(2,1)=eyx2;


    %compute t-ratios,etc.


    f=zeros(n,neq);


    f(:,1)=y2x-c(1,1)+(-2*eyx_z*inezz*zy_d'-ey2_z*inezz*zx_d')';
    f(:,2)=yx2-c(2,1)+(-ex2_z*inezz*zy_d'-2*eyx_z*inezz*zx_d')';


    ff=f'*f;

    w =    inv(ff./n);

    g=-eye(nc);
    gwg=g'*w*g;


    nww=n*c'*w*c;
    test=1-chi2cdf(nww,2);


    csave(year,1)=nww;
    csave(year,2)=test;


    year=year+1;
end;


% -------------This business makes a LaTeX table-----------

% if latex==1;
% 
%     outfmt1 = '&  %8.3f   \\\\ \n';
%     outfmt2 = '&( %8.3f ) \\\\ \n';
%     for jj=1:nyr;
%         fprintf(fid,outfmt1,csave(jj,1));
%         fprintf(fid,outfmt2,csave(jj,2));
%     end;
% 
% 
% else;

%     outfmt1 = '  %8.3f   \n';
%     outfmt2 = '( %8.3f ) \n';
%     for jj=1:nyr;
%         fprintf(fid,outfmt1,csave(jj,1));
%         fprintf(fid,outfmt2,csave(jj,2));
%     end;
% 
% 
% end;
% 
% fprintf(fid,'\n');
% fprintf(fid,'start time:  %20s \n',datestr(tttttt));
% sss=datestr(now);
% fprintf(fid,'end time:    %20s',sss);
% fclose(fid);
% toc

Stat = csave(:,1);
Stat_P = csave(:,2);
