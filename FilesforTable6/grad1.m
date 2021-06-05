%subroutine to compute gradient matrix.


function[dvdc] = grad1(a,neq);

dvdc1=zeros(neq,1);
dvdc2=zeros(neq,1);
dvdc3=zeros(neq,1);
dvdc4=zeros(neq,1);
dvdc5=zeros(neq,1);


dvdc1(1,1)=-2*a(1,1)*a(2,1);
dvdc1(2,1)=-a(2,1);
dvdc1(3,1)=0;
dvdc1(4,1)=-2*a(1,1)*a(5,1);
dvdc1(5,1)=-a(5,1);


dvdc2(1,1)=-(a(1,1).^2);
dvdc2(2,1)=-a(1,1);
dvdc2(3,1)=-1;
dvdc2(4,1)=0;
dvdc2(5,1)=0;


dvdc3(1,1)=-1;
dvdc3(2,1)=0;
dvdc3(3,1)=0;
dvdc3(4,1)=0;
dvdc3(5,1)=0;


dvdc4(1,1)=0;
dvdc4(2,1)=0;
dvdc4(3,1)=-1;
dvdc4(4,1)=0;
dvdc4(5,1)=0;


dvdc5(1,1)=0;
dvdc5(2,1)=0;
dvdc5(3,1)=0;
dvdc5(4,1)=-(a(1,1).^2);
dvdc5(5,1)=-a(1,1);




dvdc=[dvdc1,dvdc2,dvdc3,dvdc4,dvdc5];
