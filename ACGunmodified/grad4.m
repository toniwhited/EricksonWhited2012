%subroutine to compute gradient matrix.


function[dvdc] = grad4(a,neq);


dvdc1=zeros(neq,1);
dvdc2=zeros(neq,1);
dvdc3=zeros(neq,1);
dvdc4=zeros(neq,1);
dvdc5=zeros(neq,1);
dvdc6=zeros(neq,1);
dvdc7=zeros(neq,1);
dvdc8=zeros(neq,1);
dvdc9=zeros(neq,1);
dvdc10=zeros(neq,1);
dvdc11=zeros(neq,1);
dvdc12=zeros(neq,1);


dvdc1(1,1)=-2*a(1,1)*a(2,1);
dvdc1(2,1)=-a(2,1);
dvdc1(3,1)=0;
dvdc1(4,1)=-2*a(1,1)*a(5,1);
dvdc1(5,1)=-a(5,1);
dvdc1(6,1)=-3*(a(1,1).^2)*a(6,1)-3*a(2,1)*a(3,1);
dvdc1(7,1)=-2*a(1,1)*(a(6,1)+a(2,1)*a(4,1));
dvdc1(8,1)=-a(6,1)-3*a(2,1)*a(4,1);
dvdc1(9,1)=0;
dvdc1(10,1)=-3*(a(1,1).^2)*a(5,1);
dvdc1(11,1)=-4*(a(1,1).^3)*a(9,1)-12*a(1,1)*a(5,1)*a(3,1)...
               -4*a(2,1)*a(8,1);
dvdc1(12,1)=-3*(a(1,1).^2)*(a(9,1)+a(5,1)*a(4,1))...
               -3*a(5,1)*a(3,1);
dvdc1(13,1)=-2*a(1,1)*(a(9,1)+3*a(5,1)*a(4,1)+a(2,1)*a(7,1));
dvdc1(14,1)=-(a(9,1)+6*a(5,1)*a(4,1)+4*a(2,1)*a(7,1));


dvdc1(15,1)=-4*(a(1,1).^3)*a(6,1)-12*a(1,1)*a(2,1)*a(3,1);
dvdc1(16,1)=0;
dvdc1(17,1)=-5*(a(1,1).^4)*a(12,1)-30*(a(1,1).^2)*a(6,1)*a(3,1)...
                -20*a(1,1)*a(5,1)*a(8,1)-5*a(2,1)*a(10,1);
dvdc1(18,1)=-4*(a(1,1).^3)*(a(12,1)+a(6,1)*a(4,1))...
                  -12*a(1,1)*(a(6,1)*a(3,1)+a(2,1)*a(3,1)*a(4,1))...
              -4*a(5,1)*a(8,1);
dvdc1(19,1)=-3*(a(1,1).^2)*(a(12,1)+3*a(6,1)*a(4,1)+a(5,1)*a(7,1))...
                  -(3*a(6,1)*a(3,1)+9*a(2,1)*a(3,1)*a(4,1));
dvdc1(20,1)=-2*a(1,1)*(a(12,1)+6*a(6,1)*a(4,1)+4*a(5,1)*a(7,1)...
                           +a(2,1)*a(11,1));
dvdc1(21,1)=-(a(12,1)+10*a(6,1)*a(4,1)+10*a(5,1)*a(7,1)...
                    +5*a(2,1)*a(11,1));




dvdc2(1,1)=-(a(1,1).^2);
dvdc2(2,1)=-a(1,1);
dvdc2(3,1)=-1;
dvdc2(4,1)=0;
dvdc2(5,1)=0;
dvdc2(6,1)=-3*a(1,1)*a(3,1);
dvdc2(7,1)=-(a(1,1).^2)*a(4,1)-a(3,1);
dvdc2(8,1)=-3*a(1,1)*a(4,1);
dvdc2(9,1)=0;
dvdc2(10,1)=0;
dvdc2(11,1)=-4*a(1,1)*a(8,1);
dvdc2(12,1)=-a(8,1);
dvdc2(13,1)=-(a(1,1).^2)*a(7,1);
dvdc2(14,1)=-4*a(1,1)*a(7,1);


dvdc2(15,1)=-6*(a(1,1).^2)*a(3,1);
dvdc2(16,1)=-6*a(4,1);
dvdc2(17,1)=-5*a(1,1)*a(10,1);
dvdc2(18,1)=-6*(a(1,1).^2)*a(3,1)*a(4,1)-a(10,1);
dvdc2(19,1)=-9*a(1,1)*a(3,1)*a(4,1);
dvdc2(20,1)=-(a(1,1).^2)*a(11,1)-6*a(3,1)*a(4,1);
dvdc2(21,1)=-5*a(1,1)*a(11,1);






dvdc3(1,1)=-1;
dvdc3(2,1)=0;
dvdc3(3,1)=0;
dvdc3(4,1)=0;
dvdc3(5,1)=0;
dvdc3(6,1)=-3*a(1,1)*a(2,1);
dvdc3(7,1)=-(a(2,1)+a(4,1));
dvdc3(8,1)=0;
dvdc3(9,1)=0;
dvdc3(10,1)=0;
dvdc3(11,1)=-6*(a(1,1).^2)*a(5,1);
dvdc3(12,1)=-3*a(1,1)*a(5,1);
dvdc3(13,1)=-(a(5,1)+a(7,1));
dvdc3(14,1)=0;


dvdc3(15,1)=-6*(a(1,1).^2)*a(2,1);
dvdc3(16,1)=0;
dvdc3(17,1)=-10*(a(1,1).^3)*a(6,1);
dvdc3(18,1)=-6*(a(1,1).^2)*(a(6,1)+a(2,1));
dvdc3(19,1)=-a(1,1)*(3*a(6,1)+9*a(2,1)*a(4,1));
dvdc3(20,1)=-(a(6,1)+6*a(2,1)*a(4,1)+a(11,1));
dvdc3(21,1)=0;






dvdc4(1,1)=0;
dvdc4(2,1)=0;
dvdc4(3,1)=-1;
dvdc4(4,1)=0;
dvdc4(5,1)=0;
dvdc4(6,1)=0;
dvdc4(7,1)=-(a(1,1).^2)*a(2,1)-a(3,1);
dvdc4(8,1)=-3*a(1,1)*a(2,1);
dvdc4(9,1)=0;
dvdc4(10,1)=0;
dvdc4(11,1)=0;
dvdc4(12,1)=-(a(1,1).^3)*a(5,1)-a(8,1);
dvdc4(13,1)=-3*(a(1,1).^2)*a(5,1);
dvdc4(14,1)=-6*a(1,1)*a(5,1);


dvdc4(15,1)=0;
dvdc4(16,1)=-6*a(2,1);
dvdc4(17,1)=0;
dvdc4(18,1)=-(a(1,1).^4)*a(6,1)-6*(a(1,1).^2)*a(2,1)*a(3,1)-a(10,1);
dvdc4(19,1)=-3*(a(1,1).^3)*a(6,1)-9*a(1,1)*a(2,1)*a(3,1);
dvdc4(20,1)=-6*(a(1,1).^2)*a(6,1)-6*a(2,1)*a(3,1);
dvdc4(21,1)=-10*a(1,1)*a(6,1);






dvdc5(1,1)=0;
dvdc5(2,1)=0;
dvdc5(3,1)=0;
dvdc5(4,1)=-(a(1,1).^2);
dvdc5(5,1)=-a(1,1);
dvdc5(6,1)=0;
dvdc5(7,1)=0;
dvdc5(8,1)=0;
dvdc5(9,1)=-1;
dvdc5(10,1)=-(a(1,1).^3);
dvdc5(11,1)=-6*(a(1,1).^2)*a(3,1);
dvdc5(12,1)=-(a(1,1).^3)*a(4,1)-3*a(1,1)*a(3,1);
dvdc5(13,1)=-3*(a(1,1).^2)*a(4,1)-a(3,1);
dvdc5(14,1)=-6*a(1,1)*a(4,1);


dvdc5(15,1)=0;
dvdc5(16,1)=0;
dvdc5(17,1)=-10*(a(1,1).^2)*a(8,1);
dvdc5(18,1)=-4*a(1,1)*a(8,1);
dvdc5(19,1)=-(a(1,1).^3)*a(7,1)-a(8,1);
dvdc5(20,1)=-4*(a(1,1).^2)*a(7,1);
dvdc5(21,1)=-10*a(1,1)*a(7,1);






dvdc6(1,1)=0;
dvdc6(2,1)=0;
dvdc6(3,1)=0;
dvdc6(4,1)=0;
dvdc6(5,1)=0;
dvdc6(6,1)=-(a(1,1).^3);
dvdc6(7,1)=-(a(1,1).^2);
dvdc6(8,1)=-a(1,1);
dvdc6(9,1)=0;
dvdc6(10,1)=0;
dvdc6(11,1)=0;
dvdc6(12,1)=0;
dvdc6(13,1)=0;
dvdc6(14,1)=0;


dvdc6(15,1)=-(a(1,1).^4);
dvdc6(16,1)=-1;
dvdc6(17,1)=-10*(a(1,1).^3)*a(3,1);
dvdc6(18,1)=-(a(1,1).^4)*a(4,1)-6*(a(1,1).^2)*a(3,1);
dvdc6(19,1)=-3*(a(1,1).^3)*a(4,1)-3*a(1,1)*a(3,1);
dvdc6(20,1)=-6*(a(1,1).^2)*a(4,1)-a(3,1);
dvdc6(21,1)=-10*a(1,1)*a(4,1);






dvdc7(1,1)=0;
dvdc7(2,1)=0;
dvdc7(3,1)=0;
dvdc7(4,1)=0;
dvdc7(5,1)=0;
dvdc7(6,1)=0;
dvdc7(7,1)=0;
dvdc7(8,1)=0;
dvdc7(9,1)=-1;
dvdc7(10,1)=0;
dvdc7(11,1)=0;
dvdc7(12,1)=0;
dvdc7(13,1)=-(a(1,1).^2)*a(2,1)-a(3,1);
dvdc7(14,1)=-4*a(1,1)*a(2,1);


dvdc7(15,1)=0;
dvdc7(16,1)=0;
dvdc7(17,1)=0;
dvdc7(18,1)=0;
dvdc7(19,1)=-(a(1,1).^3)*a(5,1)-a(8,1);
dvdc7(20,1)=-4*(a(1,1).^2)*a(5,1);
dvdc7(21,1)=-10*a(1,1)*a(5,1);






dvdc8(1,1)=0;
dvdc8(2,1)=0;
dvdc8(3,1)=0;
dvdc8(4,1)=0;
dvdc8(5,1)=0;
dvdc8(6,1)=0;
dvdc8(7,1)=0;
dvdc8(8,1)=0;
dvdc8(9,1)=0;
dvdc8(10,1)=-1;
dvdc8(11,1)=-4*a(1,1)*a(2,1);
dvdc8(12,1)=-(a(2,1)+a(4,1));
dvdc8(13,1)=0;
dvdc8(14,1)=0;


dvdc8(15,1)=0;
dvdc8(16,1)=0;
dvdc8(17,1)=-10*(a(1,1).^2)*a(5,1);
dvdc8(18,1)=-4*a(1,1)*a(5,1);
dvdc8(19,1)=-(a(5,1)+a(7,1));
dvdc8(20,1)=0;
dvdc8(21,1)=0;






dvdc9(1,1)=0;
dvdc9(2,1)=0;
dvdc9(3,1)=0;
dvdc9(4,1)=0;
dvdc9(5,1)=0;
dvdc9(6,1)=0;
dvdc9(7,1)=0;
dvdc9(8,1)=0;
dvdc9(9,1)=0;
dvdc9(10,1)=0;
dvdc9(11,1)=-(a(1,1).^4);
dvdc9(12,1)=-(a(1,1).^3);
dvdc9(13,1)=-(a(1,1).^2);
dvdc9(14,1)=-a(1,1);


dvdc9(15,1)=0;
dvdc9(16,1)=0;
dvdc9(17,1)=0;
dvdc9(18,1)=0;
dvdc9(19,1)=0;
dvdc9(20,1)=0;
dvdc9(21,1)=0;




dvdc10(1,1)=0;
dvdc10(2,1)=0;
dvdc10(3,1)=0;
dvdc10(4,1)=0;
dvdc10(5,1)=0;
dvdc10(6,1)=0;
dvdc10(7,1)=0;
dvdc10(8,1)=0;
dvdc10(9,1)=0;
dvdc10(10,1)=0;
dvdc10(11,1)=0;
dvdc10(12,1)=0;
dvdc10(13,1)=0;
dvdc10(14,1)=0;


dvdc10(15,1)=-1;
dvdc10(16,1)=0;
dvdc10(17,1)=-5*a(1,1)*a(2,1);
dvdc10(18,1)=-(a(2,1)+a(4,1));
dvdc10(19,1)=0;
dvdc10(20,1)=0;
dvdc10(21,1)=0;




dvdc11(1,1)=0;
dvdc11(2,1)=0;
dvdc11(3,1)=0;
dvdc11(4,1)=0;
dvdc11(5,1)=0;
dvdc11(6,1)=0;
dvdc11(7,1)=0;
dvdc11(8,1)=0;
dvdc11(9,1)=0;
dvdc11(10,1)=0;
dvdc11(11,1)=0;
dvdc11(12,1)=0;
dvdc11(13,1)=0;
dvdc11(14,1)=0;


dvdc11(15,1)=0;
dvdc11(16,1)=-1;
dvdc11(17,1)=0;
dvdc11(18,1)=0;
dvdc11(19,1)=0;
dvdc11(20,1)=-(a(1,1).^2)*a(2,1)-a(3,1);
dvdc11(21,1)=-5*a(1,1)*a(2,1);




dvdc12(1,1)=0;
dvdc12(2,1)=0;
dvdc12(3,1)=0;
dvdc12(4,1)=0;
dvdc12(5,1)=0;
dvdc12(6,1)=0;
dvdc12(7,1)=0;
dvdc12(8,1)=0;
dvdc12(9,1)=0;
dvdc12(10,1)=0;
dvdc12(11,1)=0;
dvdc12(12,1)=0;
dvdc12(13,1)=0;
dvdc12(14,1)=0;


dvdc12(15,1)=0;
dvdc12(16,1)=0;
dvdc12(17,1)=-(a(1,1).^5);
dvdc12(18,1)=-(a(1,1).^4);
dvdc12(19,1)=-(a(1,1).^3);
dvdc12(20,1)=-(a(1,1).^2);
dvdc12(21,1)=-a(1,1);

dvdc=[dvdc1,dvdc2,dvdc3,dvdc4,dvdc5,dvdc6,dvdc7,dvdc8,dvdc9,dvdc10,dvdc11,dvdc12];
