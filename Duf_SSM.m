function [SSM1, SSM2]=Duf_SSM(lam1,lam2,k,kappa2,kappa3,sys)

T=[1 1; lam1 lam2];
if sys==1
 a2= (3*sqrt(-k/kappa3)*kappa3)/((lam1 - 2*lam2)*(lam1 - lam2));
 a3=(kappa3 + 12*a2*sqrt(-k/kappa3)*kappa3)/((lam1 - 3*lam2)*(lam1 - lam2));
 a4= (5*(a2*kappa3 + 3*a2^2*sqrt(-k/kappa3)*kappa3 + 3*a3*sqrt(-k/kappa3)*kappa3))/((lam1 - 4*lam2)*(lam1 - lam2));
 a5=  (-((6*a2^2*kappa3)/(lam1 - lam2)) - (3*a3*kappa3)/(lam1 - lam2) - (6*a2^3*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - ...
      (30*a2*a3*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - (12*a4*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) + (3*a2^2*kappa3)/(-lam1 + lam2) + ...
      (3*a3*kappa3)/(-lam1 + lam2) + (6*a2*a3*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) + (6*a4*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2))/...
     (-lam1 + 5*lam2);
 a6= (1/(-lam1 + 6*lam2))*(-((6*a2^3*kappa3)/(lam1 - lam2)) - (15*a2*a3*kappa3)/(lam1 - lam2) - (4*a4*kappa3)/(lam1 - lam2) - ...
      (21*a2^2*a3*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - (18*a3^2*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - ...
      (36*a2*a4*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - (15*a5*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) + (a2^3*kappa3)/(-lam1 + lam2) + ...
      (6*a2*a3*kappa3)/(-lam1 + lam2) + (3*a4*kappa3)/(-lam1 + lam2) + (3*a3^2*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) + ...
      (6*a2*a4*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) + (6*a5*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2));
 a7= (1/(-lam1 +      7*lam2))*(-((2*a2^4*kappa3)/(lam1 - lam2)) - (21*a2^2*a3*      kappa3)/(lam1 - lam2) - (9*a3^2*kappa3)/(lam1 - lam2) - ...
         (18*a2*a4*kappa3)/(lam1 - lam2) - (5*a5*kappa3)/(lam1 -       lam2) - (24*a2*a3^2*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - ...
         (24*a2^2*a4*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - (42*a3*      a4*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - ...
         (42*a2*a5*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - (18*a6* sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) + (3*a2^2*a3*kappa3)/(-lam1 + lam2) + ...
         (3*a3^2*kappa3)/(-lam1 + lam2) + (6*a2*a4*kappa3)/(-lam1 + lam2) + (3*a5*kappa3)/(-lam1 + lam2) + (6*a3*a4*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) + ...
     (6*a2*a5* sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) +  (6*a6*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2));
a8=0;
a9=0;
a10=0;
a11=0;
a12=0;
a13=0;
a14=0;
a15=0;

p2=linspace(-1,1,100);
 p1=(a2.*p2.^2+a3.*p2.^3+a4.*p2.^4+a5.*p2.^5+a6.*p2.^6+a7.*p2.^7+a8.*p2.^8+a9.*p2.^9+a10.*p2.^10+a11.*p2.^11+a12.*p2.^12+a13.*p2.^13+a14.*p2.^14+a15.*p2.^15);
 SSM2=T*[p1;p2];  

 
  a2=(3*sqrt(-k/kappa3)*kappa3)/((lam1 - lam2)*(2*lam1 - lam2));
 a3= (kappa3 + 12*a2*sqrt(-k/kappa3)*kappa3)/((lam1 - lam2)*(3*lam1 - lam2));
 a4= (5*(a2*kappa3 + 3*a2^2*sqrt(-k/kappa3)*kappa3 + 3*a3*sqrt(-k/kappa3)*kappa3))/((lam1 - lam2)*(4*lam1 - lam2));
 a5=((3*a2^2*kappa3)/(lam1 - lam2) + (3*a3*kappa3)/(lam1 - lam2) + (6*a2*a3*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) + ...
      (6*a4*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - (6*a2^2*kappa3)/(-lam1 + lam2) - (3*a3*kappa3)/(-lam1 + lam2) - ...
      (6*a2^3*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) - (30*a2*a3*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) - ...
      (12*a4*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2))/(5*lam1 - lam2);
      
  a6= (1/(6*lam1 - lam2))*((a2^3*kappa3)/(lam1 - lam2) + (6*a2*a3*kappa3)/(lam1 - lam2) + (3*a4*kappa3)/(lam1 - lam2) + ...
      (3*a3^2*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) + (6*a2*a4*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) + ...
      (6*a5*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - (6*a2^3*kappa3)/(-lam1 + lam2) - (15*a2*a3*kappa3)/(-lam1 + lam2) - ...
      (4*a4*kappa3)/(-lam1 + lam2) - (21*a2^2*a3*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) - (18*a3^2*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) - ...
      (36*a2*a4*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) - (15*a5*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2));
  
 a7=(1/(7*lam1 - lam2))*((3*a2^2*a3*kappa3)/(lam1 - lam2) + (3*a3^2*kappa3)/(lam1 - lam2) + (6*a2*a4*kappa3)/(lam1 - lam2) + ...
      (3*a5*kappa3)/(lam1 - lam2) + (6*a3*a4*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) + (6*a2*a5*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) + ...
      (6*a6*sqrt(-k/kappa3)*kappa3)/(lam1 - lam2) - (2*a2^4*kappa3)/(-lam1 + lam2) - (21*a2^2*a3*kappa3)/(-lam1 + lam2) - ...
      (9*a3^2*kappa3)/(-lam1 + lam2) - (18*a2*a4*kappa3)/(-lam1 + lam2) - (5*a5*kappa3)/(-lam1 + lam2) - ...
      (24*a2*a3^2*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) - (24*a2^2*a4*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) - ...
      (42*a3*a4*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) - (42*a2*a5*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2) - ...
      (18*a6*sqrt(-k/kappa3)*kappa3)/(-lam1 + lam2));
  
  
 a8=0;
a9=0;
a10=0;
a11=0;
a12=0;
a13=0;
a14=0;
a15=0;

 p1=linspace(-1,1,100);
 p2=(a2.*p1.^2+a3.*p1.^3+a4.*p1.^4+a5.*p1.^5+a6.*p1.^6+a7.*p1.^7++a8.*p1.^8+a9.*p1.^9+a10.*p1.^10+a11.*p1.^11+a12.*p1.^12+a13.*p1.^13+a14.*p1.^14+a15.*p1.^15);
 
 SSM1=T*[p1;p2];

else
a2= kappa2/((lam1 - 2*lam2)*(lam1 - lam2));
a3= (kappa3 + 4*kappa2*a2)/((lam1 - 3*lam2)*(lam1 - lam2));
a4=(5*(kappa3*a2 + kappa2*a2^2 + kappa2*a3))/((lam1 - 4*lam2)*(lam1 - lam2));
a5=(-((9*kappa3*a2^2)/(lam1 - lam2)) - (2*kappa2*a2^3)/(lam1 - lam2) - (6*kappa3*a3)/(lam1 - lam2) - (12*kappa2*a2*a3)/(lam1 - lam2) - ...
      (6*kappa2*a4)/(lam1 - lam2))/(-lam1 + 5*lam2);
  
a6=(-((7*kappa3*a2^3)/(lam1 - lam2)) - (21*kappa3*a2*a3)/(lam1 - lam2) - (7*kappa2*a2^2*a3)/(lam1 - lam2) - ...
      (7*kappa2*a3^2)/(lam1 - lam2) - (7*kappa3*a4)/(lam1 - lam2) - (14*kappa2*a2*a4)/(lam1 - lam2) - (7*kappa2*a5)/(lam1 - lam2))/...
     (-lam1 + 6*lam2);
     
a7= (-((2*kappa3*a2^4)/(lam1 - lam2)) - (24*kappa3*a2^2*a3)/(lam1 - lam2) - (12*kappa3*a3^2)/(lam1 - lam2) - ...
      (8*kappa2*a2*a3^2)/(lam1 - lam2) - (24*kappa3*a2*a4)/(lam1 - lam2) - (8*kappa2*a2^2*a4)/(lam1 - lam2) - ...
      (16*kappa2*a3*a4)/(lam1 - lam2) - (8*kappa3*a5)/(lam1 - lam2) - (16*kappa2*a2*a5)/(lam1 - lam2) - (8*kappa2*a6)/(lam1 - lam2))/...
     (-lam1 + 7*lam2);
     
a8=(1/(-lam1 + 8*lam2))*(-((9*kappa3*a2^3*a3)/(lam1 - lam2)) - (27*kappa3*a2*a3^2)/(lam1 - lam2) - (3*kappa2*a3^3)/(lam1 - lam2) - ...
      (27*kappa3*a2^2*a4)/(lam1 - lam2) - (27*kappa3*a3*a4)/(lam1 - lam2) - (18*kappa2*a2*a3*a4)/(lam1 - lam2) - ...
      (9*kappa2*a4^2)/(lam1 - lam2) - (27*kappa3*a2*a5)/(lam1 - lam2) - (9*kappa2*a2^2*a5)/(lam1 - lam2) - (18*kappa2*a3*a5)/(lam1 - lam2) - ...
      (9*kappa3*a6)/(lam1 - lam2) - (18*kappa2*a2*a6)/(lam1 - lam2) - (9*kappa2*a7)/(lam1 - lam2));
      
a9= (1/(-lam1 + 9*lam2))*(-((15*kappa3*a2^2*a3^2)/(lam1 - lam2)) - (10*kappa3*a3^3)/(lam1 - lam2) - (10*kappa3*a2^3*a4)/(lam1 - lam2) - ...
      (60*kappa3*a2*a3*a4)/(lam1 - lam2) - (10*kappa2*a3^2*a4)/(lam1 - lam2) - (15*kappa3*a4^2)/(lam1 - lam2) - ...
      (10*kappa2*a2*a4^2)/(lam1 - lam2) - (30*kappa3*a2^2*a5)/(lam1 - lam2) - (30*kappa3*a3*a5)/(lam1 - lam2) - ...
      (20*kappa2*a2*a3*a5)/(lam1 - lam2) - (20*kappa2*a4*a5)/(lam1 - lam2) - (30*kappa3*a2*a6)/(lam1 - lam2) - ...
      (10*kappa2*a2^2*a6)/(lam1 - lam2) - (20*kappa2*a3*a6)/(lam1 - lam2) - (10*kappa3*a7)/(lam1 - lam2) - (20*kappa2*a2*a7)/(lam1 - lam2) - ...
      (10*kappa2*a8)/(lam1 - lam2));
      
a10=(1/(-lam1 + 10*lam2))*(-((11*kappa3*a2*a3^3)/(lam1 - lam2)) - (33*kappa3*a2^2*a3*a4)/(lam1 - lam2) - ...
      (33*kappa3*a3^2*a4)/(lam1 - lam2) - (33*kappa3*a2*a4^2)/(lam1 - lam2) - (11*kappa2*a3*a4^2)/(lam1 - lam2) - ...
      (11*kappa3*a2^3*a5)/(lam1 - lam2) - (66*kappa3*a2*a3*a5)/(lam1 - lam2) - (11*kappa2*a3^2*a5)/(lam1 - lam2) - ...
      (33*kappa3*a4*a5)/(lam1 - lam2) - (22*kappa2*a2*a4*a5)/(lam1 - lam2) - (11*kappa2*a5^2)/(lam1 - lam2) - ...
      (33*kappa3*a2^2*a6)/(lam1 - lam2) - (33*kappa3*a3*a6)/(lam1 - lam2) - (22*kappa2*a2*a3*a6)/(lam1 - lam2) - ...
      (22*kappa2*a4*a6)/(lam1 - lam2) - (33*kappa3*a2*a7)/(lam1 - lam2) - (11*kappa2*a2^2*a7)/(lam1 - lam2) - ...
      (22*kappa2*a3*a7)/(lam1 - lam2) - (11*kappa3*a8)/(lam1 - lam2) - (22*kappa2*a2*a8)/(lam1 - lam2) - (11*kappa2*a9)/(lam1 - lam2));
      
      
a11=(1/(-lam1 + 11*lam2))*(-((3*kappa3*a3^4)/(lam1 - lam2)) - (36*kappa3*a2*a3^2*a4)/(lam1 - lam2) - ...
      (18*kappa3*a2^2*a4^2)/(lam1 - lam2) - (36*kappa3*a3*a4^2)/(lam1 - lam2) - (4*kappa2*a4^3)/(lam1 - lam2) - ...
      (36*kappa3*a2^2*a3*a5)/(lam1 - lam2) - (36*kappa3*a3^2*a5)/(lam1 - lam2) - (72*kappa3*a2*a4*a5)/(lam1 - lam2) - ...
      (24*kappa2*a3*a4*a5)/(lam1 - lam2) - (18*kappa3*a5^2)/(lam1 - lam2) - (12*kappa2*a2*a5^2)/(lam1 - lam2) - ...
      (12*kappa3*a2^3*a6)/(lam1 - lam2) - (72*kappa3*a2*a3*a6)/(lam1 - lam2) - (12*kappa2*a3^2*a6)/(lam1 - lam2) - ...
      (36*kappa3*a4*a6)/(lam1 - lam2) - (24*kappa2*a2*a4*a6)/(lam1 - lam2) - (24*kappa2*a5*a6)/(lam1 - lam2) - ...
      (36*kappa3*a2^2*a7)/(lam1 - lam2) - (36*kappa3*a3*a7)/(lam1 - lam2) - (24*kappa2*a2*a3*a7)/(lam1 - lam2) - ...
      (24*kappa2*a4*a7)/(lam1 - lam2) - (36*kappa3*a2*a8)/(lam1 - lam2) - (12*kappa2*a2^2*a8)/(lam1 - lam2) - ...
      (24*kappa2*a3*a8)/(lam1 - lam2) - (12*kappa3*a9)/(lam1 - lam2) - (24*kappa2*a2*a9)/(lam1 - lam2) - (12*kappa2*a10)/(lam1 - lam2));
  
a12=(1/(-lam1 + 12*lam2))*(-((13*kappa3*a3^3*a4)/(lam1 - lam2)) - (39*kappa3*a2*a3*a4^2)/(lam1 - lam2) - ...
      (13*kappa3*a4^3)/(lam1 - lam2) - (39*kappa3*a2*a3^2*a5)/(lam1 - lam2) - (39*kappa3*a2^2*a4*a5)/(lam1 - lam2) - ...
      (78*kappa3*a3*a4*a5)/(lam1 - lam2) - (13*kappa2*a4^2*a5)/(lam1 - lam2) - (39*kappa3*a2*a5^2)/(lam1 - lam2) - ...
      (13*kappa2*a3*a5^2)/(lam1 - lam2) - (39*kappa3*a2^2*a3*a6)/(lam1 - lam2) - (39*kappa3*a3^2*a6)/(lam1 - lam2) - ...
      (78*kappa3*a2*a4*a6)/(lam1 - lam2) - (26*kappa2*a3*a4*a6)/(lam1 - lam2) - (39*kappa3*a5*a6)/(lam1 - lam2) - ...
      (26*kappa2*a2*a5*a6)/(lam1 - lam2) - (13*kappa2*a6^2)/(lam1 - lam2) - (13*kappa3*a2^3*a7)/(lam1 - lam2) - ...
      (78*kappa3*a2*a3*a7)/(lam1 - lam2) - (13*kappa2*a3^2*a7)/(lam1 - lam2) - (39*kappa3*a4*a7)/(lam1 - lam2) - ...
      (26*kappa2*a2*a4*a7)/(lam1 - lam2) - (26*kappa2*a5*a7)/(lam1 - lam2) - (39*kappa3*a2^2*a8)/(lam1 - lam2) - ...
      (39*kappa3*a3*a8)/(lam1 - lam2) - (26*kappa2*a2*a3*a8)/(lam1 - lam2) - (26*kappa2*a4*a8)/(lam1 - lam2) - ...
      (39*kappa3*a2*a9)/(lam1 - lam2) - (13*kappa2*a2^2*a9)/(lam1 - lam2) - (26*kappa2*a3*a9)/(lam1 - lam2) - (13*kappa3*a10)/(lam1 - lam2) - ...
      (26*kappa2*a2*a10)/(lam1 - lam2) - (13*kappa2*a11)/(lam1 - lam2));
      
a13=(1/(-lam1 + 13*lam2))*(-((21*kappa3*a3^2*a4^2)/(lam1 - lam2)) - (14*kappa3*a2*a4^3)/(lam1 - lam2) - ...
      (14*kappa3*a3^3*a5)/(lam1 - lam2) - (84*kappa3*a2*a3*a4*a5)/(lam1 - lam2) - (42*kappa3*a4^2*a5)/(lam1 - lam2) - ...
      (21*kappa3*a2^2*a5^2)/(lam1 - lam2) - (42*kappa3*a3*a5^2)/(lam1 - lam2) - (14*kappa2*a4*a5^2)/(lam1 - lam2) - ...
      (42*kappa3*a2*a3^2*a6)/(lam1 - lam2) - (42*kappa3*a2^2*a4*a6)/(lam1 - lam2) - (84*kappa3*a3*a4*a6)/(lam1 - lam2) - ...
      (14*kappa2*a4^2*a6)/(lam1 - lam2) - (84*kappa3*a2*a5*a6)/(lam1 - lam2) - (28*kappa2*a3*a5*a6)/(lam1 - lam2) - ...
      (21*kappa3*a6^2)/(lam1 - lam2) - (14*kappa2*a2*a6^2)/(lam1 - lam2) - (42*kappa3*a2^2*a3*a7)/(lam1 - lam2) - ...
      (42*kappa3*a3^2*a7)/(lam1 - lam2) - (84*kappa3*a2*a4*a7)/(lam1 - lam2) - (28*kappa2*a3*a4*a7)/(lam1 - lam2) - ...
      (42*kappa3*a5*a7)/(lam1 - lam2) - (28*kappa2*a2*a5*a7)/(lam1 - lam2) - (28*kappa2*a6*a7)/(lam1 - lam2) - ...
      (14*kappa3*a2^3*a8)/(lam1 - lam2) - (84*kappa3*a2*a3*a8)/(lam1 - lam2) - (14*kappa2*a3^2*a8)/(lam1 - lam2) - ...
      (42*kappa3*a4*a8)/(lam1 - lam2) - (28*kappa2*a2*a4*a8)/(lam1 - lam2) - (28*kappa2*a5*a8)/(lam1 - lam2) - ...
      (42*kappa3*a2^2*a9)/(lam1 - lam2) - (42*kappa3*a3*a9)/(lam1 - lam2) - (28*kappa2*a2*a3*a9)/(lam1 - lam2) - ...
      (28*kappa2*a4*a9)/(lam1 - lam2) - (42*kappa3*a2*a10)/(lam1 - lam2) - (14*kappa2*a2^2*a10)/(lam1 - lam2) - ...
      (28*kappa2*a3*a10)/(lam1 - lam2) - (14*kappa3*a11)/(lam1 - lam2) - (28*kappa2*a2*a11)/(lam1 - lam2) - (14*kappa2*a12)/(lam1 - lam2));
      
a14=(1/(-lam1 + 14*lam2))*(-((15*kappa3*a3*a4^3)/(lam1 - lam2)) - (45*kappa3*a3^2*a4*a5)/(lam1 - lam2) - ...
      (45*kappa3*a2*a4^2*a5)/(lam1 - lam2) - (45*kappa3*a2*a3*a5^2)/(lam1 - lam2) - (45*kappa3*a4*a5^2)/(lam1 - lam2) - ...
      (5*kappa2*a5^3)/(lam1 - lam2) - (15*kappa3*a3^3*a6)/(lam1 - lam2) - (90*kappa3*a2*a3*a4*a6)/(lam1 - lam2) - ...
      (45*kappa3*a4^2*a6)/(lam1 - lam2) - (45*kappa3*a2^2*a5*a6)/(lam1 - lam2) - (90*kappa3*a3*a5*a6)/(lam1 - lam2) - ...
      (30*kappa2*a4*a5*a6)/(lam1 - lam2) - (45*kappa3*a2*a6^2)/(lam1 - lam2) - (15*kappa2*a3*a6^2)/(lam1 - lam2) - ...
      (45*kappa3*a2*a3^2*a7)/(lam1 - lam2) - (45*kappa3*a2^2*a4*a7)/(lam1 - lam2) - (90*kappa3*a3*a4*a7)/(lam1 - lam2) - ...
      (15*kappa2*a4^2*a7)/(lam1 - lam2) - (90*kappa3*a2*a5*a7)/(lam1 - lam2) - (30*kappa2*a3*a5*a7)/(lam1 - lam2) - ...
      (45*kappa3*a6*a7)/(lam1 - lam2) - (30*kappa2*a2*a6*a7)/(lam1 - lam2) - (15*kappa2*a7^2)/(lam1 - lam2) - ...
      (45*kappa3*a2^2*a3*a8)/(lam1 - lam2) - (45*kappa3*a3^2*a8)/(lam1 - lam2) - (90*kappa3*a2*a4*a8)/(lam1 - lam2) - ...
      (30*kappa2*a3*a4*a8)/(lam1 - lam2) - (45*kappa3*a5*a8)/(lam1 - lam2) - (30*kappa2*a2*a5*a8)/(lam1 - lam2) - ...
      (30*kappa2*a6*a8)/(lam1 - lam2) - (15*kappa3*a2^3*a9)/(lam1 - lam2) - (90*kappa3*a2*a3*a9)/(lam1 - lam2) - ...
      (15*kappa2*a3^2*a9)/(lam1 - lam2) - (45*kappa3*a4*a9)/(lam1 - lam2) - (30*kappa2*a2*a4*a9)/(lam1 - lam2) - ...
      (30*kappa2*a5*a9)/(lam1 - lam2) - (45*kappa3*a2^2*a10)/(lam1 - lam2) - (45*kappa3*a3*a10)/(lam1 - lam2) - ...
      (30*kappa2*a2*a3*a10)/(lam1 - lam2) - (30*kappa2*a4*a10)/(lam1 - lam2) - (45*kappa3*a2*a11)/(lam1 - lam2) - ...
      (15*kappa2*a2^2*a11)/(lam1 - lam2) - (30*kappa2*a3*a11)/(lam1 - lam2) - (15*kappa3*a12)/(lam1 - lam2) - ...
      (30*kappa2*a2*a12)/(lam1 - lam2) - (15*kappa2*a13)/(lam1 - lam2));
      
a15=(1/(-lam1 + 15*lam2))*(-((4*kappa3*a4^4)/(lam1 - lam2)) - (48*kappa3*a3*a4^2*a5)/(lam1 - lam2) - ...
      (24*kappa3*a3^2*a5^2)/(lam1 - lam2) - (48*kappa3*a2*a4*a5^2)/(lam1 - lam2) - (16*kappa3*a5^3)/(lam1 - lam2) - ...
      (48*kappa3*a3^2*a4*a6)/(lam1 - lam2) - (48*kappa3*a2*a4^2*a6)/(lam1 - lam2) - (96*kappa3*a2*a3*a5*a6)/(lam1 - lam2) - ...
      (96*kappa3*a4*a5*a6)/(lam1 - lam2) - (16*kappa2*a5^2*a6)/(lam1 - lam2) - (24*kappa3*a2^2*a6^2)/(lam1 - lam2) - ...
      (48*kappa3*a3*a6^2)/(lam1 - lam2) - (16*kappa2*a4*a6^2)/(lam1 - lam2) - (16*kappa3*a3^3*a7)/(lam1 - lam2) - ...
      (96*kappa3*a2*a3*a4*a7)/(lam1 - lam2) - (48*kappa3*a4^2*a7)/(lam1 - lam2) - (48*kappa3*a2^2*a5*a7)/(lam1 - lam2) - ...
      (96*kappa3*a3*a5*a7)/(lam1 - lam2) - (32*kappa2*a4*a5*a7)/(lam1 - lam2) - (96*kappa3*a2*a6*a7)/(lam1 - lam2) - ...
      (32*kappa2*a3*a6*a7)/(lam1 - lam2) - (24*kappa3*a7^2)/(lam1 - lam2) - (16*kappa2*a2*a7^2)/(lam1 - lam2) - ...
      (48*kappa3*a2*a3^2*a8)/(lam1 - lam2) - (48*kappa3*a2^2*a4*a8)/(lam1 - lam2) - (96*kappa3*a3*a4*a8)/(lam1 - lam2) - ...
      (16*kappa2*a4^2*a8)/(lam1 - lam2) - (96*kappa3*a2*a5*a8)/(lam1 - lam2) - (32*kappa2*a3*a5*a8)/(lam1 - lam2) - ...
      (48*kappa3*a6*a8)/(lam1 - lam2) - (32*kappa2*a2*a6*a8)/(lam1 - lam2) - (32*kappa2*a7*a8)/(lam1 - lam2) - ...
      (48*kappa3*a2^2*a3*a9)/(lam1 - lam2) - (48*kappa3*a3^2*a9)/(lam1 - lam2) - (96*kappa3*a2*a4*a9)/(lam1 - lam2) - ...
      (32*kappa2*a3*a4*a9)/(lam1 - lam2) - (48*kappa3*a5*a9)/(lam1 - lam2) - (32*kappa2*a2*a5*a9)/(lam1 - lam2) - ...
      (32*kappa2*a6*a9)/(lam1 - lam2) - (16*kappa3*a2^3*a10)/(lam1 - lam2) - (96*kappa3*a2*a3*a10)/(lam1 - lam2) - ...
      (16*kappa2*a3^2*a10)/(lam1 - lam2) - (48*kappa3*a4*a10)/(lam1 - lam2) - (32*kappa2*a2*a4*a10)/(lam1 - lam2) - ...
      (32*kappa2*a5*a10)/(lam1 - lam2) - (48*kappa3*a2^2*a11)/(lam1 - lam2) - (48*kappa3*a3*a11)/(lam1 - lam2) - ...
      (32*kappa2*a2*a3*a11)/(lam1 - lam2) - (32*kappa2*a4*a11)/(lam1 - lam2) - (48*kappa3*a2*a12)/(lam1 - lam2) - ...
      (16*kappa2*a2^2*a12)/(lam1 - lam2) - (32*kappa2*a3*a12)/(lam1 - lam2) - (16*kappa3*a13)/(lam1 - lam2) - ...
      (32*kappa2*a2*a13)/(lam1 - lam2) - (16*kappa2*a14)/(lam1 - lam2));
  
 
p2=linspace(-1,1,100);
 p1=(a2.*p2.^2+a3.*p2.^3+a4.*p2.^4+a5.*p2.^5+a6.*p2.^6+a7.*p2.^7+a8.*p2.^8+a9.*p2.^9+a10.*p2.^10+a11.*p2.^11+a12.*p2.^12+a13.*p2.^13+a14.*p2.^14+a15.*p2.^15);
 SSM2=T*[p1;p2];  
 
  a2=kappa2/((lam1 - lam2)*(2*lam1 - lam2));
 
 a3= (kappa3 + 4*kappa2*a2)/((lam1 - lam2)*(3*lam1 - lam2));
 
 a4=(5*(kappa3*a2 + kappa2*a2^2 + kappa2*a3))/((lam1 - lam2)*(4*lam1 - lam2));
 
 a5=((9*kappa3*a2^2)/(lam1 - lam2) + (2*kappa2*a2^3)/(lam1 - lam2) + (6*kappa3*a3)/(lam1 - lam2) + (12*kappa2*a2*a3)/(lam1 - lam2) + ...
      (6*kappa2*a4)/(lam1 - lam2))/(5*lam1 - lam2);
      
 a6=((7*kappa3*a2^3)/(lam1 - lam2) + (21*kappa3*a2*a3)/(lam1 - lam2) + (7*kappa2*a2^2*a3)/(lam1 - lam2) + (7*kappa2*a3^2)/(lam1 - lam2) + ...
      (7*kappa3*a4)/(lam1 - lam2) + (14*kappa2*a2*a4)/(lam1 - lam2) + (7*kappa2*a5)/(lam1 - lam2))/(6*lam1 - lam2);
      
 a7=((2*kappa3*a2^4)/(lam1 - lam2) + (24*kappa3*a2^2*a3)/(lam1 - lam2) + (12*kappa3*a3^2)/(lam1 - lam2) + ...
      (8*kappa2*a2*a3^2)/(lam1 - lam2) + (24*kappa3*a2*a4)/(lam1 - lam2) + (8*kappa2*a2^2*a4)/(lam1 - lam2) + ...
      (16*kappa2*a3*a4)/(lam1 - lam2) + (8*kappa3*a5)/(lam1 - lam2) + (16*kappa2*a2*a5)/(lam1 - lam2) + (8*kappa2*a6)/(lam1 - lam2))/...
     (7*lam1 - lam2);
     
 a8=(1/(8*lam1 - lam2))*((9*kappa3*a2^3*a3)/(lam1 - lam2) + (27*kappa3*a2*a3^2)/(lam1 - lam2) + (3*kappa2*a3^3)/(lam1 - lam2) + ...
      (27*kappa3*a2^2*a4)/(lam1 - lam2) + (27*kappa3*a3*a4)/(lam1 - lam2) + (18*kappa2*a2*a3*a4)/(lam1 - lam2) + ...
      (9*kappa2*a4^2)/(lam1 - lam2) + (27*kappa3*a2*a5)/(lam1 - lam2) + (9*kappa2*a2^2*a5)/(lam1 - lam2) + (18*kappa2*a3*a5)/(lam1 - lam2) + ...
      (9*kappa3*a6)/(lam1 - lam2) + (18*kappa2*a2*a6)/(lam1 - lam2) + (9*kappa2*a7)/(lam1 - lam2));
      
 a9=(1/(9*lam1 - lam2))*((15*kappa3*a2^2*a3^2)/(lam1 - lam2) + (10*kappa3*a3^3)/(lam1 - lam2) + (10*kappa3*a2^3*a4)/(lam1 - lam2) + ...
      (60*kappa3*a2*a3*a4)/(lam1 - lam2) + (10*kappa2*a3^2*a4)/(lam1 - lam2) + (15*kappa3*a4^2)/(lam1 - lam2) + ...
      (10*kappa2*a2*a4^2)/(lam1 - lam2) + (30*kappa3*a2^2*a5)/(lam1 - lam2) + (30*kappa3*a3*a5)/(lam1 - lam2) + ...
      (20*kappa2*a2*a3*a5)/(lam1 - lam2) + (20*kappa2*a4*a5)/(lam1 - lam2) + (30*kappa3*a2*a6)/(lam1 - lam2) + ...
      (10*kappa2*a2^2*a6)/(lam1 - lam2) + (20*kappa2*a3*a6)/(lam1 - lam2) + (10*kappa3*a7)/(lam1 - lam2) + (20*kappa2*a2*a7)/(lam1 - lam2) + ...
      (10*kappa2*a8)/(lam1 - lam2));
      
      
 a10= (1/(10*lam1 - lam2))*((11*kappa3*a2*a3^3)/(lam1 - lam2) + (33*kappa3*a2^2*a3*a4)/(lam1 - lam2) + ...
      (33*kappa3*a3^2*a4)/(lam1 - lam2) + (33*kappa3*a2*a4^2)/(lam1 - lam2) + (11*kappa2*a3*a4^2)/(lam1 - lam2) + ...
      (11*kappa3*a2^3*a5)/(lam1 - lam2) + (66*kappa3*a2*a3*a5)/(lam1 - lam2) + (11*kappa2*a3^2*a5)/(lam1 - lam2) +... 
      (33*kappa3*a4*a5)/(lam1 - lam2) + (22*kappa2*a2*a4*a5)/(lam1 - lam2) + (11*kappa2*a5^2)/(lam1 - lam2) + ...
      (33*kappa3*a2^2*a6)/(lam1 - lam2) + (33*kappa3*a3*a6)/(lam1 - lam2) + (22*kappa2*a2*a3*a6)/(lam1 - lam2) + ...
      (22*kappa2*a4*a6)/(lam1 - lam2) + (33*kappa3*a2*a7)/(lam1 - lam2) + (11*kappa2*a2^2*a7)/(lam1 - lam2) + ...
      (22*kappa2*a3*a7)/(lam1 - lam2) + (11*kappa3*a8)/(lam1 - lam2) + (22*kappa2*a2*a8)/(lam1 - lam2) + (11*kappa2*a9)/(lam1 - lam2));
      
 a11=(1/(11*lam1 - lam2))*((3*kappa3*a3^4)/(lam1 - lam2) + (36*kappa3*a2*a3^2*a4)/(lam1 - lam2) + (18*kappa3*a2^2*a4^2)/(lam1 - lam2) + ...
      (36*kappa3*a3*a4^2)/(lam1 - lam2) + (4*kappa2*a4^3)/(lam1 - lam2) + (36*kappa3*a2^2*a3*a5)/(lam1 - lam2) + ...
      (36*kappa3*a3^2*a5)/(lam1 - lam2) + (72*kappa3*a2*a4*a5)/(lam1 - lam2) + (24*kappa2*a3*a4*a5)/(lam1 - lam2) + ...
      (18*kappa3*a5^2)/(lam1 - lam2) + (12*kappa2*a2*a5^2)/(lam1 - lam2) + (12*kappa3*a2^3*a6)/(lam1 - lam2) + ...
      (72*kappa3*a2*a3*a6)/(lam1 - lam2) + (12*kappa2*a3^2*a6)/(lam1 - lam2) + (36*kappa3*a4*a6)/(lam1 - lam2) + ...
      (24*kappa2*a2*a4*a6)/(lam1 - lam2) + (24*kappa2*a5*a6)/(lam1 - lam2) + (36*kappa3*a2^2*a7)/(lam1 - lam2) + ...
      (36*kappa3*a3*a7)/(lam1 - lam2) + (24*kappa2*a2*a3*a7)/(lam1 - lam2) + (24*kappa2*a4*a7)/(lam1 - lam2) + ...
      (36*kappa3*a2*a8)/(lam1 - lam2) + (12*kappa2*a2^2*a8)/(lam1 - lam2) + (24*kappa2*a3*a8)/(lam1 - lam2) + (12*kappa3*a9)/(lam1 - lam2) + ...
      (24*kappa2*a2*a9)/(lam1 - lam2) + (12*kappa2*a10)/(lam1 - lam2));
      
 a12=(1/(12*lam1 - lam2))*((13*kappa3*a3^3*a4)/(lam1 - lam2) + (39*kappa3*a2*a3*a4^2)/(lam1 - lam2) + (13*kappa3*a4^3)/(lam1 - lam2) + ...
      (39*kappa3*a2*a3^2*a5)/(lam1 - lam2) + (39*kappa3*a2^2*a4*a5)/(lam1 - lam2) + (78*kappa3*a3*a4*a5)/(lam1 - lam2) + ...
      (13*kappa2*a4^2*a5)/(lam1 - lam2) + (39*kappa3*a2*a5^2)/(lam1 - lam2) + (13*kappa2*a3*a5^2)/(lam1 - lam2) + ...
      (39*kappa3*a2^2*a3*a6)/(lam1 - lam2) + (39*kappa3*a3^2*a6)/(lam1 - lam2) + (78*kappa3*a2*a4*a6)/(lam1 - lam2) + ...
      (26*kappa2*a3*a4*a6)/(lam1 - lam2) + (39*kappa3*a5*a6)/(lam1 - lam2) + (26*kappa2*a2*a5*a6)/(lam1 - lam2) + ...
      (13*kappa2*a6^2)/(lam1 - lam2) + (13*kappa3*a2^3*a7)/(lam1 - lam2) + (78*kappa3*a2*a3*a7)/(lam1 - lam2) + ...
      (13*kappa2*a3^2*a7)/(lam1 - lam2) + (39*kappa3*a4*a7)/(lam1 - lam2) + (26*kappa2*a2*a4*a7)/(lam1 - lam2) + ...
      (26*kappa2*a5*a7)/(lam1 - lam2) + (39*kappa3*a2^2*a8)/(lam1 - lam2) + (39*kappa3*a3*a8)/(lam1 - lam2) + ...
      (26*kappa2*a2*a3*a8)/(lam1 - lam2) + (26*kappa2*a4*a8)/(lam1 - lam2) + (39*kappa3*a2*a9)/(lam1 - lam2) +... 
      (13*kappa2*a2^2*a9)/(lam1 - lam2) + (26*kappa2*a3*a9)/(lam1 - lam2) + (13*kappa3*a10)/(lam1 - lam2) + ...
      (26*kappa2*a2*a10)/(lam1 - lam2) + (13*kappa2*a11)/(lam1 - lam2));
      
 a13= (1/(13*lam1 - lam2))*((21*kappa3*a3^2*a4^2)/(lam1 - lam2) + (14*kappa3*a2*a4^3)/(lam1 - lam2) + (14*kappa3*a3^3*a5)/(lam1 - lam2) + ...
      (84*kappa3*a2*a3*a4*a5)/(lam1 - lam2) + (42*kappa3*a4^2*a5)/(lam1 - lam2) + (21*kappa3*a2^2*a5^2)/(lam1 - lam2) + ...
      (42*kappa3*a3*a5^2)/(lam1 - lam2) + (14*kappa2*a4*a5^2)/(lam1 - lam2) + (42*kappa3*a2*a3^2*a6)/(lam1 - lam2) + ...
      (42*kappa3*a2^2*a4*a6)/(lam1 - lam2) + (84*kappa3*a3*a4*a6)/(lam1 - lam2) + (14*kappa2*a4^2*a6)/(lam1 - lam2) + ...
      (84*kappa3*a2*a5*a6)/(lam1 - lam2) + (28*kappa2*a3*a5*a6)/(lam1 - lam2) + (21*kappa3*a6^2)/(lam1 - lam2) + ...
      (14*kappa2*a2*a6^2)/(lam1 - lam2) + (42*kappa3*a2^2*a3*a7)/(lam1 - lam2) + (42*kappa3*a3^2*a7)/(lam1 - lam2) + ...
      (84*kappa3*a2*a4*a7)/(lam1 - lam2) + (28*kappa2*a3*a4*a7)/(lam1 - lam2) + (42*kappa3*a5*a7)/(lam1 - lam2) + ...
      (28*kappa2*a2*a5*a7)/(lam1 - lam2) + (28*kappa2*a6*a7)/(lam1 - lam2) + (14*kappa3*a2^3*a8)/(lam1 - lam2) + ...
      (84*kappa3*a2*a3*a8)/(lam1 - lam2) + (14*kappa2*a3^2*a8)/(lam1 - lam2) + (42*kappa3*a4*a8)/(lam1 - lam2) + ...
      (28*kappa2*a2*a4*a8)/(lam1 - lam2) + (28*kappa2*a5*a8)/(lam1 - lam2) + (42*kappa3*a2^2*a9)/(lam1 - lam2) + ...
      (42*kappa3*a3*a9)/(lam1 - lam2) + (28*kappa2*a2*a3*a9)/(lam1 - lam2) + (28*kappa2*a4*a9)/(lam1 - lam2) + ...
      (42*kappa3*a2*a10)/(lam1 - lam2) + (14*kappa2*a2^2*a10)/(lam1 - lam2) + (28*kappa2*a3*a10)/(lam1 - lam2) + ...
      (14*kappa3*a11)/(lam1 - lam2) + (28*kappa2*a2*a11)/(lam1 - lam2) + (14*kappa2*a12)/(lam1 - lam2));
      
 a14=(1/(14*lam1 - lam2))*((15*kappa3*a3*a4^3)/(lam1 - lam2) + (45*kappa3*a3^2*a4*a5)/(lam1 - lam2) + ...
      (45*kappa3*a2*a4^2*a5)/(lam1 - lam2) + (45*kappa3*a2*a3*a5^2)/(lam1 - lam2) + (45*kappa3*a4*a5^2)/(lam1 - lam2) + ...
      (5*kappa2*a5^3)/(lam1 - lam2) + (15*kappa3*a3^3*a6)/(lam1 - lam2) + (90*kappa3*a2*a3*a4*a6)/(lam1 - lam2) + ...
      (45*kappa3*a4^2*a6)/(lam1 - lam2) + (45*kappa3*a2^2*a5*a6)/(lam1 - lam2) + (90*kappa3*a3*a5*a6)/(lam1 - lam2) + ...
      (30*kappa2*a4*a5*a6)/(lam1 - lam2) + (45*kappa3*a2*a6^2)/(lam1 - lam2) + (15*kappa2*a3*a6^2)/(lam1 - lam2) + ...
      (45*kappa3*a2*a3^2*a7)/(lam1 - lam2) + (45*kappa3*a2^2*a4*a7)/(lam1 - lam2) + (90*kappa3*a3*a4*a7)/(lam1 - lam2) + ...
      (15*kappa2*a4^2*a7)/(lam1 - lam2) + (90*kappa3*a2*a5*a7)/(lam1 - lam2) + (30*kappa2*a3*a5*a7)/(lam1 - lam2) + ...
      (45*kappa3*a6*a7)/(lam1 - lam2) + (30*kappa2*a2*a6*a7)/(lam1 - lam2) + (15*kappa2*a7^2)/(lam1 - lam2) + ...
      (45*kappa3*a2^2*a3*a8)/(lam1 - lam2) + (45*kappa3*a3^2*a8)/(lam1 - lam2) + (90*kappa3*a2*a4*a8)/(lam1 - lam2) + ...
      (30*kappa2*a3*a4*a8)/(lam1 - lam2) + (45*kappa3*a5*a8)/(lam1 - lam2) + (30*kappa2*a2*a5*a8)/(lam1 - lam2) + ...
      (30*kappa2*a6*a8)/(lam1 - lam2) + (15*kappa3*a2^3*a9)/(lam1 - lam2) + (90*kappa3*a2*a3*a9)/(lam1 - lam2) + ...
      (15*kappa2*a3^2*a9)/(lam1 - lam2) + (45*kappa3*a4*a9)/(lam1 - lam2) + (30*kappa2*a2*a4*a9)/(lam1 - lam2) + ...
      (30*kappa2*a5*a9)/(lam1 - lam2) + (45*kappa3*a2^2*a10)/(lam1 - lam2) + (45*kappa3*a3*a10)/(lam1 - lam2) + ...
      (30*kappa2*a2*a3*a10)/(lam1 - lam2) + (30*kappa2*a4*a10)/(lam1 - lam2) + (45*kappa3*a2*a11)/(lam1 - lam2) + ...
      (15*kappa2*a2^2*a11)/(lam1 - lam2) + (30*kappa2*a3*a11)/(lam1 - lam2) + (15*kappa3*a12)/(lam1 - lam2) + ...
      (30*kappa2*a2*a12)/(lam1 - lam2) + (15*kappa2*a13)/(lam1 - lam2));
      
 a15=(1/(15*lam1 - lam2))*((4*kappa3*a4^4)/(lam1 - lam2) + (48*kappa3*a3*a4^2*a5)/(lam1 - lam2) + (24*kappa3*a3^2*a5^2)/(lam1 - lam2) + ...
      (48*kappa3*a2*a4*a5^2)/(lam1 - lam2) + (16*kappa3*a5^3)/(lam1 - lam2) + (48*kappa3*a3^2*a4*a6)/(lam1 - lam2) + ...
      (48*kappa3*a2*a4^2*a6)/(lam1 - lam2) + (96*kappa3*a2*a3*a5*a6)/(lam1 - lam2) + (96*kappa3*a4*a5*a6)/(lam1 - lam2) + ...
      (16*kappa2*a5^2*a6)/(lam1 - lam2) + (24*kappa3*a2^2*a6^2)/(lam1 - lam2) + (48*kappa3*a3*a6^2)/(lam1 - lam2) + ...
      (16*kappa2*a4*a6^2)/(lam1 - lam2) + (16*kappa3*a3^3*a7)/(lam1 - lam2) + (96*kappa3*a2*a3*a4*a7)/(lam1 - lam2) + ...
      (48*kappa3*a4^2*a7)/(lam1 - lam2) + (48*kappa3*a2^2*a5*a7)/(lam1 - lam2) + (96*kappa3*a3*a5*a7)/(lam1 - lam2) + ...
      (32*kappa2*a4*a5*a7)/(lam1 - lam2) + (96*kappa3*a2*a6*a7)/(lam1 - lam2) + (32*kappa2*a3*a6*a7)/(lam1 - lam2) + ...
      (24*kappa3*a7^2)/(lam1 - lam2) + (16*kappa2*a2*a7^2)/(lam1 - lam2) + (48*kappa3*a2*a3^2*a8)/(lam1 - lam2) + ...
      (48*kappa3*a2^2*a4*a8)/(lam1 - lam2) + (96*kappa3*a3*a4*a8)/(lam1 - lam2) + (16*kappa2*a4^2*a8)/(lam1 - lam2) + ...
      (96*kappa3*a2*a5*a8)/(lam1 - lam2) + (32*kappa2*a3*a5*a8)/(lam1 - lam2) + (48*kappa3*a6*a8)/(lam1 - lam2) + ...
      (32*kappa2*a2*a6*a8)/(lam1 - lam2) + (32*kappa2*a7*a8)/(lam1 - lam2) + (48*kappa3*a2^2*a3*a9)/(lam1 - lam2) + ...
      (48*kappa3*a3^2*a9)/(lam1 - lam2) + (96*kappa3*a2*a4*a9)/(lam1 - lam2) + (32*kappa2*a3*a4*a9)/(lam1 - lam2) + ...
      (48*kappa3*a5*a9)/(lam1 - lam2) + (32*kappa2*a2*a5*a9)/(lam1 - lam2) + (32*kappa2*a6*a9)/(lam1 - lam2) + ...
      (16*kappa3*a2^3*a10)/(lam1 - lam2) + (96*kappa3*a2*a3*a10)/(lam1 - lam2) + (16*kappa2*a3^2*a10)/(lam1 - lam2) + ...
      (48*kappa3*a4*a10)/(lam1 - lam2) + (32*kappa2*a2*a4*a10)/(lam1 - lam2) + (32*kappa2*a5*a10)/(lam1 - lam2) + ...
      (48*kappa3*a2^2*a11)/(lam1 - lam2) + (48*kappa3*a3*a11)/(lam1 - lam2) + (32*kappa2*a2*a3*a11)/(lam1 - lam2) + ...
      (32*kappa2*a4*a11)/(lam1 - lam2) + (48*kappa3*a2*a12)/(lam1 - lam2) + (16*kappa2*a2^2*a12)/(lam1 - lam2) + ...
      (32*kappa2*a3*a12)/(lam1 - lam2) + (16*kappa3*a13)/(lam1 - lam2) + (32*kappa2*a2*a13)/(lam1 - lam2) + (16*kappa2*a14)/(lam1 - lam2));
 
 
 
 
 p1=linspace(-1,1,100);
 p2=(a2.*p1.^2+a3.*p1.^3+a4.*p1.^4+a5.*p1.^5+a6.*p1.^6+a7.*p1.^7++a8.*p1.^8+a9.*p1.^9+a10.*p1.^10+a11.*p1.^11+a12.*p1.^12+a13.*p1.^13+a14.*p1.^14+a15.*p1.^15);
 
 SSM1=T*[p1;p2];

end

 
