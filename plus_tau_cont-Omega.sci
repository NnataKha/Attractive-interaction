clear;  

// function that determines p_i
function y=fun1(x)
  y=sin(5*x*%pi)+1;
endfunction

// function that determines r_i
function y=fun2(x)
  //y=cos(2*x*%pi)+1;
  y=(-0.2+x)^2;
endfunction

cp=intg(0,1,fun1);
cr=intg(0,1,fun2);
function y=fun11(x)
  y=fun1(x)/cp;
endfunction
function y=fun21(x)
  y=fun2(x)/cr;
endfunction

d=0.01;//the length of the intervals of  partitioning at the graphic drawn with dash line
m=1+1/d;
x=0;
for i=1:m
  p(i)=fun11(x);
  r(i)=fun21(x);
  T(i)=x;
  x=x+d;
end
p0=p;
r0=r;
po1=0.0/d;
po2=0.3/d;
ro1=0.4/d;
ro2=0.6/d;
po3=0.7/d;
po4=0.8/d;
ro3=0.9/d;
ro4=1/d;

N=1000;// number of steps (graphic drawn with dash line)
for k=1:N
//  Theta = sum(sqrt(p.*r)*d);
  

Theta1=0;
Theta2=0;
for i=1:m
    if (0<=i)&(i<po1) then
        tau(i)=0
    end    
    if (po1<=i)&(i<po2) then
        tau(i)=p(i);
        Theta1=Theta1*d+p(i)*d;
    end    
    if (po2<=i)&(i<ro1) then
        tau(i)=0
    end    
    if (ro1<=i)&(i<ro2) then
        tau(i)=r(i);
        Theta2=Theta2*d+r(i)*d;
    end    
    if (ro2<=i)&(i<po3) then
        tau(i)=0
    end    
    if (po3<=i)&(i<po4) then
        tau(i)=p(i);
        Theta1=Theta1*d+p(i)*d;
    end    
    if (po4<i)&(i<ro3) then
        tau(i)=0
    end    
    if (ro3<=i)&(i<ro4) then
        tau(i)=r(i);        
        Theta2=Theta2*d+r(i)*d;
    end    
    if ro4<=i then
        tau(i)=0;
    end 
end
    Theta=Theta1*Theta2;
  W = sum(tau.*d);
  z = 1+W+Theta;
  p = (p.*(1+Theta)+tau)./z;
  r = (r.*(1+Theta)+tau)./z;
end

if 1==1 then
    mu0(1)=p0(1)*d;
    nu0(1)=r0(1)*d;
    mu(1)=p(1)*d;
    nu(1)=r(1)*d;
    mtau(1)=tau(1)*d;
    for i=2:m
        mtau(i)=mtau(i-1)+tau(i)*d;
        mu0(i)=mu0(i-1)+p0(i)*d;
        nu0(i)=nu0(i-1)+r0(i)*d;
        mu(i)=mu(i-1)+p(i)*d;
        nu(i)=nu(i-1)+r(i)*d;
    end

    plot(T,mu0,'b--');    
    plot(T,mu,'b.-');    
    plot(T,nu0,'black--');    
    plot(T,nu,'black');    
    plot(T,mtau,'r');    
end
if 2==1 then
    plot(T,p,'b.-');    
    plot(T,p0,'b--');
    plot(T,r0,'black--');
    plot(T,r,'black');
    //plot(T,tau1,'*');
  //  legend('p(x)','r(x)','limit');
end
