clear;  

// function that determines p_i
function y=fun1(x)
  y=sin(4*x*%pi)+1;
endfunction

// function that determines r_i
function y=fun2(x)
  y=cos(2*x*%pi)+1;
endfunction

// function that determines tau
function y=fun3(x)
  y=1/(x+1);
endfunction

cp=intg(0,1,fun1);
cr=intg(0,1,fun2);
ct=intg(0,1,fun3);
function y=fun11(x)
  y=fun1(x)/cp;
endfunction
function y=fun21(x)
  y=fun2(x)/cr;
endfunction
function y=fun31(x)
  y=fun3(x)/ct;
endfunction

d=0.01;//the length of the intervals of  partitioning at the graphic drawn with dash line
m=1+1/d;
x=0;
for i=1:m
  p(i)=fun11(x);
  r(i)=fun21(x);
  //tau(i)=fun31(x);
  T(i)=x;
  x=x+d;
end
p0=p;
r0=r;

//  W = sum(tau.*d);
  //tau1=r-tau./W;
  c=1;
//W = sum(tau.*d);
    W=1;

N=10000;// number of steps (graphic drawn with dash line)
for k=1:N
  Theta = cos(%pi*k/2)+1.1;//sum(sqrt(p.*r)*d);
  c=c*(1-W/(1+W+Theta));
  z = 2+Theta;
  //z = 1+W+Theta;
//  p = (p.*(1+Theta)+tau)./z;
//  r = (r.*(1+Theta)+tau)./z;
  p = (p.*(1+Theta)+r)./z;
  r = (r.*(1+Theta)+p)./z;

end
//tau1=tau1*c+tau./W;
if 6==1 then
    mu0(1)=p0(1)*d;
    nu0(1)=r0(1)*d;
    mu(1)=p(1)*d;
    nu(1)=r(1)*d;
    for i=2:m
        mu0(i)=mu0(i-1)+p0(i)*d;
        nu0(i)=nu0(i-1)+r0(i)*d;
        mu(i)=mu(i-1)+p(i)*d;
        nu(i)=nu(i-1)+r(i)*d;
    end

    plot(T,mu0,'b--');    
    plot(T,mu,'b.-');    
    plot(T,nu0,'black--');    
    plot(T,nu,'black');    
end
if 1==1 then
    plot(T,p,'b.-');    
    plot(T,p0,'b--');
    plot(T,r0,'black--');
    plot(T,r,'black.-');
    //plot(T,tau1,'r');
  //  legend('p(x)','r(x)','limit');
end
