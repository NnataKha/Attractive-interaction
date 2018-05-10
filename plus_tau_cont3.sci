clear;  

// function p_i
function y=fun1(x)
  y=sin(4*x)+1;
endfunction

// function r_i
function y=fun2(x)
  y=cos(4*x)+1;
endfunction

cp=intg(0,1,fun1);
cr=intg(0,1,fun2);
function y=fun11(x)
  y=fun1(x)/cp;
endfunction
function y=fun21(x)
  y=fun2(x)/cr;
endfunction

d=0.001;//the length of the intervals of  partitioning at the graphic drawn with dash line
m=1+1/d;
x=0;
for i=1:m
  p(i)=fun11(x);
  r(i)=fun21(x);
  T(i)=x;
  x=x+d;
end
p1=p;
r1=r;
tau=min(p,r);
ctau=sum(tau.*d);
tau1=tau./ctau;
p_v(1,:)=p';
r_v(1,:)=r';


N=100;// number of steps (graphic drawn with dash line)
for k=1:N
    Omega_plus=list();
    Omega_minus=list();
    for i=1:m
       if p(i) > r(i) then
           Omega_plus($+1)=i;
       else Omega_minus($+1)=i
       end
    end
    op=length(Omega_plus);
    om=length(Omega_minus);
  Theta = sum(sqrt(p.*r)*d);
  tau=min(p,r);
  W = sum(tau.*d');
  z = 1+W+Theta;
  p_t = (p.*(1+Theta)+tau)./z;
  r_t = (r.*(1+Theta)+tau)./z;
  
  dif=abs(p-r);
  sd=sum(dif);
  dif=dif./sd;
  for i=1:op
      r_t(Omega_plus(i))=r_t(Omega_plus(i))*(1-dif(Omega_plus(i)));
  end
  for i=1:om
      p_t(Omega_minus(i))=p_t(Omega_minus(i))*(1-dif(Omega_minus(i)));
  end
  
  sp=sum(p_t*d);
  sr=sum(r_t*d);
  p=p_t./sp
  r=r_t./sr
  
  p_v(k+1,:)=p';
  r_v(k+1,:)=r';
end
//    clf();
    //plot(T,p1,'b');
    plot(T,r_v(1,:),'g--');
    plot(T,p_v(1,:),'r--');
    plot(T,r_v(k+1,:),'g');
    plot(T,p_v(k+1,:),'r');
//    plot(T,r_v','g');
//    plot(T,p_v','r');
    plot(T,tau1,'c--');
    //plot(T,p,'b');
    //plot(T,r,'r');
    //legend('p(x)','r(x)','tau1','p^inf(x)','r^inf(x)');
