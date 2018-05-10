clear;  

// function that determines p_i
function y=fun1(x)
  y=sin(4*x*%pi)+1;
endfunction

// function that determines r_i
function y=fun2(x)
  y=cos(4*x*%pi)+1;
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

//set B_mu
B_mu0=0.0/d;
B_mu1=0.7/d;

//set B_nu
B_nu0=0.3/d;
B_nu1=1/d;

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
alpha=0.5;
be=0.5;
//set tau
for i=1:m
    if (i>=B_mu0 & i<=B_mu1)|(i>=B_nu0 & i<=B_nu1) then
        tau(i)=alpha*p(i)+be*r(i)
        else tau(i)=0;
    end
end
ctau=sum(tau.*d);
tau1=tau./ctau;
p_v(1,:)=p';

N=100;// number of steps (graphic drawn with dash line)
for k=1:N
  Theta = sum(sqrt(p.*r)*d);
    for i=1:m
        if (i>=B_mu0 & i<=B_mu1)|(i>=B_nu0 & i<=B_nu1) then
            tau(i)=alpha*p(i)+be*r(i)
            else tau(i)=0;
        end
    end
  W = sum(tau.*d);
  z = 1+W+Theta;
  p = (p.*(1+Theta)+tau)./z;
  r = (r.*(1+Theta)+tau)./z;
  //p_v(k+1,:)=r';
end

    plot(T,p,'g');    
    plot(T,p1,'g--');
    plot(T,r1,'r--');
    plot(T,r,'r');
