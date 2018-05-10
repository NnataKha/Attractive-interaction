clear;
if 31==1 then
n=5;
for i=1:n-2
    T(i)=rand(100)//*110;
    M(i)=rand(100);
end
for i=n-1:n
    T(i)=rand(100);
    M(i)=rand(100)//*50;
end
end
M=[1;2;3;4];
T=[2;5;1;6];

n=length(M);
st = sum(T);
sm = sum(M);
T = T'./st;
M = M'./sm;
eps=1/(n*n);

        
    Omega_plus=list();
    Omega_minus=list();
    for i=1:n
       if T(i) > M(i) then
           Omega_plus($+1)=i;
       else Omega_minus($+1)=i
       end
    end
    op=length(Omega_plus);
    om=length(Omega_minus);

    alpha(1)=1;

m = 500;//number of steps
for k = 1:m
    time(k)=k;
    T_v(k,:)=T;
    M_v(k,:)=M;
    Theta(k)=0;
    for i=1:n
        Theta(k) = Theta(k) + sqrt(T(i)*M(i));
    end
    tau = min(T,M);
    W = sum(tau);
    alpha(k+1)=-alpha(k)*(Theta(k)-eps)*(Theta(k)-1+eps)/abs((Theta(k)-eps)*(Theta(k)-1+eps));
    z = 1+Theta(k)+alpha(k)*W;
    T_t = (T.*(1+Theta(k))+alpha(k)*tau)./z;
    M_t = (M.*(1+Theta(k))+alpha(k)*tau)./z;
    if alpha(k)>0 then
        s=Theta(k);
        h=Theta(k);
    else 
        s=1-Theta(k);
        h=1-Theta(k);
    end
    

//*****************************************************************************
//redistribution for M
    d=zeros(1,n);
    for i=1:op
        d(Omega_plus(i))=(M_t(Omega_plus(i))-M(Omega_plus(i)))*s;
    end
    M_t=M_t-d;
//*****************************************************************************
//redistribution for T
    d=zeros(1,n);
    for i=1:om
        d(Omega_minus(i))=(T_t(Omega_minus(i))-T(Omega_minus(i)))*h;
    end
    T_t=T_t-d;
//*****************************************************************************
    sm=sum(M_t);
    st=sum(T_t);
    M=M_t./sm
    T=T_t./st;
end
subplot(211)
plot(time,T_v);
//legend('T(1)','T(2)','T(3)','T(4)','T(5)','T(6)','T(7)','T(8)');

subplot(212)
plot(time,M_v);
//legend('M(1)','M(2)','M(3)','M(4)','M(5)','M(6)','M(7)','M(8)');
//xs2pdf(gcf(),"C:\Users\mika\Desktop\imath\article\pictures\pic.eps");

