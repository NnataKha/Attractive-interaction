clear;
if 1==1 then
n=3;
for i=1:n
    T(i)=rand(100);
    M(i)=rand(100);
end

end
//n=length(M);
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

s=sum(abs(T-M))/n;
h=s;
m = 300;//number of steps
c(1)=;
for k = 1:m
    time(k)=k;
    T_v(k,:)=T;
    M_v(k,:)=M;
    Theta = sum(sqrt(T.*M));
    delta=sum(abs(T-M))/n;
    Th(k)=delta;
    check(k+1)=check(k);

    //check the value of Theta
    if Theta<=eps then
        c(k+1)='plus';
    elseif Theta>=1-eps then
        c(k+1)='minus';
    end

    tau = min(T,M);
    W = sum(tau);
    z = 1+Theta+W;
    T_t = (T.*(1+Theta)+tau)./z;
    M_t = (M.*(1+Theta)+tau)./z;
    
    s=1;
    h=1;
    //choosing “+,+” or  “-,-”


//*****************************************************************************
//redistribution for M
    d=zeros(1,n);
    for i=1:op
        d(Omega_plus(i))=abs(T_t(Omega_plus(i))-T(Omega_plus(i)))*s;
    end
    T_t=T_t+d;
//*****************************************************************************
//redistribution for T
    d=zeros(1,n);
    for i=1:om
        d(Omega_minus(i))=abs(M_t(Omega_minus(i))-M(Omega_minus(i)))*h;
    end
    M_t=M_t+d;
//*****************************************************************************
    sm=sum(M_t);
    st=sum(T_t);
    M=M_t./sm
    T=T_t./st;
end
subplot(211)
plot(time,T_v);
plot(time,M_v,'--');
subplot(212)
plot(time,Th)
