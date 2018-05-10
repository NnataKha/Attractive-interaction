clear;
if 1==1 then
n=6;
for i=1:n
    T(i)=rand(100);
    M(i)=rand(100);
end

end
//M=[0.07669220464090441;0.32497903823139684;0.1543080013234991;0.0774208472828813;0.004630362501208455;0.36196954602010994];
//T=[0.038548173413248454;0.04549940240263871;0.2224250336377795;0.23337264019335863;0.28981549496797365;0.17033925538500108];
//M=[1;2;3;4];
//T=[31;45;1;.1];
//n=length(M);
m = 150;//number of steps


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
        
for k = 1:m

    time(k)=k;
    for i=1:op
        T_plus(k,i)=T(Omega_plus(i))
        M_plus(k,i)=M(Omega_plus(i))
    end
    
    T_v(k,:)=T;
    M_v(k,:)=M;
    Theta = sum(sqrt(T.*M));
    tau = min(T,M);
    W = sum(tau);
    z = 1+Theta+W;
    T_t = (T.*(1+Theta)+tau)./z;
    M_t = (M.*(1+Theta)+tau)./z;
    Tau_plus(k)=0;
    Tau_minus(k)=0;
    for i=1:op
        Tau_plus(k)=Tau_plus(k)+tau(Omega_plus(i))
    end
    for i=1:om
        Tau_minus(k)=Tau_minus(k)+tau(Omega_minus(i))
    end
    check(k,:)=tau//.*Tau_plus(k)/Tau_minus(k)
    s=1;
    h=s;

//*****************************************************************************
//redistribution for M
    d=zeros(1,n);
    for i=1:op
//        d(Omega_plus(i))=(T(Omega_plus(i))-T_t(Omega_plus(i)))*s;
        d(Omega_plus(i))=(M_t(Omega_plus(i))-M(Omega_plus(i)))*s;
    end
    M_t=M_t+d;
//*****************************************************************************
//redistribution for T
    d=zeros(1,n);
    for i=1:om
//        d(Omega_minus(i))=(M(Omega_minus(i))-M_t(Omega_minus(i)))*h;
        d(Omega_minus(i))=(T_t(Omega_minus(i))-T(Omega_minus(i)))*h;
    end
    T_t=T_t+d;
//*****************************************************************************
    sm=sum(M_t);
    st=sum(T_t);
    M=M_t./sm
    T=T_t./st;
end
if 1==1 then
//subplot(211)
plot(time,check,'.-')
plot(time,T_v,'--');
plot(time,M_v);
//plot(time,M_plus);
//plot(time,T_plus,'--');
//subplot(212)
//plot(time,T_v-M_v)//,'.-');

else 
    subplot(211)
    plot(time,T_v,'.-');
    plot(time,M_v,'*-');
    subplot(212)
    plot(time,T_v-M_v,'.-');
end
