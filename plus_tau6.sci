//we subtract avarage grouth on Omega_+ from each M(i)
clear;
if 1==21 then
    n=4
for i=1:n
    T(i)=rand(100);
    M(i)=rand(100);
end
end
n=6;
eps=2/(n*n);
M=[3;4;1;2;3;4];
T=[1;2;2;8;4;5];
//n=length(M);
st = sum(T);
sm = sum(M);
T = T'./st;
M = M'./sm;
check(1)=1;
time1(1)=0;
m = 400;
for k = 1:m
    time1(k+1)=k
    time(k)=k;
    T_v(k,:)=T//(6);
    M_v(k,:)=M//(6);
    Theta=sum(sqrt(T.*M));
    tau = min(T,M);
    W = sum(tau);
    Omega_plus=list();
    Omega_minus=list();
    for i=1:n
       if tau(i) == M(i) then
           Omega_plus($+1)=i;
       else Omega_minus($+1)=i
       end
    end
    z = 1+Theta+W;
    T_t = (T.*(1+Theta)+tau)./z;
    M_t = (M.*(1+Theta)+tau)./z;
    op=length(Omega_plus);
    om=length(Omega_minus);
    check(k+1)=check(k);
    delta(k)=sum(abs(T_t-M_t))/n;
//    if delta(k)<=eps then
    if Theta<=eps then
        check(k+1)=-1;
    end
//    if delta(k)>=1-eps then
    if Theta>=1-eps then
        check(k+1)=1;
    end
    wmin(k)=W;
    th(k)=Theta;
    c=1//^(check(k))
    h=1+Theta*check(k)//+delta(k);
    s=h;
//*****************************************************************************
//redistribution for M
    d=zeros(1,n);
    dif(k,:)=zeros(1,n);
    for i=1:op
//        d(Omega_plus(i))=(M(Omega_plus(i))-M_t(Omega_plus(i)))*s*check(k+1);
        d(Omega_plus(i))=(T(Omega_plus(i))-M(Omega_plus(i)))*s//*check(k+1);
    end
    dif(k,:)=d;
    T_t=T_t+d;
//*****************************************************************************
//redistribution for T
    d=zeros(1,n);
    for i=1:om
//        d(Omega_minus(i))=(T(Omega_minus(i))-T_t(Omega_minus(i)))*h*check(k+1);
        d(Omega_minus(i))=(M(Omega_minus(i))-T(Omega_minus(i)))*h//*check(k+1);
        dif(k,Omega_minus(i))=d(Omega_minus(i));
    end
    M_t=M_t+d;
//*****************************************************************************
    M=M_t./sum(M_t);
    T=T_t./sum(T_t);
end

subplot(311)
plot(time,T_v);
plot(time,M_v,'--');
subplot(312)
plot(time,th);
plot(time1,check)
subplot(313)
plot(time,dif)
