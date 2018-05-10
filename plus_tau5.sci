//we subtract avarage grouth on Omega_+ from each M(i)
clear;
if 1==1 then
n=4;
for i=1:n-2
    T(i)=rand(100)*110;
    M(i)=rand(100);
end
for i=n-1:n
    T(i)=rand(100);
    M(i)=rand(100)*50;
end
end
n=length(M);
st = sum(T);
sm = sum(M);
T = T'./st;
M = M'./sm;

m = 200;
for k = 1:m
    time(k)=k;
    T_v(k,:)=T//(6);
    M_v(k,:)=M//(6);
    Theta = sum(T.*M);
    W = 0;
    Omega_plus=list();
    Omega_minus=list();
    M_N=0;
    T_N=0;
    for i=1:n
       tau(i) = min(T(i),M(i));
       W = W + tau(i);
       if tau(i) == M(i) then
           Omega_plus($+1)=i;
       else Omega_minus($+1)=i
       end
    end

    z = 1+Theta+W;
    
    T_t = (T.*(1+Theta)+tau')./z;
    M_t = (M.*(1+Theta)+tau')./z;
    op=length(Omega_plus);
    om=length(Omega_minus);
if modulo(k,1)==0 then
//*****************************************************************************
//redistribution for M

    d=0;
    for i=1:op
        d=d+abs(M(Omega_plus(i))-M_t(Omega_plus(i)));
    end
    d=d/op;
    dd=0;
    for i=1:op
        M_t(Omega_plus(i))=M_t(Omega_plus(i))-d;
        if M_t(Omega_plus(i))<0 then
            M_t(Omega_plus(i))=d;
        end
        dd=dd+ M_t(Omega_plus(i));
     end
     s=0;
     for i=1:om
         s=s+M_t(Omega_minus(i));
     end
     d_omega=zeros(n);
     for i=1:om
         d_omega(Omega_minus(i))=M_t(Omega_minus(i))/s;
     end
     
     for i=1:om
         M_t(Omega_minus(i))=M_t(Omega_minus(i))+(1-d_omega(Omega_minus(i)))*dd;
     end

//*****************************************************************************
//redistribution for T

    d=0;
    for i=1:om
        d=d+abs(T(Omega_minus(i))-T_t(Omega_minus(i)));
    end
    d=d/om;
if 1==2 then
    d=abs(M(Omega_plus(1))-M_t(Omega_plus(1)));
    for i=1:op
        m=abs(M(Omega_plus(i))-M_t(Omega_plus(i)));
        if d<m then
            d=m;
        end
    end
end
    dd=0;
    for i=1:om
        T_t(Omega_minus(i))=T_t(Omega_minus(i))-d;
        if T_t(Omega_minus(i))<0 then
            T_t(Omega_minus(i))=d;
        end
        dd=dd+ T_t(Omega_minus(i));
     end
     s=0;
     for i=1:op
         s=s+T_t(Omega_plus(i));
     end
     d_omega=zeros(n);
     for i=1:op
         d_omega(Omega_plus(i))=T_t(Omega_plus(i))/s;
     end
     
     for i=1:op
         T_t(Omega_plus(i))=T_t(Omega_plus(i))+(1-d_omega(Omega_plus(i)))*dd;
     end
//*****************************************************************************
end
    M=M_t./sum(M_t);
    T=T_t./sum(T_t);
end

plot(time,T_v);
plot(time,M_v,'--');
