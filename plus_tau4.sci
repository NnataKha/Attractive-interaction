//we subtract maximal grouth on Omega_+ from each M(i)
clear;
//n=7;
//T = [150,140,100, 125, 105,6,4];//630
//M = [5.0,4.0,3.5, 4.5, 3, 280,330];//630
//n=6;
//T = [100,120,110,150,3,4];
//M = [5,4,5,2,310,210];
n=8;
if 1==1 then
for i=1:n-2
    T(i)=rand(100)//*110;
    M(i)=rand(100);
end
for i=n-1:n
    T(i)=rand(100);
    M(i)=rand(100)//*70;
end
end
//M=[0.0038123686788460774;0.0061186711913324265;6.414758228991696E-4;0.0037405876180425087;0.5800012548666301;0.4056856418222497];
//T=[0.1754599877636935;0.2157347341479504;0.30724683197454583;0.29831983748949625;2.089925038084202E-4;0.0030296161205055653];
//n=length(T);
st = sum(T);
sm = sum(M);
T = T'./st;
M = M'./sm;
    minM=0;
    minT=0;
    for i=1:n
       tau(i) = min(T(i),M(i));
       if tau(i) == M(i) then
           minM=M(i);
       else minT=T(i);
       end
    end
m = 50;
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

    d=abs(M(Omega_plus(1))-M_t(Omega_plus(1)));
    for i=1:op
        m=abs(M(Omega_plus(i))-M_t(Omega_plus(i)));
        if d<m then
            d=m;
        end
    end
    dd=0;
    for i=1:op
        M_t(Omega_plus(i))=M_t(Omega_plus(i))-d;
        if M_t(Omega_plus(i))<0 then
            M_t(Omega_plus(i))=minM;
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

    d=abs(T(Omega_minus(1))-T_t(Omega_minus(1)));
    for i=1:om
        m=abs(T(Omega_minus(i))-T_t(Omega_minus(i)));
        if d<m then
            d=m;
        end
    end
    dd=0;
    for i=1:om
        T_t(Omega_minus(i))=T_t(Omega_minus(i))-d;
        if T_t(Omega_minus(i))<0 then
            T_t(Omega_minus(i))=minT;
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
subplot(211)
plot(time,M_v(:,1),'g');
plot(time,T_v(:,1),'g--');
subplot(212)
plot(time,T_v(:,3),'b--');
plot(time,M_v(:,3),'b');


