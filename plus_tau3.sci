clear;
n=8;
if 1==1 then
for i=1:n-2
    T(i)=rand(100)*100;
    M(i)=rand(100);
end
for i=n-1:n
    T(i)=rand(100);
    M(i)=rand(100)*150;
end
end
//n=4;
//T = [15,14,3,32]//,11,26,2,12,54,18];
//M = [21,12,13,24]//,5,32,11,6,32,23];
    st_start = sum(T);
    sm_start = sum(M);    
    T = T'./st_start;
    M = M'./sm_start;

m = 60;
for k = 1:m
    st_start = sum(T);
    sm_start = sum(M);
    time(k)=k;
    T_v(k,:)=T;
    M_v(k,:)=M;
    summaM(k) = sum(M);
    summaT(k) = sum(T);
    st = sum(T);
    sm = sum(M);
    T = T./st;
    M = M./sm;
    Theta = sum(T.*M);
    tau = min(T,M);
    W = sum(tau);
    z = 1+Theta+W;
    
    T_conf = (T.*(1+Theta)+tau)./z;
    M_conf = (M.*(1+Theta)+tau)./z;
    if 1==1 then    
    c1=1.2;
    c2=.84;
    for i=1:n
        if T_conf(i) > max(T_conf)*0.5 then
            T(i)=T_conf(i)*c1*st_start*c2;
        else T(i) = T_conf(i)*st_start*c2
        end
        if M_conf(i) > max(M_conf)*0.8 then
             M(i)=M_conf(i)*c2*sm_start*c1;
             else M(i) = M_conf(i)*sm_start*c1
        end
    end
else
       T = T_conf.*st_star;
       M = M_conf.*sm_start;

end

end

subplot(211)
plot(time,T_v);
plot(time,M_v,'--');
subplot(212)
plot(time,summaM,'b')
plot(time,summaT,'g')
legend('sumM','sumT')
