clear;
n=2;
T = [0.8, 0.2];
M = [0.25, 0.75];
h = [0.7, 1];
s = [1, 1.3];
m = 200;
for k = 1:m
    time(k)=k;
    T_v(k,:)=T;
    M_v(k,:)=M;
    Theta = sum(T.*M);
    tau = min(T,M);
    W = sum(tau);
    z = 1+Theta+W;
    
    T_conf = (T.*(1+Theta)+tau)./z;
    M_conf = (M.*(1+Theta)+tau)./z;
    
if 1==1 then   
    T = T.*h;
    M = M.*s; 
    
    st = sum(T);
    sm = sum(M);
    
    T = T./st;
    M = M./sm;
else
    T = T_conf;
    M = M_conf;
end

end

//subplot(212)
plot(time,T_v);
//subplot(212)
plot(time,M_v,'.-');
