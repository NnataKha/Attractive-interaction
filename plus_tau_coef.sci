clear;
if 1==1 then
n=30;
for i=1:n-3
    T(i)=rand(100);
    M(n-i+1)=rand(100);
end
for i=n-2:n
    T(i)=rand(100)*100;
    M(n-i+1)=rand(100)*100;
end

end
//n=length(M);
st = sum(T);
sm = sum(M);
T = T./st;
M = M./sm;
eps=1/(n);
m = 1000;//number of steps
c(1)=1;
for k = 1:m
    time(k)=k;
    T_v(k,:)=T;
    M_v(k,:)=M;
    Theta(k) = sum(sqrt(T.*M));
    c(k+1)=c(k)*(1-eps-Theta(k))*(Theta(k)-eps)/abs((1-eps-Theta(k))*(Theta(k)-eps));

    tau = min(T,M);
    W = sum(tau);
    z = 1+Theta(k)+W;
    for i=1:n
       if T(i) > M(i) then
           M_t(i)=M(i)*(1+Theta(k))+c(k+1)*(1-Theta(k))*tau(i);
           T_t(i)=T(i)*(1+Theta(k))+tau(i);
       else            
           T_t(i)=T(i)*(1+Theta(k))+c(k+1)*(1-Theta(k))*tau(i);
           M_t(i)=M(i)*(1+Theta(k))+tau(i);
       end
    end
    sm=sum(M_t);
    st=sum(T_t);
    M=M_t'./sm
    T=T_t'./st;
    Th(k)=Theta(k)-eps;
end
subplot(211)
plot(time,T_v);
plot(time,M_v,'--');
subplot(212)
plot(time,Th)
