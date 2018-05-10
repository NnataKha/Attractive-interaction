clear;
if 1==1 then
n=15;
for i=1:n
    T(i)=rand(100)//*110;
    M(i)=rand(100);
end
end
n=length(M);
st = sum(T);
sm = sum(M);
T = T'./st;
M = M'./sm;
eps=1/(n);
m = 100;//number of steps
check(1)='plus';
for k = 1:m
    time(k)=k;
    T_v(k,:)=T;
    M_v(k,:)=M;
    Theta=0;
    for i=1:n
        Theta = Theta + sqrt(T(i)*M(i));
    end
    Th(k)=Theta;
    
    check(k+1)=check(k);

    //check the value of Theta
    if Theta<=eps then
        check(k+1)='plus';
    elseif Theta>=1-eps then
        check(k+1)='minus';
    end

    tau = min(T,M);
    W = sum(tau);
    //choosing “+,+” or  “-,-”
    if check(k+1)=='plus' then
        z = 1+Theta+W;
        T_t = (T.*(1+Theta)+tau)./z;
        M_t = (M.*(1+Theta)+tau)./z;
        elseif check(k+1)=='minus' then
            z = 1+Theta-W;
            T_t = (T.*(1+Theta)-tau)./z;
            M_t = (M.*(1+Theta)-tau)./z;
        end

    sm=sum(M_t);
    st=sum(T_t);
    M=M_t./sm;
    T=T_t./st;
end
plot(T_v,M_v,'.-');

