clear;
if 1==1 then
n=4;
for i=1:n
    T(i)=rand(100);
    M(i)=rand(100);
end

end
M=[0.3076270773539253	0.2862893346342559	0.34336478245952484	0.06271880555229399];
T=[0.4423636783115587	0.17600591303341878	0.3087916240091037	0.07283878464591878];
//M=[28.21907889849865;1.6328567228213704;3.317358400570099;21.9556898377048;4.167228835048353;4.411700166545918;30.353352995717998;5.942734143092812];
//T=[16.591281920036003;4.87892862433848;22.049985297842547;17.09988125069187;2.109881336443525;4.969067996286483;17.03596544495946;15.265008129401636];
//n=8
//n=length(M);
st = sum(T);
sm = sum(M);
T = T'./st;
M = M'./sm;
eps=1/(n*n*n);

        
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

m = 200;//number of steps
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
        s=0.6//Theta;
        h=0.6//Theta;
        elseif check(k+1)=='minus' then
            z = 1+Theta-W;
            T_t = (T.*(1+Theta)-tau)./z;
            M_t = (M.*(1+Theta)-tau)./z;
            s=0.6//1-Theta;
            h=0.6//1-Theta;
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

//subplot(211)
//plot(time,M_v(:,4),'black');
//plot(time,T_v(:,1),'black--');
//subplot(212)
//plot(time,T_v(:,4),'black--');
//legends(['T', 'M'],[1,1--])
//plot(time,M_v(:,n-1),'black');
//plot(time,T_v(:,n-1),'black--');
//subplot(313)
plot(time,T_v,'black');
plot(time,M_v,'black--');

