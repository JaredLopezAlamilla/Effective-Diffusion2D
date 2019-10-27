function [J2xxk,J2yyk,J2xyk,K0]=DiffusionTensor2D(Pkq,Akq,Jxk,Jyk,kBTx,kBTy,gx,gy,vx,vy,n,M)

    LL=(-1i*2*pi)^2*Akq;

    lx=-vx*spdiags(ones(M^2,1),0,M^2,M^2);
    ly=-vy*spdiags(ones(M^2,1),0,M^2,M^2);

    J2x=-1i*2*pi*(Jxk);
    J2y=-1i*2*pi*(Jyk);

    Chix=-1i*2*pi*(kBTx/gx)*spdiags(sort(repmat(n,[1 M])).',0,M^2,M^2)+J2x;    
    Chiy=-1i*2*pi*(kBTy/gy)*spdiags(repmat(n,[1 M]).',0,M^2,M^2)+J2y;

    PreU1xk=-(Chix+lx)*Pkq(:); PreU1xk(abs(PreU1xk)<=1e-4)=0;
    PreU1yk=-(Chiy+ly)*Pkq(:); PreU1yk(abs(PreU1yk)<=1e-4)=0;

    K0=1+(M^2-1)/2;

    LL2=[LL(1:K0-1,:);LL(K0+1:end,:)];        

    U1xk=LL2\PreU1xk([1:K0-1 K0+1:end]); %U1xk(abs(U1xk)<=1e-4)=0;        
    U12xk=U1xk+(1-U1xk(K0))*Pkq(:);   %U12xk(abs(U12xk)<=1e-4)=0;     

    U1yk=LL2\PreU1yk([1:K0-1 K0+1:end]); %U1yk(abs(U1yk)<=1e-4)=0;
    U12yk=U1yk+(1-U1yk(K0))*Pkq(:);     %U12yk(abs(U12yk)<=1e-4)=0;  

    J2xxk=sparse((Chix+lx)*U12xk); 
    J2yyk=sparse((Chiy+ly)*U12yk); 
    J2xyk=sparse((Chix+lx)*U12yk);  
end
