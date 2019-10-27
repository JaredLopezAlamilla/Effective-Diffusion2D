%-------------------------------------------------------------------------
%          Calculates effective diffusion  and entropy production
%            on the steady-state for overdamped Brownian motion 
%                 over  a tilted periodic potential in 2D
%-------------------------------------------------------------------------
clear;
x=-2:0.01:2; [xgrid,ygrid]=meshgrid(x);

Fys=0:1:8;  Fxs=0:1:8; 

DDeffs=zeros(2,2,length(Fxs),length(Fys));
Dlambda=zeros(length(Fxs),length(Fys));
vs=zeros(1,2,length(Fxs),length(Fys)); 

S_ss=zeros(length(Fxs),length(Fys));
S_cg=zeros(length(Fxs),length(Fys));
S_tury=zeros(length(Fxs),length(Fys));
S_turx=zeros(length(Fxs),length(Fys)); 

jj=1;

 for j=Fxs
     fx=j;
     kkk=1;     
    for kk=Fys
        fy=kk;
        SteadyState2D_V2P;
        [J2xxk,J2yyk,J2xyk,K0]=DiffusionTensor2D(Pkq,Akq,Jxk,Jyk,kBTx,kBTy,gx,gy,vx,vy,n,M);    
        
        D0=[kBTx/gx 0;0 kBTy/gy]; 
        diffusion=real(D0-[J2xxk(K0) J2xyk(K0);J2xyk(K0) J2yyk(K0)]); 
            
        DDeffs(:,:,jj,kkk)=diffusion(:,:);
        Dlambda(jj,kkk)=det(diffusion(:,:));
        vs(:,:,jj,kkk)=[vx vy];         
        
        S_ss(jj,kkk)=[vx vy]*[fx fy].';
        S_cg(jj,kkk)=(diffusion(:,:)\[vx vy].').'*[vx vy].';
        S_tury(jj,kkk)=vx^2/diffusion(1,1);
        S_turx(jj,kkk)=vy^2/diffusion(2,2);
        
        kkk=kkk+1;   
    end    
jj=jj+1;  
end

[GridFx,GridFy]=meshgrid(Fxs,Fys);

S_cg(abs(S_cg)<=1e-10)=0;
Rel_cg=(S_ss-S_cg)./S_ss;
Rel_cg(isnan(Rel_cg))=0;

Rel_tury=(S_ss-S_tury)./S_ss;
Rel_tury(isnan(Rel_tury))=0;

Rel_turx=(S_ss-S_turx)./S_ss;
Rel_turx(isnan(Rel_turx))=0;
