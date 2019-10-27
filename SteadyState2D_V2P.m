%-------------------------------------------------------------------------
%    This code computes the steady-state distribution in position space
%               and calculates its k-space form 
%-------------------------------------------------------------------------

% ---- k-space parameters ----
M=19;MM=1001; [n,m]=deal(-floor(M/2):floor(M/2));
% ---- potential parameters
kBTx=1; kBTy=1; gx=1; gy=1; L=1;
Vx=3;Vy=2;Vxy=5; 
% ---- Creates potential ----  
[xf,yf]=deal(L*(0:M-1)/M);
[Xf,Yf]=meshgrid(xf,yf);
V0=cos_potential2D(Xf,Yf,Vx,Vy,Vxy);

% ---- Fourier transform ---
Vkq=fftshift(fft2(V0)/(M^2));

% ---- Creates coefficient matrix ----
[Akq,Jxk,Jyk]=find_matrix(Vkq,n,M,kBTx,kBTy,gx,gy,fx,fy,L);

% ---- Solve for coefficients ---- 
[Pkq, Jxkq, Jykq, ~, ~, ~,Norm]=null_solver(Akq,Jxk,Jyk,M,L,n);

mid = round(length(Jxkq)/2);
vx = real(Jxkq(mid,mid))*L;
vy = real(Jykq(mid,mid))*L;
