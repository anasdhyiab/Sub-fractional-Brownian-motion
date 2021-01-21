%The simulation of the sample paths of sub-fractional Brownian motion 
% S^H_t= K^{-1}_H\int_R[t-s]^{H-1/2}_+ +[-t-s]^{H-1/2}_+ - 2[-s]^{H-1/2}_+ dW_s
%K_H=\sqrt{2\int_R|[1-s]^{H-1/2}_+ - [-s]^{H-1/2}_+|^2ds+\frac{1}{2H}}
% Note that: Imaginary parts of complex X and/or Y arguments ignored
%randn('state',100)
Xzero=1; M = 2000; T = 1;  dt=T/M;  H= 0.7;
tic;
R = 4; Dt = R*dt; m = M/R;            % Steps of size Dt = R*dt
S = zeros(1,m);                       % Preallocate for efficiency
St = Xzero;
for j = 1:m
    for k = 1:M
b = @(x)      (x).^(H-1/2)-(x-1).^(H-1/2)+((x-2*k).^(H-1/2)-(x-2*k+1).^(H-1/2));
dW = sqrt(dt)*randn(1,M);
Winc = sum(dW(R*(j-1)+1:R*j));         % Discretized Brownian path
K2=1/m*sqrt((2*(abs((j/m).^(H-1/2)-(j/m-1).^(H-1/2)).^2)+1/(2*H)));
St = St+K2.^(-1)*b(j/m).*Winc;
S(j) = St; 
    end 
end
toc;
subplot(2,1,2);
plot([0:Dt:T],[Xzero,S],'b')
ll=legend("H="+H);
set(ll,'Fontsize',12);
xlabel('Time vector','FontSize',12)
ylabel('sfBm','FontSize',12,'Rotation',0,'HorizontalAlignment','right');