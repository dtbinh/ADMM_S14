%% parameters
J = 6;
Nplus = [1:J];
sigma = [0.8 0.15 0.5 1 0.2 0.7];
N = 2;
M = 4;
trial = 100;
SNR = 10;
sigma_n = 10^(-SNR/10);
nj = zeros(M,1,J);
xj = zeros(M,1,J);
Hj = zeros(M,N,J);
error_zf = 0;
error_mmse = 0;
iter = 100;
K = [1:3:iter];
aMSE_zf = [];
aMSE_mmse = [];
aMSE_dice = [];


%% simulation
for k = K
    error_dice = 0;
 
    for i = 1:trial
        s = rand(N,1);
        H = [];
        n = [];
        for j=1:J
          nj(:,:,j) = sqrt(sigma(j)*sigma_n)*rand(M,1);
          Hj(:,:,j) = sqrt(sigma(j))*rand(M,N);
          xj(:,:,j) = Hj(:,:,j)*s+nj(:,:,j);
          sj_hat = inv(Hj(:,:,j)'*Hj(:,:,j))*Hj(:,:,j)'*xj(:,:,j);
          H = [H Hj(:,:,j)'];
          n = [n nj(:,:,j)'];
        end
        H = H';
        n = n';
        x = H*s+n;
        % uncomment the following line to run admm in parallel (make sure
        % to configure the local parallel profile to have J=number of cores
        %[sj]=paradmm(Hj,xj,J,N,J,k,sigma_n);
        
        % comment the following line if the previous one is uncommented
        [sj]=admm(Hj,xj,J,N,J,k,sigma_n);
        
        for j=1:J
            error_dice = error_dice + norm(s-sj(:,:,j))^2;
        end
               
        s_zf = inv(H'*H)*H'*x;
        s_mmse = inv(H'*H+sigma_n*eye(N,N))*H'*x;

        error_zf = error_zf + norm(s-s_zf)^2;
        error_mmse = error_mmse + norm(s-s_mmse)^2;

    end
aMSE_dice = [aMSE_dice error_dice/J/trial]
end
ksize = size(K);
aMSE_zf = ones(size(K))*error_zf/trial/ksize(2)
aMSE_mmse = ones(size(K))*error_mmse/trial/ksize(2)

semilogy(K,aMSE_zf,K,aMSE_mmse,K,aMSE_dice)
axis([0,iter,10^-2,10^0])
legend('Central-ZF','Central-MMSE','ADMM')
xlabel('Iterations (k)')
ylabel('MSE')
grid on    


%%
figure()
J=6;
M=4;
fz = ones(size(K))*J*(1+M);
fadmm = 2 * J * K;
plot(K,fz,'-',K,fadmm,'-');
hold on
M=150;
fz = ones(size(K))*J*(1+M);
fadmm = 2 * J * K;
plot(K,fz,'--',K,fadmm,'--');
legend('Central-ZF M=4','ADMM M=4','Central-ZF M=150','ADMM M=150')
xlabel('Iterations (k)')
ylabel('MSE')
grid on    
