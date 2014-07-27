%% Sequential implementation of ADMM

function [sj]=admm(Hj,xj,Nplus,N,J,K,sigma_n)

mu = 1; % this is rho in the final report

zj = zeros(N,1,J);
sj = zeros(N,1,J);
lamda = zeros(N,1,J,Nplus);
for  k = 1:K %number of iterations
    for j = 1:J
        mult1 = inv(Hj(:,:,j)'*Hj(:,:,j)+Nplus/mu*eye(N,N));
        mult2 = Hj(:,:,j)'*xj(:,:,j);
        for i = 1:Nplus
            mult2 = mult2 + zj(:,:,i)/mu+lamda(:,:,j,i);
        end
        sj(:,:,j) = mult1*mult2;
    end
    for j=1:J
        temp = 0;
        for i = 1:Nplus
            temp = temp+1/mu*sj(:,:,i)-lamda(:,:,i,j);
        end
        zj(:,:,j) =temp*mu/Nplus;
    end
    for j=1:J
        for i=1:Nplus
            lamda(:,:,i,j)=lamda(:,:,i,j)-1/mu*(sj(:,:,i)-zj(:,:,j));
        end
    end
end
