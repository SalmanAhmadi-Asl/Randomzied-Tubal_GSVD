clc;clear all
n=400;
M=400;
P=400;
R=50;
% A=t_prod(randn(M,R,n),randn(R,P,n));
% B=t_prod(randn(M,R,n),randn(R,P,n));
for i=1:M
for j=1:P
    for k=1:n
A(i,j,k)=1/((i^2+j^2+k^2)^(1/2));
    end
end
end
% 
for i=1:n
for j=1:M
    for k=1:P
        B(i,j,k)=1/((i^2+j^2+k^2)^(1/2));
    end
end
end

%%
tic
A_1=fft(A,[],3);
B_1=fft(B,[],3);
halfn3 = ceil((n+1)/2);
for i=1:ceil((n+1)/2)
[U(:,:,i),V(:,:,i),X(:,:,i),C(:,:,i),S(:,:,i)] = rgsvd(A_1(:,:,i),B_1(:,:,i),R);
end
for i = halfn3+1 : n
        U(:,:,i) = conj(U(:,:,n+2-i));
        V(:,:,i) = conj(V(:,:,n+2-i));
        X(:,:,i) = conj(X(:,:,n+2-i));
        C(:,:,i) = conj(C(:,:,n+2-i));
        S(:,:,i) = conj(S(:,:,n+2-i));
end

U_1=ifft(U,[],3);
V_1=ifft(V,[],3);
X_1=ifft(X,[],3);
C_1=ifft(C,[],3);
S_1=ifft(S,[],3);
toc
E_1=A-t_prod(t_prod(U_1,C_1),t_trans(X_1));
% norm(E(:))

E_2=B-t_prod(t_prod(V_1,S_1),t_trans(X_1));
% norm(E(:))
(norm(E_1(:))+norm(E_2(:)))/(norm(A(:)+norm(B(:))))

%%
tic
A_1=t_prod(A,randn(P,R,n));
B_1=t_prod(B,randn(P,R,n));
Q_1=QR_tubal(A_1);
Q_2=QR_tubal(B_1);
[U_1,V_1,X,C,S]=gsvd_tube(t_prod(t_trans(Q_1),A),t_prod(t_trans(Q_2),B));
U=t_prod(Q_1,U_1);
V=t_prod(Q_2,V_1);
toc

E_1=A-t_prod(t_prod(U,C),t_trans(X));
% norm(E(:))

E_2=B-t_prod(t_prod(V,S),t_trans(X));
% norm(E(:))

(norm(E_1(:))+norm(E_2(:)))/(norm(A(:)+norm(B(:))))