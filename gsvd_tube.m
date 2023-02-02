function [U_1,V_1,X,C,S]=gsvd_tube(A,B)
[~,~,n]=size(A);
A_1=fft(A,[],3);
B_1=fft(B,[],3);
halfn3 = ceil((n+1)/2);
tic
for i=1:ceil((n+1)/2)
[U(:,:,i),V(:,:,i),X(:,:,i),C(:,:,i),S(:,:,i)] = gsvd(A_1(:,:,i),B_1(:,:,i));
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
X=ifft(X,[],3);
C=ifft(C,[],3);
S=ifft(S,[],3);

end