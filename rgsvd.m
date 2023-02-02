function [U,V,X,C,S]=rgsvd(A,B,R)
[m,n]=size(A);
A_1=A*randn(n,R);
B_1=B*randn(n,R);
[Q_1,~]=qr(A_1,0);
[Q_2,~]=qr(B_1,0);
[U_1,V_1,X,C,S]=gsvd(Q_1'*A,Q_2'*B);
U=Q_1*U_1;
V=Q_2*V_1;
end