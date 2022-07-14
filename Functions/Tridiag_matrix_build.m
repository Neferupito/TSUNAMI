function [out,e1,e2]=Tridiag_matrix_build(a,b,c,N)

% Build a tridiagonal matrix of dimension m x m considering
% a : vector of length m 
% b : vector of length m 
% c : vector of length m 
% e1, e2 : extrem values first and last line

out=zeros(N,N);
out(2:1+N:N*N-N)=a(2:end);
out(1:1+N:N*N)=b;
out(N+1:1+N:N*N)=c(1:end-1);
e1=a(1);
e2=c(end);

end
