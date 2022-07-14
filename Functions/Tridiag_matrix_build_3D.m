function out=Tridiag_matrix_build_3D(a,b,c,N)


% Build a 3D tridiagonal matrix of dimension m x m x m considering
% a : vector of length m 
% b : vector of length m 
% c : vector of length m 


out=zeros(N,N,N);
out(2:1+N:N*N-N)=a(2:end);
out(1:1+N:N*N)=b;
out(N+1:1+N:N*N)=c(1:end-1);

end