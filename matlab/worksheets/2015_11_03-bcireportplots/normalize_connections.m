function [ B ] = normalize_connections( A)
%normalizes conections (predicted or actual) by subtracting the mean and
%dividing by the standard deviation

[m,n]=size(A);

C=reshape(A,1,m*n);

d=std(C);

mn=mean(C);

B=(A-mn)/d;

end

