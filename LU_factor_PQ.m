function [dif_Voltage]=LU_factor_PQ(J,dif_PQ) % using Crout's Algorithm for computing LU 
%A=[1 3 4 8;2 1 2 3;4 3 5 8;9 2 7 4];
%b=[1;1;1;1];
A=J;
b=dif_PQ;
n=length(A); %get the dimention of A matrix    
Q=zeros(n);
Q_1=zeros(n); 

for j=1:n
% elements of jth column in Q matrix
    for k=j:n 
        for i=1:j-1
            Q_1(k,j)=Q_1(k,j)+Q(k,i)*Q(i,j);
        end
        Q(k,j)=A(k,j)-Q_1(k,j);
     end
% elements of jth row in Q matrix
if Q(j,j)~=0     
for k=j+1:n    
    for i=1:j-1
        Q_1(j,k)=Q_1(j,k)+Q(j,i)*Q(i,k);
    end
    Q(j,k)=(A(j,k)-Q_1(j,k))/Q(j,j);
end
%else
    %a=max(A(:,j));
   % [row, col]=find(a=A(:,j));
end
  j=j+1;
end
%% get vector y by forward substitution
y_1=zeros(n,1);
y=zeros(n,1);
for k=1:n
    for j=1:k-1
        y_1(k)= y_1(k)+Q(k,j)*y(j);
    end
    y(k)=(b(k)-y_1(k))/Q(k,k);    
end
%% get the solution vector x
x_1=zeros(n,1);
x=zeros(n,1);
for k=n:-1:1
    for j=k+1:n
    x_1(k)=x_1(k)+Q(k,j)*x(j);
    end
    x(k)=y(k)-x_1(k);
end
  dif_Voltage=x;
end
 
   



    