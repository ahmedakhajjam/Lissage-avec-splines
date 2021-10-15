clear all;
clc;
X=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38];
y=[1 0.5 2.3 4.6 3.4 3.2 4.3 1 1.5 7 5 4.8 6.4 2.7 5 4.7 49 35 34 24 23 31 26 23 20 15 14 20 23 39 18 34 19 14 13 28 23 25]; 
n=length(X);
N=n-1;
%lambda=5;
sigma=1;
for j=1:N
    h(j)=X(j+1)-X(j);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% spline cubique de lissage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:N
W(i)=2*(h(i-1)+h(i));
end
for i=2:N-1
R(i)=h(i);
end
for j=1:N+1
    for i=1:N+1
        T(i,j)=0;
    end
end
T=diag(W)+diag(R,-1)+diag(R,1);
T(1,1)=1;
T(N+1,N+1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% la Marice Q
for i=1:N
z(i)=1/h(i);
end
for j=1:N+1
    for i=1:N+1
        Q(i,j)=0;
    end
end
for i=1:N+1
    for j=2:N
        if i==j
            Q(i,j)=-(z(i)+z(i-1));
        else if i==j-1
                Q(i,j)=z(i);
            else if i==j+1
                    Q(i,j)=z(i-1);
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matrice Sigma
Sigma=zeros(n,n);
vect=zeros(1,n);
for i=2:n-1
    vect(1,i)=sigma;
end
Sigma=diag(vect);
%%%%%%%%%%%%%%%%%%%%%%%%%%% la matrice Ac=v
o=1;
for lambda=0:0.0001:0.05
AA=Q'*Sigma.^2*Q+lambda*T;
v=lambda*Q'*y';
%%%%%%%%%%%%%%%%%%%%%%%%%% la valeur de c
c=LDLFact(AA,v);
c(1)=0;
c(N+1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% la valeur de a
a=y'-lambda.^(-1)*Sigma.^2*Q*c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% la valeur de d
for i=1:N
    d(i)=(c(i+1)-c(i))/3*h(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% la valeur de b
for i=1:N
    b(i)=((a(i+1)-a(i))/h(i))-c(i)*h(i)-d(i)*h(i).^2;
end

 for i=1:N
      for xx=X(i):0.01:X(i+1)
       S=d(i)*(xx-X(i)).^3+c(i)*(xx-X(i)).^2+b(i)*(xx-X(i))+a(i);
       F(i)=S;
      end
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%% fin  spline cubique de lissage
M=(Q'*Sigma^2*Q+lambda*T);
invM=inv(M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=Q*invM*Q';
%%%%%%%%%%%%%%%%%%%% score de validation croisé
ccv=0;
for i=2:n-1
    ccv=ccv+(1/n)*((((y(i)-F(i))/H(i,i)))^2); 
end
CV(o)=ccv;
Lambda(o)=lambda;
    o=o+1;
end
plot(Lambda,CV,'-b','linewidth',1);
 grid on
box on
xlabel('lambda')
ylabel('CV(lambda) ')
  hold on;
 legend('Le score de validation croisé')

