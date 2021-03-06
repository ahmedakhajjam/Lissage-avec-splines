%PFE MASTER-AHMED AKHAJJAM
%2020/2021
clear all;
clc;
for i=1:58
X(i)=i;
end
yy = load('data-age_0_19.csv');
%yy = load('data-age_20_39.csv');
%yy = load('data-age_40_59.csv');
%yy = load('data-age_60+.csv');
y=yy'; 
n=length(X);
N=n-1;
lambda=0.002;
sigma=1;
for j=1:N
    h(j)=X(j+1)-X(j);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% spline cubique de lissage
for k=1:N+1
    alpha(k)=y(k);
end
for i=1:N
    B(i)=0;
end
for i=2:N
    B(i)=(3/h(i))*(alpha(i+1)-alpha(i))-(3/h(i-1))*(alpha(i)-alpha(i-1));
end
B(1)=0;
B(N+1)=0;
for j=1:N+1
    for i=1:N+1
        A(i,j)=0;
    end
end
for j=1:N+1
    for i=2:N
        if i==j
            A(i,j)=2*(h(i-1)+h(i));
        else if i==j+1
                A(i,j)=h(i-1);
            else if i==j-1
                    A(i,j)=h(i);
                end
            end
        end
    end
end
A(1,1)=1;
A(N+1,N+1)=1;
gamma=A\B';
for i=1:N
    beta(i)=(1/h(i))*(alpha(i+1)-alpha(i))-(h(i)/3)*(2*gamma(i)+gamma(i+1));
    delta(i)=(gamma(i+1)-gamma(i))/(3*h(i));
end
 %%%%%%% fin  spline cubique  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% spline cubique de lissage
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
%%%%%%%%%%%%%%%%%%%%%% la Matrice Q
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
AA=Q'*Sigma.^2*Q+lambda*T;
v=lambda*Q'*y';
%%%%%%%%%%%%%%%%%%%%%%%%%% le valeur de c
c=LDLFact(AA,v);
c(1)=0;
c(N+1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% le valeur de a
a=y'-lambda.^(-1)*Sigma.^2*Q*c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% le valeur de d
for i=1:N
    d(i)=(c(i+1)-c(i))/3*h(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% le valeur de b
for i=1:N
    b(i)=((a(i+1)-a(i))/h(i))-c(i)*h(i)-d(i)*h(i).^2;
end
 %%%%%%% la fin de spline cubique de lissage
%%%%%%%%%%%%%%%%%%%%%%%%%% trace les figures
 for i=1:N
    xx=X(i):0.01:X(i+1);
    S=d(i)*(xx-X(i)).^3+c(i)*(xx-X(i)).^2+b(i)*(xx-X(i))+a(i);
    s=alpha(i)+beta(i)*(xx-X(i))+gamma(i)*(xx-X(i)).^2+delta(i)*(xx-X(i)).^3; 
    % plot(xx,S,'-r','linewidth',1);%% lissage
    plot(xx,s,'-r','linewidth',1);
    plot(X,y,'r+')%,'linewidth',1.1);
     hold on;
 end
%grid on
box on
xlabel('Date ')
ylabel('Nouveaux cas quotidiens ')
%legend('Nombre infect?s','spline cubique de lissage pour lambda=0.001','spline cubique de lissage pour lambda=0.01')
