%PFE MASTER-AHMED AKHAJJAM
%2020/2021
function liss=LISSAGESPL(X,y,lambda,sigma)
%y = load('data2.csv');%disp(yy);
%for i=1:44
%X(i)=i;
%end
X=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]% 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37];
y=[10 24 23 39 24 32 49 60 45 59 48 64 27 50 47 49 35 34 24 23]% 31 26 23 20 15 14 20 23 39 18 34 19 14 13 28 23 25];
n=length(X);
N=n-1;
lambda=0.1; %%%%%% paramètre de lissage
sigma=1;%%%%%%%%%%%%% valeur de sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% la Marice Q
for j=1:N
    h(j)=X(j+1)-X(j);
end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% la Marice Q
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matrice Sigma
Sigma=zeros(n,n);
vect=zeros(1,n);
for i=2:n-1
    vect(1,i)=sigma;
end
Sigma=diag(vect);
%%%%%%%%%%%%%%%%%%%%%%%%%%% la matrice Ac=v
A=Q'*Sigma.^2*Q+lambda*T;
v=lambda*Q'*y';
%%%%%%%%%%%%%%%%%%%%%%%%%% la valeur de c
c=LDLFact(A,v);
c(1)=0;
c(N+1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% la valeur de a
a=y'-lambda.^(-1)*Sigma.^2*Q*c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% la valeur de d
for i=1:N
    d(i)=(c(i+1)-c(i))/3*h(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% la valeur de b
for i=1:N
    b(i)=((a(i+1)-a(i))/h(i))-c(i)*h(i)-d(i)*h(i).^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trace le figure
 for i=1:N
       xx=X(i):0.1:X(i+1);
       S=d(i)*(xx-X(i)).^3+c(i)*(xx-X(i)).^2+b(i)*(xx-X(i))+a(i);
       plot(xx,S,'-r','linewidth',1);
       hold on;
 end
 grid on
box on
xlabel('x ')
ylabel('y ')
  hold on;
 SPLINCUBIQUE(X,y);
 legend('spline cubique de lissage')
end
