%PFE MASTER-AHMED AKHAJJAM
%2020/2021
function SP= SPLINCUBIQUE(X,y)
X=[1 2 3 4 5 6 7 8 9 10 11 12];
y=[2 5 6 8 7 4 5 6 4 2 5 8];
n=length(X);
N=n-1;
for j=1:N
    h(j)=X(j+1)-X(j);
end
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
%%%%%%%%%%%%%%%
for j=1:N
   xx=[X(j):0.01:X(j+1)];
   s=alpha(j)+beta(j)*(xx-X(j))+gamma(j)*(xx-X(j)).^2+delta(j)*(xx-X(j)).^3; 
   plot(xx,s,'-b','linewidth',1.5);
   hold on; 
end
legend('spline cubique');
end

