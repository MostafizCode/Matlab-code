# Matlab-code
% Summation, Multiplication, Division, Reminder Method.

❤❤%input section❤❤
x=input("Enter the value of x[...]");
y=input("Enter the value of y[...]");

%Formula section
Summation=x+y;
Multiplication=x.*y;
Division=x/y;
Reminder=rem(x,y);






❤❤% Equation Solve❤❤

%input Section
syms x y
eq1=input("Enter the first equation");
eq2=input("Enter the second equation");

%Formula section
sol=solve([eq1,eq2],[x,y]);
disp(sol.x);
disp(sol.y);





❤❤% Equation Plot❤❤

%Input section
syms y
x=-12:0.2:12;
y=input("Enter the equation");

%Formula section
plot(x,y);

❤❤% Curve Fitting Method❤❤

%Input section
x=input("enter_matrix",'s');
y=input("enter_matrix",'s');

%Formula Section
S1=sum(x,"all");
S2=sum(y,"all");
M1=(x.*y);
S3=sum(M1,"all");
M2=(x.*x);
S4=sum(M2,"all");

%Equation section
syms a b
eq1=S3*a+b*S1==S2;
eq2=S1*a+b*S4==S3;

%Solution section
sol=solve([eq1,eq2],[a,b]);
disp(sol.a);
disp(sol.b);

%Plot section
syms a. b
a=0:30;
eq1=S3*a+b*S1==S2;
plot(a,eq1);





❤❤%Newton Forward Method❤❤

%Input section
x=input("Enter first matrix");
y=input("Enter second matrix");

%Formula section
n=length(x);
p=3;
d(:,1)=y';
for j=2:n
    for i=j:n
        d(i,j)=(d(i-1,j-1)-d(i,j-1))/(x(i-j+1)-x(i));
    end
end
a=diag(d)';
df(1)=1;
c(1)=a(1);
for j=2:n
    df(j)=(p-x(j-1)).*df(j-1)
    c(j)=a(j).*df(j)
end
fp=sum(c);





❤❤%Newton backward method❤❤

%Input section
x=input("Enter first matrix value");
fx=input("Enter second matrix value");
dt=zeros(6,7);

%Loop section
for i=1:6
    dt(i,1)=x(i);
    dt(i,2)=fx(i);
end
n=5;
for j=3:7
    for i=1:n
        dt(i,j)=dt(i+1,j-1)-dt(i,j-1)
    end
    n=n-1
end
h=x(2)-x(1)
xp=27;
for i=1:6
    q=(xp-x(i))/h;
    if(q>0&&q<1)
        p=q;
    
    end
end

l=xp-(p*h)
for i=1:6
    if(l==x(i))
        r=i
    end
end

%Formula section
f0=fx(r);
f01=dt((r-1),3);
f02=dt((r-2),(3+1));
f03=dt((r-3),(3+2));
f04=dt((r-3),(3+2));
fp=(f0)+((p*f01)+(p*(p+1)*f02)/(2))+((p*(p+1)*(p+2)*f03)/(6));






❤❤%Gauss elimination method❤❤

%Input section
a=input("Enter matrix[..]");

% [m,n)=size(a);
[m,n]=size(a);

%Start loop
for j=1:m-1
    for z=2:m
    if a(j,j)==0
t=a(j,:);a(j,:)=a(z,:);
a(z,:)=t;
    end
end
    for i=j+1:m
a(i,:)=a(i,:)-a(j,:)*(a(i,j)/a(j,j));
    end
end
x=zeros(1,m);
for s=m:-1:1
c=0;
    for k=2:m
    c=c+a(s,k)*x(k);
    end
    x(s)=(a(s,n)-c)/a(s,s);
end
disp('Gauss elimination method:');
a
x'





❤❤% Gauss-Jordan method❤❤

%Input section
a=input("Enter the value of matrix[...]");

[m,n]=size(a);
%Start loop
for j=1:m-1
    for z=2:m
        if a(j,j)==0
            t=a(1,:);a(1,:)=a(z,:);
            a(z,:)=t;
        end
    end
    for i=j+1:m
        a(i,:)=a(i,:)-a(j,:)*(a(i,j)/a(j,j));
    end
end

for j=m:-1:2
    for i=j-1:-1:1
        a(i,:)=a(i,:)-a(j,:)*(a(i,j)/a(j,j));
    end
end

for s=1:m
    a(s,:)=a(s,:)/a(s,s);
    x(s)=a(s,n);
end

%Display method
disp('Gauss-Jordan method:');
a
x'





❤❤%Gauss Seidel MEthod❤❤

%input section
A=input('Enter first metrix');
b=input('Enter second matrix');
p=[0; 0; 0];
n=input('enter iteration number');
toll=0.0001;
N=length(b);
x=zeros(N,1);
y=zeros(N,1);

%Loop section
for j=1:n
    for i=1:N
        x(i)=(b(i)/A(i,i))-(A(i,[1:i-1,i+1:N])*p([1:i-1,i+1:N]))/A(i,i);
        p(i)=x(i);
    end
    fprintf('iteration number %d\n',j)
    x
    if abs(x-y)<toll
        break
    end
    y=x
end





❤❤%Bisection Method❤❤

syms x
y=input('Enter nono-linear equation');
a=input('Enter first gauss');
b=input('Enter second gauss');
e=input('Tolerable error');

fa=eval(subs(y,x,a));
fb=eval(subs(y,x,b));

if fa*fb>0
    disp('Given initial values');
else
    c=(a+b)/2;
    fc=eval(subs(y,x,c));
    fprintf('\n\na\t\t\tb\t\t\tc\t\t\tf(c)\n');
    while abs(fc)>e
        fprintf('%f\t%f\t%f\t%f\n',a,b,c,fc);
        if fa*fc<0
            b=c;
        else
            a=c;
        end
        c=(a+b)/2;
        fc=eval(subs(y,x,c));
    end
    fprintf('\nRoot is\n',c);
end



❤❤%Newton Raphson Method❤❤


syms x;

y=input('Enter no linear equation');
a=input('Enter initial gauss');
e=input('Tolrable error');
N=input('Enter maximum number of steps');
step=1;

g=diff(y,x);

fa=eval(subs(y,x,a));

while abs(fa)>e
    fa=eval(subs(y,x,a));
    ga=eval(subs(g,x,a));
    if ga==0
        disp('Division by zero');
        break
    end
    b=a-fa/ga;
    frintf('Step=%d\ta=%f\tf(a)=%f\n',step,a,fa);
    a=b;
    if step>N
        disp('Not covergent');
        break
    end
    step=ste+1;
end
fprintf('Root is %f\n',a);
