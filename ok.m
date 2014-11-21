
close all;
xite = 0.25;
alfa = 0.05;
s = 2;
IN = 4;H = 5;Out = 3;
if s == 1
Wi = 0.50*rand(H,IN);
Wi_1 = Wi;Wi_2=Wi;Wi_3=Wi;
Wo=0.50*rand(Out,H);
Wo_1 = Wo;Wo_2=Wo;Wo_3=Wo;
end

if s == 2
Wi = 0.50*rand(H,IN);
Wi_1 = Wi;Wi_2=Wi;Wi_3=Wi;
Wo=0.50*rand(Out,H);
Wo_1 = Wo;Wo_2=Wo;Wo_3=Wo;
end

x=[0,0,0];
u_1=0;u_2=0;u_3=0;u_4=0;u_5=0;
y_1=0;y_2=0;y_3=0;

Oh = zeros(H,1);
l=Oh;
error_2=0;error_1=0;
 
ts=0.001;
for k=1:1:6000
    time(k) = k*ts;
    
    if(s==1)
        rin(k)=1.0;
    else
        rin(k) = sin(1*2*pi*k*ts);
    end
    
    a(k)=1.2*(1-0.8*exp(-0.1*k));
    yout(k) = a(k)*y_1/(1+y_1^2)+u_1;
    error(k) = rin(k)-yout(k);
    
    xi=[rin(k),yout(k),error(k),1];
    
    x(1) = error(k)-error_1;
    x(2) = error(k);
    x(3) = error(k)-2*error_1+error_2;
    
    epid = [x(1);x(2);x(3)];
    I = xi*Wi';
    for j=1:1:H
        Oh(j) = (exp(I(j))-exp(-I(j))/(exp(I(j)))+exp(-I(j)));
    end
    K = Wo*Oh;
    for ll = 1:1:Out
        K(ll)=exp(K(ll))/(exp(K(ll))+exp(-K(ll)));
    end
    kp(k) = K(1);ki(k)=K(2);kd(k)=K(3);
    Kpid = [kp(k),ki(k),kd(k)];
    
    du(k) = Kpid*epid;
    u(k) = u_1+du(k);
    u(k) = max(min(10,u(k)),-10);
    
    dyu(k)=sign((yout(k)-y_1)/(u(k)-u_1+0.0000001));
    
    for j=1:1:Out
        dK(j) = 2/(exp(K(j))+exp(-K(j)))^2;
    end
    for j=1:1:Out
        delta3(j)=error(k)*dyu(k)*epid(j)*dK(j);
    end
    for j=1:1:Out
        for i = 1:1:H
            d_Wo=xite*delta3(j)*Oh(i)+alfa*(Wo_1-Wo_2);
        end
    end
    Wo = Wo_1+d_Wo+alfa*(Wo_1-Wo_2);
    for i = 1:1:H
        dO(i)=4/(exp(I(i))+exp(-I(i)))^2;
    end
    segma=delta3*Wo;
    for i=1:1:H
        delta2(i)=dO(i)*segma(i);
    end
    
    d_Wi=xite*delta2'*xi;
    Wi=Wi_1+d_Wi+alfa*(Wi_1-Wi_2);
    
    u_5=u_4;u_4=u_3;u_3=u_2;u_2=u_1;u_1=u(k);
    y_2=y_1;y_1=yout(k);
    
    Wo_3=Wo_2;
    Wo_2=Wo_1;
    Wo_1=Wo;
    
    Wi_3=Wi_2;
    Wi_2=Wi_1;
    Wi_1=Wi;
    
    error_2=error_1;
    error_1=error(k);
end
figure(1);
plot(time,rin,'r',time,yout,'b');
xlabel('time(s)');ylabel('rin,yout');
figure(2);
plot(time,error,'r');
xlabel('time(s)');ylabel('error');
figure(3);
plot(time,u,'r');
xlabel('time(s)');ylabel('u');
figure(4);
subplot(311);
plot(time,kp,'r');
xlabel('time(s)');ylabel('kp');
subplot(312);
plot(time,ki,'g');
xlabel('time(s)');ylabel('ki');
subplot(313);
plot(time,kd,'r');
xlabel('time(s)');ylabel('kd');