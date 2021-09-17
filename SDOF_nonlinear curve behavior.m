clear all;
clc;
%% problem information
%units------>kg & N & m
solution1=1;
solution2=0;
beta=1/6;                                                     %linear acceleratin
gama=1/2;
deltaT=0.01;
ultimateT=16;
m=20000;
k0=(55000*5)*pi;
k1=0.0008;
kesi=0.05;
fy=55000;
c=2*m*((k0/m)^0.5)*kesi;
%% TIME increment
t=zeros((ultimateT/deltaT)+1,1);
for i=1:1:(ultimateT/deltaT)+1
    t(i)=(i-1)*deltaT;
end

%% force
p=zeros((ultimateT/deltaT)+1,1);
for i=1:1:(ultimateT/deltaT)+1
   
   if t(i)<=2
       p(i)=40000*(1+sin((pi/2)*(t(i)+3)));
   
   elseif 2<t(i) && t(i)<4
       p(i)=70000+10000*sin(4*pi*(t(i)+(17/8)));
   elseif 4<t(i) && t(i)<8
       p(i)=0;
   
  
   elseif 8<t(i) && t(i)<10
       p(i)=-40000*(1+sin((pi/2)*(t(i)+3)));
 
   elseif 10<=t(i) && t(i)<12
       p(i)=-70000-10000*sin(4*pi*(t(i)+(17/8)));
 
   elseif t(i)>12
       p(i)=0;
 
   end
end

%% responses
u=zeros((ultimateT/deltaT)+1,1);
ud=zeros((ultimateT/deltaT)+1,1);
udd=zeros((ultimateT/deltaT)+1,1);
deltaU=zeros((ultimateT/deltaT),1);
deltaUd=zeros((ultimateT/deltaT),1);
deltaUdd=zeros((ultimateT/deltaT),1);
fs=zeros((ultimateT/deltaT)+1,1);
k=zeros((ultimateT/deltaT)+1,1);
khat=zeros((ultimateT/deltaT)+1,1);
udd(1)=(p(1)-c*ud(1)-k0*u(1))/m;

a1=(m/(beta*deltaT))+gama*c/beta;
a2=(gama*c*deltaT/(2*beta))+(0.5*m/beta)-c*deltaT;
if solution1==1

for i=1:1:(ultimateT/deltaT)     
    if abs(fs(i))<fy
        k(i)=(55000*5)*pi*cos((5*pi)*u(i));
    else 
        k(i)=0.00008;
    end

    khat(i)=k(i)+m/(beta*deltaT^2)+(c*gama/(beta*deltaT));
    deltaU(i)=(p(i+1)-p(i)+a1*ud(i)+a2*udd(i))/khat(i);
    deltaUd(i)=(1-0.5*gama/beta)*deltaT*udd(i)+(gama/(beta*deltaT))*deltaU(i)-gama*ud(i)/beta;
    deltaUdd(i)=(1/(beta*(deltaT^2)))*deltaU(i)-(1/(beta*deltaT))*ud(i)-0.5*udd(i)/beta;
    u(i+1)=u(i)+deltaU(i);
    ud(i+1)=deltaUd(i)+ud(i);
    fs(i+1)=fs(i)+k(i)*deltaU(i);
   if fs(i+1)>fs(i)
       if fs(i+1)>=fy
           fs(i+1)=fy;
       else
           fs(i+1)=55000*sin(5*pi*u(i+1));
       end
   end    
   if fs(i+1)<fs(i)
       if abs(fs(i+1))>=fy
       fs(i+1)=-fy;
       else
       fs(i+1)=55000*sin(5*pi*u(i+1));
       end
   end
    udd(i+1)=(p(i+1)-c*ud(i+1)-fs(i+1))/m;
end

else
    for i=1:1:(ultimateT/deltaT)
    if abs(fs(i))<fy
        k(i)=550000;
    else
        k(i)=0.0008 ;
    end
    khat(i)=k(i)+m/(beta*deltaT^2)+(c*gama/(beta*deltaT));
    deltaU(i)=(p(i+1)-p(i)+a1*ud(i)+a2*udd(i))/khat(i);
    deltaUd(i)=(1-0.5*gama/beta)*deltaT*udd(i)+(gama/(beta*deltaT))*deltaU(i)-gama*ud(i)/beta;
    deltaUdd(i)=(1/(beta*(deltaT^2)))*deltaU(i)-(1/(beta*deltaT))*ud(i)-0.5*udd(i)/beta;
    u(i+1)=u(i)+deltaU(i);
    ud(i+1)=deltaUd(i)+ud(i);
    fs(i+1)=fs(i)+k(i)*deltaU(i);
   if fs(i+1)>fy
       fs(i+1)=fy;
   end
   if fs(i+1)<-fy
       fs(i+1)=-fy;
   end
    
  
    udd(i+1)=(p(i+1)-c*ud(i+1)-fs(i+1))/m;
    end
end
%% plotting
figure(1)
plot(t,u,'b','LineWidth',2);
xlabel('Time(sec)','FontSize',20);
ylabel('Displacement(m)','FontSize',20);
figure(2)
plot(t,ud,'b','LineWidth',2);
xlabel('Time(sec)','FontSize',20);
ylabel('velocity(m/sec)','FontSize',20);
figure(3)
plot(t,udd,'b','LineWidth',2);
xlabel('Time(sec)','FontSize',20);
ylabel('acceleration(m/sec^2)','FontSize',20);
figure(4)
plot(u,fs,'b','LineWidth',2);
xlabel('u(m)','FontSize',20);
ylabel('fs(N)','FontSize',20);
figure(5)
plot(t,p,'b','LineWidth',2);
xlabel('Time(sec)','FontSize',20);
ylabel('Load(N)','FontSize',20);
% plot(t,p);