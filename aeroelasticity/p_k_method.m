% identify parameters
clear;
clc;
epsilon=10^-4;
b=2;l=1;
omega_alpha=52;
miu=5;
e=0.4;
x_alpha=0.25;
omega_halpha=0.5;
r=sqrt(0.5);
a1=50;a2=200;a3=0.01;
% record gamma,omega,V data to draw the curve
X=zeros((a2-a1)/a3+1,4);
Z=zeros((a2-a1)/a3+1,1);
for V=a1:a3:a2
y=(omega_alpha*b/V)^2;
k=1;
% calculate the initial k1 and k2
[L_h,L_alpha,M_h,M_alpha]=parameter(k);
A=miu^2*(r^2-x_alpha^2);
B=miu^2*r^2*y*(1+omega_halpha^2)+miu*k^2*x_alpha*(M_h-2*e*L_h+L_alpha)-....
miu*k^2*r^2*L_h-miu*k^2*(M_alpha-e*(L_alpha+M_h)+e^2*L_h);
C=miu^2*r^2*y^2*omega_halpha^2-miu*k^2*omega_halpha^2*y*(M_alpha-e*(M_h+L_alpha)+e^2*L_h)....
-miu*k^2*r^2*L_h*y+k^4*(M_alpha*L_h-L_alpha*M_h);
P_a=(-B+sqrt(B^2-4*A*C))/2/A;
P_b=(-B-sqrt(B^2-4*A*C))/2/A;
p1=sqrt(P_a); 
p2=sqrt(P_b);
k1=abs(imag(p1));
k2=abs(imag(p2));
% calculate the first branch
for m=1:1000
    [L_h,L_alpha,M_h,M_alpha]=parameter(k1);
    A=miu^2*(r^2-x_alpha^2);
   B=miu^2*r^2*y*(1+omega_halpha^2)+miu*k1^2*x_alpha*(M_h-2*e*L_h+L_alpha)...
   -miu*k1^2*r^2*L_h-miu*k1^2*(M_alpha-e*(L_alpha+M_h)+e^2*L_h);
   C=miu^2*r^2*y^2*omega_halpha^2-miu*k1^2*omega_halpha^2*y*(M_alpha-e*(M_h+L_alpha)+e^2*L_h)...
-miu*k1^2*r^2*L_h*y+k1^4*(M_alpha*L_h-L_alpha*M_h);
   P_a=(-B+sqrt(B^2-4*A*C))/2/A;
   P_b=(-B-sqrt(B^2-4*A*C))/2/A;
   p3=sqrt(P_a);
  p4=sqrt(P_b);
% select the positive root
  if(imag(p3)<0)
      p3=-p3;
  end
 if(imag(p4)<0)
      p4=-p4;
  end
  k3=imag(p3);
  k4=imag(p4);
  error1=abs(k1-k3);
  error2=abs(k1-k4);
  if(error1<error2)
      k1=k3; 
      error=error1;
       gamma1=real(p3)/imag(p3);
  else 
      k1=k4;
       gamma1=real(p4)/imag(p4);
      error=error2;
  end
  if(error<epsilon)
      omega1=k1*V/b;
      break
  end
end
% calculate the second branch 
for m=1:1000
    [L_h,L_alpha,M_h,M_alpha]=parameter(k2);
    A=miu^2*(r^2-x_alpha^2);
   B=miu^2*r^2*y*(1+omega_halpha^2)+miu*k2^2*x_alpha*(M_h-2*e*L_h+L_alpha)...
   -miu*k2^2*r^2*L_h-miu*k2^2*(M_alpha-e*(L_alpha+M_h)+e^2*L_h);
   C=miu^2*r^2*y^2*omega_halpha^2-miu*k2^2*omega_halpha^2*y*(M_alpha-e*(M_h+L_alpha)+e^2*L_h)...
-miu*k2^2*r^2*L_h*y+k2^4*(M_alpha*L_h-L_alpha*M_h);
   P_a=(-B+sqrt(B^2-4*A*C))/2/A;
   P_b=(-B-sqrt(B^2-4*A*C))/2/A;
   p3=sqrt(P_a);
  p4=sqrt(P_b);
    if(imag(p3)<0)
      p3=-p3;
    end
    if(imag(p4)<0)
      p4=-p4;
    end
  k3=imag(p3);
  k4=imag(p4);
  error1=abs(k2-k3);
  error2=abs(k2-k4);
  if(error1<error2)
      k2=k3;
      gamma2=real(p3)/imag(p3);
      error=error1;
  else 
      k2=k4;
      gamma2=real(p4)/imag(p4);
      error=error2;
  end
  if(error<epsilon)
      omega2=k2*V/b;
      break
  end
end
X(l,1)=gamma1;X(l,2)=gamma2;
X(l,3)=omega1;X(l,4)=omega2;
Z(l,1)=V;
l=l+1;
end
% draw the V-gamma curve 
figure(1);
plot(Z,X(:,1),'b',Z,X(:,2),'r','LineWidth',2);
% draw the reference line 
line([a1,a2],[0,0],'LineWidth',2,'color','k','Linestyle','--');
% draw the V-omega curve 
figure(2);
plot(Z,X(:,3),'b',Z,X(:,4),'r','LineWidth',2);