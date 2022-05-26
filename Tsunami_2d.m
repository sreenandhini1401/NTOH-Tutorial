clc;
clear all;
x=100;
dx=1;
g=9.81;
n=101;
h=zeros(n,n);
q=0;
for j=n:-1:1
        q=q+1;
        h(q,:)= round(50-49*tanh(-(dx*(j-1)-20)/8));
end
t=150;
dt=1/100;
a= dt/(2*dx);
c= dt/dx;
eta= zeros(n,n);
m=zeros(n,n);
o=zeros(n,n);
D=zeros(n,n);
final=zeros(n,n,t);
x1=0:dx:x;
g=9.81;

for i=2:n-1
    for j=2:n-1
        xl= dx*(i-1);
        yl= dx*(j-1);
        eta(i,j)=5*exp(-(xl-30)^2/10-(yl-50)^2/20);
    end
end
   eta(1,:)=eta(2,:); eta(n,:)=eta(n-1,:);
   eta(:,1)=eta(:,2); eta(:,n)=eta(:,n-1);
    D= h+eta;
    final(:,:,2)=eta;
    
    [un1,un2,un3,F1,F2,F3,G1,G2,G3]=dependent(eta,m,o,g,D);  
    for i=2:n-1
        for j= 2:n-1
           U2(i,j) = un2(i,j) - a * (F2(i,j+1) - F2(i,j-1)) - a*(G2(i+1,j) - G2(i-1,j))-a*g*D(i,j)*(eta(i,j+1)-eta(i,j-1));
           U3(i,j) = un3(i,j) - a * (F3(i,j+1) - F3(i,j-1)) - a*(G3(i+1,j) - G3(i-1,j))-a*g*D(i,j)*(eta(i+1,j)-eta(i-1,j));
        end
    end
                U2(1,:)=U2(2,:);U2(n,:)=U2(n-1,:);U2(:,1)=U2(:,2);U2(:,n)=U2(:,n-1);
                U3(1,:)=U3(2,:);U3(n,:)=U3(n-1,:);U3(:,1)=U3(:,2);U3(:,n)=U3(:,n-1);
                   
    for i=2:n-1
        for j= 2:n-1
              U1(i,j) = un1(i,j) - a * (U2(i,j+1) - U2(i,j-1)) - a*(U3(i+1,j) - U3(i-1,j)); 
        end
    end   
           U1(1,:)=U1(2,:); U1(n,:)=U1(n-1,:);U1(:,1)=U1(:,2); U1(:,n)=U1(:,n-1);
           final(:,:,3)=U1;
           D1= h+U1;
           [un11,un21,un31,F11,F21,F31,G11,G21,G31]=dependent(U1,U2,U3,g,D1);
 
 for k= 3:t
    for i= 2:n-1
        for j= 2:n-1
           u2(i,j) = un2(i,j) - c * (F21(i,j+1) - F21(i,j-1)) - c*(G21(i+1,j) - G21(i-1,j))-c*g*D1(i,j)*(U1(i,j+1)-U1(i,j-1));
           u3(i,j) = un3(i,j) - c * (F31(i,j+1) - F31(i,j-1)) - c*(G31(i+1,j) - G31(i-1,j))-c*g*D1(i,j)*(U1(i+1,j)-U1(i-1,j));  
        end
    end
        u2(1,:)=u2(2,:);u2(n,:)=u2(n-1,:);u2(:,1)=u2(:,2);u2(:,n)=u2(:,n-1);
        u3(1,:)=u3(2,:);u3(n,:)=u3(n-1,:);u3(:,1)=u3(:,2);u3(:,n)=u3(:,n-1);
    for i= 2:n-1
        for j= 2:n-1
            u1(i,j) = un1(i,j) - c * (u2(i,j+1) - u2(i,j-1)) - c*(u3(i+1,j) - u3(i-1,j)); 
          end
    end
        u1(1,:)=u1(2,:); u1(n,:)=u1(n-1,:);u1(:,1)=u1(:,2); u1(:,n)=u1(:,n-1); 
        final(:,:,k+1)= u1;

     ao=u1;b=u2;kl=u3;
     u1=U1;u2=U2;u3=U3;
     U1=ao;U2=b;U3=kl; 

    D=h+U1;
    D1= h+u1;
    [un1,un2,un3,F1,F2,F3,G1,G2,G3]=dependent(U1,U2,U3,g,D);
    [un11,un21,un31,F11,F21,F31,G11,G21,G31]=dependent(u1,u2,u3,g,D1);

 end
% %  h=-h/100;
% %  ui= plot(h);
% %  colormap('default');
% %  hold on;
%  for index=1:t
%     g= final(:,:,index);
%     ci=mesh(g);
%     axis ([0 100 0 100 -5 5]);
%     ci.FaceColor='interp';
%     colorbar;
%     xlabel('X Domain [m]');
%     ylabel('Y Domain [m]');
%     zlabel('Height [m]');
%     drawnow;
%     refreshdata(ci);
%  end
    
  for index=1:t
     g= final(50,:,index);
     ci=plot(g);
     axis ([0 100 0 10 ]);%-100 100]);
     ci.FaceColor='interp';
     colorbar;
     xlabel('X Domain [m]');
     ylabel('Y Domain [m]');
     zlabel('Height [m]');
     drawnow;
     refreshdata(ci);
     pause(0.2);
  end



 
function [un1, un2, un3,F1, F2,F3, G1,G2,G3] = dependent(h,u,v,g,D)
un1 = h;                 
un2 = u; 
un3 = v;

F1 = u;
F2 = (u.^2) ./D; % + (g.*h.^2)/2;
F3 = (u.*v)./D;

G1= v;
G2= (u.*v)./D;
G3 = (v.^2) ./D; %+ (g.*h.^2)/2;

end

 



