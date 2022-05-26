function oe21d032_plot(n,Re,caseno)
%n=51;
length=1;
dx=length/50;
xe=0:dx:1;
ye=0:dx:1;
if caseno==1
         filename1=  'z1.csv';
          filename2 =  's1.csv';
          filename3=  'u1.csv';
          filename4=   'v1.csv';
           filename5=  'P1.csv';

elseif (caseno==2)

             filename1=  'z2.csv';
          filename2 =  's2.csv';
          filename3=  'u2.csv';
          filename4=   'v2.csv';
           filename5=  'P2.csv';

elseif (caseno==3) 

          filename1=  'z3.csv';
          filename2 =  's3.csv';
          filename3=  'u3.csv';
          filename4=   'v3.csv';
           filename5=  'P3.csv';

elseif (caseno==4) 

          filename1=  'z4.csv';
          filename2 =  's4.csv';
          filename3=  'u4.csv';
          filename4=   'v4.csv';
          filename5=  'P4.csv';
end 

c=readmatrix(filename1);
d=readmatrix(filename2);
e=readmatrix(filename3);
f=readmatrix(filename4);
g=readmatrix(filename5);

c=flip(flip(c),2);
d=flip(flip(d),2);
e=flip(flip(e),2);
f=flip(flip(f),2);
g=flip(flip(g),2);

figure(1);

contourf(c,[-15:4:15]),xlabel('nx'),...
ylabel('ny');axis('square','tight');colorbar
title( {['Vorticity for Re=',num2str(Re)]} );

figure(2) 
contourf(d,[-0.2:0.01:0.2],'LineColor','none'),xlabel('nx'),...
ylabel('ny'),title({['Stream function for Re=',num2str(Re)]} );axis('square','tight');colorbar

figure(3)
contourf(g,[-550:20:500],'LineColor','none'),xlabel('nx'),...
ylabel('ny'),title({['Pressure distribution for Re=',num2str(Re)]});axis('square','tight');colorbar


figure(4)
contourf(e,[-1:0.1:1],'LineColor','none'),xlabel('nx'),...
ylabel('ny'),title({['u-velocity for Re=',num2str(Re)]});axis('square','tight');colorbar


figure(5)
contourf(f,[-1:0.1:1],'LineColor','none'),xlabel('nx'),...
ylabel('ny'),title({['v-velocity for Re=',num2str(Re)]});axis('square','tight');colorbar



