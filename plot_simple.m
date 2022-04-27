l=1.0;
n=50;
dx= l/(n-1);
x=0:dx:1;
y=0:dx:1;
u=readmatrix('outputu.txt');
u_=flip(u,2);
v=readmatrix('outputv.txt');
v_=flip(v,2);
p=readmatrix('outputp.txt')
p_=flip(p,2);
figure(1) 
contourf(x,y,u_(:,:)',50, 'edgecolor','none');colormap jet
colorbar;
axis([0 1 0 1]); 
title('steady Ux for Re(2000)'); 
xlabel('Length');
ylabel('Height');
%subplot(2,2,2)
 figure(2)
contourf(x,y,v_(:,:)',50, 'edgecolor','none');colormap jet
colorbar;
axis([0 1 0 1]); 
title('steady Vy for Re(2000)'); 
xlabel('Length');
ylabel('Height');
%subplot(2,2,3)
 figure(3)
contourf(x,y,p_(:,:)',50, 'edgecolor','none');colormap jet
colorbar;
axis([0 1 0 1]); 
title('steady P for Re(2000)'); 
xlabel('Length');
ylabel('Height');
%subplot(2,2,4)
figure(4)
e=u';f=v';
e=flip(e);f=flip(f,2);
streamslice(e,f)
axis([1 51 1 51]); 
title('Streamlines for Re(2000)'); 
xlabel('nx');
ylabel('ny');


