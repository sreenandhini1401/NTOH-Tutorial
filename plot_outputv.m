function plot_outputv
    a=readmatrix('vvel 1s.txt');
    b=readmatrix('vvel 3s.txt');
    c=readmatrix('vvel 5s.txt');
    d=readmatrix('vvel 7s.txt');
%     e=readmatrix('eta4.0s.txt');
%     f=readmatrix('eta5.0s.txt');
%     g=readmatrix('eta6.0s.txt');
%     d=readmatrix('eta7.0s.txt');
    subplot(2,2,1)
    h=surf(a,'EdgeColor','interp');
    title(['Flow at 1s after dam break']);
    subplot(2,2,2)
    h=surf(b,'EdgeColor','interp');
    title(['Flow at 3s after dam break']);
    subplot(2,2,3)
    h=surf(c,'EdgeColor','interp');
    title(['Flow at 5s after dam break']);
    subplot(2,2,4)
    h=surf(d,'EdgeColor','interp');
    title(['Flow at 7s after dam break']);
%     subplot(3,2,5)
%     h=surf(f,'EdgeColor','interp');
%     subplot(3,2,6)
%     h=surf(g,'EdgeColor','interp');
    
%     subplot(2,2,2)
%     j=surf(b,'EdgeColor','interp');
%    
%     subplot(2,2,3)
%     i=surf(c,'EdgeColor','interp');
    
end