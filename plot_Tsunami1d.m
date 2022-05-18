n=801;
for i=1:801
    h(1,i)=round(50-45*tanh(((del_x*(i-1))-70)/8));
end
dt=0.001
h=-h/100;
del_x=0.125;
c=readmatrix('eta_1d.txt');
c=c';
x= 0: 0.125:100;
[a b]= size(c);
% time=[1,2,3];
% for i= 1:a
%          e= c(i,:);
%          p=plot(x,e,x,h);
%          title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(dt*i)]});
%          drawnow;
%          refreshdata(p);
% end
s=c(500,:);g=c(1000,:);o=c(1500,:);l=c(2000,:);
             figure(1);
             plot(x,s,x,h);
             title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(0.5)]});
             figure(2);
             plot(x,g,x,h);
             title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(1)]});
             figure(3);
             plot(x,o,x,h);
             title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(1.5)]});
             figure(4);
             plot(x,l,x,h);
             title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(2)]});
       
