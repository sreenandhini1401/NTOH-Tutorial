n=801;
del_x=0.125;
for i=1:801
    h(1,i)=round(50-45*tanh(((del_x*(i-1))-70)/8));
end
dt=0.001
h=-h/100;

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
s=c(500,:);g=c(1000,:);o=c(1500,:);l=c(2000,:);b=c(3000,:);
             figure(1);
             plot(x,s,x,h);
             xlabel('Distance from offshore(m)');
             ylabel('Height(m)');
             legend('Height','Bottom profile');
             title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(0.5)]});
             figure(2);
             plot(x,g,x,h);
             xlabel('Distance from offshore(m)');
             ylabel('Height(m)');
             legend('Height','Bottom profile');
             title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(1)]});
             figure(3);
             plot(x,o,x,h);
             xlabel('Distance from offshore(m)');
             ylabel('Height(m)');
             legend('Height','Bottom profile');
             title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(1.5)]});
             figure(4);
             plot(x,l,x,h);
             xlabel('Distance from offshore(m)');
             ylabel('Height(m)');
             legend('Height','Bottom profile');
             title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(2)]});
             figure(5);
             plot(x,b,x,h);
             xlabel('Distance from offshore(m)');
             ylabel('Height(m)');
             legend('Height','Bottom profile');
             title({['Tsunami propagation in seconds'];['time(\itt) = ',num2str(3)]});
       
