%% ESM final project- stommel thermohaline model under different time steps, forcing regimes, numerical methods
%Hudson Moss
clear all
close all

%no climate forcing
%delta t= .01
for i=1:1
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;

f= @(t,T,S) T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all

figure(2)
%delta t= .1
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;

f= @(t,T,S) T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all

figure(3)
%delta t= 1
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;

f= @(t,T,S) T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all

figure(4)
%delta t= 10
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;

f= @(t,T,S) T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
end

%+.5 C per time step
%delta t= .01
for i=1:1
figure(5)
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=.5+T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=.5+T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=.5+T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;

f= @(t,T,S) .5+T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=.5+T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all

figure(6)
%delta t= .1
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=.5+T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=.5+T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=.5+T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;

f= @(t,T,S) .5+T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=.5+T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all

figure(7)
%delta t= 1
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=.5+T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=.5+T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=.5+T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;

f= @(t,T,S) .5+T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=.5+T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all

figure(8)
%delta t= 10
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=.5+T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=.5+T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=.5+T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;

f= @(t,T,S) .5+T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=.5+T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
end

%+1 C per time step
%delta t= .01
for i=1:1
figure(9)
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=1+T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=1+T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=1+T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;

f= @(t,T,S) 1+T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.01;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=1+T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=.01)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all

figure(10)
%delta t= .1
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=1+T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=1+T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=.5+T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;

f= @(t,T,S) 1+T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=.1;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=1+T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=.1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all

figure(11)
%delta t= 1
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=1+T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=1+T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=1+T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;

f= @(t,T,S) .5+T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=1;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=1+T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=1)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all

figure(12)
%delta t= 10
for i=1:1 %forward euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;

for i=1:(length(t)-1)
    t(i+1)=t(i)+delta_t;
    T(i+1)=1+T(i)+((c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t);
    S(i+1)=S(i)+((d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t);      
end
    subplot(2,2,1)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Forward Euler Scheme, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %center euler
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;
  
for j=1:(length(t)-1)
    t(j+1)=t(j)+delta_t;

    T(j+1)=1+T(j)+((c.*(tau-T(j))-abs(2.*q).*T(j)).*delta_t);
    S(j+1)=S(j)+((d.*(s-S(j))-abs(2.*q).*S(j)).*delta_t);
     
    T(j+2)=1+T(j+1)-T(j)/((2*delta_t));
    S(j+2)=S(j+1)-S(j)/((2*delta_t));

end
    subplot(2,2,2)
    plot(t,S(1:length(S)-1),'r')
    hold on
    plot(t,T(1:length(T)-1),'b')
    title('Stommel Thermohaline Model, Center Euler Scheme, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
for i=1:1 %%runge-kutta 

c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;

f= @(t,T,S) 1+T+(c.*(tau-T)-abs(2.*q)).*T;
g= @(t,T,S) S+(d.*(s-S)-abs(2.*q)).*S;

for i=1:(length(t)-1)
    
    k1=f(t(i),T(i),S(i));
    l1=g(t(i),T(i),S(i));
    
    k2=f(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    l2=g(t(i)+delta_t/2,(T(i)+0.5*k1*delta_t),(S(i)+(0.5*l1*delta_t)));
    
    k3=f(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    l3=g(t(i)+delta_t/2,(T(i)+0.5*k2*delta_t),(S(i)+(0.5*l2*delta_t)));
    
    k4=f(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    l4=g(t(i)+delta_t,(T(i)+k3*delta_t),(S(i)+l3*delta_t));
    
    T(i+1)=T(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    S(i+1)=S(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    
end

    subplot(2,2,3)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Runge-Kutta Scheme, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')

end
clear all
for i=1:1 %predictor-corrector method
c=1;    %temperature transfer coefficient 
d=1;    %salinity transfer coefficient 
T(1)=9;   %temperature of system
S(1)=35.1;   %salinity of system
tau=7; %external temperature (constant)
s=36;    %external salinity (constant)
q=1;    %flow 

delta_t=10;
t=0:delta_t:100;

for i=1:(length(t)-1)

    k1_T=(c.*(tau-T(i))-abs(2.*q).*T(i)).*delta_t;
    k1_S=(d.*(s-S(i))-abs(2.*q).*S(i)).*delta_t;
   
    k2_T=(c.*(tau-(T(i)+k1_T))-abs(2.*q).*(T(i)+k1_T)).*delta_t;
    k2_S=(d.*(s-(S(i)+k1_S))-abs(2.*q).*(S(i)+k1_S)).*delta_t;
 
    T(i+1)=1+T(i)+0.5.*(k1_T+k2_T); 
    S(i+1)=S(i)+0.5.*(k1_S+k2_S);
    
end

%III. Plot it all 

    subplot(2,2,4)
    plot(t,S,'r')
    hold on
    plot(t,T,'b')
    title('Stommel Thermohaline Model, Predictor Corrector, delta_t=10)')
    legend({'Salinity','Temperature'})
    xlabel('t')
    ylabel('Degrees C & g/kg Respectively')
end 
clear all
end
