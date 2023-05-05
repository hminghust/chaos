% Largest Lyapunov Exponents
Z=[];  
a=10; 
c=8/3; 
d0=1e-7; 
bs = linspace(0,300,301);
transient = 50;
for b=bs
    params = [a,b,c];
    lsum=0; 
    x=1;y=1;z=1;   % initial  value
    x1=1;y1=1;z1=1+d0; 
    for i=1:100  
        [T1,Y1]=ode45(@(t,X) Lorenz(t,X,params),[0 1],[x;y;z]); 
        [T2,Y2]=ode45(@(t,X) Lorenz(t,X,params),[0 1],[x1;y1;z1]); 
        n1=length(Y1);n2=length(Y2);     
        x=Y1(n1,1);y=Y1(n1,2);z=Y1(n1,3);   
        x1=Y2(n2,1);y1=Y2(n2,2);z1=Y2(n2,3);   
        d1=sqrt((x-x1)^2+(y-y1)^2+(z-z1)^2);     
        % Update Location
        x1=x+(d0/d1)*(x1-x);   
        y1=y+(d0/d1)*(y1-y);    
        z1=z+(d0/d1)*(z1-z); 
        if i> transient
            lsum=lsum+log(d1/d0);  
        end  
    end  
    Z=[Z lsum/(i-transient)]; 
end 

