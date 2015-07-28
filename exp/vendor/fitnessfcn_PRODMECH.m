function [ Prod ] = fitnessfcn_PRODMECH ( co, varargin  )

tsim = varargin{1}{1};
Tf = varargin{1}{2};


Xw = co(1);
WFEO = co(2);
af = co(3);
b = co(4);

Tref = 573;
param = [0.626 6.208 2.068 6.042 0.018 0.016 0.230 0.024 0.361 1.919 23.170 16.9 13.740 14.3 35.5 31.2 5.8 17.870 2.138];
paramAc = [0.024 0.361 17.87 2.138 1.919];
x0 = [1,0,0,0];
a = [];
Xo = [];
% spmd
for k=1:size(tsim,2)
    
    %******MECH******
    if k == 1
       %Modelo de conocimiento actividad
       T = min(feval(Tf,1,af,b),673.0);
        kd = paramAc(1)*exp(-paramAc(3)*(1/T-1/Tref)/1.987e-03);
        kdw = paramAc(2)*exp(-paramAc(4)*(1/T-1/Tref)/1.987e-03);
        fda=@(t,a) -kd*exp(-kdw*(Xw+0.643))*a^paramAc(5);
        [ts,ys]=ode45(fda,[0,0.0001],1);
        a=[a ; ys(end)];
        T = min(feval(Tf,a(k),af,b),673.0);
        
        %Modelo de Conocimiento reacción principal                          
        k1 = param(1)*exp(-param(11)*(1/T-1/Tref)/1.987e-03);
        k2 = param(2)*exp(-param(12)*(1/T-1/Tref)/1.987e-03);
        k3 = param(3)*exp(-param(13)*(1/T-1/Tref)/1.987e-03);
        k4 = param(4)*exp(-param(14)*(1/T-1/Tref)/1.987e-03);
        k5 = param(5)*exp(-param(15)*(1/T-1/Tref)/1.987e-03);
        k6 = param(6)*exp(-param(16)*(1/T-1/Tref)/1.987e-03);
        kw = param(7)*exp(-param(17)*(1/T-1/Tref)/1.987e-03);                     
        tita = exp(-kw*(Xw+0.643)); 

        fdx = @(tao,x) [-((k1+k5+(k2+k4)*x(2))*x(1)-k6*x(3))*tita*a(k);(k1*x(1)+(k2-k4)*x(2)*x(1)-k3*x(2)+k6*x(3))*tita*a(k);(k3*x(2)+2*k4*x(1)*x(2)-3*k6*x(3))*tita*a(k);(k5*x(1)+k6*x(3))*tita*a(k)];           
        taospan = [0,WFEO];
        [taos,xs] = ode45(fdx,taospan,x0);
        Xo=[Xo ; xs(end,2)];
    else
        %Modelo de conocimiento actividad
        T = min(feval(Tf,a(k-1),af,b),673.0);
        kd = paramAc(1)*exp(-paramAc(3)*(1/T-1/Tref)/1.987e-03);
        kdw = paramAc(2)*exp(-paramAc(4)*(1/T-1/Tref)/1.987e-03);
        fda=@(t,a) -kd*exp(-kdw*(Xw+0.643))*a^paramAc(5);
        [ts,ys]=ode45(fda,[tsim(k-1),tsim(k)],a(k-1));
        a=[a ; ys(end)];
        T = min(feval(Tf,a(k),af,b),673.0);
        
        %Modelo de Conocimiento reacción principal                          
        k1 = param(1)*exp(-param(11)*(1/T-1/Tref)/1.987e-03);
        k2 = param(2)*exp(-param(12)*(1/T-1/Tref)/1.987e-03);
        k3 = param(3)*exp(-param(13)*(1/T-1/Tref)/1.987e-03);
        k4 = param(4)*exp(-param(14)*(1/T-1/Tref)/1.987e-03);
        k5 = param(5)*exp(-param(15)*(1/T-1/Tref)/1.987e-03);
        k6 = param(6)*exp(-param(16)*(1/T-1/Tref)/1.987e-03);
        kw = param(7)*exp(-param(17)*(1/T-1/Tref)/1.987e-03);                     
        tita = exp(-kw*(Xw+0.643)); 

        fdx = @(tao,x) [-((k1+k5+(k2+k4)*x(2))*x(1)-k6*x(3))*tita*a(k);(k1*x(1)+(k2-k4)*x(2)*x(1)-k3*x(2)+k6*x(3))*tita*a(k);(k3*x(2)+2*k4*x(1)*x(2)-3*k6*x(3))*tita*a(k);(k5*x(1)+k6*x(3))*tita*a(k)];           
        taospan = [0,WFEO];
        [taos,xs] = ode45(fdx,taospan,x0);
        Xo=[Xo ; xs(end,2)];
    end
    if a(k) < 0.10 || Xo(k) < 0.10
        break;
    end
end
% end
%Cálculo del área bajo la curva:
if k > 2
    XoInt = trapz(tsim(1:k), Xo);
else
    XoInt = 0;
end
Prod = (-1)*(XoInt / WFEO);
end