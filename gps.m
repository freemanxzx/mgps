function [t, y]=gps(ode, x, h, tspan, type)
    if strcmp(type, "eeuler")
        disp("eeuler");
        [t, y]=eeuler(ode, x, h, tspan);
    elseif strcmp(type, "ieuler")
        disp("ieuler");
        [t, y]=ieuler(ode, x, h, tspan);
    elseif strcmp(type, "gps")
        disp("gps");
        [t, y]=gps1(ode, x, h, tspan);
    elseif strcmp(type, "enhanced-gps")
        disp("enhanced-gps");
        [t, y]=enhanced_gps(ode, x, h, tspan);
    elseif strcmp(type, "midpoint-gps")
        disp("midpoint-gps");
        [t, y]=midpoint_gps(ode, x, h, tspan);
    elseif strcmp(type, "opt_midpoint_gps")
        disp("opt_midpoint_gps");
        [t, y]=opt_midpoint_gps(ode, x, h, tspan);
    elseif strcmp(type, "gps3")
        disp("gps3");
        [t, y]=gps3(ode, x, h, tspan, false);
    elseif strcmp(type, "cone-gps3")
        [t, y]=gps3(ode, x, h, tspan, true);
        disp("cone-gps3");
    elseif strcmp(type, "gps4")
         disp("gps4");
        [t, y]=gps4(ode, x, h, tspan, false);
    elseif strcmp(type, "cone-gps4")
         disp("cone-gps4");
        [t, y]=gps4(ode, x, h, tspan, true);
    elseif strcmp(type, "gps4-exp")
         disp("gps4-exp");
        [t, y]=gps4_exp(ode, x, h, tspan);
%     elseif strcmp(type, "gps4_exp_harmonic_mean")
%          disp("gps4-exp");
%         [t, y]=gps4_exp_harmonic_mean(ode, x, h, tspan);    
    end
        
function [t, y]=gps1(ode, x, h, tspan)
    xn = x;
    y = [];
    t=tspan(1):h:tspan(2);
    for tn = t
        y = [y; xn]; 
        fn = ode(tn, xn);
%         fn_norm = sqrt(fn' * fn);
%         xn_norm = sqrt(xn * xn');
%         temp = (h * (fn_norm))/(xn_norm);
%         an = cosh(temp);
%         bn = sinh(temp);
%         xn = xn + (((an - 1)*(xn*fn) + bn*xn_norm*fn_norm)/(fn' * fn))*fn';
        %disp((fn')/xn)
        xn = xn + gps_h(xn, fn, h)*fn';
    end    

 
function [t, y]=enhanced_gps(ode, x,h,tspan)
    xn = x;
    y = [];
    t=tspan(1):h:tspan(2);
    for tn = t
        y = [y; xn]; 
        
        fn = ode(tn, xn);
        xn_ = xn + (h*fn')/2;
        fn_ = ode(tn+h/2, xn_);
        
        fn_norm_ = sqrt(fn_' * fn_);
        xn_norm = sqrt(xn * xn');
        xn_norm_ = sqrt(xn_ * xn_');
        
        temp = (h * (fn_norm_))/(xn_norm_);
        an = cosh(temp);
        bn = sinh(temp);
        xn = xn + (((an - 1)*(xn*fn_) + bn*xn_norm*fn_norm_)/(fn_' * fn_))*fn_';
    end

function [t, y]=midpoint_gps(ode, x,h,tspan)
    xn = x;
    y = [];
    t=tspan(1):h:tspan(2);
    for tn = t
        y = [y; xn]; 
        fn = ode(tn, xn);
        xn_ = xn + gps_h(xn, fn, h)*fn';
        xn_ = 0.5*(xn + xn_);
        fn1 = ode(tn+h/2, xn_');
        xn = xn + gps_h(xn, fn, h)*fn1';
    end

function [t, y]=opt_midpoint_gps(ode, x,h,tspan)
    xn = x;
    y = [];
    t=tspan(1):h:tspan(2);
    for tn = t
        y = [y; xn]; 
        fn = ode(tn, xn);
        % case 1
         %xn_ = xn + (h*fn')/2;
        % case 2
        xn_ = xn + gps_h(xn, fn, h/2)*fn';
        fn_ = ode(tn+h/2, xn_);
        fn1 = ode(tn+h/2, xn + gps_h(xn, fn, h/2)*fn_');
        xn = xn + gps_h(xn, fn, h)*fn1';
        
        %fn1 = ode(tn+h/2, xn + gps_h(xn, fn, h/2)*fn');
        %fn1 = ode(tn+h/2, xn + gps_h(xn_, fn_, h/2)*fn_');
        %fn1 = ode(tn+h/2, xn + (h/2)*fn');
        % case 3 midpoint
%         xn_ = xn + gps_h(xn, fn, h)*fn';
%         xn_ = 0.5*(xn + xn_);
%         fn1 = ode(tn+h/2, xn_');
%         xn = xn + gps_h(xn, fn, h)*fn1';
    end


function gh=gps_h(xn, fn, h)
    fn_norm = sqrt(fn' * fn);
    xn_norm = sqrt(xn * xn');
    temp = (h * (fn_norm))/(xn_norm);
    an = cosh(temp);
    bn = sinh(temp);
    
    gh = ((an - 1)*(xn*fn) + bn*xn_norm*fn_norm)/(fn' * fn);
    
    
function [t, y]=gps3(ode, x, h, tspan, is_cone)

tau=h/2;

xn = x;
y = [];
t = tspan(1):h:tspan(2);
for tn = t
    y = [y;xn];
    
    xn1 = xn;   
    f1 = ode(tn,xn1);
    xn_norm = sqrt(xn*xn');
    xn1_norm = sqrt(xn1*xn1');
    c1 = (2*tau^2*f1*xn'+4*tau*xn1_norm*xn_norm)/(4*(xn1*xn1')-tau^2*(f1*f1'));
    xn2 = xn + c1*f1;
    f2 = ode(tn,xn2);
    c2 = (3*tau^2*f2*xn'+6*tau*sqrt(xn2*xn2')*xn_norm)/(4*(xn2*xn2')-tau^2*(f2*f2'));
    f3 = xn + c2*f2;
    
    Fn = (2*f1+3*f2+4*f3)/9;
    if is_cone
        xn_square = xn*xn';
        c3 = (4*h*xn_square+2*h^2*Fn*xn')/(4*xn_square-h^2*(Fn*Fn'));
        xn_1 = xn + c3*Fn;
    else
        xn_1 = xn + h*Fn;
    end
    
%     y = [y;xn_1];
    
    xn = xn_1;
end


function [t,y]=gps4(ode, x, h, tspan, is_cone)
tau=h/2;
xn = x;
y = [];
t=tspan(1):h:tspan(2);

for tn=t
    %disp(tn);
    y = [y;xn];
    
    %¡–æÿ’Û ‰≥ˆ
    f1 = ode(tn, xn);
    xn1 = xn;
    xn_norm = (sqrt(xn*xn'));
    c1 = ((2*tau^2*xn*f1 + 4*tau*sqrt(xn1*xn1')*xn_norm)/(4*(xn1*xn1')-tau^2*(f1'*f1)));
    xn2 = (xn+c1*f1');
    
    f2 = ode(tn, xn2);
    c2 = ((2*tau^2*xn*f2 + 4*tau*sqrt(xn2*xn2')*xn_norm)/(4*(xn2*xn2')-tau^2*(f2'*f2)));
    xn3 = (xn+c2*f2');
    
    f3 = ode(tn, xn3);
    c3 = ((2*tau^2*xn*f3 + 8*tau*sqrt(xn3*xn3')*xn_norm)/(4*(xn3*xn3')-tau^2*(f3'*f3)));
    xn4 = (xn+c3*f3');
    
    f4 = ode(tn, xn4);
    
    Fn = ((f1 + 2*f2 + 2*f3 + f4)/6);
    if is_cone
        xn_square = (xn*xn');
        c4 = ((4*h*(xn_square) + 2*h^2*xn*Fn)/(4*xn_square - h^2*(Fn'*Fn)));
        xn = (xn + c4*Fn');
    else
        xn = (xn + h*Fn');
    end
%     disp("y")
%     disp(y);
%     break;
    % y = [y;xn];
end

function [t, y]=gps4_exp(ode, x, h, tspan)
tau=h/2;
xn = x;
y = [];
t=tspan(1):h:tspan(2);
for tn = t
    y = [y;xn];
    
    f1 = ode(tn, xn);
    
    f1_norm = sqrt(f1'*f1);
    xn_norm = sqrt(xn*xn');
    temp = (tau*f1_norm)/xn_norm;
    a1 = cosh(temp);
    b1 = sinh(temp);
    c1 = ((a1-1)*xn*f1 + b1*f1_norm*xn_norm)/(f1'*f1);
    xn2 = xn + c1*f1';
    f2 = ode(tn + tau, xn);
    f2_norm = sqrt(f2'*f2);
    xn2_norm = sqrt(xn2*xn2');
    temp = (tau*f2_norm)/xn2_norm;
    a2 = cosh(temp);
    b2 = sinh(temp);
    c2 = ((a2-1)*xn*f2 + b2*f2_norm*xn_norm)/(f2'*f2);
    xn3 = xn + c2*f2';
    f3 = ode(tn+tau, xn3);
    f3_norm = sqrt(f3'*f3);
    xn3_norm = sqrt(xn3*xn3');
    temp = (tau*f3_norm)/xn3_norm;
    a3 = cosh(temp);
    b3 = sinh(temp);
    c3 = (2*(a3-1)*xn*f3 + 2*b3*f3_norm*xn_norm)/(f3'*f3);
    xn4 = xn + c3*f3';
    f4 = ode(tn, xn4);
    Fn = (f1 + 2*f2 + 2*f3 + f4)/6;
    Fn_norm = sqrt(Fn'*Fn);
    temp = h*Fn_norm/xn_norm;
    an = cosh(temp);
    bn = sinh(temp);
    xn = xn + (((an-1)*xn*Fn + bn*xn_norm*Fn_norm)/(Fn'*Fn))*Fn';
    
end

%explicit euler method
function [t, y]=eeuler(ode, x, h, tspan)
xn = x;
y = [];
t=tspan(1):h:tspan(2);
for tn = t
    y = [y;xn];
    f1 = ode(tn, xn);
    xn = xn + h*f1';
end

%implicit euler method
function [t, y]= ieuler(ode, x, h, tspan)
xn = x;
y = [];
t=tspan(1):h:tspan(2);
for tn = t
    y = [y;xn];
    f1 = ode(tn, xn);
    xn_1 = xn + h*f1';
    f2 = ode(tn+h, xn_1);
    xn = xn + h*f2'
end

% function [t, y]=gps4_exp_harmonic_mean(ode, x, h, tspan)
% tau=h/2;
% xn = x;
% y = [];
% t=tspan(1):h:tspan(2);
% for tn = t
%     y = [y;xn];
%     
%     f1 = ode(tn, xn);
%     
%     f1_norm = sqrt(f1'*f1);
%     xn_norm = sqrt(xn*xn');
%     temp = (tau*f1_norm)/xn_norm;
%     a1 = cosh(temp);
%     b1 = sinh(temp);
%     c1 = ((a1-1)*xn*f1 + b1*f1_norm*xn_norm)/(f1'*f1);
%     xn2 = xn + c1*f1';
%     f2 = ode(tn + tau, xn);
%     f2_norm = sqrt(f2'*f2);
%     xn2_norm = sqrt(xn2*xn2');
%     temp = (tau*f2_norm)/xn2_norm;
%     a2 = cosh(temp);
%     b2 = sinh(temp);
%     c2 = ((a2-1)*xn*f2 + b2*f2_norm*xn_norm)/(f2'*f2);
%     xn3 = xn + c2*f2';
%     f3 = ode(tn+tau, xn3);
%     f3_norm = sqrt(f3'*f3);
%     xn3_norm = sqrt(xn3*xn3');
%     temp = (tau*f3_norm)/xn3_norm;
%     a3 = cosh(temp);
%     b3 = sinh(temp);
%     c3 = (2*(a3-1)*xn*f3 + 2*b3*f3_norm*xn_norm)/(f3'*f3);
%     xn4 = xn + c3*f3';
%     f4 = ode(tn, xn4);
%     % Fn = (f1 + 2*f2 + 2*f3 + f4)/6;
% %     disp(f1)
% %     disp("Fn")
% %     disp(f1.*f2)
% %     disp(f1+f2)
% %     disp((f1.*f2)/(f1+f2))
% %     Fn = (f1.*f2)/(f1+f2) + (f2.*f3)/(f2+f3) + (f3.*f4)/(f3+f4);
% %     disp(Fn)
%     Fn = (sqrt(abs(f1.*f2)) + sqrt(abs(f2.*f3)) + sqrt(abs(f3.*f4)))/3
%     Fn_norm = sqrt(Fn'*Fn);
%     temp = h*Fn_norm/xn_norm;
%     an = cosh(temp);
%     bn = sinh(temp);
%     xn = xn + (((an-1)*xn*Fn + bn*xn_norm*Fn_norm)/(Fn'*Fn))*Fn';
%     
% end

