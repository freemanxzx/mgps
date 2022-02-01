function [t,y]=gr_slex(x, h, tspan)
    xn = x;
    y = [];
    t=tspan(1):h:tspan(2);
    for tn = t
        y = [y; xn];
        xnplus1 = xnplus1_predictor(xn(1), xn(2), h);
        bar_xn = (xnplus1+xn(1))/2;
        omega = sqrt(cos(bar_xn));
        delta = (2/omega)*tan((h*omega)/2);
        xn_21 = -delta*((-cos(xnplus1) + cos(xn(1)))/(xnplus1-xn(1))) + xn(2);
        xn_11 = (delta/2)*(xn_21 + xn(2)) + xn(1);
        
        xn = [xn_11 xn_21];       
    end
end

function xnplus1=xnplus1_predictor(xn,pn,h)
    Vp = sin(xn);
    Vpp = cos(xn);
    omegan = sqrt(abs(Vpp));
    tmp = omegan*h;
    if Vpp > 0
        xnplus1 = xn + sin(tmp)*pn/omegan - (1-cos(tmp))*Vp/omegan^2;
    elseif Vpp == 0
        xnplus1 = xn + h*pn - h^2*Vp/omegan^2;
    else
        xnplus1 = xn + sin(tmp)*pn/omegan - (1-cos(tmp))*Vp/omegan^2;
    end
    
end

function omegan=compute_omgean(bar_xn, type)
     omegan = sqrt(cos(bar_xn));
end