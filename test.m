test_energy_conservation();

function dy=ode(t, x)
    dy1 = x(2);
    dy2 = -sin(x(1));
    dy = [dy1;dy2];
end

function energy=energy_conservation(y)
        energy = 0.5*y(:,1).^2-cos(y(:,2));
end

function test_energy_conservation()

    tspan=[0 200];
    x=[0 1.8];

    %h=0.25;
    h=0.001;
    
    [t,y]=gps(@ode, x, h, tspan, "midpoint-gps");
    energy=energy_conservation(y);
    delta_energy=energy - energy(1,:);
    plot(t, delta_energy, '-r', 'DisplayName', 'MGPS', 'linewidth', 2);
    hold on;

    [t,y]=gr_slex(x, h, tspan);
    
    energy=energy_conservation(y);
    delta_energy=energy - energy(1,:);
    plot(t, delta_energy, '-.k', 'DisplayName', 'GR_SLEX','linewidth', 1);
    hold on;
    
    xlabel('$t$', 'Interpreter', 'LaTex');
    ylabel('$H(x)-H(0)$', 'Interpreter', 'LaTex');
    legend('MGPS', 'GR\_SLEX');
end