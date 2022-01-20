load('CookAssignemnt1UnknownCurrent.mat')

%  model parameters
ba = .1;
v0a = -35;
gBar = .01; % mS
Er = 0; % mV
taua = 8; %msec

% plot steady state activation
figure(1);
clf;

subplot(6,2,1);
v = -80 : 20;
plot(v, twoParamSig(v, [ba v0a]), 'LineWidth', 1.5);
hold on;
plot(v, twoParamSigNeg(v, [ba v0a]), 'LineWidth', 1.5);
xlabel('mV');
legend('xa(t = inf)', 'xi(t = inf)', 'Location', 'E');
title('Outward current');

% plot parameters
subplot(6,2,2);
title('Parameters');
set(gca, 'Visible', 'off');
text(0,0,['ba = ' num2str(ba)]);
text(0,1,['v0a = ' num2str(v0a) ' mV']);
text(0,2,['taua = ' num2str(taua) ' msec']);
text(2,0,['Er = ' num2str(Er) ' mV']);
text(2,1,['gBar = ' num2str(gBar) ' mS']);

ylim([0 6]);
xlim([0 4]);    

colors = ['#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'];

% run step

for p = 1:5
    tStop = 110;
    dt = .25;

    ea = 1 - exp(-dt/taua);

    tt = 0;
    j = 1;
    while tt <= tStop

        % update voltage
        v(j) = vStep(j, p);

        % update variables xa, xi
        xaInf = twoParamSig(v(j), [ba v0a]);
        if j > 1
            xa(j) = xa(j-1) + (xaInf - xa(j-1)) * ea;
            xi(j) = 1 - xa(j);
        else
            xa(j) = xaInf;
            xi(j) = 1 - xa(j);
        end

        % update conductance, current, and time
        g(j) = gBar * xa(j)^2 * xi(j)^2;
        i(j) = g(j) * (v(j) - Er);
        t(j) = tt;

        tt = tt + dt;
        j = j + 1;
    end

    % plot time course
    subplot(6,1,2);
    plot(t,v, 'LineWidth', 1.5);
    hold on;
    axis tight;
    ylabel('V step (mV)');

    subplot(6,1,3);
    plot(t,xa, 'LineWidth', 1.5);
    hold on;
    axis tight;
    ylabel('xa');
    
    subplot(6,1,4);
    plot(t,xi, 'LineWidth', 1.5);
    hold on;
    axis tight;
    ylabel('xi');
    ylim([0 inf]);

    subplot(6,1,5);
    plot(t,g, 'LineWidth', 1.5);
    hold on;
    axis tight;
    ylabel('g (mS)');
    ylim([0 inf]);

    subplot(6,1,6);
    plot(t,i, 'LineWidth', 1.5, 'Color', colors(p,:));
    hold on;
    plot(t, iUnknownCurrent(:,p), '--', 'LineWidth', 0.75, 'Color', colors(p,:));
    axis tight;
    xlabel('msec');
    ylabel('i (mA)');
end



function y = twoParamSig(x, params)

    b = params(1);
	x0 = params(2);
	
	y = 1 ./ (1 + exp(-b .* (x - x0)));
end

function y = twoParamSigNeg(x, params)

    b = params(1);
	x0 = params(2);
	
	y = 1 ./ (1 + exp(b .* (x - x0)));
end