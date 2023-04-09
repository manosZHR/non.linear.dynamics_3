clear; close all; clc

ct_1=[1 2 1 3];
ct_2=[3 2 4 3];
syms k1 k2 k3 k4 x y

xdot_1 = ct_1(1)*x-ct_1(1)*x^2-ct_1(2)*x*y;
ydot_1 = ct_1(3)*y-ct_1(3)*y^2-ct_1(4)*x*y;

xdot_2 = ct_2(1)*x-ct_2(1)*x^2-ct_2(2)*x*y;
ydot_2 = ct_2(3)*y-ct_2(3)*y^2-ct_2(4)*x*y;

[solx_1,soly_1] = solve(xdot_1 == 0, ydot_1 == 0);
[solx_2,soly_2] = solve(xdot_2 == 0, ydot_2 == 0);


J1 = jacobian([xdot_1,ydot_1],[x,y]);
J2 = jacobian([xdot_2,ydot_2],[x,y]);

J1 = matlabFunction(J1);
J2 = matlabFunction(J2);

A1 = J1(double(solx_1(1)),double(soly_1(1)));
A2 = J1(double(solx_1(2)),double(soly_1(2)));
A3 = J1(double(solx_1(3)),double(soly_1(3)));
A4 = J1(double(solx_1(4)),double(soly_1(4)));


B1 = J2(double(solx_2(1)),double(soly_2(1)));
B2 = J2(double(solx_2(2)),double(soly_2(2)));
B3 = J2(double(solx_2(3)),double(soly_2(3)));
B4 = J2(double(solx_2(4)),double(soly_2(4)));

lambda_A = [eig(A1) eig(A2) eig(A3) eig(A4)];
lambda_B = [eig(B1) eig(B2) eig(B3) eig(B4)];


tspan = linspace(0,10,100);
x0 = linspace(0,1.5,20);
y0 = linspace(0,1.5,20);


figure(1)
grid on; box on; axis tight
scale=1;
for i = 1:length(x0)
    for j = 1:length(y0)
        in_cond = [x0(i) y0(j)];
        [tout,yout] = ode89(@(t,x) dyn_s_o(t,x,ct_1),tspan,in_cond);
        subplot(1,2,1);
        p1=plot(yout(:,1),yout(:,2),'k');
        hold on
        x=yout(:,1); y=yout(:,2); xx=gradient(x); yy=gradient(y);
        subplot(1,2,2);
        p2=quiver(x,y,xx,yy,scale,'k');
        hold on
    end
end
subplot(1,2,1);
p3=plot(solx_1,soly_1,'or','LineWidth',3);
legend(p3,{'static points'});
xlabel('x'); ylabel('y'); title('phase field')


subplot(1,2,2);
p4=plot(solx_1,soly_1,'or','LineWidth',3);
legend(p4,{'static points'});
axis( [x0(1) x0(end) y0(1) y0(end)] );
xlabel('x'); ylabel('y'); title('phase field (quiver plot)')
hold off

figure(2)
grid on; box on; axis tight
scale=1;
for i = 1:length(x0)
    for j = 1:length(y0)
        in_cond = [x0(i) y0(j)];
        [tout,yout] = ode89(@(t,x) dyn_s_o(t,x,ct_2),tspan,in_cond);
        subplot(1,2,1);
        p5=plot(yout(:,1),yout(:,2),'k');
        hold on
        x=yout(:,1); y=yout(:,2); xx=gradient(x); yy=gradient(y);
        subplot(1,2,2);
        p6=quiver(x,y,xx,yy,scale,'k');
        hold on
    end
end
subplot(1,2,1);
p7=plot(solx_2,soly_2,'or','LineWidth',3);
legend(p7,{'static points'});
xlabel('x'); ylabel('y'); title('phase field')


subplot(1,2,2);
p8=plot(solx_2,soly_2,'or','LineWidth',3);
legend(p8,{'static points'});
axis( [x0(1) x0(end) y0(1) y0(end)] );
xlabel('x'); ylabel('y'); title('phase field (quiver plot)')
hold off





