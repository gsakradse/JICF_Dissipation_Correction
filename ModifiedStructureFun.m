function [diss, eta, plotstr] = ModifiedStructureFun(Sdiss)

al1 = 1.2670;
al2 = -0.02795;
be1 = -0.2737;
be2 = -0.1545;

del = 9.525e-4; %resolution of Ycol in m

nu = 1.562e-5; % m^2/s value at 25 C.

x = linspace(0,50,1000);

eq1 = al1*exp(al2*x) + be1*exp(be2*x);

eq2 = ((Sdiss * del^4)/ nu^3)*x.^(-4);

% Find interesction

s = eq1 - eq2;
ix = find(s > -0.025 & s < 0.025);
xInt = x(ix);
y1Int = eq1(ix);
y2Int = eq2(ix);

mXin = mean(xInt);
my1in = mean(y1Int);
my2in = mean(y2Int);

myInt = mean([my1in, my2in]);

diss = Sdiss/myInt;
eta = del/mXin;

plotstr.x = x;
plotstr.eq1 = eq1;
plotstr.eq2 = eq2;
plotstr.myInt = myInt;
plotstr.myXin = mXin;

%%
%{
inds = [1 5 10 15 16 17 18 20 25 30];
cols = pc3(length(inds));

fpos = FigPosition([700 500]);


figure('pos', fpos)
for ii = inds
    plot(x,eq1(:,ii),'--')
    hold on
    plot(x,eq2(:,ii),'-')
    plot(mXin,myInt,'o', 'MarkerSize', 8)
    
end
hold off
grid on

ylim([0 1.2])
xlim([0 60])
%}
