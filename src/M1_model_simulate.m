global k1 k2 k3 kd kt ka k1f k1b kp S40 S60;

k1 = 0.001;
k2 = 1e-5; % rate of 60S binding, corresponds to binding of a 60S in 0.1 sec
k3 = 0.01;
kd = 8e-4; % mRNA half-life = 20min
kt = kd*1e3; % number of mRNAs~1000, protein/mRNA ratio is ~100
ka = 0.1; % scanning step is very fast
%k1f = 1e-6; % rate of 40S binding, corresponds to binding of a 40S in 1 sec
k1f = k2/10;
k1b = k1f/100;
kp = 1e-3; % protein half-life = 30 min
S40 = 1e5;
S60 = 5e5;

%M = x(1);
%M40S = x(2);
%F = x(3);
%R = x(4);
%P = x(5);

MS = kt/kd*(k2*S60+kd)/(k1f*k2/(k1b+ka+kd)*S60*S40+kd+k2*S60+k1f*(ka+kd)/(k1b+ka+kd)*S40);
disp(sprintf('MS=%f',MS));
M40SS = kt*k1f/kd*S40*(k2*S60+kd)/(k1f*k2*S40*S60+kd*(k1b+ka+kd)+k2*(k1b+ka+kd)*S60+k1f*(ka+kd)*S40);
disp(sprintf('M40SS=%f',M40SS));
FS = kt*k1f*ka/kd*S40/(k1f*k2*S40*S60+kd*(k1b+ka+kd)+k2*(k1b+ka+kd)*S60+k1f*(ka+kd)*S40);
disp(sprintf('FS=%f',FS));
RS = kt/kd*ka/(k3+kd)*S40*S60/(S40*S60+kd*(k1b+ka+kd)/k1f*k2+(k1b+ka+kd)/k1f*S60+(ka+kd)/k2*S40);
disp(sprintf('RS=%f',RS));
PS = kt/kd*k3/(k3+kd)*ka/kp*S40*S60/(S40*S60+kd*(k1b+ka+kd)/k1f*k2+(k1b+ka+kd)/k1f*S60+(ka+kd)/k2*S40);
disp(sprintf('PS=%f',PS));

tmax = 10000;

[t,x] = ode45(@M1_model,[0 tmax],[0 0 0 0 0]); 
t(1)=1;
loglog(t,x); hold on;
loglog([min(t) max(t)],[MS MS],'k--');
loglog([min(t) max(t)],[M40SS M40SS],'k--');
loglog([min(t) max(t)],[FS FS],'k--');
loglog([min(t) max(t)],[RS RS],'k--');
loglog([min(t) max(t)],[PS PS],'k--');

n = size(x,1);

disp(sprintf('Ms=%f,M40Ss=%f,Fs=%f,Rs=%f,Ps=%f',x(n,1),x(n,2),x(n,3),x(n,4),x(n,5)));

legend('M','M40S','F','R','P');

set(gcf,'Color','w');
set(gca,'FontSize',16);
xlabel('seconds','FontSize',20);
ylabel('Amount','FontSize',20);

