global k2 kd kt ka k1f ks;

% kd = 2e-3;
% kt = kd*1e3;
% k2 = 1e-5;
% k1f = 1.2e-7;
% ka = 1.1e-5;
% kr = 1.3e-3;

% numbers are taken from http://book.bionumbers.org

kd = 8e-4; % mRNA half-life = 20min
%kt = kd*1e3; % number of mRNAs~1000, protein/mRNA ratio is ~100
kt = kd*1e4; % number of mRNAs~1000, protein/mRNA ratio is ~10
k2 = 2e-5; % rate of 60S binding, corresponds to binding of a 60S in 0.1 sec
%k1f = 1e-6; % rate of 40S binding, corresponds to binding of a 40S in 1 sec
k1f = k2/10;
ka = 0.1; % scanning step is very fast
ks = 5e-4; % protein half-life = 30 min


MS = kt/kd;
a = kt/kd*ka;
m = ka/k1f;
n = kd/k2;
p = ka/k2;


disp(sprintf('MS=%e, MS*ka/ks=%e',MS,MS*ka/ks));
disp(sprintf('m=%e, n=%e, p=%e',m,n,p));
disp(sprintf('sqrt(m*n)=%e',sqrt(m*n)));

if 1

rhigh =  (MS*ka/ks - (p+m) + sqrt((MS*ka/ks - (p+m))^2-4*m*n))/2;
rlow = (MS*ka/ks - (p+m) - sqrt((MS*ka/ks - (p+m))^2-4*m*n))/2;
disp(sprintf('rhigh=%e rlow=%e',rhigh,rlow));

tmax = 50000;

deltalog = log10(rhigh*10)-log10(rlow/100);
dr = deltalog/25;
Ri = log10(rlow/100):dr:log10(rhigh*10);
Ri = 10.^Ri;

for i=1:length(Ri)
[t,R] = ode45(@one_component_model,[0 tmax],[Ri(i)]); 
t(1) = 1;
loglog(t,R,'k-'); hold on;
end
ylim([0.5 rhigh*10]);
xlim([100 tmax]);

plot([1 tmax],[rhigh rhigh],'b--','LineWidth',3);
plot([1 tmax],[rlow rlow],'r--','LineWidth',3);
plot([1 tmax],[1 1],'k--','LineWidth',2);

set(gcf,'Color','w');
set(gca,'FontSize',16);
xlabel('seconds','FontSize',20);
ylabel('Amount of ribosomes','FontSize',20);
%legend('Rf','Rp','Rr');

figure;
end


%k2v = [1e-2,1e-3,1e-4,1e-5,1e-6,1e-7];
%k2v = [1e-4,2e-5,1e-5,1e-6];
%k2v = [1e-5];
k2v = [5e-4,1e-4,5e-5,1e-5,5.5e-6];
%k2v = [1e-6,1e-7];

for i=1:length(k2v)
k2 = k2v(i); % rate of 60S binding, corresponds to binding of a 60S in 0.1 sec
k1f = k2/10;

kaa = -4:0.01:1;
kaa = 10.^kaa;

m = kaa/k1f;
n = kd/k2;
p = 0;% kaa/k2;

rhighs =  (MS*kaa/ks - (p+m) + sqrt((MS*kaa/ks - (p+m)).^2-4*m*n))/2;
rlows = (MS*kaa/ks - (p+m) - sqrt((MS*kaa/ks - (p+m)).^2-4*m*n))/2;
MS1 = MS*kaa/ks - kaa/k1f;
%MS2 = (MS*kaa/ks - (p+m)).^2-4*m*n;

MS2 = (MS*1/ks - 1/k1f)*kaa;
MS3 = sqrt(4*kaa/k1f*kd/k2);
MS3_ = m*n./MS2;

ii = find(MS3.*MS3>MS2.*MS2);
rlows(ii) = NaN;
rhighs(ii) = NaN;

loglog(kaa,rhighs,'b-','LineWidth',2); hold on;
loglog(kaa,rlows,'r-','LineWidth',2); hold on;
%loglog(kaa,MS1,'g--','LineWidth',2); hold on;
%loglog(kaa,MS2,'g--','LineWidth',2); hold on;
%loglog(kaa,MS3,'m--','LineWidth',2); hold on;
%loglog(kaa,MS3_,'b--','LineWidth',2); hold on;
set(gcf,'Color','w');
set(gca,'FontSize',16);
xlabel('Scanning rate, k_a','FontSize',20);
ylabel('Amount of ribosomes','FontSize',20);
legend('S_{High}','S_{Low}');
plot([min(kaa) max(kaa)],[1 1],'k-');

plot([min(kaa) max(kaa)],[1e5 1e5],'k--','Color',[0.5 0.5 0.5]);
plot([min(kaa) max(kaa)],[1e6 1e6],'k--','Color',[0.5 0.5 0.5]);

ylim([1e-2 2e6])
end

set(gca,'YTick',[1e-2 1e-1 1 10 100 1000 10000 100000 1000000]);