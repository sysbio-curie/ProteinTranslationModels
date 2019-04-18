function dydt = M1_model(t,x)

global k1 k2 k3 kd kt ka k1f k1b kp S40 S60;

M = x(1);
M40S = x(2);
F = x(3);
R = x(4);
P = x(5);

RT = kt;
R1 = k1f*S40*M-k1b*M40S;
RMA = ka*M40S;
R2 = k2*S60*F;
R3 = k3*R;
RMd = kd*M;
RM40Sd = kd*M40S;
RFd = kd*F;
RRd = kd*R;
Rp = kp*P;

dydt = [...
RT-R1-RMd+R2;...
R1-RMA-RM40Sd;...
RMA-R2-RFd;...
R2-R3-RRd;...
R3-Rp...
];
end