function dRdt = one_component_model(t,R)

global k2 kd kt ka k1f ks;

a = kt/kd*ka;
m = ka/k1f;
n = kd/k2;
p = ka/k2;

dRdt = a*R*R/(R*R+(p+m)*R+m*n)-ks*R;
end