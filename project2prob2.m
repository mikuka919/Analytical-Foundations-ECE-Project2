clear all
beta = 0.5 * 10^(-5);%damping ratio
J = 2 * 10^(-4);%moment on inertia
K = 0.029;%motor constant
L = 0.01;%inductor
R = 3.38;%ohms
y = logspace(-2,4,10000);
s = j * y;
H = K ./ (K^2 - beta .*(R + s.*L)- J.*s.*(R+s.*L));
Y = 20 * log10(abs(H));
semilogx(y,Y);
xlabel('frequency')
ylabel('log scale mag of transfer function')



%dc term = 35.2
%analyze bode plot to find freq where output is 1% of dc term
% find freq where output <= .35
omega_switch = 120;
f_switch = omega_switch / (2*pi);


%dc freq = 0
s_dc = 0;
H_dc = K / (K^2 + beta *(R + s_dc*L)+ J*s_dc*(R+s_dc*L));%33.8035
q = 20 * log10(H_dc/100);%target for switching freq = -9.42
%dc output = v_in * H_dc, v_in = dc_out/H_dc
%target output = 324
T = 2*pi/omega_switch;%0.0546
V_in_324 = 324 / H_dc% = 9.2072
D = V_in_324 / 12 %0.7987

%part f
D_f = 0.5;
a_f = -j * omega_switch;
a_f_pos = j*omega_switch;
alpha_f = 12/(a_f*T) * (exp(a_f*D_f*T)-1);
H_switch =  K / (K^2 + beta *(R + a_f_pos*L)+ J*a_f_pos*(R+a_f_pos*L));
f_out = alpha_f * H_switch;
c1_mag = 4 * abs(f_out)%2.5652, speed variation
c1_phase = angle(f_out);%left in rad form, = -8.562 deg
%f_signal = c1_mag * cos(omega_switch)

%part g
delta_t = T/10000;
t = [0:T/10000:10];
Va = 6*square(omega_switch*t,D*100) + 6;%Va(t) for 10 seconds
%plot(t,Va)
i_a = zeros(size(t));
ohm = zeros(size(t));
alpha_g = 12/(T*a_f) * (exp(a_f*D*T)-1);
g_out = alpha_g * H_switch;

for i = 1:1:length(t)-1
    if i == 1
        i_a(1) = 0;
        ohm(1) = 0;
        %numeric sols to ia and ohm diff eqs from earlier
        i_a(i+1) = delta_t * (Va(i)/L - K*ohm(i)/L - R/L * i_a(i)) + i_a(i);
        ohm(i+1) = delta_t * (K/J * i_a(i) - beta/J * ohm(i)) + ohm(i);
    else
        i_a(i+1) = delta_t * (Va(i)/L - K*ohm(i)/L - R/L * i_a(i)) + i_a(i);
        ohm(i+1) = delta_t * (K/J * i_a(i) - beta/J * ohm(i)) + ohm(i);
    end
end
plot(t,ohm)
xlabel('Time(s)')
ylabel('ohm(t)')

%part h
goal = 20*log10(H_dc/10000);%-49.42, find freq at this val on bode plot
f_switch2 = 2060;
a_h = -j*f_switch2;
T2 = 2*pi/f_switch2;
delta_t2 = T2/10000;
t2 = [0:T2/10000:10];
Va2 = 6*square(f_switch2*t2,D*100) + 6;%Va(t) for 10 seconds
%plot(t,Va)
i_a2 = zeros(size(t));
ohm2 = zeros(size(t));

for k = 1:1:length(t2)-1
    if k == 1
        i_a2(1) = 0;
        ohm2(1) = 0;
        i_a2(k+1) = delta_t2 * (Va2(k)/L - K*ohm2(k)/L - R/L * i_a2(k)) + i_a2(k);
        ohm2(k+1) = delta_t2 * (K/J * i_a2(k) - beta/J * ohm2(k)) + ohm2(k);
    else
        i_a2(k+1) = delta_t2 * (Va2(k)/L - K*ohm2(k)/L - R/L * i_a2(k)) + i_a2(k);
        ohm2(k+1) = delta_t2 * (K/J * i_a2(k) - beta/J * ohm2(k)) + ohm2(k);
    end
end
%plot(t2,ohm2)
xlabel('Time(s)')
ylabel('ohm(t)')


