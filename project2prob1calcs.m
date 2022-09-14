clear all
%freq, period, time vector
n = 2;%which harmonic
omega_0 = 2*pi*10000;%fund harmonic
T = 10^(-4);%period
t = [-T/2:T/10000:T/2];%period interval


%model of n(t)
x = sawtooth(omega_0*t+pi,1/2);%creates n(t)
%plot(t,x)
%grid on

%fourier transform of n(t)
a = -j * n * omega_0;%change of variable
X = x.*exp(a.*t);%function to be integrated
alpha_n = (1/T) * trapz(t,X);%complex coeff of fourier transform
abs(alpha_n)


%find vals for c_n, phi_n
c_mag = 2 * abs(alpha_n);%magnitude
c_phase = angle(alpha_n)*180/pi;%phase shift
c_freq = n * omega_0/(2*pi);


%transfer functions
omega_plot = [0:0.01:omega_0];
omega_c = 15000;
s = j*omega_plot;
H2 = (omega_c^2)./(s.^2 + (2/sqrt(2))*s*omega_c + omega_c^2);
H4 = (omega_c^4)./((s.^2 + s.*0.7654*omega_c + omega_c^2).*(s.^2 + s.*1.8478*omega_c + omega_c^2));

%test mag and phase of filters w/ different omega_c
figure(2);
subplot(2,1,1);
semilogx(omega_plot,abs(H4))%check mag at 2*pi*10000
xlabel('Omega');
ylabel('Mag of H4');
subplot(2,1,2)
semilogx(omega_plot,angle(H4)*180/pi)%check phase at 2*pi*60
xlabel('Omega');
ylabel('Angle of H4');



%H4 vector for different omega values
omega = [2*pi*60 omega_0 2*omega_0 3*omega_0];
H4_vec = [];
V_in = [];
for k = 1:1:length(omega)
    a2 = -j * omega(k);
    s2 = j * omega(k);
    if omega(k) == 2*pi*60
        V_in(k) = exp(-j*0);
    else
        T = 10^(-4);
        t = [-T/2:T/10000:T/2];
        X2 = x.*exp(a2.*t);
        alpha_n2 = (1/T) * trapz(t,X2);%complex coeff of fourier transform
        V_in(k) = alpha_n2;
    end
    H4_vec = [H4_vec (omega_c^4)/((s2^2 + s2*0.7654*omega_c + omega_c^2)*(s2^2 + s2*1.8478*omega_c + omega_c^2))];
end
check = abs(V_in(1))


%V_out Vector
V_out = [];
for q = 1:1:length(omega)
    V_out = [V_out H4_vec(q)*V_in(q)];
end

V_out_60Hz = 2 * abs(V_out(1))*cos(omega(1)*t + angle(V_out(1)));
V_outHarm1 = 2 * abs(V_out(2))*cos(omega(2)*t + angle(V_out(2)));
V_outHarm2 = 2 * abs(V_out(3))*cos(omega(3)*t + angle(V_out(3)));
V_outHarm3 = 2 * abs(V_out(4))*cos(omega(4)*t + angle(V_out(4)));
%plot(t2,V_out_60Hz)
2 * abs(V_out(1));

abs(V_out(2))








