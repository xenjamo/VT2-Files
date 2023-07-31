



w_gps = 1*2*pi;
w_rs = w_gps;

s = tf('s');
G_tp = w_gps/(s+w_gps);
G_hp = 1-G_tp;

G_tp_d = c2d(G_tp,Ts,'tustin');
G_hp_d = c2d(G_hp,Ts,'tustin');

bode(G_tp, G_hp)


