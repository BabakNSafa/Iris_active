x = 0:.1:6;
a_s = 1;
r_s = 3.4;
scale = 5;
figure
plot(x,(1-tanh((x-(r_s+a_s))*scale)).*(x>r_s))