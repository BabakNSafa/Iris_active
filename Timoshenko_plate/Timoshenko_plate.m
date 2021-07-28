function Timoshenko_plate()
syms A B C s_r s_t u r E v a b T_s
s_r     = A/(r^2)     + B*(1+2*log(r)) + 2*C;
s_t     = -A/(r^2)    + B*(3+2*log(r)) + 2*C;

%bc1; sr (@r=a) = 0
s_r_bc1 = (subs(s_r,r,a) == 0);

%bc2; st (@r=a) = T_s
s_t_bc2 = (subs(s_t,r,a) == T_s);

%bc3; u (@r=b) = 0
u       = r/E * (s_t - v * s_r);
u_bc3   = (subs(u,r,b) == 0);

y   = solve([s_r_bc1, s_t_bc2, u_bc3],[A, B, C]);

u   = subs(u,   [A, B, C],  [y.A, y.B, y.C]);
s_r = subs(s_r, [A, B, C],  [y.A, y.B, y.C]);
s_t = subs(s_t, [A, B, C],  [y.A, y.B, y.C]);

figure('Units','pixels','WindowStyle','normal','Position',[100,200,1500,600]);

nu = 0:.1:.4;
for i=1:length(nu)
    subplot(2,3,1)
    u_plot = subs(u,[T_s, a, b, E, v],[-.1, 3.4, 6, 1, nu(i)]);
    hold on
    fplot(u_plot,[3.4,6])
    ylabel('$u (r)$','interpreter','latex')
    xlabel('$r$','interpreter','latex')

    subplot(2,3,2)
    fplot(diff(u_plot),[3.4,6])
    ylabel('$\epsilon_r (r)$','interpreter','latex')
    xlabel('$r$','interpreter','latex')
    
    subplot(2,3,3)
    fplot(u_plot/r,[3.4,6])
    ylabel('$\epsilon_{\theta} (r)$','interpreter','latex')
    xlabel('$r$','interpreter','latex')
    
    subplot(2,3,4)
    s_r_plot = subs(s_r,[T_s, a, b, E, v],[-.1, 3.4, 6, 1, nu(i)]);
    hold on
    fplot(s_r_plot,[3.4,6],'Color','k')
    ylabel('$\sigma_r (r)$','interpreter','latex')
    xlabel('$r$','interpreter','latex')
    
    subplot(2,3,5)
    s_t_plot = subs(s_t,[T_s, a, b, E, v],[-.1, 3.4, 6, 1, nu(i)]);
    hold on
    fplot(s_t_plot,[3.4,6],'Color','k')
    ylabel('$\sigma_{\theta} (r)$','interpreter','latex')  
    xlabel('$r$','interpreter','latex')
    
    subplot(2,3,6)
    hold on
    fplot(0,'Color','k')
    ylabel('$\tau_{r\theta} (r)$','interpreter','latex')  
    xlabel('$r$','interpreter','latex')
end
for i=1:6
    subplot(2,3,i)
    hold on
    set(gca,'fontsize',18,'xlim',[3.4,6])
end

end