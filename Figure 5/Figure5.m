% Generating Fig.5 from 2410.18186
% (dNeff vs. f/C for tau decays to muon and axion) 

% Import the dNeff data
dNeff_data = importdata('dNeff_fa_tau_dec.dat',',',1);
fa_x = dNeff_data.data(:,1);
dNeff_GG_y = dNeff_data.data(:,2); % dNeff from nBE
dNeff_full_y = dNeff_data.data(:,3); % dNeff from fBE


figure1=figure();
axes1 = axes('Parent',figure1);
semilogx(fa_x,dNeff_full_y,'LineWidth',2); 
hold on;
plot(fa_x,dNeff_GG_y,'LineWidth',2);
yline(0.3,'--','LineWidth',2,'Color',[0.7176 0.2745 1.0]); % Planck
yline(0.1,'--','LineWidth',2,'Color',[0.929 0.694 0.125]); % Simons obs.
yline(0.054,'--','LineWidth',2,'Color',[0.466 0.674 0.188]); % CMB-S4
grid on;
xlabel('f_a/C_{\tau\mu} [{\rm GeV}]','Interpreter','tex');
ylabel('\Delta N_{\rm eff}','Interpreter','tex')
set(axes1,'FontSize',25,'XMinorTick','on','XScale','log','YMinorTick','on',...
    'YScale','linear');
% Create textboxes
annotation(figure1,'textbox',...
    [0.704,0.817,0.35,0.0577],...
    'Color',[0.7176 0.2745 1.0],...
    'String','Planck (2σ)',...
    'FontSize',25,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(figure1,'textbox',...
    [0.1484 0.3661 0.3523 0.0577],...
    'Color',[0.929 0.694 0.125],...
    'String','Simons Obs. (2σ)',...
    'FontSize',25,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(figure1,'textbox',...
    [0.1514 0.2565 0.3485 0.0577],...
    'Color',[0.466 0.674 0.188],...
    'String','CMB-S4 (2σ)',...
    'FontSize',25,...
    'FitBoxToText','off',...
    'EdgeColor','none');
hold off;

