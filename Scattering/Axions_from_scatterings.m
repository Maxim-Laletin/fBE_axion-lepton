% Axions produced via diagonal interactions with leptons
% (on the example of muon scattering processes)

%% Constants

% Planck mass
MPl = 1e+19; % GeV
% Zeta function of 3
Zeta3 = 1.20205690315959; % (~ 1.20206)
% Fine structure constant
a_el = 1/137; % 1
% Electric charge
e_g = sqrt(a_el*4*pi); % 1


% Mass of mu
mmu = 0.10565; % GeV
% Mu dof
g_mu = 2;
% Photon degrees of freedom
g_ph = 2;
% Pion decay constant
fpi = 92.3e-3; % GeV
% Pion mass
mpi = 135e-3; % GeV

% Diagonal coupling of muons to axions
C_mu = 1.0;

% Critical density of the Universe
%rho_crit = 1.878e-29*(0.673^2); % g/cm^3
rho_crit = 1.878e-29*(0.673^2)*5.6e+23/(5e+13)^3; % GeV^4
% Present-day entropy density of the Universe
s_ent0 = 2.9e-3/(5e+13)^3; % GeV^3

% Electroweak phase transition temperature
T_EW = 160; % GeV

% PQ breaking scale (QCD axion)
%fa = mpi*fpi/2/maxi; % GeV

fa = 10^6; % GeV

% Number of axion degrees of freedom
g_ax = 1;
% Mass of axion
%maxi = 0.6e-8; % GeV 
maxi = mpi*fpi/2/fa; % GeV
% FLV mixing parameter
A_FLV = 1.0;

%% /////////////// Defining a grid of x and q values ///////////////
% x is defined through the mass of the heaviest lepton 

%xi_in = mmu/T_EW;
xi_in = 0.01; % [1]
xi_end = 20.0; % [1]

qi_in = 1e-2; % [1]
qi_end = 20; % [1]

xi = logspace(log10(xi_in),log10(xi_end),1000);
qi = logspace(log10(qi_in),log10(qi_end),200);

%/////////////////////////////////////////////////////////////////////////
%% Effective number of density and entropy degrees of freedom

% Imported interpolation table from the corresponding file
gdof_array = importdata('Rel_dof_from_1606_07494.csv',' ',2);

gdof_x = mtau./(10.^(gdof_array.data(:,1))*1e-3); % x [1]
gdof_y = gdof_array.data(:,2); % g_rho
gdof_z = gdof_y./gdof_array.data(:,3); % g_s

% Constructing splines
gs_sp = spline(gdof_x,gdof_z); % g_s(x) as a cubic spline
grho_sp = spline(gdof_x,gdof_y); %g_rho(x) as a cubic spline

% Entropy degrees of freedom as a function of x
gs = @(x) ppval(gs_sp,x);
% Energy density degrees of freedom as a function of x 
grho = @(x) ppval(grho_sp,x);

% ** Calculate the derivative of gs **
% Coefficients and break-points of the spline
[breaks,coefs,l,k,d] = unmkpp(gs_sp);

% Create a piece-wise polynomial for the derivative
pp2 = mkpp(breaks,(-1/3)*repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
% g_tilda as a function of x
gtilda = @(x) (x./gs(x)).*ppval(pp2,x);

%% Cosmology

% Hubble parameter (at RD)
H = @(x) mmu^2./x.^2./(MPl/sqrt(8*pi^3.*grho(x)/90)); % GeV
H_t = @(x) H(x)./(1.0 + gtilda(x)); % GeV

% FOR CONSTANT DOF UNCOMMENT THIS BELOW
% ========================================================================
% % Hubble rate with constant dof
% const_dof = 20;
% H_t = @(x) mmu^2./x.^2./(MPl/sqrt(8*pi^3.*const_dof/90)); % GeV
% gs = @(x) const_dof;
% gtilda = @(x) 0.0; 
% ========================================================================

% Entropy density
s_ent = @(x) gs(x).*4*pi^2*(mmu./x).^3/90; % GeV^3

% Equilibrium number density of muons
% relativistic
nmeq = @(x) g_mu*mmu^3*integral(@(t) sqrt(t.^2 - x.^2).*t./(exp(t)+1.0),x,Inf)/2/pi^2./x.^3; % GeV^3
% Maxwell-Boltzmann
nmeq_MB = @(x) g_mu*mmu^3*besselk(2,x)/2/pi^2./x; % GeV^3

% Equilibrium number density of photons
% relativistic
npheq = @(x) g_ph*Zeta3*(mmu./x).^3/pi^2; % GeV^3
% Maxwell-Boltzmann
npheq_MB = @(x) g_ph*(mmu./x).^3/pi^2; % GeV^3

% Equilibrium number density of axions (relativistic)
naeq = @(x) g_ax*Zeta3*(mmu./x).^3/pi^2; % GeV^3


% Averaged decay rate
GammaD = @(x) (2*mmu^3/2/pi^2./x)*Gamma.*(besselk(1,x)./naeq(x)); % GeV

% Equilibrium comoving density of axions (relativistic)
Yax_eq = @(x) g_ax*90*Zeta3./(gs(x)*4*pi^4); % 1

% Arrays of equilibrium density values
Yaxeqi = arrayfun(@(x) Yax_eq(x),xi);

%% //////// Initial value of Y and distribution function //////// 
%  Null condition
Y_0 = 0.0;
distr_ini = @boltic; % (see the functions at the end of the file)

% Equilibrium condition
%Y_0 = Yaxeqi(1);
%distr_ini = @faeq;
%/////////////////////////////////////////////////////////////////////////
%% Axion distribution

% Equilibrium distribution of axions
f_ax_eq = @(q) (exp(q)-1.0).^(-1);

% Array of equilibrium distribution values
f_ax_eq_array = arrayfun(f_ax_eq,qi);

%% Velocity-averaged cross sections for processes
% (assuming MB distributions)

% Cross-section for mu mu -> a gamma
sigma_ann = @(s) (C_mu*e_g*mmu/fa)^2*atanh(sqrt(1-4*mmu^2./s))./(s-4*mmu^2)/4/pi; % GeV^-2
% Cross-section for Primakoff (mu gamma -> a mu)
sigma_prim = @(s) (C_mu*e_g*mmu/fa)^2*(2*s.^2.*log(s/mmu^2) - 3*s.^2 + 4*mmu^2*s - mmu^4)./s.^2./(s-mmu^2)/32/pi; % GeV^-2

% Auxiliary function lambda
lam_f = @(x,y,z) (x - (y+z).^2).*(x - (y-z).^2); % GeV^4

% Velocity-averaged cross-section for annihilation 
sigmaV_ann = @(x) x.*integral(@(s) lam_f(s,mmu,mmu).*sigma_ann(s).*besselk(1,sqrt(s)*x/mmu)./sqrt(s),4*mmu^2,4*20*mmu^2.*(1 + 1./x.^2) )./besselk(2,x).^2/8/mmu^5; % GeV^-2

% Velocity-averaged cross-section for Primakoff
sigmaV_prim = @(x) x.^3.*integral(@(s) lam_f(s,mmu,0).*sigma_prim(s).*besselk(1,sqrt(s)*x/mmu)./sqrt(s),mmu^2,20*mmu^2.*(1 + 1./x.^2) ,'AbsTol',1e-35)./besselk(2,x)/16/mmu^5; % GeV^-2

%% Solving the ODE for the number density (assuming kinetic equilibrium)

%qintode = @(x) integral(@(q) sqrt(q.^2-x.^2)./(exp(q)+1.0),x,inf);

%Rate of annihilation
Rate_ann = @(x) nmeq(x).^2.*sigmaV_ann(x); % GeV^4
Rate_ann_MB = @(x) nmeq_MB(x).^2.*sigmaV_ann(x); % GeV^4 (Maxwell-Boltzmann)
%Rate of Primakoff
Rate_prim = @(x) 2*nmeq(x).*npheq(x).*sigmaV_prim(x); % GeV^4
Rate_prim_MB = @(x) 2*nmeq_MB(x).*npheq(x).*sigmaV_prim(x); % GeV^4 (Maxwell-Boltzmann)
% factor of 2 takes into account particles and antiparticles

% CosmoRatePrim = @(x) Rate_prim(x)./(x.*H_t(x).*s_ent(x));
% CosmoRateAnn = @(x) Rate_ann(x)./(x.*H_t(x).*s_ent(x));
% 
% CosmoRatePrim_MB = @(x) Rate_prim_MB(x)./(x.*H_t(x).*s_ent(x));
% CosmoRateAnn_MB = @(x) Rate_ann_MB(x)./(x.*H_t(x).*s_ent(x));
% 
% CosmoRateFull = @(x) CosmoRatePrim(x) + CosmoRateAnn(x);

% dY/dX = ...
[xode, Yode] = ode45(@(x,y) (Rate_ann(x) + Rate_prim(x)).*(1 - y./Yax_eq(x)).*(x.*H_t(x).*s_ent(x)).^(-1),xi,Y_0,odeset('AbsTol',1e-6));
% xode - x points
% Yode - y points

%% Full integration for collision terms for i+j -> k+a processes

% Import the values of x and q for the interpolation from the table 
x_cc = importdata('x.dat',',');
q_cc = importdata('q.dat',',');

% Obtain the values from files provided by CollCalc (table)
CTAnn = importdata('CollTermAnn.dat',' ');
CTPrim = importdata('CollTermPrim.dat',' ');

%% Solving the Boltzmann equation for the distribution function 

% Parameters to pass to the PDE solver
par = [g_mu^2*(e_g*C_mu/fa)^2; mmu];

% Interpolating the logs of the corresponding collision term functions
log_interp_ann = @(x,q) interp2(x_cc,q_cc,log(CTAnn'),x,q,'linear');
log_interp_prim = @(x,q) interp2(x_cc,q_cc,log(CTPrim'),x,q,'linear');

% Collision terms based on interpolation
CollTermAnn_interp = @(x,q) (0.5./q/g_ax).*(x./mmu).*(par(1)*par(2)^2)*exp(log_interp_ann(x,q));
CollTermPrim_interp = @(x,q) (0.5./q/g_ax).*(x./mmu).*(2*par(1)*par(2)^2*(g_ph/g_mu))*exp(log_interp_prim(x,q));
% the factor of 2 difference with annihilation is because Primakoff
% scattering concerns both muons and antimuons

% Simplified collision terms (assuming Maxwell-Boltzmann distributions)
CollTermAnn_simp = @(x,q) CollTermAnn_simp_func(x,q,(0.5./g_ax).*(par(1)*par(2)^2)*mmu);
CollTermPrim_simp =  @(x,q) CollTermPrim_simp_func(x,q,(0.5./g_ax).*(2*par(1)*par(2)^2)*mmu);
% (see functions at the end of the file)

% Total
CollTermFull_interp = @(x,q) CollTermPrim_interp(x,q) + CollTermAnn_interp(x,q);
CollTermFull_simp = @(x,q) CollTermAnn_simp(x,q) + CollTermPrim_simp(x,q);

%% Integrated rates

% % Integrated rates (over the axion momentum space)
% x_probe = logspace(-3,1,100);
% IntegratedRatePrim_interp = arrayfun(@(x) (1/2/pi^2)*(mmu./x).^3.*integral(@(q) q.^2.*CollTermPrim_interp(x,q),0.0001,25,'ArrayValued',true),x_probe);
% IntegratedRateAnn_interp = arrayfun(@(x) (1/2/pi^2)*(mmu./x).^3.*integral(@(q) q.^2.*CollTermAnn_interp(x,q),0.0001,25,'ArrayValued',true),x_probe);
% IntegratedRatePrim_simp = arrayfun(@(x) (1/2/pi^2)*(mmu./x).^3.*integral(@(q) q.^2.*CollTermPrim_simp(x,q),0.0001,25,'ArrayValued',true),x_probe);
% IntegratedRateAnn_simp = arrayfun(@(x) (1/2/pi^2)*(mmu./x).^3.*integral(@(q) q.^2.*CollTermAnn_simp(x,q),0.0001,25,'ArrayValued',true),x_probe);
% Rate_prim_ar = arrayfun(@(x) Rate_prim(x),x_probe);
% Rate_ann_ar = arrayfun(@(x) Rate_ann(x),x_probe);
% 
% % Partial rates
% figure01 = figure();
% axes01 = axes('Parent',figure01);
% loglog(x_probe,fa^2*IntegratedRateAnn_interp,'LineWidth',2,'Color',[0.392156862745098 0.831372549019608 0.0745098039215686]);
% hold on;
% loglog(x_probe,fa^2*IntegratedRateAnn_simp,'--','LineWidth',2,'Color',[0 0 1]);
% loglog(x_probe,fa^2*IntegratedRatePrim_interp,'LineWidth',2);
% loglog(x_probe,fa^2*IntegratedRatePrim_simp,'--','LineWidth',2);
% grid on;
% set(axes01,'FontSize',25);
% xlim([0.01 10]);
% ylim([1e-22 1e-1]);
% xlabel('x');
% ylabel('f_a^2 \cdot Rate [GeV^6]');
% legend(' Annihilation',' Annihilation (simplified)',' Primakoff',' Primakoff (simplified)');
% hold off;
% % Step 3: Create the inset plot (smaller axes on top)
% inset_ax = axes('Position', [0.18806161745828 0.213664596273292 0.2997432605905 0.327950310559006]);  % [left, bottom, width, height]
% plot(x_probe,fa^2*IntegratedRateAnn_interp,'LineWidth',2,'Color',[0.392156862745098 0.831372549019608 0.0745098039215686]);
% hold on;
% plot(x_probe,fa^2*IntegratedRateAnn_simp,'--','LineWidth',2,'Color',[0 0 1]);
% plot(x_probe,fa^2*IntegratedRatePrim_interp,'LineWidth',2);
% plot(x_probe,fa^2*IntegratedRatePrim_simp,'--','LineWidth',2);
% xlim([0.01 0.02]);
% set(inset_ax, 'Box', 'on','FontSize',14,'XTick',[0.01 0.015 0.02],'YTick',...
%     [0.02 0.04]);  % Add a box around the inset plot
% % Step 4: Optionally, adjust the inset plot aesthetics
% %set(inset_ax, 'XTick', [], 'YTick', []);  % Remove ticks for clarity
% hold off;
%
% % Total rates
% figure02 = figure();
% axes02 = axes('Parent',figure02);
% %loglog(x_probe,IntegratedRateAnn_simp + IntegratedRatePrim_simp,'LineWidth',2); 
% loglog(x_probe,IntegratedRateAnn_interp + IntegratedRatePrim_interp,'LineWidth',2);
% hold on;
% plot(x_probe,Rate_ann_ar + Rate_prim_ar,'r--','LineWidth',2);
% grid on;
% set(axes02,'FontSize',25);
% xlabel('x');
% ylabel('Rate');
% title('Integrated total rates');
% legend(axes02,'From collision term','nBE rate','Location','southwest');
% hold off;

%% Solving the PDE

boltfun = @(q,x,u,dudq) boltfunC(q,x,u,dudq,par,gtilda,H_t,CollTermFull_interp);
boltfun_simp = @(q,x,u,dudq) boltfunC(q,x,u,dudq,par,gtilda,H_t,CollTermFull_simp);

% /////////// Solving the PDE ////////////
sol = pdepe(0,boltfun_simp,distr_ini,@boltbc,qi,xi,odeset('RelTol',1e-5));
%sol = pdepe(0,boltfun,distr_ini,@boltbc,qi,xi,odeset('RelTol',1e-5));

% Obtaining the comoving density from the distrubiton
sol_Y = g_ax*(1/2/pi^2)*(trapz(qi,sol,2)').*(mmu./xi).^3./arrayfun(s_ent,xi);

% Obtaining the temperature of the axion gas (ratio to the SM temperature)
T_ratio = (30*Zeta3/pi^4)*trapz(qi,qi.*sol,2)./trapz(qi,sol,2);

%% Plot the result and compare it to the ODE solution for the number density
% and to the equilibrium density of axions

figure1 = figure; 
axes1 = axes('Parent',figure1);
hold(axes1,'on');
loglog(xi,sol_Y,'LineWidth',2); % fBE solution
plot(xode,Yode,'LineWidth',2); % nBE solution
plot(xi,Yaxeqi,'LineWidth',2); % Equilibrium
grid(axes1,'on');
xlim([xi(1) xi(end)]);
ylabel('Y_a');
xlabel('x');
title(['f_a = ' num2str(fa,2) ' GeV']);
set(axes1,'FontSize',25,'XMinorTick','on','XScale','log','YMinorTick','on',...
    'YScale','log');
legend(axes1,'fBE','nBE','Equilibrium','Location','southeast');
hold off;


% Compare distributions at a given point in x (xis)
xis = 200;
figure2 = figure();
axes2 = axes('Parent',figure2);
hold(axes2,'on');
plot(qi,sol(xis,:), 'LineWidth',2);
plot(qi,(Yode(xis)./Yaxeqi(xis)).*qi.^2.*f_ax_eq_array, 'LineWidth',2); % normalized to the solution of nBE
grid on;
xlim([qi(1) qi(end)]);
ylabel({'q^2 \cdot f_a(q)',''});
xlabel('q');
title(['f_a = ' num2str(fa,2) ' GeV, x = ' num2str(xi(xis)) ]);
legend(axes2,'fBE','nBE','Location','northeast');
set(axes2,'FontSize',25);
hold off;

% Ratio of the axion energy density to the SM plasma density (radiation)
rho_ratio = trapz(qi,qi.*sol,2)/(pi^4/15)./grho(xi)';

%% Creating a movie with distribution functions
% 
% fmov = figure;
% plot(qi,(mmu/xi(end)).^3.*sol(end,:)/s_ent(xi(end)),'LineWidth',2);
% ax = gca;
% ax.YLabel.String = 'q^2 \cdot f(q)  T^3/ s';
% ax.XLabel.String = 'q';
% title('fa = ',num2str(fa,'%.1e') );
% ax.XLim = [0 15];
% hold(ax,'on');
% grid(ax,'on');
% dim = [.25 .5 .3 .3];
% annot = annotation('textbox',dim,'String',['x = ',num2str(xi(1)), ', Y_{nbe} = ',num2str(Yode(1)),', Y_{fbe} = ',num2str(sol_Y(1))],'FitBoxToText','on');
% % Plot of the equilibrium distr. with the same NORMALIZATION as the actual
% % solution
% 
% plot(qi,(mmu/xi(end)).^3.*( Yode(end)./Yaxeqi(end) ).*qi.^2.*f_ax_eq_array,'LineWidth',2);
% hold(ax,'off');
% delete(annot);
% 
% axis tight manual % squeeze the y-axis to fit the curves
% ax.NextPlot = 'replaceChildren'; % replace curves in the axes
% 
% loops = 200;
% M(loops) = struct('cdata',[],'colormap',[]);
% 
% fmov.Visible = 'on';
% 
% % Write a movie to export
% video_name = ['DistrEvolutionMuScattering' num2str(fa,2) '.avi'];
% v = VideoWriter(video_name,'Motion JPEG AVI');
% v.Quality = 100;
% v.FrameRate = 10;
% open(v);
% 
% video_step = 5;
% frameIdx = 1;
% for j = 1:video_step:length(xi)
% 
%     plot(qi,(mmu/xi(j)).^3.*sol(j,:)/s_ent(xi(j)),'LineWidth',2);
%     ax.NextPlot = 'add'; % add the next plot to the axes
%     plot(qi,(mmu/xi(j)).^3.*(Yode(j)./Yaxeqi(j)).*qi.^2.*f_ax_eq_array/s_ent(xi(j)),'LineWidth',2);
%     str = ['x = ',num2str(xi(j))];
%     title('fa = ',num2str(fa,'%.1e') );
%     annot = annotation('textbox',dim,'String',['x = ',num2str(xi(j)), ', Y_{nbe} = ',num2str(Yode(j),'%.2e'),', Y_{fbe} = ',num2str(sol_Y(j),'%.2e')],'FitBoxToText','on');
%     drawnow;
%     M(frameIdx) = getframe(gcf);
%     writeVideo(v,getframe(gcf));
% 
%        frameIdx = frameIdx + 1;
% 
%     %ax.NextPlot = 'replaceChildren'; %clear the axes for the next curve(s)
%     delete(annot);
%     cla(ax);
% 
% end
% 
% close(v);

%% Calculation of dNeff (Delta N Effective)

% Simplified formula for dNeff
dNeff_GG = 74.85*Yode(end)^(4/3);

% Corrected by the A^(1/3) factor
dNeff_simplified_corrected = 74.85*Yode(xis)^(4/3)/(Yode(xis)/Yaxeqi(xis))^(1/3);

% Simplified formula that takes into account the axion contribution to the
% radiation bath
dNeff_with_axion_impact = (4/7)*((11/4)*(2*pi^4*(43/11)*Yode(end)/45/Zeta3)/(1 - 2*pi^4*Yode(end)/45/Zeta3))^(4/3);

% Omega of axions
Omega_ax = (pi^2/30)*(pi^2*Yode(end)*s_ent0/Zeta3)^(4/3)/rho_crit; % 1

% Full formula to compute at later times in x
axion_distrib = sol(end,:)./qi.^2; % axion dist. function at the end of evolution
rho_gamma = @(T) 2*pi^2*T.^4/30; % photon energy density [GeV^4]
dNeff_general = @(T,gT)  (8/7)*(11/4)^(4/3)*(g_ax*(1/2/pi^2)*(gs(xi_end)./gT).^(-1).*trapz(sqrt(qi.^2 + (maxi./T).^2),qi.*sqrt(qi.^2 + (maxi./T).^2).*sqrt((gT/gs(xi_end)).^(2/3).*qi.^2 + (maxi./T).^2).*axion_distrib,2).*T.^4)./rho_gamma(T);
% gT is the number of entropy degrees of freedom at T

dNeff_full_lateUni = dNeff_general(100E-9,43/11);
% at T = 100 eV, so that the axions are definetely relativistic

%% ****** Functions that are used for the solution of ODEs and PDEs *******

% Describes the PDE (with a possibility to pass external params)
function [c,f,s] = boltfunC(q,x,u,dudq,params,gtild,ht,ctann)
    c = 1; % factor in front of df/dx
    f = gtild(x).*u.*q./x; % factor in front of 
    % ==== Source function ====
     s = (q.^2./ht(x)./x).*ctann(x,q).*(1 - u./faeq(q)) - 3*gtild(x).*u./x; 
end

% Another possible realization
% function [c,f,s] = boltfunC(q,x,u,dudq,params,gtild,ht,ctann)
%     c = x./q.^3; % factor in front of df/dx
%     f = gtild(x).*u./q.^2; % factor in front of 
%     % ==== Source function ====
%      %s = q.*ctann(x,q).*(1 - u./faeq(q)); 
% end


% Zero initial distribution function
function fx0 = boltic(q)
    fx0 = 0;
end

% Boundary conditions on the dist. function
function [pl,ql,pr,qr] = boltbc(xl,ul,xr,ur,x)
    pl = ul;
    ql = 0;
    pr = ur;
    qr = 0;
end

% Equilibrium distribution function
function equil_func = faeq(q)
    equil_func = q.^2./(exp(q)-1);
end


function y = CollTermAnn_simp_func(x,q,norm)

    if x < 0.001 
        ek_max = 20.0;
    else
        ek_max = Inf;
    end

    y = norm*(q.*x)^(-1).*(exp(-q)./q/2/(2*pi)^3 ).*integral(@(ek) ((2*ek.*q - x.^2).*atanh(sqrt(1 - x.^2./ek./q)) - ek.*q.*sqrt(1 - x.^2./ek./q))./(exp(ek)),x.^2./q,ek_max);
end

function y = CollTermPrim_simp_func(x,q,norm)

    Integral_s = @(s,x) 2*s.*log(s./x.^2) + 4*x.^2.*log(s) - 5*s + x.^4./s;
       
    if x < 0.001
        ek_max = 20.0;
    else
        ek_max = Inf;
    end

    y = norm*(q.*x)^(-1).*(exp(-q)./q/32/(2*pi)^3 ).*integral(@(ek) (Integral_s(x.^2 + 2*ek.*q + 2*sqrt(ek.^2 - x.^2).*q,x) - Integral_s(x.^2 + 2*ek.*q - 2*sqrt(ek.^2 - x.^2).*q,x))./(exp(ek)+1.0),x,ek_max);   
end


