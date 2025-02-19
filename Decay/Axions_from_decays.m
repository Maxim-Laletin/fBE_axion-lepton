% Lepton flavour-violating decays to axions
% (on the example of tau decay to muon and axion)

%% Constants 

% Planck mass
MPl = 1e+19; % GeV
% Zeta function of 3
Zeta3 = 1.20205690315959; % (~ 1.20206)

% Mass of tau
mtau = 1.776; % GeV
% Tau dof
g_tau = 2;
% Mass of mu
mmu = 0.10565; % GeV
% Mu dof
g_mu = 2;

% Pion decay constant
fpi = 92.3e-3; % GeV
% Pion mass
mpi = 135e-3; % GeV

% Critical density of the Universe
%rho_crit = 1.878e-29*(0.673^2); % g/cm^3
rho_crit = 1.878e-29*(0.673^2)*5.6e+23/(5e+13)^3; % GeV^4
% Present-day entropy density of the Universe
s_ent0 = 2.9e-3/(5e+13)^3; % GeV^3

% ** Possible initial conditions ** 
% Electroweak phase transition temperature
T_EW = 160; % GeV

%% Axion physics

% PQ breaking scale (QCD axion)
fa = 10^8; % GeV
%fa = mpi*fpi/2/maxi; % GeV

% Number of axion degrees of freedom
g_ax = 1;
% Mass of axion 
%maxi = 0.6e-8; % GeV 
maxi = mpi*fpi/2/fa; % GeV (calculated from fa for QCD axion)
% FLV mixing parameter
A_FLV = 1.0;

% Rate of ( tau -> mu + a ) decay
Gamma = A_FLV^2*mtau^3*(1 - (mmu/mtau)^2)^3/64/pi/fa^2; % GeV

%% /////////////// The grid of x and q values ///////////////
% x = m/T is defined through the mass of the heaviest lepton 
% q = p/T

%xi_in = mtau/T_EW;
xi_in = 0.1; % [1]
xi_end = 30.0; % [1]

qi_in = 1e-2; % [1]
qi_end = 20; % [1]
 
%xi = linspace(xi_in,xi_end,1000);
xi = logspace(log10(xi_in),log10(xi_end),1000);
%qi = linspace(qi_in,qi_end,200);
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

% Hubble parameter (at RD stage)
H = @(x) mtau^2./x.^2./(MPl/sqrt(8*pi^3.*grho(x)/90)); % GeV
H_t = @(x) H(x)./(1.0 + gtilda(x)); % GeV

% FOR CONSTANT DOF TEST UNCOMMENT BELOW
% ========================================================================
% % Hubble rate with constant dof
% const_dof = 75;
% H_t = @(x) mtau^2./x.^2./(MPl/sqrt(8*pi^3.*const_dof/90)); % GeV
% % h_interpol = @(x) const_dof; 
% gs = @(x) const_dof;
% gtilda = @(x) 0.0; 
% ========================================================================

% Entropy density
s_ent = @(x) gs(x).*4*pi^2*(mtau./x).^3/90; % GeV^3

% Equilibrium number density of tau (Maxwell-Boltzmann)
nteq = @(x) g_tau*mtau^3*besselk(2,x)/2/pi^2./x; % GeV^3
% Equilibrium number density of axions (relativistic)
naeq = @(x) g_ax*Zeta3*(mtau./x).^3/pi^2; % GeV^3

% Averaged decay rate
GammaD = @(x) (2*mtau^3/2/pi^2./x)*Gamma.*(besselk(1,x)./naeq(x)); % GeV

% Equilibrium comoving density of axions (relativistic)
Yax_eq = @(x) g_ax*90*Zeta3./(gs(x)*4*pi^4); % 1

% Array of equilibrium density values
Yaxeqi = arrayfun(@(x) Yax_eq(x),xi);



%% //////// Initial value of Y and distribution function //////// 
% Null condition
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

%% Solving the ODE for the number density (assuming kinetic equilibrium)

% Relativistic correction to the density of tau leptons
qintode = @(x) integral(@(q) sqrt(q.^2-x.^2)./(exp(q)+1.0),x,inf);
% For Maxwell-Boltzmann
%qintode = @(x) besselk(1,x).*x;

% dY/dX = ...
[xode, Yode] = ode45(@(x,y) g_tau*(2*Gamma)*mtau^3*qintode(x).*(1 - y./Yax_eq(x)).*(2*pi^2*H_t(x).*s_ent(x).*x.^3).^(-1),xi,Y_0,odeset('AbsTol',1e-8,'RelTol',1e-3));
% xode - x points
% Yode - y points

%% Solving the Boltzmann equation for the distribution function 

% (Constant) M^2
Msquared = A_FLV^2*mtau^4*(1 - (mmu/mtau)^2)^2/4/fa^2; % GeV^2

% Parameters to pass to the PDE solver
par = [(2*g_tau*Msquared/mtau/16/pi/g_ax); mtau^2*sqrt(8*pi^3/90)/MPl; mmu/mtau]; 
% [factor in front of A(x,q) and B(x,q); factor in H(x); ratio in A(x,q) and B (x,q)]

% Packing the parameters into a function that is passed to the PDE solver
boltfun = @(q,x,u,dudq) boltfunC(q,x,u,dudq,par,gtilda,H_t);

% /////////// Solving the PDE ////////////
sol=pdepe(0,boltfun,distr_ini,@boltbc,qi,xi);
% distr_ini - initial distribution 
% boltbc - boundary conditions

%% Results

% Obtaining the comoving density from the distrubiton
sol_Y = g_ax*(1/2/pi^2)*trapz(qi,sol,2).*(mtau./xi').^3./arrayfun(s_ent,xi)';

% Obtaining the temperature of the axion gas (ratio to the SM temperature)
T_ratio = (30*Zeta3/pi^4)*trapz(qi,qi.*sol,2)./trapz(qi,sol,2);

% Plot the result and compare it to the ODE solution for the number density
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


%% Creating a movie with distribution functions

% fmov = figure;
% plot(qi,(mtau/xi(end)).^3.*sol(end,:)/s_ent(xi(end)),'LineWidth',2);
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
% plot(qi,(mtau/xi(end)).^3.*( Yode(end)./Yaxeqi(end) ).*qi.^2.*f_ax_eq_array,'LineWidth',2);
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
% video_name = ['DistrEvolutionTauDecay' num2str(fa,2) '.avi'];
% v = VideoWriter(video_name,'Motion JPEG AVI');
% v.Quality = 100;
% v.FrameRate = 10;
% open(v);
% 
% video_step = 5;
% frameIdx = 1;
% for j = 1:video_step:length(xi)
% 
%     plot(qi,(mtau/xi(j)).^3.*sol(j,:)/s_ent(xi(j)),'LineWidth',2);
%     ax.NextPlot = 'add'; % add the next plot to the axes
%     plot(qi,(mtau/xi(j)).^3.*(Yode(j)./Yaxeqi(j)).*qi.^2.*f_ax_eq_array/s_ent(xi(j)),'LineWidth',2);
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

% Full dNeff in late Universe, but when axions are still relativistic
dNeff_full_lateUni = dNeff_general(100E-9,43/11);
% at T = 100 eV, so that the axions are definetely relativistic

%% ****** Functions that are used for the solution of ODEs and PDEs *******

% Describes the PDE (with a possibility to pass external params)
function [c,f,s] = boltfunC(q,x,u,dudq,params,gtild,ht)
    c = 1; % factor in front of df/dx
    f = gtild(x).*u.*q./x; % factor in front of
    % ==== Source function ====
    s = (B(x,q,params,ht) - A(x,q,params,ht)).*(faeq(q) - u) - 3*gtild(x).*u./x; 

end

% Boundary conditions on the dist. function
function [pl,ql,pr,qr] = boltbc(xl,ul,xr,ur,x)
    pl = ul;
    ql = 0;
    pr = ur;
    qr = 0;
end

% Zero initial distribution functions
function fx0 = boltic(q)
    fx0 = 0;
end

% Equilibrium distribution function
function equil_func = faeq(q)
    equil_func = q.^2./(exp(q)-1);
    %equil_func = q.^2./exp(q); % Maxwell-Boltzmann
end

function A_func = A(x,q,pars,ht)
    A_func = pars(1)*(ht(x)*q.^2).^(-1).*log(1+exp(-max([x (x*pars(3)+q) x.^2.*(1 + 4*q.^2./(x.^2.*(1-pars(3)^2)) - pars(3)^2)/4./q])));
    % For Maxwell-Boltzmann dist
    %A_func = pars(1)*(ht(x).*q.^2).^(-1).*exp(-max([x (x*pars(3)+q) x.^2.*(1 + 4*q.^2./(x.^2.*(1-pars(3)^2)) - pars(3)^2)/4./q])); % takes into account condition |costheta| <= 1
   

end

function B_func = B(x,q,pars,ht)
    B_func = pars(1)*(ht(x)*q.^2).^(-1).*log(1+exp(-max([(x-q) x*pars(3) x.^2/4./q])));
    % For Maxwell-Boltzmann dist for muons
    %B_func = pars(1)*(ht(x)*q.^2).^(-1).*exp(-max([(x-q) x*pars(3) x.^2/4./q]));
    
end
