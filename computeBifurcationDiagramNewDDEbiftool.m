% Compute bifurcation diagram for full system
% 10.12.2024
% Kyle Wedgwood

% Initialise a few things
close all; clear; clc;
run('setDefaultParameters.m');
addpath('~/Dropbox/ddebiftool-git/ddebiftool/');
addpath('~/Dropbox/ddebiftool-git/ddebiftool_utilities',...
    '~/Dropbox/ddebiftool-git/ddebiftool_extra_nmfm',...
    '~/Dropbox/ddebiftool-git/ddebiftool_extra_psol') 
addpath(genpath('HelperFunctions'));

% Set up ddebiftool parameters
ind_tau = 3;
pars = [model_pars.p_e, model_pars.p_i, 0.0];

% Find oscillations by pushing through Hopf
funcs = set_funcs(...
  'sys_rhs', @(xx,pars) WC_rhs(0, xx, pars, model_pars), ...
  'sys_tau', ind_tau);

% Compute steady state
pars(2) = -8.0;
u0 = fsolve(@(u) WC_rhs(0, u, pars, model_pars), [0.4;0.5], ...
  optimset('Display', 'iter'));

% Run single cell continuation
pt0 = dde_stst_create('x', u0 ,'parameter', pars);

[branch,suc] = SetupStst(funcs, 'point', pt0, 'step',0.1,...
    'contpar', 1); 

branch = br_contn(funcs, branch, 30);
branch = br_rvers(branch);
branch = br_contn(funcs, branch, 100);
branch = br_stabl(funcs, branch, 0, 0);

%% Branch at hopf to find steady state
[branch,eqs_tests,ind_bifeqs,bifeqs_types] = LocateSpecialPoints(funcs, branch);
ind_hopf = ind_bifeqs(strcmp(bifeqs_types,'hopf'));

[po_branch,suc] = SetupPsol(funcs, branch, ind_hopf, 'contpar', 1, ...,
  'hopfcorrection', 1);

% Compute branch of periodics
po_branch = br_contn(funcs, po_branch, 30);
po_branch = br_stabl(funcs, po_branch, 0, 0);

%% Plot all branches and stability
fig = figure;
ax = axes(fig);
hold(ax, 'on');
xlabel(ax, 'p_e');
ylabel(ax, 'E', 'Rotation', 0);
xlim([-12,12]);
Plot2dBranch(branch, 'ax', ax,'funcs', funcs);
Plot2dBranch(po_branch, 'ax', ax);
set(ax, 'Fontsize', 20);

%% Set up coupled system
figure;

% Redefine parameters
gain = 0.2;
delay = 2.0;
pars = [gain, delay];  % THIS NEEDS TO BE A ROW VECTOR
ind_gain = 1;
ind_tau  = 2;

% Create interpolant for orbit near homoclinic
ind = length(po_branch.point)-6;
Tp = po_branch.point(ind).period;
orbit = @(t) interp1(po_branch(ind).mesh'*Tp, po_branch.point(ind), t);
mesh = Tp*po_branch.point(ind).mesh';
profile = po_branch.point(ind).profile';
model_pars.p_e = po_branch.point(ind).parameter(1);

% Set up ddebiftool
funcs = set_funcs(...
  'sys_rhs', @(xx,pars) WC_rhs_coupled(xx, pars, model_pars), ... % from below
  'sys_tau', ind_tau);

max_tries = 100;

%% Set up branch using new approach
[pfuncs,sync_branch] = SetupBranchSym('in-phase', funcs, pars, ind_tau, mesh, profile, Tp);
sync_branch = br_contn(pfuncs, sync_branch, max_tries);
sync_branch = br_rvers(sync_branch);
sync_branch = br_contn(pfuncs, sync_branch, max_tries);
sync_branch = br_stabl(pfuncs, sync_branch, 0, 0);

[sync_branch,sync_nunst,sync_bifs,sync_bifind] = MonitorChange(pfuncs, ...
  sync_branch, 'range', 2:length(sync_branch.point), ...
  'printlevel', 1, 'print_residual_info', 0, 'min_iterations', 5);

[pfuncs,anti_branch] = SetupBranchSym('anti-phase', funcs, pars, ind_tau, mesh, profile, Tp);
anti_branch = br_contn(pfuncs, anti_branch, max_tries);
anti_branch = br_rvers(anti_branch);
anti_branch = br_contn(pfuncs, anti_branch, max_tries);
anti_branch = br_stabl(pfuncs, anti_branch, 0, 0);

[anti_branch,anti_nunst,anti_bifs,anti_bifind] = MonitorChange(pfuncs, ...
  anti_branch, 'range', 2:length(anti_branch.point), ...
  'printlevel', 1, 'print_residual_info', 0, 'min_iterations', 5);

% Branch off antiphase to out-of-phase branch
[nfuncs,oop_branch,suc] = SetupBranchFromPitchfork('anti-phase', ...
  funcs, anti_branch, anti_bifind(1), ind_tau);

% Continue branch
oop_branch = br_contn(nfuncs, oop_branch, max_tries);
oop_branch = br_stabl(nfuncs, oop_branch, 0, 0);

%% Plot branches computed thus far
fig = figure;
ax = axes(fig);
hold(ax, 'on');
xlabel(ax, '\tau');
ylabel(ax, '{\Delta\phi} (x 2\pi)');
xlim([pi/2,pi]);
ylim([-0.01,0.51]);
Plot2dBranch(sync_branch, 'ax', ax, 'y', @(p) p_delta_phi_remesh(p));
Plot2dBranch(anti_branch, 'ax', ax, 'y', @(p) p_delta_phi_remesh(p));
Plot2dBranch(oop_branch, 'ax', ax, 'y', @(p) p_delta_phi_remesh(p));

set(ax, 'Fontsize', 20);

%% Continue pitchforks - from sync branch
% [nunst] = GetStability(sync_branch);
% ind = find(diff(nunst) == 1);
% ind = ind(end);
% 
% sync_branch.point = rmfield(sync_branch.point, 'd_phi');
% sync_branch.point = rmfield(sync_branch.point, 'd_peak');
% 
% %% Run continuation
% [pfuncs,sync_fold_branch,suc] = SetupPOfold(funcs, sync_branch, ind, ...
%   'contpar', [1,2], 'dir', 1);
% sync_fold_branch.parameter.max_step = [ind_gain, 0.01];
% 
% sync_fold_branch = br_contn(pfuncs, sync_fold_branch, 100);