% SAVESYSTEMINFO_FFM.m saves the system information obtained from the
% discretizer for a fixed-fixed membrane. The saved system can be used for further numerical
% simulation using various stepping schemes.
path = uigetdir(cd);
% ------------------------------------------------------------------------
%                      General System Setup
% ------------------------------------------------------------------------
System.System = 'Non-uniform Fixed-fixed Membrane';
System.Date = string(datetime);
System.SimulationMethod = Method;
System.SimulationSpace = SimSpace;
System.DampingType = DampingType;
% ------------------------------------------------------------------------
%                Forcing Setup for the Simulation
% ------------------------------------------------------------------------
System.ForcingType = Forcing;
System.ExciteForm = ExciteForm;
System.ForceWeight = ForceWeight;
System.ExciteLocation = ExciteLoc;
% ------------------------------------------------------------------------
%                 Membrane Mode Shapes Setup   
% ------------------------------------------------------------------------
System.DiscretizationMethod = 'AMM';     % Discretization Method Used
System.SpatialResolution.dx = dx;        % Spatial Resolution in x 
System.SpatialResolution.dy = dy;        % Spatial Resolution in y
System.Dimensions.a0 = a0;               % Membrane Dimension in x
System.Dimensions.b0 = b0;               % Membrane Dimension in y
System.UniformDensity = rho0;            % Uniform Material Density
System.Density = rho_xy_str;             % Nonuniform Density (definition)
System.UniformTension = T0;              % Uniform Tension
System.Tension = T_xy_str;               % Nonuniform Tension (definition)
System.Force = F_xy_str;                 % Forcing (definition)
% ------------------------------------------------------------------------
%                    Stepping Specification
% ------------------------------------------------------------------------
System.fs = fs;                          % Sampling Frequency (Hz)
System.dt = dt;                          % Sampling Time
System.rr = rr;                          % Resampling Rate
System.tend = tend;                      % Ending Time for Simulation
% ------------------------------------------------------------------------
% Discretization Setup
% ------------------------------------------------------------------------
System.nm = nm;                          % Total Number of Modes
System.mm = mm;                          % Number of Modes in x direction
System.nn = nn;                          % Number of Modes in y direction

% Discretization Results
System.M = M_am;                         % Mass Matrix
System.K = K_am;                         % Stiffness Matrix
System.T = T;                            % Configuration Space Modal Matrix
System.Lambda = Lambda;                  % Configuration Space Eigenvalues
System.phi_am = phi_am;                  % Approximated mode shapes
System.phi_tilde = phi_tilde;            % Assumed Modes
System.weight = weight;                  % Forcing Weigth Vector
System.omegan = omega_am;                % Natural Angular Frequency

save([path,'\System.mat'],'System')      % Save the System Matrix

disp('SYSTEM INFORMATION SAVED!')