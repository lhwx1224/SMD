% SAVESYSTEMINFO_PPB.m saves the system information obtained from the
% discretizer for a pin-pined beam. The saved system can be used for further numerical
% simulation using various stepping schemes.
path = uigetdir(cd);
% ------------------------------------------------------------------------
%                      General System Setup
% ------------------------------------------------------------------------
System.System = 'Non-uniform Pin-Pined Beam';
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
%                 Beam Mode Shapes Setup   
% ------------------------------------------------------------------------
System.DiscretizationMethod = 'AMM';
System.SpatialResolution = dx;
System.Dimensions = L0;
System.UniformDensity = m0;
System.Density = mx;
System.UniformRigidity = EI0;
System.FlexuralRigidity = EIx;
System.Force = Fx;
% ------------------------------------------------------------------------
%                    Stepping Specification
% ------------------------------------------------------------------------
System.fs = fs;
System.dt = dt;
System.rr = rr;
System.tend = tend;
% ------------------------------------------------------------------------
% Discretization Setup
% ------------------------------------------------------------------------
System.nm = nm;
% Discretization Results
System.M = M_am;
System.K = K_am;
System.T = T;
System.Lambda = Lambda;
System.phi_am = phi_am;
System.phi_tilde = phi_tilde;
System.weight = weight;
System.omegan = omega_am;

save([path,'\System.mat'],'System')

disp('SYSTEM INFORMATION SAVED!')