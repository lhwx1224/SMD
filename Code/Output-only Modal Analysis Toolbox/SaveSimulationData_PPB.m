% ------------------------------------------------------------------------
%                      General System Setup
% ------------------------------------------------------------------------
Data.System = 'Non-uniform Pin-Pined Beam';
Data.Date = string(datetime);
Data.SimulationMethod = Method;
Data.SimulationSpace = SimSpace;
Data.DampingType = DampingType;
% ------------------------------------------------------------------------
%                Forcing Setup for the Simulation
% ------------------------------------------------------------------------
Data.ForcingType = Forcing;
Data.ExciteForm = ExciteForm;
Data.ForceWeight = ForceWeight;
Data.ExciteLocation = ExciteLoc;
% ------------------------------------------------------------------------
%                 Beam Mode Shapes Setup   
% ------------------------------------------------------------------------
Data.DiscretizationMethod = 'AMM';
Data.SpatialResolution = dx;
Data.Dimensions = L0;
Data.UniformDensity = m0;
Data.Density = mx;
Data.UniformRigidity = EI0;
Data.FlexuralRigidity = EIx;
Data.Force = Fx;
% ------------------------------------------------------------------------
%                    Stepping Specification
% ------------------------------------------------------------------------
Data.fs = fs;
Data.dt = dt;
Data.rr = rr;
Data.tend = tend;
Data.Forcing = forcing;
% Initial Conditions
Data.IW = IW;
Data.ModeIndices = ModeIndices;
% ------------------------------------------------------------------------
% Discretization Setup
% ------------------------------------------------------------------------
Data.nm = nm;
% Discretization Results
Data.M = M_am;
Data.K = K_am;
Data.T = T;
Data.Lambda = Lambda;
Data.phi_am = phi_am;
Data.phi_tilde = phi_tilde;
Data.C = C_am;
Data.weight = weight;
Data.omegan = omega_am;
Data.omegad = omegad_am;
% Mode Shapes
Data.Phi_tilde = Phi_tilde;
Data.Phi_x = Phi_x;
Data.Phi_v = Phi_v;
% Response Data
Data.Eta = eta_am;
Data.Y = y_modal;
% Response Size
Data.ResponseMatSize = size(y_modal);
save([path,'Data.mat'],'Data','-v7.3')

disp('DATA SET SAVED!')