package NewHTHP
  extends Modelica.Icons.Package;
  import Modelica.Units.SI;
  import Refrigerant = Fluids.R601;
  import Fluid = Fluids.Water;

  package Machines "This package contains blocks that convert mechanical shaft work into 
                                            a change in fluid enthalpy"
    extends Modelica.Icons.Package;

    class Compressor
      import NewHTHP.Utilities.InitialisationMethod;

      function set_hs
        /* Off-design compressor map function */
        input SI.ThermodynamicTemperature T_in "[K] inlet temperature";
        input SI.Pressure p_in "[Pa] inlet pressure";
        input SI.Pressure p_out "[Pa] outlet pressure";
        input SI.SpecificEntropy s_in "[J/kgK] inlet entropy";
        output SI.SpecificEnthalpy hs "[J/kg] isentropic outlet enthalpy";
      protected
        Fluid.SaturationProperties sat "fluid saturated state";
        Fluid.ThermodynamicState state_l "saturated liquid thermodynamic state";
        Fluid.ThermodynamicState state_v "saturated vapour thermodynamic state";
        Real x(unit = "-") "vapour quality";
        Real k(unit = "-") "isentropic exponent";
        SI.ThermodynamicTemperature T_s "[K] isentropic outlet temperature";
      algorithm
        sat := Fluid.setSat_p(p_out);
        state_l := Fluid.setBubbleState(sat, 1);
        state_v := Fluid.setDewState(sat, 1);
        x := (s_in - Fluid.specificEntropy(state_l))/(Fluid.specificEntropy(state_v) - Fluid.specificEntropy(state_l));
        if x > 1.01 then
          k := Fluid.specificHeatCapacityCp(Fluid.setState_pT(p_in, T_in))/Fluid.specificHeatCapacityCv(Fluid.setState_pT(p_in, T_in));
          T_s := T_in*(p_out/p_in)^((k - 1)/k);
          hs := Fluid.specificEnthalpy_pT(p_out, T_s);
        else
          hs := x*Fluid.specificEnthalpy(state_v) + (1 - x)*Fluid.specificEnthalpy(state_l);
        end if;
      end set_hs;

      function offdesign_compressor_map
        /* Off-design compressor map function */
        input SI.MassFlowRate m_dot "[kg/s] inlet mass flow rate";
        input SI.ThermodynamicTemperature T_in "[K] inlet temperature";
        input SI.Pressure p_in "[Pa] inlet pressure";
        input SI.AngularVelocity omega "[rad/s] shaft speed";
        output Real PI "[-] compression ratio";
        output Real eta "[-] isentropic efficiency";
      protected
        parameter Real[6] design_values = {16.9, 70.11 + 273.15, 116700, 30000*2*Modelica.Constants.pi/60, 10.26, 0.8};
        parameter Real p = 1.8;
        parameter Real m = 1.8;
        parameter Real c4 = 0.3;
        Real G_dot;
        Real n_dot;
        Real[3] c;
      algorithm
        G_dot := (m_dot*sqrt(T_in)/p_in)/(design_values[1]*sqrt(design_values[2])/design_values[3]);
        n_dot := (omega/sqrt(T_in))/(design_values[4]/sqrt(design_values[2]));
        c[1] := n_dot/(p*(1 - m/n_dot) + n_dot*(n_dot - m)^2);
        c[2] := (p - 2*m*n_dot^2)/(p*(1 - m/n_dot) + n_dot*(n_dot - m)^2);
        c[3] := -(p*m*n_dot - m^2*n_dot^3)/(p*(1 - m/n_dot) + n_dot*(n_dot - m)^2);
        PI := design_values[5]*(c[1]*G_dot^2 + c[2]*G_dot + c[3]);
        eta := design_values[6]*((1 - c4*(1 - n_dot))*(n_dot/G_dot)*(2 - n_dot/G_dot));
// The simplest - not physical
//PI := design_values[5];
//eta := design_values[6];
// Simple - physical
// PI := if omega <= (30000*2*Modelica.Constants.pi/60) then 0 + ((design_values[5])/(30000*2*Modelica.Constants.pi/60))*omega else design_values[5];
// eta := ((-design_values[6]+0.5)/(30000*2*Modelica.Constants.pi/60)^2)*(omega - (30000*2*Modelica.Constants.pi/60))^2 + design_values[6];
      end offdesign_compressor_map;

      import offdesign_compressor_map;
      import set_hs;
      // Input/Output
      parameter Boolean omega_fromInput = false "Direct omega input" annotation(
        Dialog(tab = "Input/Output", group = "From Input"));
      parameter Boolean enableOutput = false "Output connection for selectable quantity" annotation(
        Dialog(tab = "Input/Output", group = "Connections"));
      parameter Boolean enableHeatPort = false "Include access heatport" annotation(
        Dialog(tab = "Input/Output", group = "Connections"));
      // Model parameters
      parameter SI.MomentOfInertia J_p = 5e-4 "Moment of inertia" annotation(
        Dialog(group = "Parameters", enable = not omega_fromInput));
      // Regularisation parameters
      parameter SI.Density rho_min = 1e-3 "Minimum allowed density" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      parameter SI.Pressure p_min = 5000 "Minimum allowed pressure" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      parameter SI.MassFlowRate epsilon_mdot = 0.001 "Flow regularisation parameter" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      parameter SI.AngularVelocity epsilon_omega = 1 "Angular velocity regularisation parameter" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      // Advanced
      parameter StateSelect omegaStateSelect = if omega_fromInput then StateSelect.default else StateSelect.prefer "State select for omega" annotation(
        Dialog(tab = "Advanced"));
      // Initialisation
      parameter Boolean initPhi = true "If true phi is initialised" annotation(
        Dialog(tab = "Initialisation", group = "Angular", enable = not omega_fromInput));
      parameter SI.Angle phi_0 = 0 "Initial value for phi" annotation(
        Dialog(tab = "Initialisation", group = "Angular", enable = initPhi));
      parameter Boolean initTau = true "If true tau is initialised" annotation(
        Dialog(tab = "Initialisation", group = "Rotational", enable = not omega_fromInput));
      parameter SI.Torque tau_0 = 1900000*60*0.97/(30000*2*Modelica.Constants.pi) "Initial value for tau" annotation(
        Dialog(tab = "Initialisation", group = "Rotational", enable = initTau));
      parameter InitialisationMethod initOmega = NewHTHP.Utilities.InitialisationMethod "Initialisation method for omega" annotation(
        Dialog(tab = "Initialisation", group = "Rotational", enable = not omega_fromInput));
      parameter SI.AngularVelocity omega_0 = 25000*(2*Modelica.Constants.pi/60) "Initial value for omega" annotation(
        Dialog(tab = "Initialisation", group = "Rotational", enable = (initOmega == InitialisationMethod.state)));
      parameter SI.AngularAcceleration omega_prime_0 = 0 "Initial value for der(omega)" annotation(
        Dialog(tab = "Initialisation", group = "Rotational", enable = (initOmega == InitialisationMethod.derivative)));
      parameter InitialisationMethod initMdot = NewHTHP.Utilities.InitialisationMethod.none "Initialisation method for m_flow" annotation(
        Dialog(tab = "Initialisation", group = "Flow"));
      parameter SI.MassFlowRate mdot_0 = 16.9 "Initial value for mass flow rate" annotation(
        Dialog(tab = "Initialisation", group = "Flow", enable = (initMdot == InitialisationMethod.state)));
      parameter Boolean initEnthalpy = false "If true hs is initialised" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = not omega_fromInput));
      parameter SI.SpecificEnthalpy hs_0 = 500000 "Initial value for outlet specific enthalpy" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = initEnthalpy));
      //parameter Real[3] coeff = {0, 0, 1/28.561}  "Caracteristic curve coefficients";
      //Connections
      Connections.FluidInlet inlet annotation(
        Placement(transformation(extent = {{-115, 15}, {-85, -15}})));
      Connections.FluidOutlet outlet annotation(
        Placement(transformation(extent = {{85, 15}, {115, -15}})));
      Modelica.Mechanics.Rotational.Interfaces.Flange_a flange(phi = phi, tau = tau) if not omega_fromInput "Flange to receive mechanical power" annotation(
        Placement(transformation(extent = {{-20, -20}, {20, 20}}, origin = {0, -100}, rotation = -90), iconTransformation(extent = {{-20, -20}, {20, 20}}, origin = {0, -100}, rotation = -90)));
      Modelica.Blocks.Interfaces.RealInput omega_input(unit = "rad/s") = omega if omega_fromInput "Input to directly set pump speed [rad/s]" annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {80, -100}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatport(Q_flow = Q_t) if enableHeatPort "Access-heat dumping port" annotation(
        Placement(transformation(extent = {{-20, -20}, {20, 20}}, origin = {0, 100}, rotation = 90), iconTransformation(extent = {{-20, -20}, {20, 20}}, origin = {0, 100}, rotation = 90)));
      //Modelica.Blocks.Interfaces.RealOutput output_val(unit = Sensors.Internal.getFlowUnit(outputQuantity)) = getQuantity(inlet.State, mdot, outputQuantity, rho_min) if enableOutput "Measured value [variable]" annotation(
      //   Placement(transformation(extent = {{-20, -20}, {20, 20}}, origin = {-80, -100}, rotation = 270)));
      Modelica.Blocks.Interfaces.RealOutput output_val(unit = "rad/s") if enableOutput "Measured value [variable]" annotation(
        Placement(transformation(extent = {{-20, -20}, {20, 20}}, origin = {-80, -100}, rotation = 270)));
      Real PI(unit = "-");
      Real eta(unit = "-");
      SI.SpecificEnthalpy h_out;
      // Guess values
      SI.MassFlowRate mdot;
      SI.SpecificEnthalpy hs;
      SI.SpecificEnthalpy dh;
      SI.Angle phi;
      SI.Torque tau;
      SI.Torque tau_s;
      SI.Torque tau_normalized "moment after zero massflow normalization";
      SI.AngularVelocity omega;
      SI.Power Q_t(start = 0) "work that could not be performed on the fluid, and is dumped to heat port";
      SI.Power W_t "technichal work performed on fluid";
    initial equation
      if initPhi then
        phi = phi_0;
      end if;
      if initTau then
        tau = tau_0;
      end if;
      if initOmega == InitialisationMethod.state then
        omega = omega_0;
      elseif initOmega == InitialisationMethod.derivative then
        der(omega) = omega_prime_0;
      elseif initOmega == InitialisationMethod.steadyState then
        der(omega) = 0;
      end if;
      if initMdot == InitialisationMethod.state then
        mdot = mdot_0;
      elseif initMdot == InitialisationMethod.steadyState then
        der(mdot) = 0;
      end if;
      if initEnthalpy then
        hs = hs_0;
      end if;
    equation
      inlet.mdot + outlet.mdot = 0;
      mdot = inlet.mdot;
      (PI, eta) = offdesign_compressor_map(mdot, inlet.State.T, inlet.p, omega);
      outlet.p = max(PI*inlet.p, p_min);
      hs = set_hs(inlet.State.T, inlet.p, outlet.p, inlet.State.s);
      h_out = inlet.State.h + (hs - inlet.State.h)/eta;
//outlet.State = Fluid.setState_ph(outlet.p, h_out);
      outlet.State = Fluid.setState_ph(outlet.p, inlet.State.h + dh);
      tau_s = (h_out - inlet.State.h)*mdot*omega/(omega^2 + (epsilon_omega)^2);
//A characteristic function of mdot on rotational speed - CANCEL IF APPROPRIATE
//mdot = epsilon_mdot + (16.9 - epsilon_mdot)*omega/(30000*(2*Modelica.Constants.pi)/60);
// (outlet.p - inlet.p) = coeff[1] + coeff[2]*mdot + coeff[3]*mdot^2;
      W_t = tau_s*omega;
      dh = (mdot*W_t)/(mdot^2 + (epsilon_mdot)^2);
      Q_t = W_t - mdot*dh;
      tau_normalized = tau_s;
      if noEvent(omega_fromInput) then
        tau = 0;
        phi = 0;
      else
        J_p*der(omega) = tau - tau_normalized;
        omega = der(phi);
      end if;
      if enableOutput then
        output_val = omega;
      end if;
//Annotation
      annotation(
        Icon(graphics = {Line(points = {{-75, 0}, {-100, 0}}, thickness = 1.5), Line(points = {{75, 0}, {100, 0}}, thickness = 1.5), Ellipse(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-75, 75}, {75, -75}}), Line(points = {{-53.03, 53.03}, {64.95, 37.5}}, thickness = 1), Line(points = {{-53.03, -53.03}, {64.95, -37.5}}, thickness = 1)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true)),
        Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true)),
        version = "",
        uses);
    end Compressor;

    package Test "Simple examples to test the models in this package"
      extends Modelica.Icons.ExamplesPackage;

      model Compressor_test
        NewHTHP.Boundaries.Source source(T0_par = 70.11 + 273.15, p0_par(displayUnit = "Pa") = 116700, pressureFromInput = false, setEnthalpy = false) annotation(
          Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        NewHTHP.Boundaries.Sink sink(p0_par(displayUnit = "Pa") = 300000, pressureFromInput = false) annotation(
          Placement(visible = true, transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        Modelica.Blocks.Sources.Step step(height = 30000*(2*Modelica.Constants.pi/60), offset = 25000*(2*Modelica.Constants.pi/60), startTime = 300) annotation(
          Placement(visible = true, transformation(origin = {70, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Continuous.LimPID pid(Ti = 10, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.Init.InitialOutput, k = 0.1, yMax = 60000) annotation(
          Placement(visible = true, transformation(origin = {32, -18}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Sources.Torque torque(useSupport = false) annotation(
          Placement(visible = true, transformation(origin = {-32, 16}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Gain shaft_eta(k = 0.97) annotation(
          Placement(visible = true, transformation(origin = {-2, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Compressor compressor(J_p = 5, enableHeatPort = false, enableOutput = true, initEnthalpy = false, initOmega = NewHTHP.Utilities.InitialisationMethod.state, initPhi = true, mdot_0 = 0, omega_fromInput = false) annotation(
          Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      equation
        connect(step.y, pid.u_s) annotation(
          Line(points = {{60, 0}, {51, 0}, {51, -18}, {44, -18}}, color = {0, 0, 127}));
        connect(pid.y, shaft_eta.u) annotation(
          Line(points = {{21, -18}, {15, -18}, {15, 0}, {10, 0}}, color = {0, 0, 127}));
        connect(shaft_eta.y, torque.tau) annotation(
          Line(points = {{-13, 0}, {-17.5, 0}, {-17.5, 16}, {-20, 16}}, color = {0, 0, 127}));
        connect(compressor.output_val, pid.u_m) annotation(
          Line(points = {{-60, -8}, {-29, -8}, {-29, -78}, {32.5, -78}, {32.5, -30}, {32, -30}}, color = {0, 0, 127}));
        connect(torque.flange, compressor.flange) annotation(
          Line(points = {{-42, 16}, {-52, 16}, {-52, 0}, {-60, 0}}));
        connect(source.outlet, compressor.inlet) annotation(
          Line(points = {{-70, -22}, {-70, -10}}, color = {85, 0, 255}, thickness = 1));
        connect(compressor.outlet, sink.inlet) annotation(
          Line(points = {{-70, 10}, {-70, 22}}, color = {85, 0, 255}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
      end Compressor_test;

      model Compressor_withLoop_test
        Modelica.Blocks.Sources.Step step(height = 27000*(2*Modelica.Constants.pi/60), offset = 25000*(2*Modelica.Constants.pi/60), startTime = 300) annotation(
          Placement(visible = true, transformation(origin = {70, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Sources.Torque torque(useSupport = false) annotation(
          Placement(visible = true, transformation(origin = {-32, 16}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Gain shaft_eta(k = 0.97) annotation(
          Placement(visible = true, transformation(origin = {-2, 16}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        NewHTHP.Machines.Compressor compressor(J_p = 5, enableHeatPort = false, enableOutput = true, initEnthalpy = false, initMdot = NewHTHP.Utilities.InitialisationMethod.none, initOmega = NewHTHP.Utilities.InitialisationMethod.state, initPhi = true, mdot_0 = 16.9, omega_fromInput = false) annotation(
          Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        Modelica.Blocks.Sources.Step step1(height = 0, offset = 1814, startTime = 0) annotation(
          Placement(visible = true, transformation(origin = {76, -80}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Continuous.LimPID pid(Ti = 10, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.Init.InitialOutput, k = 0.1, yMax = 60000, y_start = 600) annotation(
          Placement(visible = true, transformation(origin = {32, -18}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Valves.JTValve jTValve annotation(
          Placement(visible = true, transformation(origin = {-94, 0}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
        Volumes.MixedVolume mixedVolume(T_0 = 114 + 273.15, initPressure = true, p_0 = 116700000000) annotation(
          Placement(visible = true, transformation(origin = {-80, 44}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      equation
        connect(torque.flange, compressor.flange) annotation(
          Line(points = {{-42, 16}, {-52, 16}, {-52, 0}, {-60, 0}}));
        connect(pid.y, shaft_eta.u) annotation(
          Line(points = {{21, -18}, {15, -18}, {15, 16}, {10, 16}}, color = {0, 0, 127}));
        connect(step.y, pid.u_s) annotation(
          Line(points = {{60, 0}, {51, 0}, {51, -18}, {44, -18}}, color = {0, 0, 127}));
        connect(compressor.output_val, pid.u_m) annotation(
          Line(points = {{-60, -8}, {-29, -8}, {-29, -78}, {32.5, -78}, {32.5, -30}, {32, -30}}, color = {0, 0, 127}));
        connect(compressor.inlet, jTValve.outlet) annotation(
          Line(points = {{-70, -10}, {-70, -34}, {-94, -34}, {-94, -6}}, color = {85, 0, 255}));
        connect(mixedVolume.outlet, jTValve.inlet) annotation(
          Line(points = {{-90, 44}, {-94, 44}, {-94, 6}}, color = {85, 0, 255}));
        connect(compressor.outlet, mixedVolume.inlet) annotation(
          Line(points = {{-70, 10}, {-70, 44}}, color = {85, 0, 255}));
        connect(shaft_eta.y, torque.tau) annotation(
          Line(points = {{-13, 16}, {-20, 16}}, color = {0, 0, 127}));
        annotation(
          Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
      end Compressor_withLoop_test;
    equation

    end Test;
    annotation(
      Icon(graphics = {Line(points = {{-75, 0}, {-100, 0}}, thickness = 1.5), Line(points = {{75, 0}, {100, 0}}, thickness = 1.5), Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1.5, extent = {{-75, 75}, {75, -75}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true)),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true)),
      version = "",
      uses);
  end Machines;

  package HeatExchangers
    extends Modelica.Icons.Package;

    class CoFlow_Evaporator_Upwind
      function correct_Shah
        input Real Co "[-] convective number";
        input Real Bo "[-] boiling number";
        output Real F_e "[-] correction factor";
      protected
        Real h_b[4](each unit = "[-]") "coefficients for individual boiling regimes";
      algorithm
        h_b[1] := 1.8*Co^(-0.8);
        if Bo >= 0.3e-4 then
          h_b[2] := 230*Bo^0.5;
        else
          h_b[2] := 1 + 46*Bo^0.5;
        end if;
        if Bo >= 11e-4 then
          h_b[3] := 14.70*(Bo^0.5)*Modelica.Math.exp(2.74*(Co)^(-0.1));
        else
          h_b[3] := 15.43*(Bo^0.5)*Modelica.Math.exp(2.74*(Co)^(-0.1));
        end if;
        if Bo >= 11e-4 then
          h_b[4] := 14.70*(Bo^0.5)*Modelica.Math.exp(2.47*(Co)^(-0.15));
        else
          h_b[4] := 15.43*(Bo^0.5)*Modelica.Math.exp(2.47*(Co)^(-0.15));
        end if;
        if Co <= 0.1 then
          F_e := max(h_b[1], h_b[4]);
        elseif Co > 1 then
          F_e := max(h_b[1], h_b[2]);
        else
          F_e := max(h_b[1], h_b[3]);
        end if;
      end correct_Shah;

      function alpha_Shah
        import correct_Shah;
        input Fluid.ThermodynamicState state "fluid thermodyanmic state";
        input Real d_h "[m] inner channel hydraulic diameter";
        input Real m_dot "[kg/s] mass flow rate through single channel";
        output Real alpha "[W/m2K] convective HTC";
      protected
        Fluid.SaturationProperties sat "fluid saturated state";
        Fluid.ThermodynamicState state_l "saturated liquid thermodynamic state";
        Fluid.ThermodynamicState state_v "saturated vapour thermodynamic state";
        Real rho[3](each unit = "kg/m3") "fluid density";
        Real mu[3](each unit = "Pas") "fluid viscosity";
        Real lambda[3](each unit = "W/mK") "fluid thermal conductivity";
        Real Pr[3](each unit = "-") "Prandtl number";
        Real alpha_t[2](each unit = "W/m2K") "temporary convective HTC";
        Real Re(unit = "-") "Reynolds number";
        Real Nu(unit = "-") "Nusselt number";
        Real Co(unit = "-") "Convection number";
        Real Bo(unit = "-") "Boiling number";
        Real F_e(unit = "-") "Correction factor";
        Real f(unit = "-") "Moody's friction factor, i.e. 4C_f frinction coefficient";
        Real x(unit = "-") "vapour quality";
        Real A_c(unit = "m2") "inner pipe cross-section";
      algorithm
        sat := Fluid.setSat_p(state.p);
        state_l := Fluid.setBubbleState(sat, 1);
        state_v := Fluid.setDewState(sat, 1);
        rho[1] := Fluid.density(state_l);
        rho[2] := Fluid.density(state);
        rho[3] := Fluid.density(state_v);
        mu[1] := Fluid.dynamicViscosity(state_l);
        mu[2] := Fluid.dynamicViscosity(state);
        mu[3] := Fluid.dynamicViscosity(state_v);
        lambda[1] := Fluid.thermalConductivity(state_l);
        lambda[2] := Fluid.thermalConductivity(state);
        lambda[3] := Fluid.thermalConductivity(state_v);
        Pr[1] := Fluid.prandtlNumber(state_l);
        Pr[2] := Fluid.prandtlNumber(state);
        Pr[3] := Fluid.prandtlNumber(state_v);
        x := (Fluid.specificEnthalpy(state) - Fluid.specificEnthalpy(state_l))/(Fluid.specificEnthalpy(state_v) - Fluid.specificEnthalpy(state_l));
        A_c := Modelica.Constants.pi*((d_h)^2)/4;
        if x <= 0 or x >= 1 then
          Re := (rho[2]*d_h/mu[2])*(m_dot/(A_c*rho[2]));
          f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
          Nu := ((f/8)*(Re - 1000)*Pr[2])/(1 + 12.7*((f/8)^0.5)*(Pr[2]^(2/3) - 1));
          alpha := Nu*lambda[2]/d_h;
        else
          Re := (rho[1]*d_h/mu[1])*(m_dot*(1 - x)/(A_c*rho[1]));
          if Re >= 2300 then
            if x >= 0.05 then
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := ((f/8)*(Re - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
              alpha := Nu*lambda[1]/d_h;
              Co := (((1 - x)/x)^0.8)*(rho[3]/rho[1])^0.5;
              Bo := (((34*(state.p/1e5)^(0.18))/(1 - 0.0045*(state.p/1e5)))^3)*A_c/(m_dot*(Fluid.specificEnthalpy(state_v) - Fluid.specificEnthalpy(state_l)));
              F_e := correct_Shah(Co, Bo);
              alpha := F_e*alpha;
            else
              Re := (rho[1]*d_h/mu[1])*(m_dot/(A_c*rho[1]));
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := ((f/8)*(Re - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
              alpha_t[1] := Nu*lambda[1]/d_h;
              Re := (rho[1]*d_h/mu[1])*(m_dot*0.95/(A_c*rho[1]));
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := ((f/8)*(Re - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
              alpha_t[2] := Nu*lambda[1]/d_h;
              Co := (((1 - 0.05)/0.05)^0.8)*(rho[3]/rho[1])^0.5;
              Bo := (((34*(state.p/1e5)^(0.18))/(1 - 0.0045*(state.p/1e5)))^3)*A_c/(m_dot*(Fluid.specificEnthalpy(state_v) - Fluid.specificEnthalpy(state_l)));
              F_e := correct_Shah(Co, Bo);
              alpha_t[2] := F_e*alpha_t[2];
              alpha := alpha_t[2] + ((alpha_t[1] - alpha_t[2])/(0 - 0.05))*(x - 0.05);
            end if;
          else
            Re := (rho[3]*d_h/mu[3])*(m_dot/(A_c*rho[3]));
            f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
            Nu := ((f/8)*(Re - 1000)*Pr[3])/(1 + 12.7*((f/8)^0.5)*(Pr[3]^(2/3) - 1));
            alpha_t[1] := Nu*lambda[3]/d_h;
            f := (0.790*Modelica.Math.log(2300) - 1.64)^(-2);
            Nu := ((f/8)*(2300 - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
            alpha_t[2] := Nu*lambda[1]/d_h;
            Co := (((1 - (1 - 2300*mu[1]*A_c/(d_h*m_dot)))/(1 - 2300*mu[1]*A_c/(d_h*m_dot)))^0.8)*(rho[3]/rho[1])^0.5;
            Bo := (((34*(state.p/1e5)^(0.18))/(1 - 0.0045*(state.p/1e5)))^3)*A_c/(m_dot*(Fluid.specificEnthalpy(state_v) - Fluid.specificEnthalpy(state_l)));
            F_e := correct_Shah(Co, Bo);
            alpha_t[2] := F_e*alpha_t[2];
            alpha := alpha_t[2] + ((alpha_t[1] - alpha_t[2])/(1 - (1 - 2300*mu[1]*A_c/(d_h*m_dot))))*(x - (1 - 2300*mu[1]*A_c/(d_h*m_dot)));
          end if;
        end if;
      end alpha_Shah;

      import alpha_Shah;
      import correct_Shah;
      //Model parameters
      parameter Integer N = 12 "# of HX nodes" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length d_h = 0.0254 "[mm] pipe A hydraulic diameter" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real N_p(unit = "-") = 4 "# of pipes A - so that u_max < 1 m/s" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Area A = 7000 "[m2] heat transfer surface" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length L = 10 "[m] linear dimension heat exchanger" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real csi(unit = "-") = 0.2 "material share" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Density rho_w = 7800 "[kg/m3] Density HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.SpecificHeatCapacity cp_w = 477 "[J/kgK] Density HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.ThermalConductivity lambda_w = 14.9 "[W/mK] Thermal conductivity HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.Density rho_HTF = 1000 "[kg/m3] Density HTF" annotation(
        Dialog(tab = "Parameters", group = "Heat transfer"));
      parameter SI.SpecificHeatCapacity cp_HTF = 4186 "[J/kgK] Density HTF" annotation(
        Dialog(tab = "Parameters", group = "Heat transfer"));
      parameter SI.CoefficientOfHeatTransfer alpha_HTF = 1500 "[W/m2K] Heat trasnfer coefficient HTF" annotation(
        Dialog(tab = "Parameters", group = "Heat transfer"));
      parameter SI.Pressure Delta_p = 20000 "[Pa] Pressure drop across HX" annotation(
        Dialog(tab = "Parameters", group = "Fluid"));
      // Initialisation
      parameter SI.ThermodynamicTemperature T_0 = 20 + 273.15 "intial temperature" annotation(
        Dialog(tab = "Initialisation", group = "Thermal"));
      // Connections
      NewHTHP.Connections.FluidInlet inletA annotation(
        Placement(visible = true, transformation(origin = {0, -20}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, -20}, extent = {{-115, 15}, {-85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidOutlet outletA annotation(
        Placement(visible = true, transformation(origin = {0, -20}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {0, -20}, extent = {{85, 15}, {115, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidInlet_Set inletB annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, 20}, extent = {{-115, 15}, {-85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidOutlet_Set outletB annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {0, 20}, extent = {{85, 15}, {115, -15}}, rotation = 0)));
      SI.ThermodynamicTemperature T[2, N](each start = T_0);
      SI.ThermodynamicTemperature T_w[N](each start = T_0);
      Fluid.ThermodynamicState stateA[N] "fluid thermodynamic state";
      SI.CoefficientOfHeatTransfer alpha[N](each start = 10000);
      SI.SpecificEnthalpy h[N](each start = Fluid.specificEnthalpy_pT(100000, T_0));
      SI.Density rho[N](each start = Fluid.density_pT(100000, T_0));
    protected
      parameter SI.Volume V = A/450 "[m3] volume";
    equation
      T[1, 1] = inletA.State.T;
      stateA[1] = inletA.State;
      h[1] = stateA[1].h;
      rho[1] = stateA[1].d;
      T[2, 1] = inletB.T;
      alpha[1:N] = alpha_Shah(stateA[1:N], d_h, inletA.mdot/N_p);
//System of ODEs
      (V*(1 - csi)/2)/(N - 1)*rho[2:N - 1].*der(h[2:N - 1]) = inletA.mdot*(h[1:N - 2] - h[2:N - 1]) - A/(N - 1)*alpha[2:N - 1].*(T[1, 2:N - 1] - T_w[2:N - 1]);
      (V*(1 - csi)/2)/(N - 1)/2*rho[N].*der(h[N]) = inletA.mdot*(h[N - 1] - h[N]) - A/(N - 1)/2*alpha[N].*(T[1, N] - T_w[N]);
      for ii in 2:N loop
        stateA[ii] = Fluid.setState_ph(inletA.p - ((Delta_p)/(N - 1))*(ii - 1), h[ii]);
      end for;
      T[1, 2:N] = stateA[2:N].T;
      rho[2:N] = stateA[2:N].d;
      V*(csi)/(N - 1)/2*rho_w*cp_w*der(T_w[1]) = lambda_w*(T_w[3] - 2*T_w[2] + T_w[1])/((L/(N - 1))^2) + A/(N - 1)/2*alpha[1].*(T[1, 1] - T_w[1]) + A/(N - 1)/2*alpha_HTF*(T[2, 1] - T_w[1]);
      V*(csi)/(N - 1)*rho_w*cp_w*der(T_w[2:N - 1]) = lambda_w*(T_w[3:N] - 2*T_w[2:N - 1] + T_w[1:N - 2])/((L/(N - 1))^2) + A/(N - 1)*alpha[2:N - 1].*(T[1, 2:N - 1] - T_w[2:N - 1]) + A/(N - 1)*alpha_HTF*(T[2, 2:N - 1] - T_w[2:N - 1]);
      V*(csi)/(N - 1)/2*rho_w*cp_w*der(T_w[N]) = lambda_w*(T_w[N] - 2*T_w[N - 1] + T_w[N - 2])/((L/(N - 1))^2) + A/(N - 1)/2*alpha[N].*(T[1, N] - T_w[N]) + A/(N - 1)/2*alpha_HTF*(T[2, N] - T_w[N]);
      (V*(1 - csi)/2)/(N - 1)*cp_HTF*rho_HTF*der(T[2, 2:N - 1]) = inletB.mdot*cp_HTF*(T[2, 1:N - 2] - T[2, 2:N - 1]) - A/(N - 1)*alpha_HTF.*(T[2, 2:N - 1] - T_w[2:N - 1]);
      (V*(1 - csi)/2)/(N - 1)/2*cp_HTF*rho_HTF*der(T[2, N]) = inletB.mdot*cp_HTF*(T[2, N - 1] - T[2, N]) - A/(N - 1)/2*alpha_HTF*(T[2, N] - T_w[N]);
      outletA.p = stateA[N].p;
      outletA.State = stateA[N];
      outletB.T = T[2, N];
      outletB.p = inletB.p;
      outletA.mdot + inletA.mdot = 0;
      outletB.mdot + inletB.mdot = 0;
// Annotation
      annotation(
        Icon(graphics = {Rectangle(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-75, 50}, {75, -50}}), Line(points = {{-100, 20}, {-75, 20}}, thickness = 1.5), Line(points = {{-100, -20}, {-75, -20}}, thickness = 1.5), Line(points = {{75, 20}, {100, 20}}, thickness = 1.5), Line(points = {{100, -20}, {75, -20}}, thickness = 1.5), Bitmap(extent = {{-25, -25}, {25, 25}}, imageSource = "iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAMAAADVRocKAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAzUExURQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAKMFRskAAAAQdFJOUwAQIDBAUGBwgI+fr7/P3+8jGoKKAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAEvElEQVRoQ+2Z2YKDIAxFh7oAopL//9q5gYjgQrXLW8/LoG2RhHCTOH8/fvx4iurMYFq5+ALtTMyg5PrTaKLR9IbIy40P8/Ck+a+ayIUbn8aQjYMH0Vec5KiT0Uhf2eiZHjKy0VefZl33assRqtXWdtecqLp+WTWvW/YWu73e3dFMIZbJyHUNNeKLvVz8KRkrt+z2EdrTZLV2RJPcqWBowheTsR2R052ZaT63H9bZ8GkzXwhmhyVPWcT00fix4qDVOhhc26mApWkgn02neuuqWqTIJ+t6GmR0Cs7stc1KtDTKiM/jLKMKvW5kdA2d7T+skdEH6bKdba7E0V3yVZvne/ACq9IiYL+hWFDaGBWN+4oB4axMpu0HHOhrcnQbHOFAPNDXaDrdH/uzNePk/eRs/nFnxuFOeLcijzuFU1o+YfyrKUJZT36yo8ccpRX84Nm2jVINKyDNT6XnCGiz11EeMV8+BXR5zK65pLmlKoKhOXnT5MIHFYwPTnBVI8PrND76xc0D5nbrDJp8Wr6DC3kVDXJBvHOdRd7hHqxepbzf+cxdYav5e+2FFLAhZd925INpxc3Q8TxoHh2XezC1uyugj/UHYWhEafTumDZj+Op0c6NzxZ0wdyebMB/o2MgPb26asM0ZOgrYodLHHHz06Ar4UZaROTbDnqcqtSDov70ZSEPhU0RRiJLjyi7s8OLEqyDwMtWy8uspv5kI/rydiC3NyUk9+TjxWgYXcBjcz/Q4RbF4VAMtwT+fx+L9B2AbyDtjJ4jqUqbiDEBBkSP2hrxSSvQxRfmsYmwd1JuZ7eYZ7YWKdE/bG50V/CEDIYsNnCQ2YXkcwfdAIvBaFt5xmZ4bIYH8Diz72ZS4lNhiuverOSQCiSWle54ZXcpa5I4vZbUc1GvLFOhHaMIjcL6XjcXTNyp7myETaj0gvmBOOu5tnoZeA3JcqCUSDa4llWL+dx20j0I4BU4KomteSMk7pp0PZnbPADlBzVGUGUULfBkiGSRC18+G9VPx7LIFvspB29WxX7ixUWX8aJrzFvgiR1LGjsh7P2HbAl/jTIwPOtZdC3yJswfkzaVwvwUOnOQy9seOuy1wYCntSlT19cstcJIPpjI3S4kaWYGdgETcjZZz0KbmZS/z5ltI1eg2Py6oXcoNRbVReX/0DBWze97hocrIgx6t06Z3A+dt6YY+NYBZ4uVeYDGisTHpFPShIyFyT0MLOWsM1mOZmQiFV9huMHbATN5t/IO2cVnVbrs2ILzTN8rOQrNUMt5tl4+Q8rHEUWxduHfG0s8wiJ7S4A61Ut/sdheLWl/ON3M9BYVcsjA8szcyyKLczLZhWRW5KJXNXFIu9Ahx/bzPWBHa7HB5SCn++cuyc1xaxgMVGSK11lKVFhwo8gG56mle0rE6CoU2c15/Stb2ygLzQNmR3oaD9YV4rVoovZo3vocg5JZCATWCeKhaLewznjS+J+Akx2nRESyCVq0Wtm0ve7a6dVybx/YpVez1amHcnpZn/YL8tyGr2erVQjoHAuw9MTahOqNzRXhSLYQ3FgmExgvldrVawC6sTofsXjk890CzIAmCXyd+4+Usb9vSWG/e6n0KpLzAWHHlm8RkIRc/fvzI+fv7BzDINdD0LCjWAAAAAElFTkSuQmCC")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end CoFlow_Evaporator_Upwind;

    class CtrFlow_Condenser_Upwind
      function alpha_Shah
        input Fluid.ThermodynamicState state "fluid thermodyanmic state";
        input SI.Length d_h "[m] inner channel hydraulic diameter";
        input SI.MassFlowRate m_dot "[kg/s] mass flow rate through single channel";
        output SI.CoefficientOfHeatTransfer alpha "[W/m2K] convective HTC";
      protected
        Fluid.SaturationProperties sat "fluid saturated state";
        Fluid.ThermodynamicState state_l "saturated liquid thermodynamic state";
        Fluid.ThermodynamicState state_v "saturated vapour thermodynamic state";
        Real rho[3](each unit = "kg/m3") "fluid density";
        Real mu[3](each unit = "Pas") "fluid viscosity";
        Real lambda[3](each unit = "W/mK") "fluid thermal conductivity";
        Real Pr[3](each unit = "-") "Prandtl number";
        Real alpha_t[2](each unit = "W/m2K") "temporary convective HTC";
        Real Re(unit = "-") "Reynolds number";
        Real Nu(unit = "-") "Nusselt number";
        Real Z(unit = "-") "A correction parameter";
        Real f(unit = "-") "Moody's friction factor, i.e. 4C_f frinction coefficient";
        Real x(unit = "-") "vapour quality";
        Real A_c(unit = "m2") "inner pipe cross-section";
      algorithm
        sat := Fluid.setSat_p(state.p);
        state_l := Fluid.setBubbleState(sat, 1);
        state_v := Fluid.setDewState(sat, 1);
        rho[1] := Fluid.density(state_l);
        rho[2] := Fluid.density(state);
        rho[3] := Fluid.density(state_v);
        mu[1] := Fluid.dynamicViscosity(state_l);
        mu[2] := Fluid.dynamicViscosity(state);
        mu[3] := Fluid.dynamicViscosity(state_v);
        lambda[1] := Fluid.thermalConductivity(state_l);
        lambda[2] := Fluid.thermalConductivity(state);
        lambda[3] := Fluid.thermalConductivity(state_v);
        Pr[1] := Fluid.prandtlNumber(state_l);
        Pr[2] := Fluid.prandtlNumber(state);
        Pr[3] := Fluid.prandtlNumber(state_v);
        x := (Fluid.specificEnthalpy(state) - Fluid.specificEnthalpy(state_l))/(Fluid.specificEnthalpy(state_v) - Fluid.specificEnthalpy(state_l));
        A_c := Modelica.Constants.pi*((d_h)^2)/4;
        if x <= 0 or x >= 1 then
          Re := (rho[2]*d_h/mu[2])*(m_dot/(A_c*rho[2]));
          f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
          Nu := ((f/8)*(Re - 1000)*Pr[2])/(1 + 12.7*((f/8)^0.5)*(Pr[2]^(2/3) - 1));
          alpha := Nu*lambda[2]/d_h;
        else
          Re := (rho[1]*d_h/mu[1])*(m_dot*(1 - x)/(A_c*rho[1]));
          if Re >= 100 then
            if x >= 0.05 then
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := 0.023*(Re^0.8)*(Pr[1]^0.4);
              alpha := Nu*lambda[1]/d_h;
              Z := (((1 - x)/x)^0.8)*(state.p/Fluid.getCriticalPressure())^0.4;
              alpha := alpha*(1 + 3.8/(Z^0.95));
            else
              Re := (rho[1]*d_h/mu[1])*(m_dot/(A_c*rho[1]));
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := ((f/8)*(Re - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
              alpha_t[1] := Nu*lambda[1]/d_h;
              Re := (rho[1]*d_h/mu[1])*(m_dot*0.95/(A_c*rho[1]));
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := 0.023*(Re^0.8)*(Pr[1]^0.4);
              alpha_t[2] := Nu*lambda[1]/d_h;
              Z := (((1 - 0.05)/0.05)^0.8)*(state.p/Fluid.getCriticalPressure())^0.4;
              alpha_t[2] := alpha_t[2]*(1 + 3.8/(Z^0.95));
              alpha := alpha_t[2] + ((alpha_t[1] - alpha_t[2])/(0 - 0.05))*(x - 0.05);
            end if;
          else
            Re := (rho[3]*d_h/mu[3])*(m_dot/(A_c*rho[3]));
            f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
            Nu := ((f/8)*(Re - 1000)*Pr[3])/(1 + 12.7*((f/8)^0.5)*(Pr[3]^(2/3) - 1));
            alpha_t[1] := Nu*lambda[3]/d_h;
            f := (0.790*Modelica.Math.log(100) - 1.64)^(-2);
            Nu := 0.023*(100^0.8)*(Pr[1]^0.4);
            alpha_t[2] := Nu*lambda[1]/d_h;
            Z := (((1 - (1 - 100*mu[1]*A_c/(d_h*m_dot)))/(1 - 100*mu[1]*A_c/(d_h*m_dot)))^0.8)*(state.p/Fluid.getCriticalPressure())^0.4;
            alpha_t[2] := alpha_t[2]*(1 + 3.8/(Z^0.95));
            alpha := alpha_t[2] + ((alpha_t[1] - alpha_t[2])/(1 - (1 - 100*mu[1]*A_c/(d_h*m_dot))))*(x - (1 - 100*mu[1]*A_c/(d_h*m_dot)));
          end if;
        end if;
      end alpha_Shah;

      import alpha_Shah;
      //Model parameters
      parameter Integer N = 12 "# of HX nodes" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length d_h = 0.0254 "[m] pipe A hydraulic diameter" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real N_p(unit = "-") = 4 "# of pipes A - so that u_max < 1 m/s" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Area A = 7000 "[m2] heat transfer surface" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length L = 10 "[m] linear dimension heat exchanger" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real csi(unit = "-") = 0.2 "material share" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Density rho_w = 7800 "[kg/m3] Density HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.SpecificHeatCapacity cp_w = 477 "[J/kgK] Density HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.ThermalConductivity lambda_w = 14.9 "[W/mK] Thermal conductivity HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.Density rho_HTF = 1000 "[kg/m3] Density HTF" annotation(
        Dialog(tab = "Parameters", group = "Heat transfer"));
      parameter SI.SpecificHeatCapacity cp_HTF = 4186 "[J/kgK] Density HTF" annotation(
        Dialog(tab = "Parameters", group = "Heat transfer"));
      parameter SI.CoefficientOfHeatTransfer alpha_HTF = 1500 "[W/m2K] Heat trasnfer coefficient HTF" annotation(
        Dialog(tab = "Parameters", group = "Heat transfer"));
      parameter SI.Pressure Delta_p = 20000 "[Pa] Pressure drop across HX" annotation(
        Dialog(tab = "Parameters", group = "Fluid"));
      // Initialisation
      parameter SI.ThermodynamicTemperature T_0 = 20 + 273.15 "[K] intial temperature" annotation(
        Dialog(tab = "Initialisation", group = "Thermal"));
      // Connections
      NewHTHP.Connections.FluidInlet inletA annotation(
        Placement(visible = true, transformation(origin = {0, -20}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, -20}, extent = {{115, 15}, {85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidOutlet outletA annotation(
        Placement(visible = true, transformation(origin = {0, -20}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {0, -20}, extent = {{-85, 15}, {-115, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidInlet_Set inletB annotation(
        Placement(visible = true, transformation(origin = {0, 20}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, 20}, extent = {{-115, 15}, {-85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidOutlet_Set outletB annotation(
        Placement(visible = true, transformation(origin = {0, 20}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {0, 20}, extent = {{85, 15}, {115, -15}}, rotation = 0)));
      SI.ThermodynamicTemperature T[2, N](each start = T_0);
      SI.ThermodynamicTemperature T_w[N](each start = T_0);
      Fluid.ThermodynamicState stateA[N] "fluid thermodynamic state";
      SI.CoefficientOfHeatTransfer alpha[N](each start = 10000);
      SI.SpecificEnthalpy h[N](each start = Fluid.specificEnthalpy_pT(100000, T_0));
      SI.Density rho[N](each start = Fluid.density_pT(100000, T_0));
    protected
      parameter SI.Volume V = A/450 "[m3] volume";
    equation
      T[1, 1] = inletA.State.T;
      stateA[1] = inletA.State;
      h[1] = stateA[1].h;
      rho[1] = stateA[1].d;
      T[2, N] = inletB.T;
      alpha[1:N] = alpha_Shah(stateA[1:N], d_h, inletA.mdot/N_p);
//System of ODEs
      (V*(1 - csi)/2)/(N - 1)*rho[2:N - 1].*der(h[2:N - 1]) = inletA.mdot*(h[1:N - 2] - h[2:N - 1]) - A/(N - 1)*alpha[2:N - 1].*(T[1, 2:N - 1] - T_w[2:N - 1]);
      (V*(1 - csi)/2)/(N - 1)/2*rho[N].*der(h[N]) = inletA.mdot*(h[N - 1] - h[N]) - A/(N - 1)/2*alpha[N].*(T[1, N] - T_w[N]);
      for ii in 2:N loop
        stateA[ii] = Fluid.setState_ph(inletA.p - ((Delta_p)/(N - 1))*(ii - 1), h[ii]);
      end for;
      T[1, 2:N] = stateA[2:N].T;
      rho[2:N] = stateA[2:N].d;
      V*(csi)/(N - 1)/2*rho_w*cp_w*der(T_w[1]) = lambda_w*(T_w[3] - 2*T_w[2] + T_w[1])/((L/(N - 1))^2) + A/(N - 1)/2*alpha[1].*(T[1, 1] - T_w[1]) + A/(N - 1)/2*alpha_HTF*(T[2, 1] - T_w[1]);
      V*(csi)/(N - 1)*rho_w*cp_w*der(T_w[2:N - 1]) = lambda_w*(T_w[3:N] - 2*T_w[2:N - 1] + T_w[1:N - 2])/((L/(N - 1))^2) + A/(N - 1)*alpha[2:N - 1].*(T[1, 2:N - 1] - T_w[2:N - 1]) + A/(N - 1)*alpha_HTF*(T[2, 2:N - 1] - T_w[2:N - 1]);
      V*(csi)/(N - 1)/2*rho_w*cp_w*der(T_w[N]) = lambda_w*(T_w[N] - 2*T_w[N - 1] + T_w[N - 2])/((L/(N - 1))^2) + A/(N - 1)/2*alpha[N].*(T[1, N] - T_w[N]) + A/(N - 1)/2*alpha_HTF*(T[2, N] - T_w[N]);
      (V*(1 - csi)/2)/(N - 1)/2*cp_HTF*rho_HTF*der(T[2, 1]) = inletB.mdot*cp_HTF*(T[2, 2] - T[2, 1]) - A/(N - 1)/2*alpha_HTF*(T[2, 1] - T_w[1]);
      (V*(1 - csi)/2)/(N - 1)*cp_HTF*rho_HTF*der(T[2, 2:N - 1]) = inletB.mdot*cp_HTF*(T[2, 3:N] - T[2, 2:N - 1]) - A/(N - 1)*alpha_HTF.*(T[2, 2:N - 1] - T_w[2:N - 1]);
      outletA.p = stateA[N].p;
      outletA.State = stateA[N];
      outletB.T = T[2, 1];
      outletB.p = inletB.p;
      inletA.mdot + outletA.mdot = 0;
      inletB.mdot + outletB.mdot = 0;
// Annotation
      annotation(
        Icon(graphics = {Rectangle(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-75, 50}, {75, -50}}), Line(points = {{-100, 20}, {-75, 20}}, thickness = 1.5), Line(points = {{-100, -20}, {-75, -20}}, thickness = 1.5), Line(points = {{75, 20}, {100, 20}}, thickness = 1.5), Line(points = {{100, -20}, {75, -20}}, thickness = 1.5), Bitmap(extent = {{-25, -25}, {25, 25}}, imageSource = "iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAMAAADVRocKAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAzUExURQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAKMFRskAAAAQdFJOUwAQIDBAUGBwgI+fr7/P3+8jGoKKAAAACXBIWXMAAA7DAAAOwwHHb6hkAAACYUlEQVRoQ+2Yy5qDIAxGBa8IQt7/aScqrXXUmgRYzDeeRRs3hFz+oFQPDw//kLaNRiFUCDqaZbAALppF6CEE6OJDAXSAvoWCSbJg559iScIEKayzL5UkhQma/4slaZwTNFMoSbjxerUUQAG5qQmGaM6+sBiZMeC2Rd0rW/loQmiiiWCSPp6y4MBEa2EAH61M9OD3Wfdry+YCJfBLXE1eMZhjUZexkQvcbpTARlYxuE0CG12+Oh8qvOLP3Eo4VnglW52Hq9lmYYxWEniMXahWA+QIwVzv8z3AU8AhdLlNrM6hfdmctuiLDK16M/p9stq+BjCHEKIlpD3X2EbqK8ZNAMkh3AaQGoK7P1ZwD9ES0NwHkNZIljIvO5iixaamvf54sZxH2sAfpBMJxyhpWOLhKRuq5J2NwnPBU9/etExsLb07JpHYSD260knKjIcJuXRYZko/72FkSDaQGBmS5YilT8Xvo5o3YfjjouOJx7A/F0Ze2Vp2ESxvymv2zPa8AcavMtNBxXYQmNoEiAaVidd3mn30O96XNlM2CLOxmbJBmA76/UUAAaZyLHuc8k5yzuHxgrUn1uERYY141uER4eSI+ga1h/EJeXJXQgBDIIpZh+NdCYXLL/zfjGwRrFCveNcLYQkN6eqvSbiPx73dJhc7SJagBQv+xkPtYRQmaGGCz+vSI5ifKWX9OYZvY7UPn/fBMgYAeyHTGr2b1PWXXgJz4kKZAEHcPzsMurD7pVQz4vJnfkXoEV2AHdpGV0rXbe8CPjvRfLhA9W728Sa4IefyK52xzkPwkzV9emkfHh7+DlX1A+MlFydyV9TxAAAAAElFTkSuQmCC")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end CtrFlow_Condenser_Upwind;

    class CtrFlow_IHX_Upwind
      function alpha_Gnielinski
        input Fluid.ThermodynamicState state "fluid thermodyanmic state";
        input SI.Length d_h "[m] inner channel hydraulic diameter";
        input SI.MassFlowRate m_dot "[kg/s] mass flow rate through single channel";
        output SI.CoefficientOfHeatTransfer alpha "[W/m2K] convective HTC";
      protected
        SI.Density rho "fluid density";
        SI.DynamicViscosity mu "fluid viscosity";
        SI.ThermalConductivity lambda "fluid thermal conductivity";
        Real Pr(unit = "-") "Prandtl number";
        Real Re(unit = "-") "Reynolds number";
        Real Nu(unit = "-") "Nusselt number";
        Real f(unit = "-") "Moody's friction factor, i.e. 4C_f frinction coefficient";
        SI.Area A_c "inner pipe cross-section";
      algorithm
        rho := Fluid.density(state);
        mu := Fluid.dynamicViscosity(state);
        lambda := Fluid.thermalConductivity(state);
        Pr := Fluid.prandtlNumber(state);
        A_c := Modelica.Constants.pi*((d_h)^2)/4;
        Re := (rho*d_h/mu)*(m_dot/(A_c*rho));
        f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
        Nu := ((f/8)*(Re - 1000)*Pr)/(1 + 12.7*((f/8)^0.5)*(Pr^(2/3) - 1));
        alpha := Nu*lambda/d_h;
      end alpha_Gnielinski;

      import alpha_Gnielinski;
      parameter Integer N = 12 "# of HX nodes" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length d_hA = 0.0254 "[m] pipe A hydraulic diameter" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real N_pA(unit = "-") = 4 "# of pipes A - so that u_max < 1 m/s" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length d_hB = 0.0254 "[m] pipe B hydraulic diameter" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real N_pB(unit = "-") = 5 "# of pipes B - so that u_max < 1 m/s" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Area A = 7000 "[m2] heat transfer surface" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length L = 10 "[m] linear dimension heat exchanger" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real csi(unit = "-") = 0.2 "material share" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Density rho_w = 7800 "[kg/m3] Density HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.SpecificHeatCapacity cp_w = 477 "[J/kgK] Density HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.ThermalConductivity lambda_w = 14.9 "[W/mK] Thermal conductivity HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.Pressure Delta_p = 20000 "[Pa] Pressure drop across HX" annotation(
        Dialog(tab = "Parameters", group = "Fluid"));
      // Initialisation
      parameter SI.ThermodynamicTemperature T_0 = 20 + 273.15 "[K] intial temperature" annotation(
        Dialog(tab = "Initialisation", group = "Thermal"));
      // Connections
      NewHTHP.Connections.FluidInlet inletA annotation(
        Placement(visible = true, transformation(origin = {0, 38}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, 38}, extent = {{-115, 15}, {-85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidOutlet outletA annotation(
        Placement(visible = true, transformation(origin = {-200, -38}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {0, -38}, extent = {{-85, 15}, {-115, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidInlet inletB annotation(
        Placement(visible = true, transformation(origin = {100, -100}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 90)));
      NewHTHP.Connections.FluidOutlet outletB annotation(
        Placement(visible = true, transformation(origin = {-100, 100}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{85, -15}, {115, 15}}, rotation = 90)));
      SI.ThermodynamicTemperature T[2, N](each start = T_0);
      SI.ThermodynamicTemperature T_w[N](each start = T_0);
      Fluid.ThermodynamicState stateA[N] "fluid thermodynamic state";
      Fluid.ThermodynamicState stateB[N] "fluid thermodynamic state";
      SI.CoefficientOfHeatTransfer alpha[2, N](each start = 1500);
      SI.SpecificEnthalpy h[2, N](each start = Fluid.specificEnthalpy_pT(100000, T_0));
      SI.Density rho[2, N](each start = Fluid.density_pT(100000, T_0));
    protected
      parameter SI.Volume V = A/450 "[m3] volume";
    equation
      T[1, 1] = inletA.State.T;
      stateA[1] = inletA.State;
      h[1, 1] = stateA[1].h;
      rho[1, 1] = stateA[1].d;
      T[2, N] = inletB.State.T;
      stateB[N] = inletB.State;
      h[2, N] = stateB[N].h;
      rho[2, N] = stateB[N].d;
      alpha[1, 1:N] = alpha_Gnielinski(stateA[1:N], d_hA, inletA.mdot/N_pA);
      alpha[2, 1:N] = alpha_Gnielinski(stateB[1:N], d_hB, inletB.mdot/N_pB);
//System of ODEs
      (V*(1 - csi)/2)/(N - 1)*rho[1, 2:N - 1].*der(h[1, 2:N - 1]) = inletA.mdot*(h[1, 1:N - 2] - h[1, 2:N - 1]) - A/(N - 1)*alpha[1, 2:N - 1].*(T[1, 2:N - 1] - T_w[2:N - 1]);
      (V*(1 - csi)/2)/(N - 1)/2*rho[1, N].*der(h[1, N]) = inletA.mdot*(h[1, N - 1] - h[1, N]) - A/(N - 1)/2*alpha[1, N].*(T[1, N] - T_w[N]);
      for ii in 2:N loop
        stateA[ii] = Fluid.setState_ph(inletA.p - ((Delta_p)/(N - 1))*(ii - 1), h[1, ii]);
      end for;
      T[1, 2:N] = stateA[2:N].T;
      rho[1, 2:N] = stateA[2:N].d;
      V*(csi)/(N - 1)/2*rho_w*cp_w*der(T_w[1]) = lambda_w*(T_w[3] - 2*T_w[2] + T_w[1])/((L/(N - 1))^2) + A/(N - 1)/2*alpha[1, 1].*(T[1, 1] - T_w[1]) + A/(N - 1)/2*alpha[2, 1]*(T[2, 1] - T_w[1]);
      V*(csi)/(N - 1)*rho_w*cp_w*der(T_w[2:N - 1]) = lambda_w*(T_w[3:N] - 2*T_w[2:N - 1] + T_w[1:N - 2])/((L/(N - 1))^2) + A/(N - 1)*alpha[1, 2:N - 1].*(T[1, 2:N - 1] - T_w[2:N - 1]) + A/(N - 1)*alpha[2, 2:N - 1].*(T[2, 2:N - 1] - T_w[2:N - 1]);
      V*(csi)/(N - 1)/2*rho_w*cp_w*der(T_w[N]) = lambda_w*(T_w[N] - 2*T_w[N - 1] + T_w[N - 2])/((L/(N - 1))^2) + A/(N - 1)/2*alpha[1, N].*(T[1, N] - T_w[N]) + A/(N - 1)/2*alpha[2, N].*(T[2, N] - T_w[N]);
      (V*(1 - csi)/2)/(N - 1)/2*rho[2, 1]*der(h[2, 1]) = inletB.mdot*(h[2, 2] - h[2, 1]) - A/(N - 1)/2*alpha[2, 1]*(T[2, 1] - T_w[1]);
      (V*(1 - csi)/2)/(N - 1)*rho[2, 2:N - 1].*der(h[2, 2:N - 1]) = inletB.mdot*(h[2, 3:N] - h[2, 2:N - 1]) - A/(N - 1)*alpha[2, 2:N - 1].*(T[2, 2:N - 1] - T_w[2:N - 1]);
      for ii in 1:N - 1 loop
        stateB[ii] = Fluid.setState_ph(inletB.p - ((Delta_p)/(N - 1))*((N - 1) - (ii - 1)), h[2, ii]);
      end for;
// stateB[1:N - 1] = Fluid.setState_ph(inletB.p, h[2, 1:N - 1]);
      T[2, 1:N - 1] = stateB[1:N - 1].T;
      rho[2, 1:N - 1] = stateB[2:N].d;
      outletA.p = stateA[N].p;
      outletA.State = stateA[N];
      outletB.p = stateB[1].p;
      outletB.State = stateB[1];
      inletA.mdot + outletA.mdot = 0;
      inletB.mdot + outletB.mdot = 0;
// Annotation
      annotation(
        Icon(graphics = {Line(points = {{0, 75}, {0, 100}}, thickness = 1.5), Line(points = {{0, -75}, {0, -100}}, thickness = 1.5), Ellipse(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-75, 75}, {75, -75}}), Line(points = {{20, -37.5}, {-100, -37.5}}, thickness = 1.5), Line(points = {{20, 37.5}, {-100, 37.5}}, thickness = 1.5), Line(points = {{20, 37.5}, {-30, 0}, {20, -37.5}}, thickness = 1.5)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true)),
        Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true)),
        version = "",
        uses);
    end CtrFlow_IHX_Upwind;

    class CtrFlow_IHX_Upwind_Refrigerant
      function alpha_Shah_cond
        input Fluid.ThermodynamicState state "fluid thermodyanmic state";
        input SI.Length d_h "[m] inner channel hydraulic diameter";
        input SI.MassFlowRate m_dot "[kg/s] mass flow rate through single channel";
        output SI.CoefficientOfHeatTransfer alpha "[W/m2K] convective HTC";
      protected
        Refrigerant.SaturationProperties sat "fluid saturated state";
        Refrigerant.ThermodynamicState state_l "saturated liquid thermodynamic state";
        Refrigerant.ThermodynamicState state_v "saturated vapour thermodynamic state";
        Real rho[3](each unit = "kg/m3") "fluid density";
        Real mu[3](each unit = "Pas") "fluid viscosity";
        Real lambda[3](each unit = "W/mK") "fluid thermal conductivity";
        Real Pr[3](each unit = "-") "Prandtl number";
        Real alpha_t[2](each unit = "W/m2K") "temporary convective HTC";
        Real Re(unit = "-") "Reynolds number";
        Real Nu(unit = "-") "Nusselt number";
        Real Z(unit = "-") "A correction parameter";
        Real f(unit = "-") "Moody's friction factor, i.e. 4C_f frinction coefficient";
        Real x(unit = "-") "vapour quality";
        Real A_c(unit = "m2") "inner pipe cross-section";
      algorithm
        sat := Refrigerant.setSat_p(state.p);
        state_l := Refrigerant.setBubbleState(sat, 1);
        state_v := Refrigerant.setDewState(sat, 1);
        rho[1] := Refrigerant.density(state_l);
        rho[2] := Refrigerant.density(state);
        rho[3] := Refrigerant.density(state_v);
        mu[1] := Refrigerant.dynamicViscosity(state_l);
        mu[2] := Refrigerant.dynamicViscosity(state);
        mu[3] := Refrigerant.dynamicViscosity(state_v);
        lambda[1] := Refrigerant.thermalConductivity(state_l);
        lambda[2] := Refrigerant.thermalConductivity(state);
        lambda[3] := Refrigerant.thermalConductivity(state_v);
        Pr[1] := Refrigerant.prandtlNumber(state_l);
        Pr[2] := Refrigerant.prandtlNumber(state);
        Pr[3] := Refrigerant.prandtlNumber(state_v);
        x := (Refrigerant.specificEnthalpy(state) - Refrigerant.specificEnthalpy(state_l))/(Refrigerant.specificEnthalpy(state_v) - Refrigerant.specificEnthalpy(state_l));
        A_c := Modelica.Constants.pi*((d_h)^2)/4;
        if x <= 0 or x >= 1 then
          Re := (rho[2]*d_h/mu[2])*(m_dot/(A_c*rho[2]));
          f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
          Nu := ((f/8)*(Re - 1000)*Pr[2])/(1 + 12.7*((f/8)^0.5)*(Pr[2]^(2/3) - 1));
          alpha := Nu*lambda[2]/d_h;
        else
          Re := (rho[1]*d_h/mu[1])*(m_dot*(1 - x)/(A_c*rho[1]));
          if Re >= 100 then
            if x >= 0.05 then
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := 0.023*(Re^0.8)*(Pr[1]^0.4);
              alpha := Nu*lambda[1]/d_h;
              Z := (((1 - x)/x)^0.8)*(state.p/Refrigerant.getCriticalPressure())^0.4;
              alpha := alpha*(1 + 3.8/(Z^0.95));
            else
              Re := (rho[1]*d_h/mu[1])*(m_dot/(A_c*rho[1]));
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := ((f/8)*(Re - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
              alpha_t[1] := Nu*lambda[1]/d_h;
              Re := (rho[1]*d_h/mu[1])*(m_dot*0.95/(A_c*rho[1]));
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := 0.023*(Re^0.8)*(Pr[1]^0.4);
              alpha_t[2] := Nu*lambda[1]/d_h;
              Z := (((1 - 0.05)/0.05)^0.8)*(state.p/Refrigerant.getCriticalPressure())^0.4;
              alpha_t[2] := alpha_t[2]*(1 + 3.8/(Z^0.95));
              alpha := alpha_t[2] + ((alpha_t[1] - alpha_t[2])/(0 - 0.05))*(x - 0.05);
            end if;
          else
            Re := (rho[3]*d_h/mu[3])*(m_dot/(A_c*rho[3]));
            f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
            Nu := ((f/8)*(Re - 1000)*Pr[3])/(1 + 12.7*((f/8)^0.5)*(Pr[3]^(2/3) - 1));
            alpha_t[1] := Nu*lambda[3]/d_h;
            f := (0.790*Modelica.Math.log(100) - 1.64)^(-2);
            Nu := 0.023*(100^0.8)*(Pr[1]^0.4);
            alpha_t[2] := Nu*lambda[1]/d_h;
            Z := (((1 - (1 - 100*mu[1]*A_c/(d_h*m_dot)))/(1 - 100*mu[1]*A_c/(d_h*m_dot)))^0.8)*(state.p/Fluid.getCriticalPressure())^0.4;
            alpha_t[2] := alpha_t[2]*(1 + 3.8/(Z^0.95));
            alpha := alpha_t[2] + ((alpha_t[1] - alpha_t[2])/(1 - (1 - 100*mu[1]*A_c/(d_h*m_dot))))*(x - (1 - 100*mu[1]*A_c/(d_h*m_dot)));
          end if;
        end if;
      end alpha_Shah_cond;

      function correct_Shah
        input Real Co "[-] convective number";
        input Real Bo "[-] boiling number";
        output Real F_e "[-] correction factor";
      protected
        Real h_b[4](each unit = "[-]") "coefficients for individual boiling regimes";
      algorithm
        h_b[1] := 1.8*Co^(-0.8);
        if Bo >= 0.3e-4 then
          h_b[2] := 230*Bo^0.5;
        else
          h_b[2] := 1 + 46*Bo^0.5;
        end if;
        if Bo >= 11e-4 then
          h_b[3] := 14.70*(Bo^0.5)*Modelica.Math.exp(2.74*(Co)^(-0.1));
        else
          h_b[3] := 15.43*(Bo^0.5)*Modelica.Math.exp(2.74*(Co)^(-0.1));
        end if;
        if Bo >= 11e-4 then
          h_b[4] := 14.70*(Bo^0.5)*Modelica.Math.exp(2.47*(Co)^(-0.15));
        else
          h_b[4] := 15.43*(Bo^0.5)*Modelica.Math.exp(2.47*(Co)^(-0.15));
        end if;
        if Co <= 0.1 then
          F_e := max(h_b[1], h_b[4]);
        elseif Co > 1 then
          F_e := max(h_b[1], h_b[2]);
        else
          F_e := max(h_b[1], h_b[3]);
        end if;
      end correct_Shah;

      function alpha_Shah_eva
        import correct_Shah;
        input Fluid.ThermodynamicState state "fluid thermodyanmic state";
        input Real d_h "[m] inner channel hydraulic diameter";
        input Real m_dot "[kg/s] mass flow rate through single channel";
        output Real alpha "[W/m2K] convective HTC";
        //output Real x_out;
        //output Real hh[3];
      protected
        Fluid.SaturationProperties sat "fluid saturated state";
        Fluid.ThermodynamicState state_l "saturated liquid thermodynamic state";
        Fluid.ThermodynamicState state_v "saturated vapour thermodynamic state";
        Real rho[3](each unit = "kg/m3") "fluid density";
        Real mu[3](each unit = "Pas") "fluid viscosity";
        Real lambda[3](each unit = "W/mK") "fluid thermal conductivity";
        Real Pr[3](each unit = "-") "Prandtl number";
        Real alpha_t[2](each unit = "W/m2K") "temporary convective HTC";
        Real Re(unit = "-") "Reynolds number";
        Real Nu(unit = "-") "Nusselt number";
        Real Co(unit = "-") "Convection number";
        Real Bo(unit = "-") "Boiling number";
        Real F_e(unit = "-") "Correction factor";
        Real f(unit = "-") "Moody's friction factor, i.e. 4C_f frinction coefficient";
        Real x(unit = "-") "vapour quality";
        Real A_c(unit = "m2") "inner pipe cross-section";
      algorithm
        sat := Fluid.setSat_p(state.p);
        state_l := Fluid.setBubbleState(sat, 1);
        state_v := Fluid.setDewState(sat, 1);
        rho[1] := Fluid.density(state_l);
        rho[2] := Fluid.density(state);
        rho[3] := Fluid.density(state_v);
        mu[1] := Fluid.dynamicViscosity(state_l);
        mu[2] := Fluid.dynamicViscosity(state);
        mu[3] := Fluid.dynamicViscosity(state_v);
        lambda[1] := Fluid.thermalConductivity(state_l);
        lambda[2] := Fluid.thermalConductivity(state);
        lambda[3] := Fluid.thermalConductivity(state_v);
        Pr[1] := Fluid.prandtlNumber(state_l);
        Pr[2] := Fluid.prandtlNumber(state);
        Pr[3] := Fluid.prandtlNumber(state_v);
        x := (Fluid.specificEnthalpy(state) - Fluid.specificEnthalpy(state_l))/(Fluid.specificEnthalpy(state_v) - Fluid.specificEnthalpy(state_l));
        A_c := Modelica.Constants.pi*((d_h)^2)/4;
        if x <= 0 or x >= 1 then
          Re := (rho[2]*d_h/mu[2])*(m_dot/(A_c*rho[2]));
          f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
          Nu := ((f/8)*(Re - 1000)*Pr[2])/(1 + 12.7*((f/8)^0.5)*(Pr[2]^(2/3) - 1));
          alpha := Nu*lambda[2]/d_h;
        else
          Re := (rho[1]*d_h/mu[1])*(m_dot*(1 - x)/(A_c*rho[1]));
          if Re >= 2300 then
            if x >= 0.05 then
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := ((f/8)*(Re - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
              alpha := Nu*lambda[1]/d_h;
              Co := (((1 - x)/x)^0.8)*(rho[3]/rho[1])^0.5;
              Bo := (((34*(state.p/1e5)^(0.18))/(1 - 0.0045*(state.p/1e5)))^3)*A_c/(m_dot*(Fluid.specificEnthalpy(state_v) - Fluid.specificEnthalpy(state_l)));
              F_e := correct_Shah(Co, Bo);
              alpha := F_e*alpha;
            else
              Re := (rho[1]*d_h/mu[1])*(m_dot/(A_c*rho[1]));
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := ((f/8)*(Re - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
              alpha_t[1] := Nu*lambda[1]/d_h;
              Re := (rho[1]*d_h/mu[1])*(m_dot*0.95/(A_c*rho[1]));
              f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
              Nu := ((f/8)*(Re - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
              alpha_t[2] := Nu*lambda[1]/d_h;
              Co := (((1 - 0.05)/0.05)^0.8)*(rho[3]/rho[1])^0.5;
              Bo := (((34*(state.p/1e5)^(0.18))/(1 - 0.0045*(state.p/1e5)))^3)*A_c/(m_dot*(Fluid.specificEnthalpy(state_v) - Fluid.specificEnthalpy(state_l)));
              F_e := correct_Shah(Co, Bo);
              alpha_t[2] := F_e*alpha_t[2];
              alpha := alpha_t[2] + ((alpha_t[1] - alpha_t[2])/(0 - 0.05))*(x - 0.05);
            end if;
          else
            Re := (rho[3]*d_h/mu[3])*(m_dot/(A_c*rho[3]));
            f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
            Nu := ((f/8)*(Re - 1000)*Pr[3])/(1 + 12.7*((f/8)^0.5)*(Pr[3]^(2/3) - 1));
            alpha_t[1] := Nu*lambda[3]/d_h;
            f := (0.790*Modelica.Math.log(2300) - 1.64)^(-2);
            Nu := ((f/8)*(2300 - 1000)*Pr[1])/(1 + 12.7*((f/8)^0.5)*(Pr[1]^(2/3) - 1));
            alpha_t[2] := Nu*lambda[1]/d_h;
            Co := (((1 - (1 - 2300*mu[1]*A_c/(d_h*m_dot)))/(1 - 2300*mu[1]*A_c/(d_h*m_dot)))^0.8)*(rho[3]/rho[1])^0.5;
            Bo := (((34*(state.p/1e5)^(0.18))/(1 - 0.0045*(state.p/1e5)))^3)*A_c/(m_dot*(Fluid.specificEnthalpy(state_v) - Fluid.specificEnthalpy(state_l)));
            F_e := correct_Shah(Co, Bo);
            alpha_t[2] := F_e*alpha_t[2];
            alpha := alpha_t[2] + ((alpha_t[1] - alpha_t[2])/(1 - (1 - 2300*mu[1]*A_c/(d_h*m_dot))))*(x - (1 - 2300*mu[1]*A_c/(d_h*m_dot)));
          end if;
        end if;
//x_out := x;
//hh[1] := Re;
//hh[2] := f;
//hh[3] := Pr[2];
      end alpha_Shah_eva;

      //  function alpha_Gnielinski_A
      //    input Refrigerant.ThermodynamicState state "fluid thermodyanmic state";
      //    input SI.Length d_h "[m] inner channel hydraulic diameter";
      //    input SI.MassFlowRate m_dot "[kg/s] mass flow rate through single channel";
      //    output SI.CoefficientOfHeatTransfer alpha "[W/m2K] convective HTC";
      //  protected
      //    SI.Density rho "fluid density";
      //    SI.DynamicViscosity mu "fluid viscosity";
      //    SI.ThermalConductivity lambda "fluid thermal conductivity";
      //    Real Pr(unit = "-") "Prandtl number";
      //    Real Re(unit = "-") "Reynolds number";
      //    Real Nu(unit = "-") "Nusselt number";
      //    Real f(unit = "-") "Moody's friction factor, i.e. 4C_f frinction coefficient";
      //    SI.Area A_c "inner pipe cross-section";
      //  algorithm
      //    rho := Refrigerant.density(state);
      //    mu := Refrigerant.dynamicViscosity(state);
      //    lambda := Refrigerant.thermalConductivity(state);
      //    Pr := Refrigerant.prandtlNumber(state);
      //    A_c := Modelica.Constants.pi*((d_h)^2)/4;
      //   Re := (rho*d_h/mu)*(m_dot/(A_c*rho));
      //    f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
      //    Nu := ((f/8)*(Re - 1000)*Pr)/(1 + 12.7*((f/8)^0.5)*(Pr^(2/3) - 1));
      //    alpha := Nu*lambda/d_h;
      //  end alpha_Gnielinski_A;
      //  function alpha_Gnielinski_B
      //    input Fluid.ThermodynamicState state "fluid thermodyanmic state";
      //    input SI.Length d_h "[m] inner channel hydraulic diameter";
      //    input SI.MassFlowRate m_dot "[kg/s] mass flow rate through single channel";
      //    output SI.CoefficientOfHeatTransfer alpha "[W/m2K] convective HTC";
      //  protected
      //    SI.Density rho "fluid density";
      //    SI.DynamicViscosity mu "fluid viscosity";
      //    SI.ThermalConductivity lambda "fluid thermal conductivity";
      //    Real Pr(unit = "-") "Prandtl number";
      //    Real Re(unit = "-") "Reynolds number";
      //   Real Nu(unit = "-") "Nusselt number";
      //    Real f(unit = "-") "Moody's friction factor, i.e. 4C_f frinction coefficient";
      //    SI.Area A_c "inner pipe cross-section";
      //  algorithm
      //    rho := Fluid.density(state);
      //    mu := Fluid.dynamicViscosity(state);
      //    lambda := Fluid.thermalConductivity(state);
      //    Pr := Fluid.prandtlNumber(state);
      //    A_c := Modelica.Constants.pi*((d_h)^2)/4;
      //   Re := (rho*d_h/mu)*(m_dot/(A_c*rho));
      //    f := (0.790*Modelica.Math.log(Re) - 1.64)^(-2);
      //    Nu := ((f/8)*(Re - 1000)*Pr)/(1 + 12.7*((f/8)^0.5)*(Pr^(2/3) - 1));
      //    alpha := Nu*lambda/d_h;
      //  end alpha_Gnielinski_B;
      //  import alpha_Gnielinski_A;
      //  import alpha_Gnielinski_B;
      import alpha_Shah_cond;
      import alpha_Shah_eva;
      import correct_Shah;
      parameter Integer N = 12 "# of HX nodes" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length d_hA = 0.0254 "[m] pipe A hydraulic diameter" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real N_pA(unit = "-") = 4 "# of pipes A - so that u_max < 1 m/s" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length d_hB = 0.0254 "[m] pipe B hydraulic diameter" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real N_pB(unit = "-") = 5 "# of pipes B - so that u_max < 1 m/s" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Area A = 7000 "[m2] heat transfer surface" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Length L = 10 "[m] linear dimension heat exchanger" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter Real csi(unit = "-") = 0.2 "material share" annotation(
        Dialog(tab = "Parameters", group = "Geometry"));
      parameter SI.Density rho_w = 7800 "[kg/m3] Density HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.SpecificHeatCapacity cp_w = 477 "[J/kgK] Density HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.ThermalConductivity lambda_w = 14.9 "[W/mK] Thermal conductivity HX material" annotation(
        Dialog(tab = "Parameters", group = "Thermal"));
      parameter SI.Pressure Delta_p = 20000 "[Pa] Pressure drop across HX" annotation(
        Dialog(tab = "Parameters", group = "Fluid"));
      // Initialisation
      parameter SI.ThermodynamicTemperature T_0 = 20 + 273.15 "[K] intial temperature" annotation(
        Dialog(tab = "Initialisation", group = "Thermal"));
      // Connections
      NewHTHP.Connections.FluidInlet_refrigerant inletA annotation(
        Placement(visible = true, transformation(origin = {0, 38}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, 38}, extent = {{-115, 15}, {-85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidOutlet_refrigerant outletA annotation(
        Placement(visible = true, transformation(origin = {-200, -38}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {0, -38}, extent = {{-85, 15}, {-115, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidInlet inletB annotation(
        Placement(visible = true, transformation(origin = {100, -100}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 90)));
      NewHTHP.Connections.FluidOutlet outletB annotation(
        Placement(visible = true, transformation(origin = {-100, 100}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{85, -15}, {115, 15}}, rotation = 90)));
      SI.ThermodynamicTemperature T[2, N](each start = T_0);
      SI.ThermodynamicTemperature T_w[N](each start = T_0);
      Refrigerant.ThermodynamicState stateA[N] "fluid thermodynamic state";
      Fluid.ThermodynamicState stateB[N] "fluid thermodynamic state";
      SI.CoefficientOfHeatTransfer alpha[2, N](each start = 1500);
      //SI.SpecificEnthalpy h[2, N];
      //SI.Density rho[2, N];
      SI.SpecificEnthalpy h[2, N](each start = Refrigerant.specificEnthalpy_pT(1713000, T_0));
      SI.Density rho[2, N](each start = Refrigerant.density_pT(1713000, T_0));
    protected
      parameter SI.Volume V = A/450 "[m3] volume";
      //initial equation
      //h[1, 1] = Refrigerant.specificEnthalpy_pT(1713000, T_0);
      //h[2, N] = Fluid.specificEnthalpy_pT(400000, T_0);
      //rho[1, 1] = Refrigerant.density_pT(1713000, T_0);
      //rho[2, N] = Fluid.density_pT(400000, T_0);
    equation
      T[1, 1] = inletA.State.T;
      stateA[1] = inletA.State;
      h[1, 1] = stateA[1].h;
      rho[1, 1] = stateA[1].d;
      T[2, N] = inletB.State.T;
      stateB[N] = inletB.State;
      h[2, N] = stateB[N].h;
      rho[2, N] = stateB[N].d;
      alpha[1, 1:N] = alpha_Shah_cond(stateA[1:N], d_hA, inletA.mdot/N_pA);
      alpha[2, 1:N] = alpha_Shah_eva(stateB[1:N], d_hB, inletB.mdot/N_pB);
//System of ODEs
      (V*(1 - csi)/2)/(N - 1)*rho[1, 2:N - 1].*der(h[1, 2:N - 1]) = inletA.mdot*(h[1, 1:N - 2] - h[1, 2:N - 1]) - A/(N - 1)*alpha[1, 2:N - 1].*(T[1, 2:N - 1] - T_w[2:N - 1]);
      (V*(1 - csi)/2)/(N - 1)/2*rho[1, N].*der(h[1, N]) = inletA.mdot*(h[1, N - 1] - h[1, N]) - A/(N - 1)/2*alpha[1, N].*(T[1, N] - T_w[N]);
      for ii in 2:N loop
        stateA[ii] = Refrigerant.setState_ph(inletA.p - ((Delta_p)/(N - 1))*(ii - 1), h[1, ii]);
      end for;
      T[1, 2:N] = stateA[2:N].T;
      rho[1, 2:N] = stateA[2:N].d;
      V*(csi)/(N - 1)/2*rho_w*cp_w*der(T_w[1]) = lambda_w*(T_w[3] - 2*T_w[2] + T_w[1])/((L/(N - 1))^2) + A/(N - 1)/2*alpha[1, 1].*(T[1, 1] - T_w[1]) + A/(N - 1)/2*alpha[2, 1]*(T[2, 1] - T_w[1]);
      V*(csi)/(N - 1)*rho_w*cp_w*der(T_w[2:N - 1]) = lambda_w*(T_w[3:N] - 2*T_w[2:N - 1] + T_w[1:N - 2])/((L/(N - 1))^2) + A/(N - 1)*alpha[1, 2:N - 1].*(T[1, 2:N - 1] - T_w[2:N - 1]) + A/(N - 1)*alpha[2, 2:N - 1].*(T[2, 2:N - 1] - T_w[2:N - 1]);
      V*(csi)/(N - 1)/2*rho_w*cp_w*der(T_w[N]) = lambda_w*(T_w[N] - 2*T_w[N - 1] + T_w[N - 2])/((L/(N - 1))^2) + A/(N - 1)/2*alpha[1, N].*(T[1, N] - T_w[N]) + A/(N - 1)/2*alpha[2, N].*(T[2, N] - T_w[N]);
      (V*(1 - csi)/2)/(N - 1)/2*rho[2, 1]*der(h[2, 1]) = inletB.mdot*(h[2, 2] - h[2, 1]) - A/(N - 1)/2*alpha[2, 1]*(T[2, 1] - T_w[1]);
      (V*(1 - csi)/2)/(N - 1)*rho[2, 2:N - 1].*der(h[2, 2:N - 1]) = inletB.mdot*(h[2, 3:N] - h[2, 2:N - 1]) - A/(N - 1)*alpha[2, 2:N - 1].*(T[2, 2:N - 1] - T_w[2:N - 1]);
      for ii in 1:N - 1 loop
        stateB[ii] = Fluid.setState_ph(inletB.p, h[2, ii]);
      end for;
// stateB[1:N - 1] = Fluid.setState_ph(inletB.p, h[2, 1:N - 1]);
      T[2, 1:N - 1] = stateB[1:N - 1].T;
      rho[2, 1:N - 1] = stateB[2:N].d;
      outletA.p = stateA[N].p;
      outletA.State = stateA[N];
      outletB.p = stateB[1].p;
      outletB.State = stateB[1];
      inletA.mdot + outletA.mdot = 0;
      inletB.mdot + outletB.mdot = 0;
// Annotation
      annotation(
        Icon(graphics = {Line(points = {{0, 75}, {0, 100}}, thickness = 1.5), Line(points = {{0, -75}, {0, -100}}, thickness = 1.5), Ellipse(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-75, 75}, {75, -75}}), Line(points = {{20, -37.5}, {-100, -37.5}}, thickness = 1.5), Line(points = {{20, 37.5}, {-100, 37.5}}, thickness = 1.5), Line(points = {{20, 37.5}, {-30, 0}, {20, -37.5}}, thickness = 1.5)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true)),
        Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true)),
        version = "",
        uses);
    end CtrFlow_IHX_Upwind_Refrigerant;
    annotation(
      Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1.5, extent = {{-75, 50}, {75, -50}}), Line(points = {{-100, 20}, {-75, 20}}, thickness = 1.5), Line(points = {{-100, -20}, {-75, -20}}, thickness = 1.5), Line(points = {{75, 20}, {100, 20}}, thickness = 1.5), Line(points = {{100, -20}, {75, -20}}, thickness = 1.5)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
  end HeatExchangers;

  package Valves
    extends Modelica.Icons.Package;

    class JTValve
      function flow_coefficient_map
        /* Characteristic map of flow coeffcient to % opening
                                                To be made selectable from component - just like Compressor model */
        input Real valve_open "[-] valve opening, from 0 to 1";
        output Real k_V "Flow factor, i.e. flow [m3/h] of H2O [1 kg/m3] over unit pressure drop [bar]";
      protected
        parameter Real[3] coeff = {868.6, 0, 0} "Characteristic valve coefficients";
      algorithm
        k_V := coeff[1] + coeff[2]*valve_open + coeff[3]*valve_open^2 "Generic 2nd-order equation";
      end flow_coefficient_map;

      import flow_coefficient_map;
      // Input/Output
      parameter Boolean enableInput = false "Output connection for selectable quantity" annotation(
        Dialog(tab = "Input/Output", group = "Connections"));
      parameter Boolean enableOutput = true "Output connection for selectable quantity" annotation(
        Dialog(tab = "Input/Output", group = "Connections"));
      // Regularisation parameters
      parameter SI.Density rho_min = 1e-10 "Minimum allowed density" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      parameter SI.Pressure p_min = 100 "Minimum allowed pressure" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      // Connections
      NewHTHP.Connections.FluidInlet inlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {44, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidOutlet outlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {-44, 0}, extent = {{85, 15}, {115, -15}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput output_val(unit = "0 to 1") if enableOutput "Measured value [variable]" annotation(
        Placement(visible = true, transformation(origin = {-80, -100}, extent = {{-20, -20}, {20, 20}}, rotation = 270), iconTransformation(origin = {40, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 270)));
      Modelica.Blocks.Interfaces.RealInput input_val(unit = "0 to 1") if enableInput "Measured value [variable]" annotation(
        Placement(visible = true, transformation(origin = {80, 100}, extent = {{-20, -20}, {20, 20}}, rotation = 270), iconTransformation(origin = {-40, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));
      parameter Real valve_open(min = 0, max = 1) = 0.5 "[-] valve opening, from 0 to 1";
      Real k_V;
      Real Dp;
      //Fluid.SaturationProperties sat "fluid saturated state";
      //Fluid.ThermodynamicState state_l "saturated liquid thermodynamic state";
      //Fluid.ThermodynamicState state_v "saturated vapour thermodynamic state";
      //Real x(unit = "-") "vapour quality";
      // Guess values
      SI.Density rho_l(start = 527);
    equation
      inlet.mdot + outlet.mdot = 0;
      k_V = flow_coefficient_map(valve_open);
      rho_l = max(Fluid.density(Fluid.setBubbleState(Fluid.setSat_p(inlet.State.p), 1)), rho_min);
      Dp = (((inlet.mdot*3600)^2)/(rho_l))/(k_V^2);
      outlet.p = max(inlet.p - 1e5*Dp, p_min);
      outlet.State = Fluid.setState_ph(outlet.p, inlet.State.h);
      if enableOutput then
        output_val = valve_open;
      elseif enableInput then
        input_val = valve_open;
      end if;
// Annotation
      annotation(
        Icon(graphics = {Line(points = {{-40, 0}, {-50, 0}}, thickness = 1.5), Line(points = {{40, 0}, {50, 0}}, thickness = 1.5), Polygon(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, points = {{0, 0}, {-40, 20}, {-40, -20}, {0, 0}}), Polygon(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, points = {{0, 0}, {40, 20}, {40, -20}, {0, 0}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end JTValve;

    package Test "Simple examples to test the models in this package"
      extends Modelica.Icons.ExamplesPackage;

      model JTValve_test
        NewHTHP.Valves.JTValve jTValve(enableInput = false, enableOutput = true, valve_open = 0.5) annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        NewHTHP.Boundaries.Source source(T0_par = 108.4 + 273.15, enthalpyFromInput = false, p0_par(displayUnit = "Pa") = 1158000, pressureFromInput = false, setEnthalpy = false, temperatureFromInput = false) annotation(
          Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        NewHTHP.Boundaries.Sink sink(p0_par(displayUnit = "Pa") = 127, pressureFromInput = false) annotation(
          Placement(visible = true, transformation(origin = {0, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      equation
        connect(source.outlet, jTValve.inlet) annotation(
          Line(points = {{0, -32}, {0, -6}}, color = {85, 0, 255}, thickness = 1));
        connect(sink.inlet, jTValve.outlet) annotation(
          Line(points = {{0, 28}, {0, 6}}, color = {85, 0, 255}, thickness = 1));
      end JTValve_test;
    equation

    end Test;
    annotation(
      Icon(graphics = {Polygon(fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1.5, points = {{0, 0}, {-60, -30}, {-60, 30}, {0, 0}}), Polygon(fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1.5, points = {{0, 0}, {-60, -30}, {-60, 30}, {0, 0}}), Polygon(rotation = 180, fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1.5, points = {{0, 0}, {-60, -30}, {-60, 30}, {0, 0}}), Line(points = {{0, 0}, {0, 40}}, thickness = 1.5), Ellipse(origin = {0, 40}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1.5, extent = {{22, 10}, {-22, -10}})}));
  end Valves;

  package Volumes
    extends Modelica.Icons.Package;

    class MixedVolume
      // Input/Output
      parameter Boolean enableInlet = true "Include inlet" annotation(
        Dialog(tab = "Input/Output", group = "From Input"));
      parameter Boolean enableOutlet = true "Include outlet" annotation(
        Dialog(tab = "Input/Output", group = "From Input"));
      parameter Boolean enableHeatPort = false "Include access heatport" annotation(
        Dialog(tab = "Input/Output", group = "From Input"));
      // Initialisation
      parameter Boolean initPressure = true "If true pressure is initialised" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state"));
      parameter SI.Pressure p_0 = 136700 "Initial Pressure" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = initPressure));
      parameter Boolean initEnergy = true "If true specific internal energy is initialised" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state"));
      parameter Boolean use_h0 = false "If true initial specific enthalpy is supplied instead of temperature" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = initEnergy));
      parameter SI.ThermodynamicTemperature T_0 = 25 + 273.15 "Initial Temperature" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = initEnergy and (not use_h0)));
      parameter SI.SpecificEnthalpy h_0 = 500000 "Initial specific enthalpy" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = initEnergy and use_h0));
      // Model parameters
      parameter Real k_volume_damping(unit = "-", min = 0) = 0.5 "Damping factor multiplicator" annotation(
        Dialog(tab = "General", group = "Parameters"));
      parameter SI.Volume V = 1 "[m3] Vessel volume" annotation(
        Dialog(tab = "General", group = "Parameters"));
      // Regularisation parameters
      parameter SI.Density rho_min = 1e-3 "Minimum allowed density" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      parameter SI.Pressure p_min = 5000 "Minimum allowed pressure" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      // Advanced
      parameter Boolean usePreferredMediumStates = false "Use medium states instead of the ones differentiated in this component" annotation(
        Dialog(tab = "Advanced"));
      // Connections
      Connections.FluidInlet inlet if enableInlet annotation(
        Placement(transformation(extent = {{-115, 15}, {-85, -15}})));
      Connections.FluidOutlet outlet if enableOutlet annotation(
        Placement(transformation(extent = {{85, 15}, {115, -15}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatport(Q_flow = Q_t, T = T) if enableHeatPort annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -90}, {10, -70}}, rotation = 0), iconTransformation(origin = {0, -10}, extent = {{-10, -90}, {10, -70}}, rotation = 0)));
      SI.Power Q_t;
      SI.ThermodynamicTemperature T(start = T_0);
      SI.Pressure p(start = p_0);
      SI.SpecificEnthalpy h(start = Fluid.specificEnthalpy(Fluid.setState_pT(p_0, T_0)));
      SI.SpecificEnthalpy h_in = if noEvent(enableInlet) then inlet.State.h else 0;
      SI.SpecificEnthalpy h_out = if noEvent(enableOutlet) then h else 0;
      SI.Mass M(start = V*Fluid.density(Fluid.setState_pT(p_0, T_0)));
      SI.MassFlowRate mdot_in(min = 0.0) = if noEvent(enableInlet) then inlet.mdot else 0;
      SI.MassFlowRate mdot_out(max = 0.0) = if noEvent(enableOutlet) then outlet.mdot else 0;
      SI.Energy U(start = V*Fluid.density(Fluid.setState_pT(p_0, T_0))*Fluid.specificInternalEnergy(Fluid.setState_pT(p_0, T_0)));
    initial equation
      if initPressure then
        p = p_0;
      end if;
      if initEnergy then
        if use_h0 then
          h = h_0;
        else
          T = T_0;
        end if;
      end if;
    equation
      der(M) = mdot_in + mdot_out;
      U = M*Fluid.specificInternalEnergy(Fluid.setState_ph(p, h));
      der(U) = Q_t + mdot_in*h_in + mdot_out*h_out;
      T = Fluid.temperature(Fluid.setState_ph(p, h));
      Q_t = 0;
      M = V*Fluid.density(Fluid.setState_ph(p, h));
//p*V = M/(0.07215)*Modelica.Constants.R*T;
/*
      der(M) = mdot_in + mdot_out;
      der(M*Fluid.specificInternalEnergy(Fluid.setState_ph(max(p,p_min),h))) = Q_t + mdot_in*h_in + mdot_out*h_out;
      T = Fluid.temperature(Fluid.setState_ph(max(p,p_min),h));
      
      Q_t = 0;
      
      M = V*max(Fluid.density(Fluid.setState_ph(max(p,p_min),h)),rho_min);
      //p*V = M/(0.07215)*Modelica.Constants.R*T;
      */
  if enableInlet then
// inlet.State.p = p;
      end if;
      if enableOutlet then
        outlet.State = Fluid.setState_ph(p, h);
        outlet.p = p;
//outlet.mdot = -(0.01 + (10/1e5)*p);
      end if;
      annotation(
        Icon(graphics = {Polygon(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, points = {{-60, 60}, {0, 80}, {60, 60}, {60, 40}, {60, -40}, {60, -60}, {0, -80}, {-60, -60}, {-60, -40}, {-60, 40}, {-60, 60}}, smooth = Smooth.Bezier), Line(points = {{60, 0}, {100, 0}}, thickness = 1.5), Line(points = {{-60, 0}, {-100, 0}}, thickness = 1.5)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end MixedVolume;
    
    class SteamAccumulator_pRulez
      // Input/Output
      parameter Boolean enableInlet = true "Include inlet" annotation(
        Dialog(tab = "Input/Output", group = "From Input"));
      parameter Boolean enableOutlet = true "Include outlet" annotation(
        Dialog(tab = "Input/Output", group = "From Input"));
      parameter Boolean enableHeatPort = false "Include access heatport" annotation(
        Dialog(tab = "Input/Output", group = "From Input"));
      // Initialisation
      parameter Boolean initPressure = true "If true pressure is initialised" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state"));
      parameter SI.Pressure p_0 = 400000 "Initial Pressure" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = initPressure));
      parameter Boolean initEnergy = true "If true specific internal energy is initialised" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state"));
      parameter Boolean use_h0 = false "If true initial specific enthalpy is supplied instead of temperature" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = initEnergy));
      parameter SI.ThermodynamicTemperature T_0 = 25 + 273.15 "Initial Temperature" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = initEnergy and (not use_h0)));
      parameter SI.SpecificEnthalpy h_0 = 500000 "Initial specific enthalpy" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state", enable = initEnergy and use_h0));
      parameter Real phi_0 = 0.6 "Initial level" annotation(
        Dialog(tab = "Initialisation", group = "Thermodynamic state"));
      // Model parameters
      parameter Real k_volume_damping(unit = "-", min = 0) = 0.5 "Damping factor multiplicator" annotation(
        Dialog(tab = "General", group = "Parameters"));
      parameter SI.Volume V = 1 "[m3] Vessel volume" annotation(
        Dialog(tab = "General", group = "Parameters"));
      parameter Real k_flow_in (unit = "kg/sPa", min = 0) = 0 annotation(
        Dialog(tab = "General", group = "Parameters"));
      parameter Real k_flow_out (unit = "kg/sPa", min = 0) = 0 annotation(
        Dialog(tab = "General", group = "Parameters"));
      // Regularisation parameters
      parameter SI.Density rho_min = 1e-3 "Minimum allowed density" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      parameter SI.Pressure p_min = 5000 "Minimum allowed pressure" annotation(
        Dialog(tab = "Regularisation", group = "Parameters"));
      // Advanced
      parameter Boolean usePreferredMediumStates = false "Use medium states instead of the ones differentiated in this component" annotation(
        Dialog(tab = "Advanced"));
      // Connections
      NewHTHP.Connections.FluidInlet inlet if enableInlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = -90)));
      Connections.FluidOutlet outlet if enableOutlet annotation(
        Placement(transformation(extent = {{85, 15}, {115, -15}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatport(Q_flow = Q_t, T = T) if enableHeatPort annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -90}, {10, -70}}, rotation = 0), iconTransformation(origin = {0, -10}, extent = {{-10, -90}, {10, -70}}, rotation = 0)));
      SI.Power Q_t;
      SI.ThermodynamicTemperature T;
      SI.Pressure p(start = p_0);
      SI.SpecificEnthalpy h;
      SI.Density rho "[kgm-3] Density";
      Fluid.SaturationProperties sat "fluid saturated state";
      Real phi(start = phi_0, fixed = true) "[-] Fluid level inside the tank";
      Fluid.ThermodynamicState state "thermodynamic state";
      SI.SpecificEnthalpy h_in = if noEvent(enableInlet) then inlet.State.h else 0;
      SI.SpecificEnthalpy h_out = if noEvent(enableOutlet) then state_out.h else 0;
      Fluid.ThermodynamicState state_out "Thermodynamic state for steady mass flow rate";
      SI.Mass M;
      SI.MassFlowRate mdot_in(min = 0.0) = if noEvent(enableInlet) then inlet.mdot else 0;
      SI.MassFlowRate mdot_out(max = 0.0) = if noEvent(enableOutlet) then outlet.mdot else 0;
    initial equation
      p = p_0;
      phi = phi_0;
      if initPressure then
        p = p_0;
      end if;
      if initEnergy then
        if use_h0 then
          h = h_0;
        else
          T = T_0;
        end if;
      end if;
    equation
// Mass conservation
      if enableOutlet then
        -mdot_out = (if time <= 500 then 0 else if time > 750 then k_flow_out*(p - outlet.p) else k_flow_out*(p - outlet.p)*(time-500)/250);
        outlet.State = Fluid.setState_ph(outlet.p, state_out.h);
      end if;
      
      if enableInlet then
        mdot_in = (if time <= 500 then 0 else if time > 510 then k_flow_in*(inlet.p - p) else k_flow_in*(inlet.p - p)*(time-500)/10);
      end if;
      
      
      sat = Fluid.setSat_p(p);
      rho = phi*sat.dl + (1 - phi)*sat.dv;
      V*(der(phi)*(sat.dl - sat.dv) + der(p)*(sat.ddldp*phi + (1 - phi)*sat.ddvdp)) = mdot_in + mdot_out;
      M = rho*V;
      if phi >= 0.9 then
        state_out = Fluid.setBubbleState(sat);
      else
        state_out = Fluid.setDewState(sat);
      end if;
// Energy conservation
      Q_t = 0;
      rho*h = phi*sat.dl*sat.hl + (1 - phi)*sat.dv*sat.hv;
      V*(der(phi)*(sat.dl*sat.hl - sat.dv*sat.hv) + (phi*(sat.hl*sat.ddldp + sat.dl*sat.dhldp) + (1 - phi)*(sat.hv*sat.ddvdp + sat.dv*sat.dhvdp))*der(p)) = mdot_in*h_in + mdot_out*h_out + V*der(p) + Q_t;
      T = Fluid.temperature(Fluid.setState_ph(p, h));
      state = Fluid.setState_ph(p, h);
      annotation(
        Icon(graphics = {Polygon(rotation = -90, fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, points = {{-60, 60}, {0, 80}, {60, 60}, {60, 40}, {60, -40}, {60, -60}, {0, -80}, {-60, -60}, {-60, -40}, {-60, 40}, {-60, 60}}, smooth = Smooth.Bezier), Line(points = {{77, 0}, {100, 0}}, thickness = 1.5), Line(points = {{0, 60}, {0, 100}}, thickness = 1.5), Line(points = {{-75, 0}, {75, 0}}, thickness = 1.5)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end SteamAccumulator_pRulez;

    package Test "Simple examples to test the models in this package"
      extends Modelica.Icons.ExamplesPackage;

      model MixedVolume_test
        Modelica.Blocks.Sources.Step step1(height = 10, offset = 6.9, startTime = 50) annotation(
          Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        NewHTHP.Boundaries.Source source(T0_par = 70.11 + 273.15, mdotFromInput = true, p0_par(displayUnit = "Pa") = 116700, pressureFromInput = false, setEnthalpy = false) annotation(
          Placement(visible = true, transformation(origin = {-10, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        MixedVolume mixedVolume(T_0 = 70 + 273.15, V = 1000, enableInlet = true, enableOutlet = true, p_0(displayUnit = "Pa") = 116700, rho_min(displayUnit = "kg/m3") = 0.001) annotation(
          Placement(visible = true, transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Boundaries.Sink sink(p0_par(displayUnit = "Pa") = 116700, pressureFromInput = false) annotation(
          Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(step1.y, source.mdot0_var) annotation(
          Line(points = {{-38, 0}, {-20, 0}}, color = {0, 0, 127}));
        connect(source.outlet, mixedVolume.inlet) annotation(
          Line(points = {{-2, 0}, {20, 0}}, color = {85, 0, 255}, thickness = 1));
        connect(mixedVolume.outlet, sink.inlet) annotation(
          Line(points = {{40, 0}, {58, 0}}, color = {85, 0, 255}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
      end MixedVolume_test;

      model Loop_withMixedVolume_test
        Modelica.Blocks.Sources.Step step1(height = 10, offset = 25, startTime = 75) annotation(
          Placement(visible = true, transformation(origin = {76, -80}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Step step(height = 30000*(2*Modelica.Constants.pi/60), offset = 15000*(2*Modelica.Constants.pi/60), startTime = 50) annotation(
          Placement(visible = true, transformation(origin = {86, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Continuous.LimPID pid(Ti = 10, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.Init.InitialOutput, k = 0.1, yMax = 6000000) annotation(
          Placement(visible = true, transformation(origin = {48, -18}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Sources.Torque torque(useSupport = false) annotation(
          Placement(visible = true, transformation(origin = {-16, 16}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Gain shaft_eta(k = 0.97) annotation(
          Placement(visible = true, transformation(origin = {14, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        NewHTHP.Machines.Compressor compressor(J_p = 0.005, enableHeatPort = false, enableOutput = true, initEnthalpy = false, initOmega = NewHTHP.Utilities.InitialisationMethod.state, initPhi = true, mdot_0 = 0, omega_fromInput = false) annotation(
          Placement(visible = true, transformation(origin = {-54, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        NewHTHP.Volumes.MixedVolume mixedVolume(T_0(displayUnit = "K") = 300, V = 10000, enableHeatPort = false, enableInlet = true, enableOutlet = true, initEnergy = true, p_min = 999.9999999999999, rho_min(displayUnit = "kg/m3") = 0.001, usePreferredMediumStates = false, use_h0 = false) annotation(
          Placement(visible = true, transformation(origin = {-74, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        NewHTHP.Valves.JTValve jTValve annotation(
          Placement(visible = true, transformation(origin = {-94, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      equation
        connect(step.y, pid.u_s) annotation(
          Line(points = {{75, 0}, {66, 0}, {66, -18}, {59, -18}}, color = {0, 0, 127}));
        connect(pid.y, shaft_eta.u) annotation(
          Line(points = {{37, -18}, {31, -18}, {31, 0}, {26, 0}}, color = {0, 0, 127}));
        connect(shaft_eta.y, torque.tau) annotation(
          Line(points = {{3, 0}, {-1.5, 0}, {-1.5, 16}, {-4, 16}}, color = {0, 0, 127}));
        connect(compressor.output_val, pid.u_m) annotation(
          Line(points = {{-44, -8}, {-13, -8}, {-13, -78}, {48.5, -78}, {48.5, -30}, {48, -30}}, color = {0, 0, 127}));
        connect(torque.flange, compressor.flange) annotation(
          Line(points = {{-26, 16}, {-36, 16}, {-36, 0}, {-44, 0}}));
        connect(compressor.inlet, mixedVolume.outlet) annotation(
          Line(points = {{-54, -10}, {-53, -10}, {-53, -60}, {-64, -60}, {-64, -60}}, color = {85, 0, 255}, thickness = 1));
        connect(compressor.outlet, jTValve.inlet) annotation(
          Line(points = {{-54, 10}, {-52, 10}, {-52, 80}, {-94, 80}, {-94, 6}}, color = {85, 0, 255}));
        connect(mixedVolume.inlet, jTValve.outlet) annotation(
          Line(points = {{-84, -60}, {-84, -14}, {-94, -14}, {-94, -6}}, color = {85, 0, 255}));
        annotation(
          Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
      end Loop_withMixedVolume_test;

      model Loop_withMixedVolume_supply_test
        Modelica.Blocks.Sources.Step step1(height = 10, offset = 25, startTime = 75) annotation(
          Placement(visible = true, transformation(origin = {76, -80}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Step step(height = 30000*(2*Modelica.Constants.pi/60), offset = 15000*(2*Modelica.Constants.pi/60), startTime = 50) annotation(
          Placement(visible = true, transformation(origin = {86, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Continuous.LimPID pid(Ti = 10, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.Init.InitialOutput, k = 0.1, yMax = 6000000) annotation(
          Placement(visible = true, transformation(origin = {48, -18}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Sources.Torque torque(useSupport = false) annotation(
          Placement(visible = true, transformation(origin = {-16, 16}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Gain shaft_eta(k = 0.97) annotation(
          Placement(visible = true, transformation(origin = {14, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        NewHTHP.Machines.Compressor compressor(J_p = 0.005, enableHeatPort = false, enableOutput = true, initEnthalpy = false, initOmega = NewHTHP.Utilities.InitialisationMethod.state, initPhi = true, mdot_0 = 0, omega_fromInput = false) annotation(
          Placement(visible = true, transformation(origin = {-55, -1}, extent = {{-11, -11}, {11, 11}}, rotation = 90)));
        NewHTHP.Volumes.MixedVolume mixedVolume(T_0(displayUnit = "K") = 141 + 273.15, V = 10000, enableHeatPort = false, enableInlet = true, enableOutlet = true, h_0 = 550000, initEnergy = true, p_0 = 1198000, p_min = 999.9999999999999, rho_min(displayUnit = "kg/m3") = 0.001, usePreferredMediumStates = false, use_h0 = false) annotation(
          Placement(visible = true, transformation(origin = {-80, 40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        NewHTHP.Valves.JTValve jTValve annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      equation
        connect(step.y, pid.u_s) annotation(
          Line(points = {{75, 0}, {66, 0}, {66, -18}, {59, -18}}, color = {0, 0, 127}));
        connect(pid.y, shaft_eta.u) annotation(
          Line(points = {{37, -18}, {31, -18}, {31, 0}, {26, 0}}, color = {0, 0, 127}));
        connect(shaft_eta.y, torque.tau) annotation(
          Line(points = {{3, 0}, {-1.5, 0}, {-1.5, 16}, {-4, 16}}, color = {0, 0, 127}));
        connect(compressor.output_val, pid.u_m) annotation(
          Line(points = {{-44, -10}, {-13, -10}, {-13, -78}, {48.5, -78}, {48.5, -30}, {48, -30}}, color = {0, 0, 127}));
        connect(torque.flange, compressor.flange) annotation(
          Line(points = {{-26, 16}, {-36, 16}, {-36, -1}, {-44, -1}}));
        connect(compressor.outlet, mixedVolume.inlet) annotation(
          Line(points = {{-55, 10}, {-54.5, 10}, {-54.5, 40}, {-60, 40}, {-60, 41}, {-70, 41}, {-70, 40}}, color = {85, 0, 255}));
        connect(mixedVolume.outlet, jTValve.inlet) annotation(
          Line(points = {{-90, 40}, {-100, 40}, {-100, 6}}, color = {85, 0, 255}));
        connect(compressor.inlet, jTValve.outlet) annotation(
          Line(points = {{-55, -12}, {-55, -40}, {-101, -40}, {-101, -6}, {-100, -6}}, color = {85, 0, 255}));
        annotation(
          Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
      end Loop_withMixedVolume_supply_test;
    equation

    end Test;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Polygon(fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1.5, points = {{-60, 60}, {0, 80}, {60, 60}, {60, 40}, {60, -40}, {60, -60}, {0, -80}, {-60, -60}, {-60, -40}, {-60, 40}, {-60, 60}}, smooth = Smooth.Bezier), Line(origin = {-80, 0}, points = {{-20, 0}, {20, 0}}, thickness = 1.5), Line(origin = {80, 0}, points = {{-20, 0}, {20, 0}}, thickness = 1.5)}));
  end Volumes;

  package Boundaries
    extends Modelica.Icons.SourcesPackage;

    class Source
      parameter Boolean setEnthalpy = false "Prescribe specific enthalpy instead of temperature?";
      parameter Boolean pressureFromInput = false "Use input connector for pressure?";
      parameter Boolean temperatureFromInput = false "Use input connector for temperature?" annotation(
        Dialog(enable = not setEnthalpy));
      parameter Boolean enthalpyFromInput = false "Use input connector for specific enthalpy" annotation(
        Dialog(enable = setEnthalpy));
      parameter Boolean mdotFromInput = false "Use input connector for mass flow rate" annotation(
        Dialog(enable = setMassFlowRate));
      parameter SI.Temperature T0_par = Medium.T_default "Temperature set value" annotation(
        Dialog(enable = not setEnthalpy and not temperatureFromInput));
      parameter SI.Pressure p0_par = Medium.p_default "Pressure set value" annotation(
        Dialog(enable = not pressureFromInput));
      parameter SI.SpecificEnthalpy h0_par = 2000 "Specific enthalpy set value" annotation(
        Dialog(enable = setEnthalpy and not enthalpyFromInput));
      parameter SI.MassFlowRate mdot0_par = 16.9 "Mass flow rate set value" annotation(
        Dialog(enable = mdotFromInput));
      Modelica.Blocks.Interfaces.RealInput p0_var(unit = "Pa") if pressureFromInput "Pressure input connector [Pa]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, 40}, {-20, 80}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-40, 40}, {0, 80}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput T0_var(unit = "K") if not setEnthalpy and temperatureFromInput "Temperature input connector [K]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-40, -20}, {0, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput h0_var(unit = "J/kg") if setEnthalpy and enthalpyFromInput "Enthalpy input connector [J/kg]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-20, -60}, extent = {{-40, -20}, {0, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput mdot0_var(unit = "kg/s") if mdotFromInput "Mass flow rate connector [kg/s]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{40, -20}, {0, 20}}, rotation = -180)));
      Connections.FluidOutlet outlet annotation(
        Placement(transformation(extent = {{60, -20}, {100, 20}})));
    protected
      Modelica.Blocks.Interfaces.RealInput p0(unit = "Pa") "Internal pressure connector";
      Modelica.Blocks.Interfaces.RealInput T0(unit = "K") "Internal temperature connector";
      Modelica.Blocks.Interfaces.RealInput h0(unit = "J/kg") "Internal enthalpy connector";
      Modelica.Blocks.Interfaces.RealInput mdot0(unit = "kg/s") "Mass flow rate connector";
    equation
      connect(T0_var, T0);
      if not temperatureFromInput or setEnthalpy then
        T0 = T0_par;
      end if;
      connect(p0_var, p0);
      if not pressureFromInput then
        p0 = p0_par;
      end if;
      connect(h0_var, h0);
      if not enthalpyFromInput or not setEnthalpy then
        h0 = h0_par;
      end if;
      connect(mdot0_var, mdot0);
      if not mdotFromInput then
        mdot0 = mdot0_par;
      end if;
      outlet.mdot = -mdot0;
      outlet.p = p0;
      outlet.State = if not setEnthalpy then Fluid.setState_pT(p0, T0) else Fluid.setState_ph(p0, h0);
      annotation(
        Icon(graphics = {Rectangle(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-20, 80}, {20, -80}}), Line(points = {{20, 0}, {80, 0}}, thickness = 1), Line(points = {{-20, 80}, {20, 0}, {-20, -80}}, thickness = 1)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Source;

    class Source_Set
      parameter Boolean setEnthalpy = false "Prescribe specific enthalpy instead of temperature?";
      parameter Boolean pressureFromInput = false "Use input connector for pressure?";
      parameter Boolean temperatureFromInput = false "Use input connector for temperature?" annotation(
        Dialog(enable = not setEnthalpy));
      parameter Boolean enthalpyFromInput = false "Use input connector for specific enthalpy" annotation(
        Dialog(enable = setEnthalpy));
      parameter Boolean mdotFromInput = false "Use input connector for mass flow rate" annotation(
        Dialog(enable = setMassFlowRate));
      parameter SI.Temperature T0_par = Medium.T_default "Temperature set value" annotation(
        Dialog(enable = not setEnthalpy and not temperatureFromInput));
      parameter SI.Pressure p0_par = Medium.p_default "Pressure set value" annotation(
        Dialog(enable = not pressureFromInput));
      parameter SI.SpecificEnthalpy h0_par = 2000 "Specific enthalpy set value" annotation(
        Dialog(enable = setEnthalpy and not enthalpyFromInput));
      parameter SI.MassFlowRate mdot0_par = 16.9 "Mass flow rate set value" annotation(
        Dialog(enable = mdotFromInput));
      Modelica.Blocks.Interfaces.RealInput p0_var(unit = "Pa") if pressureFromInput "Pressure input connector [Pa]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, 40}, {-20, 80}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-40, 40}, {0, 80}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput T0_var(unit = "K") if not setEnthalpy and temperatureFromInput "Temperature input connector [K]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-40, -20}, {0, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput h0_var(unit = "J/kg") if setEnthalpy and enthalpyFromInput "Enthalpy input connector [J/kg]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-20, -60}, extent = {{-40, -20}, {0, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput mdot0_var(unit = "kg/s") if mdotFromInput "Mass flow rate connector [kg/s]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{40, -20}, {0, 20}}, rotation = -180)));
      Connections.FluidOutlet_Set outlet annotation(
        Placement(transformation(extent = {{60, -20}, {100, 20}})));
    protected
      Modelica.Blocks.Interfaces.RealInput p0(unit = "Pa") "Internal pressure connector";
      Modelica.Blocks.Interfaces.RealInput T0(unit = "K") "Internal temperature connector";
      Modelica.Blocks.Interfaces.RealInput h0(unit = "J/kg") "Internal enthalpy connector";
      Modelica.Blocks.Interfaces.RealInput mdot0(unit = "kg/s") "Mass flow rate connector";
    equation
      connect(T0_var, T0);
      if not temperatureFromInput or setEnthalpy then
        T0 = T0_par;
      end if;
      connect(p0_var, p0);
      if not pressureFromInput then
        p0 = p0_par;
      end if;
      connect(h0_var, h0);
      if not enthalpyFromInput or not setEnthalpy then
        h0 = h0_par;
      end if;
      connect(mdot0_var, mdot0);
      if not mdotFromInput then
        mdot0 = mdot0_par;
      end if;
      outlet.mdot = -mdot0;
      outlet.p = p0;
      outlet.T = if not setEnthalpy then T0 else Fluid.temperature(Fluid.setState_ph(p0, h0));
      annotation(
        Icon(graphics = {Rectangle(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-20, 80}, {20, -80}}), Line(points = {{20, 0}, {80, 0}}, thickness = 1), Line(points = {{-20, 80}, {20, 0}, {-20, -80}}, thickness = 1)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Source_Set;

    class Sink
      parameter Boolean pressureFromInput = false "If true pressure comes from real input";
      parameter Boolean mdotFromInput = false "Use input connector for mass flow rate" annotation(
        Dialog(enable = setMassFlowRate));
      parameter Boolean pressureSet = false "If true pressure at sink value";
      parameter SI.Pressure p0_par = Medium.p_default "Pressure setpoint of Sink" annotation(
        Dialog(enable = not pressureFromInput));
      Modelica.Blocks.Interfaces.RealInput p0_var(unit = "Pa") if pressureFromInput "Pressure setpoint [Pa]" annotation(
        Placement(visible = true, transformation(origin = {20, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 180), iconTransformation(origin = {40, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 180)));
      Modelica.Blocks.Interfaces.RealInput mdot0_var(unit = "kg/s") if mdotFromInput "Mass flow rate connector [kg/s]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {82, 0}, extent = {{-40, -20}, {0, 20}}, rotation = 180)));
      NewHTHP.Connections.FluidInlet inlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{60, -20}, {100, 20}}, rotation = 0), iconTransformation(origin = {-160, 0}, extent = {{60, -20}, {100, 20}}, rotation = 0)));
    protected
      Modelica.Blocks.Interfaces.RealInput p0(unit = "Pa") "Internal pressure connector";
      Modelica.Blocks.Interfaces.RealInput mdot0(unit = "kg/s") "Internal pressure connector";
    equation
      connect(p0_var, p0);
      if not pressureFromInput then
        p0 = p0_par;
      end if;
      if pressureSet then
        inlet.p = p0;
      end if;
      connect(mdot0_var, mdot0);
      if mdotFromInput then
        inlet.mdot = mdot0;
      else
        mdot0 = -1;
      end if;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-20, 80}, {20, -80}}), Line(points = {{-80, 0}, {-20, 0}}, thickness = 1), Line(points = {{-20, 80}, {20, 0}, {-20, -80}}, thickness = 1)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Sink;

    class Sink_Set
      parameter Boolean pressureFromInput = false "If true pressure comes from real input";
      parameter Boolean pressureSet = false "If true pressure at sink value";
      parameter SI.Pressure p0_par = Medium.p_default "Pressure setpoint of Sink" annotation(
        Dialog(enable = not pressureFromInput));
      Modelica.Blocks.Interfaces.RealInput p0_var(unit = "Pa") if pressureFromInput "Pressure setpoint [Pa]" annotation(
        Placement(visible = true, transformation(origin = {20, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 180), iconTransformation(origin = {40, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 180)));
      NewHTHP.Connections.FluidInlet_Set inlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{60, -20}, {100, 20}}, rotation = 0), iconTransformation(origin = {-160, 0}, extent = {{60, -20}, {100, 20}}, rotation = 0)));
    protected
      Modelica.Blocks.Interfaces.RealInput p0(unit = "Pa") "Internal pressure connector";
    equation
      connect(p0_var, p0);
      if not pressureFromInput then
        p0 = p0_par;
      end if;
      if pressureSet then
        inlet.p = p0;
      end if;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-20, 80}, {20, -80}}), Line(points = {{-80, 0}, {-20, 0}}, thickness = 1), Line(points = {{-20, 80}, {20, 0}, {-20, -80}}, thickness = 1)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Sink_Set;

    class Source_Set_refrigerant
      parameter Boolean setEnthalpy = false "Prescribe specific enthalpy instead of temperature?";
      parameter Boolean pressureFromInput = false "Use input connector for pressure?";
      parameter Boolean temperatureFromInput = false "Use input connector for temperature?" annotation(
        Dialog(enable = not setEnthalpy));
      parameter Boolean enthalpyFromInput = false "Use input connector for specific enthalpy" annotation(
        Dialog(enable = setEnthalpy));
      parameter Boolean mdotFromInput = false "Use input connector for mass flow rate" annotation(
        Dialog(enable = setMassFlowRate));
      parameter SI.Temperature T0_par = Medium.T_default "Temperature set value" annotation(
        Dialog(enable = not setEnthalpy and not temperatureFromInput));
      parameter SI.Pressure p0_par = Medium.p_default "Pressure set value" annotation(
        Dialog(enable = not pressureFromInput));
      parameter SI.SpecificEnthalpy h0_par = 2000 "Specific enthalpy set value" annotation(
        Dialog(enable = setEnthalpy and not enthalpyFromInput));
      parameter SI.MassFlowRate mdot0_par = 16.9 "Mass flow rate set value" annotation(
        Dialog(enable = mdotFromInput));
      Modelica.Blocks.Interfaces.RealInput p0_var(unit = "Pa") if pressureFromInput "Pressure input connector [Pa]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, 40}, {-20, 80}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-40, 40}, {0, 80}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput T0_var(unit = "K") if not setEnthalpy and temperatureFromInput "Temperature input connector [K]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-40, -20}, {0, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput h0_var(unit = "J/kg") if setEnthalpy and enthalpyFromInput "Enthalpy input connector [J/kg]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-20, -60}, extent = {{-40, -20}, {0, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput mdot0_var(unit = "kg/s") if mdotFromInput "Mass flow rate connector [kg/s]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{40, -20}, {0, 20}}, rotation = -180)));
      Connections.FluidOutlet_Set outlet annotation(
        Placement(transformation(extent = {{60, -20}, {100, 20}})));
    protected
      Modelica.Blocks.Interfaces.RealInput p0(unit = "Pa") "Internal pressure connector";
      Modelica.Blocks.Interfaces.RealInput T0(unit = "K") "Internal temperature connector";
      Modelica.Blocks.Interfaces.RealInput h0(unit = "J/kg") "Internal enthalpy connector";
      Modelica.Blocks.Interfaces.RealInput mdot0(unit = "kg/s") "Mass flow rate connector";
    equation
      connect(T0_var, T0);
      if not temperatureFromInput or setEnthalpy then
        T0 = T0_par;
      end if;
      connect(p0_var, p0);
      if not pressureFromInput then
        p0 = p0_par;
      end if;
      connect(h0_var, h0);
      if not enthalpyFromInput or not setEnthalpy then
        h0 = h0_par;
      end if;
      connect(mdot0_var, mdot0);
      if not mdotFromInput then
        mdot0 = mdot0_par;
      end if;
      outlet.mdot = -mdot0;
      outlet.p = p0;
      outlet.T = if not setEnthalpy then T0 else Refrigerant.temperature(Refrigerant.setState_ph(p0, h0));
      annotation(
        Icon(graphics = {Rectangle(fillColor = {220, 38, 127}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-20, 80}, {20, -80}}), Line(points = {{20, 0}, {80, 0}}, thickness = 1), Line(points = {{-20, 80}, {20, 0}, {-20, -80}}, thickness = 1)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Source_Set_refrigerant;

    class Source_refrigerant
      parameter Boolean setEnthalpy = false "Prescribe specific enthalpy instead of temperature?";
      parameter Boolean pressureFromInput = false "Use input connector for pressure?";
      parameter Boolean temperatureFromInput = false "Use input connector for temperature?" annotation(
        Dialog(enable = not setEnthalpy));
      parameter Boolean enthalpyFromInput = false "Use input connector for specific enthalpy" annotation(
        Dialog(enable = setEnthalpy));
      parameter Boolean mdotFromInput = false "Use input connector for mass flow rate" annotation(
        Dialog(enable = setMassFlowRate));
      parameter SI.Temperature T0_par = Medium.T_default "Temperature set value" annotation(
        Dialog(enable = not setEnthalpy and not temperatureFromInput));
      parameter SI.Pressure p0_par = Medium.p_default "Pressure set value" annotation(
        Dialog(enable = not pressureFromInput));
      parameter SI.SpecificEnthalpy h0_par = 2000 "Specific enthalpy set value" annotation(
        Dialog(enable = setEnthalpy and not enthalpyFromInput));
      parameter SI.MassFlowRate mdot0_par = 16.9 "Mass flow rate set value" annotation(
        Dialog(enable = mdotFromInput));
      Modelica.Blocks.Interfaces.RealInput p0_var(unit = "Pa") if pressureFromInput "Pressure input connector [Pa]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, 40}, {-20, 80}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-40, 40}, {0, 80}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput T0_var(unit = "K") if not setEnthalpy and temperatureFromInput "Temperature input connector [K]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-40, -20}, {0, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput h0_var(unit = "J/kg") if setEnthalpy and enthalpyFromInput "Enthalpy input connector [J/kg]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-20, -60}, extent = {{-40, -20}, {0, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput mdot0_var(unit = "kg/s") if mdotFromInput "Mass flow rate connector [kg/s]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{40, -20}, {0, 20}}, rotation = -180)));
      Connections.FluidOutlet_refrigerant outlet annotation(
        Placement(transformation(extent = {{60, -20}, {100, 20}})));
    protected
      Modelica.Blocks.Interfaces.RealInput p0(unit = "Pa") "Internal pressure connector";
      Modelica.Blocks.Interfaces.RealInput T0(unit = "K") "Internal temperature connector";
      Modelica.Blocks.Interfaces.RealInput h0(unit = "J/kg") "Internal enthalpy connector";
      Modelica.Blocks.Interfaces.RealInput mdot0(unit = "kg/s") "Mass flow rate connector";
    equation
      connect(T0_var, T0);
      if not temperatureFromInput or setEnthalpy then
        T0 = T0_par;
      end if;
      connect(p0_var, p0);
      if not pressureFromInput then
        p0 = p0_par;
      end if;
      connect(h0_var, h0);
      if not enthalpyFromInput or not setEnthalpy then
        h0 = h0_par;
      end if;
      connect(mdot0_var, mdot0);
      if not mdotFromInput then
        mdot0 = mdot0_par;
      end if;
      outlet.mdot = -mdot0;
      outlet.p = p0;
      outlet.State = if not setEnthalpy then Refrigerant.setState_pT(p0, T0) else Refrigerant.setState_ph(p0, h0);
      annotation(
        Icon(graphics = {Rectangle(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-20, 80}, {20, -80}}), Line(points = {{20, 0}, {80, 0}}, thickness = 1), Line(points = {{-20, 80}, {20, 0}, {-20, -80}}, thickness = 1)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Source_refrigerant;

    class Sink_refrigerant
      parameter Boolean pressureFromInput = false "If true pressure comes from real input";
      parameter Boolean pressureSet = false "If true pressure at sink value";
      parameter SI.Pressure p0_par = Medium.p_default "Pressure setpoint of Sink" annotation(
        Dialog(enable = not pressureFromInput));
      Modelica.Blocks.Interfaces.RealInput p0_var(unit = "Pa") if pressureFromInput "Pressure setpoint [Pa]" annotation(
        Placement(visible = true, transformation(origin = {20, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 180), iconTransformation(origin = {40, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 180)));
      NewHTHP.Connections.FluidInlet_refrigerant inlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{60, -20}, {100, 20}}, rotation = 0), iconTransformation(origin = {-160, 0}, extent = {{60, -20}, {100, 20}}, rotation = 0)));
    protected
      Modelica.Blocks.Interfaces.RealInput p0(unit = "Pa") "Internal pressure connector";
    equation
      connect(p0_var, p0);
      if not pressureFromInput then
        p0 = p0_par;
      end if;
      if pressureSet then
        inlet.p = p0;
      end if;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-20, 80}, {20, -80}}), Line(points = {{-80, 0}, {-20, 0}}, thickness = 1), Line(points = {{-20, 80}, {20, 0}, {-20, -80}}, thickness = 1)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Sink_refrigerant;
    annotation(
      Documentation(revisions = "<html>
  <p><img src=\"modelica:/ThermofluidStream/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</p>
  
  </html>", info = "<html>
  <p>The boundaries are Sorces and Sinks, as well as Volumes, that are conceptually a source and a sink with extra equations and act as loop breakers in closes cycles, and therefore are also boundaries.</p>
  </html>"));
  end Boundaries;

  package Networks "This package contains blocks that tranport and transform flow variables"
    extends Modelica.Icons.Package;

    model Tee
      NewHTHP.Connections.FluidInlet inlet1 annotation(
        Placement(visible = true, transformation(origin = {60, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {66, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidInlet inlet2 annotation(
        Placement(visible = true, transformation(origin = {100, -60}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, 66}, extent = {{-115, 15}, {-85, -15}}, rotation = 90)));
      NewHTHP.Connections.FluidOutlet outlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {-40, 0}, extent = {{85, 15}, {115, -15}}, rotation = 0)));
      SI.SpecificEnthalpy h;
      
    equation
      inlet1.mdot + inlet2.mdot + outlet.mdot = 0;
      inlet1.State.h*inlet1.mdot + inlet2.State.h*inlet2.mdot + h*outlet.mdot = 0;
      inlet1.p = outlet.p;
      
      outlet.State = Fluid.setState_ph(inlet1.p, h);
      annotation(
        Icon(graphics = {Line(points = {{20, 0}, {60, 0}}, color = {85, 0, 255}, thickness = 1.5), Ellipse(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-20, 20}, {20, -20}}), Line(points = {{20, 0}, {60, 0}}, color = {85, 0, 255}, thickness = 1.5)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Tee;
    
    model Tee_pMatch
      NewHTHP.Connections.FluidInlet inlet1 annotation(
        Placement(visible = true, transformation(origin = {60, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {66, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidInlet inlet2 annotation(
        Placement(visible = true, transformation(origin = {100, -60}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {0, 66}, extent = {{-115, 15}, {-85, -15}}, rotation = 90)));
      NewHTHP.Connections.FluidOutlet outlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {-40, 0}, extent = {{85, 15}, {115, -15}}, rotation = 0)));
      SI.SpecificEnthalpy h;
      
    equation
      inlet1.mdot + inlet2.mdot + outlet.mdot = 0;
      inlet1.State.h*inlet1.mdot + inlet2.State.h*inlet2.mdot + h*outlet.mdot = 0;
      inlet1.p = outlet.p;
      inlet1.p = inlet2.p;
      
      outlet.State = Fluid.setState_ph(inlet1.p, h);
      annotation(
        Icon(graphics = {Line(points = {{20, 0}, {60, 0}}, color = {85, 0, 255}, thickness = 1.5), Ellipse(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-20, 20}, {20, -20}}), Line(points = {{20, 0}, {60, 0}}, color = {85, 0, 255}, thickness = 1.5)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Tee_pMatch;
    annotation(
      Icon(graphics = {Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1, extent = {{40, 10}, {60, -10}}), Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1, extent = {{-40, -40}, {-60, -60}}), Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward, lineThickness = 1, extent = {{-60, 60}, {-40, 40}}), Line(points = {{-40, 50}, {40, 0}, {-40, -50}}, thickness = 1)}));
  end Networks;

  package Utilities
    extends Modelica.Icons.UtilitiesPackage;
    
    model PRV
      parameter Boolean pressureFromInput = false "If true pressure comes from real input";
      parameter SI.Pressure p0_par = 200000 "Pressure setpoint of Sink";
      NewHTHP.Connections.FluidInlet inlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 0), iconTransformation(origin = {40, 0}, extent = {{-115, 15}, {-85, -15}}, rotation = 0)));
      NewHTHP.Connections.FluidOutlet outlet annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{85, 15}, {115, -15}}, rotation = 0), iconTransformation(origin = {-40, 0}, extent = {{85, 15}, {115, -15}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput p0_var(unit = "Pa") if pressureFromInput "Imposed downstream pressure [kg/s]" annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-60, -40}, {-20, 0}}, rotation = 0), iconTransformation(origin = {0, 20}, extent = {{-40, -20}, {0, 20}}, rotation = 270))); 
      
    equation
      if pressureFromInput then
        outlet.p = p0_var;
      else
        outlet.p = p0_par;
      end if;
      
      inlet.mdot + outlet.mdot = 0;
      outlet.State = Fluid.setState_ph(outlet.p,inlet.State.h);
      
      annotation(
        Icon(graphics = {Line(points = {{-60, 0}, {60, 0}}, color = {85, 0, 255}, thickness = 1.5), Polygon(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, points = {{0, 0}, {-40, 20}, {-40, -20}, {0, 0}}), Polygon(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1.5, points = {{0, 0}, {40, 20}, {40, -20}, {0, 0}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end PRV;
    type InitialisationMethod = enumeration(none "No initialization", steadyState "Steady state initialization (derivatives of states are zero)", state "Initialization with initial states", derivative "Initialization with initial derivatives of states") "Choices for initialization of a state.";
  end Utilities;

  package Connections
    extends Modelica.Icons.InterfacesPackage;

    connector FluidInlet "Inlet port for a fluid"
      SI.Pressure p "Absolute pressure";
      input Fluid.ThermodynamicState State "Thermodynamic state";
      flow SI.MassFlowRate mdot "Mass flow rate";
      //Annotation
      annotation(
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Polygon(lineColor = {85, 0, 255}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, points = {{-40, 0}, {-80, 80}, {80, 0}, {-80, -80}, {-40, 0}})}));
    end FluidInlet;

    connector FluidInlet_Set "Inlet port for a fluid"
      SI.Pressure p "Absolute pressure";
      SI.ThermodynamicTemperature T "Temperature";
      flow SI.MassFlowRate mdot "Mass flow rate";
      //Annotation
      annotation(
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Polygon(fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, lineThickness = 1, points = {{-40, 0}, {-80, 80}, {80, 0}, {-80, -80}, {-40, 0}})}));
    end FluidInlet_Set;

    connector FluidOutlet "Outlet port for a fluid"
      SI.Pressure p "Absolute pressure";
      output Fluid.ThermodynamicState State "Thermodynamic state";
      flow SI.MassFlowRate mdot "Mass flow rate";
      //Annotation
      annotation(
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Polygon(lineColor = {85, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, points = {{-40, 0}, {-80, 80}, {80, 0}, {-80, -80}, {-40, 0}})}));
    end FluidOutlet;

    connector FluidOutlet_Set "Outlet port for a fluid"
      SI.Pressure p "Absolute pressure";
      SI.ThermodynamicTemperature T "Temperature";
      flow SI.MassFlowRate mdot "Mass flow rate";
      //Annotation
      annotation(
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Polygon(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, points = {{-40, 0}, {-80, 80}, {80, 0}, {-80, -80}, {-40, 0}})}));
    end FluidOutlet_Set;

    connector FluidInlet_refrigerant "Inlet port for a fluid"
      SI.Pressure p "Absolute pressure";
      input Refrigerant.ThermodynamicState State "Thermodynamic state";
      flow SI.MassFlowRate mdot "Mass flow rate";
      //Annotation
      annotation(
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Polygon(lineColor = {220, 38, 127}, fillColor = {254, 97, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, points = {{-40, 0}, {-80, 80}, {80, 0}, {-80, -80}, {-40, 0}})}));
    end FluidInlet_refrigerant;

    connector FluidOutlet_refrigerant "Outlet port for a fluid"
      SI.Pressure p "Absolute pressure";
      output Refrigerant.ThermodynamicState State "Thermodynamic state";
      flow SI.MassFlowRate mdot "Mass flow rate";
      //Annotation
      annotation(
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Polygon(lineColor = {220, 38, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, points = {{-40, 0}, {-80, 80}, {80, 0}, {-80, -80}, {-40, 0}})}));
    end FluidOutlet_refrigerant;
  end Connections;

  package Examples "Examples for this Library"
    extends Modelica.Icons.ExamplesPackage;

    model Energies_Condenser
      NewHTHP.Boundaries.Sink_refrigerant sink_refrigerant(p0_par = 1713000) annotation(
        Placement(visible = true, transformation(origin = {-30, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source_refrigerant source_refrigerant(T0_par(displayUnit = "K") = 434.4, mdot0_par = 19.75, mdotFromInput = true, p0_par(displayUnit = "Pa") = 1713000, temperatureFromInput = true) annotation(
        Placement(visible = true, transformation(origin = {50, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source source(T0_par(displayUnit = "K") = 416.75, mdot0_par = 2.33, mdotFromInput = true, p0_par(displayUnit = "Pa") = 399999.9999999999) annotation(
        Placement(visible = true, transformation(origin = {-30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      NewHTHP.Boundaries.Sink sink(mdotFromInput = false, p0_par(displayUnit = "Pa") = 399999.9999999999) annotation(
        Placement(visible = true, transformation(origin = {50, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      NewHTHP.HeatExchangers.CtrFlow_IHX_Upwind_Refrigerant ctrFlow_IHX_Upwind_Refrigerant(A = 910, N_pA = 800, N_pB = 50, T_0 = 145 + 273.15) annotation(
        Placement(visible = true, transformation(origin = {10, 10}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Ramp Feed_flow(duration = 100, height = -(2.33 - 2.83), offset = 2.33, startTime = 500)  annotation(
        Placement(visible = true, transformation(origin = {-70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_flow(duration = 100, height = -(19.75 - 23.95), offset = 19.75, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, -50}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_Ts(duration = 100, height = -(434.4 - 434), offset = 434.4, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(source.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletB) annotation(
        Line(points = {{-30, 42}, {0, 42}, {0, 10}}, color = {85, 0, 255}));
      connect(sink.inlet, ctrFlow_IHX_Upwind_Refrigerant.outletB) annotation(
        Line(points = {{50, 42}, {20, 42}, {20, 10}}, color = {85, 0, 255}));
      connect(source_refrigerant.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletA) annotation(
        Line(points = {{42, -30}, {14, -30}, {14, 0}}, color = {220, 38, 127}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletA, sink_refrigerant.inlet) annotation(
        Line(points = {{6, 0}, {-22, 0}, {-22, -30}}, color = {220, 38, 127}));
      connect(source_refrigerant.mdot0_var, WF_flow.y) annotation(
        Line(points = {{60, -30}, {60, -50}, {79, -50}}, color = {0, 0, 127}));
      connect(source_refrigerant.T0_var, WF_Ts.y) annotation(
        Line(points = {{54, -30}, {54, -10}, {79, -10}}, color = {0, 0, 127}));
      connect(source.mdot0_var, Feed_flow.y) annotation(
        Line(points = {{-30, 60}, {-30, 70}, {-59, 70}}, color = {0, 0, 127}));
      annotation(
        uses(Modelica(version = "4.0.0")),
        Diagram);
    end Energies_Condenser;

    model Energies_Condenser_and_SA_TU
    
      NewHTHP.Boundaries.Sink_refrigerant sink_refrigerant(p0_par (displayUnit = "Pa")= 1713000) annotation(
        Placement(visible = true, transformation(origin = {-30, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source_refrigerant source_refrigerant(T0_par(displayUnit = "K") = 434.4, mdot0_par = 19.75, mdotFromInput = true, p0_par(displayUnit = "Pa") = 1713000, temperatureFromInput = true) annotation(
        Placement(visible = true, transformation(origin = {50, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source source(T0_par(displayUnit = "K") = 416.75, mdot0_par = 2.33, mdotFromInput = true, p0_par(displayUnit = "bar") = 399999.9999999999, pressureFromInput = false, setEnthalpy = false) annotation(
        Placement(visible = true, transformation(origin = {-50, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      NewHTHP.HeatExchangers.CtrFlow_IHX_Upwind_Refrigerant ctrFlow_IHX_Upwind_Refrigerant(A = 910, N_pA = 800, N_pB = 50, T_0 = 145 + 273.15) annotation(
        Placement(visible = true, transformation(origin = {12, 34}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
      NewHTHP.Boundaries.Sink sink(p0_par(displayUnit = "Pa") = 1) annotation(
        Placement(visible = true, transformation(origin = {106, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      NewHTHP.Volumes.SteamAccumulator_pRulez steamAccumulator_pRulez(V = 116, enableOutlet = false, initEnergy = false, initPressure = true, k_flow_in = (2.83 - 2.33)/200000, p_0 = 200000, phi_0 = 0.748)  annotation(
        Placement(visible = true, transformation(origin = {32, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp Feed_flow(duration = 100, height = -(2.33 - 2.83), offset = 2.33, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {-90, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_Ts(duration = 100, height = -(434.4 - 434), offset = 434.4, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_flow(duration = 100, height = -(19.75 - 23.95), offset = 19.75, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(source.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletB) annotation(
        Line(points = {{-50, 62}, {2, 62}, {2, 34}}, color = {85, 0, 255}));
      connect(source_refrigerant.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletA) annotation(
        Line(points = {{42, -30}, {42, -28}, {16, -28}, {16, 24}}, color = {220, 38, 127}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletA, sink_refrigerant.inlet) annotation(
        Line(points = {{8, 24}, {8, 0}, {-22, 0}, {-22, -30}}, color = {220, 38, 127}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletB, steamAccumulator_pRulez.inlet) annotation(
        Line(points = {{22, 34}, {32, 34}, {32, 26}}, color = {85, 0, 255}));
      connect(sink.inlet, ctrFlow_IHX_Upwind_Refrigerant.outletB) annotation(
        Line(points = {{98, 64}, {22, 64}, {22, 34}}, color = {85, 0, 255}));
      connect(Feed_flow.y, source.mdot0_var) annotation(
        Line(points = {{-78, 50}, {-50, 50}, {-50, 80}}, color = {0, 0, 127}));
      connect(WF_Ts.y, source_refrigerant.T0_var) annotation(
        Line(points = {{80, 30}, {54, 30}, {54, -30}}, color = {0, 0, 127}));
      connect(WF_flow.y, source_refrigerant.mdot0_var) annotation(
        Line(points = {{79, -30}, {60, -30}}, color = {0, 0, 127}));
      annotation(
        uses(Modelica(version = "4.0.0")),
        Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));

    end Energies_Condenser_and_SA_TU;

    model Energies_Condenser_and_SA_TD
    
      NewHTHP.Boundaries.Sink_refrigerant sink_refrigerant(p0_par = 1713000) annotation(
        Placement(visible = true, transformation(origin = {-30, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source_refrigerant source_refrigerant(T0_par(displayUnit = "K") = 434.4, mdot0_par = 19.75, mdotFromInput = true, p0_par(displayUnit = "Pa") = 1713000, temperatureFromInput = true) annotation(
        Placement(visible = true, transformation(origin = {50, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source source(T0_par(displayUnit = "K") = 416.75, mdot0_par = 2.33, mdotFromInput = true, p0_par(displayUnit = "bar") = 399999.9999999999) annotation(
        Placement(visible = true, transformation(origin = {-30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      NewHTHP.HeatExchangers.CtrFlow_IHX_Upwind_Refrigerant ctrFlow_IHX_Upwind_Refrigerant(A = 910, N_pA = 800, N_pB = 50, T_0 = 145 + 273.15) annotation(
        Placement(visible = true, transformation(origin = {12, 34}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
      NewHTHP.Boundaries.Sink sink(p0_par(displayUnit = "Pa") = 1) annotation(
        Placement(visible = true, transformation(origin = {106, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      NewHTHP.Volumes.SteamAccumulator_pRulez steamAccumulator_pRulez(V = 116,enableInlet = false, enableOutlet = true, initEnergy = false, k_flow_out = (2.33 - 1.21)/200000, p_0 = 399999.9999999999, phi_0 = 0.8)  annotation(
        Placement(visible = true, transformation(origin = {32, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      NewHTHP.Networks.Tee_pMatch tee_pMatch annotation(
        Placement(visible = true, transformation(origin = {74, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      NewHTHP.Utilities.PRV prv(p0_par = 200000, pressureFromInput = false)  annotation(
        Placement(visible = true, transformation(origin = {52, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp Feed_flow(duration = 250, height = -(2.33 - 1.21), offset = 2.33, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {-90, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_Ts(duration = 250, height = -(434.4 - 433), offset = 434.4, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_flow(duration = 250, height = -(19.75 - 10.28), offset = 19.75, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(source.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletB) annotation(
        Line(points = {{-30, 42}, {-30, 44}, {2, 44}, {2, 34}}, color = {85, 0, 255}));
      connect(source_refrigerant.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletA) annotation(
        Line(points = {{42, -30}, {42, -28}, {16, -28}, {16, 24}}, color = {220, 38, 127}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletA, sink_refrigerant.inlet) annotation(
        Line(points = {{8, 24}, {8, 0}, {-22, 0}, {-22, -30}}, color = {220, 38, 127}));
      connect(steamAccumulator_pRulez.outlet, tee_pMatch.inlet2) annotation(
        Line(points = {{42, 14}, {74, 14}, {74, 52}}, color = {85, 0, 255}));
      connect(tee_pMatch.outlet, sink.inlet) annotation(
        Line(points = {{80, 56}, {98, 56}, {98, 64}}, color = {85, 0, 255}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletB, prv.inlet) annotation(
        Line(points = {{22, 34}, {46, 34}, {46, 56}}, color = {85, 0, 255}));
      connect(prv.outlet, tee_pMatch.inlet1) annotation(
        Line(points = {{58, 56}, {70, 56}}, color = {85, 0, 255}));
      connect(source_refrigerant.mdot0_var, WF_flow.y) annotation(
        Line(points = {{60, -30}, {79, -30}}, color = {0, 0, 127}));
      connect(WF_Ts.y, source_refrigerant.T0_var) annotation(
        Line(points = {{80, 30}, {54, 30}, {54, -30}}, color = {0, 0, 127}));
      connect(Feed_flow.y, source.mdot0_var) annotation(
        Line(points = {{-78, 50}, {-30, 50}, {-30, 60}}, color = {0, 0, 127}));
      annotation(
        uses(Modelica(version = "4.0.0")),
        Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));

    end Energies_Condenser_and_SA_TD;

    model Energies_Condenser_and_8barSA_TU
    
      NewHTHP.Boundaries.Sink_refrigerant sink_refrigerant(p0_par (displayUnit = "Pa")= 1713000) annotation(
        Placement(visible = true, transformation(origin = {-30, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source_refrigerant source_refrigerant(T0_par(displayUnit = "K") = 434.4, mdot0_par = 19.75, mdotFromInput = true, p0_par(displayUnit = "Pa") = 1713000, temperatureFromInput = true) annotation(
        Placement(visible = true, transformation(origin = {50, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source source(T0_par(displayUnit = "K") = 416.75, mdot0_par = 2.33, mdotFromInput = true, p0_par(displayUnit = "bar") = 399999.9999999999, pressureFromInput = true, setEnthalpy = false) annotation(
        Placement(visible = true, transformation(origin = {-50, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      NewHTHP.HeatExchangers.CtrFlow_IHX_Upwind_Refrigerant ctrFlow_IHX_Upwind_Refrigerant(A = 910, N_pA = 800, N_pB = 50, T_0 = 145 + 273.15) annotation(
        Placement(visible = true, transformation(origin = {12, 34}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.TimeTable Feed_p_table(table = [0, 400000; 500, 400000; 500.01, 800000; 4101, 800000]) annotation(
        Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Sink sink(p0_par(displayUnit = "Pa") = 1) annotation(
        Placement(visible = true, transformation(origin = {106, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      NewHTHP.Volumes.SteamAccumulator_pRulez steamAccumulator_pRulez(V = 103, enableOutlet = false, initEnergy = false, initPressure = true, k_flow_in = (2.56 - 2.33)/400000, p_0 = 399999.9999999999, phi_0 = 0.736)  annotation(
        Placement(visible = true, transformation(origin = {32, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp Feed_flow(duration = 10, height = -(2.33 - 2.56), offset = 2.33, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {-90, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_Ts(duration = 10, height = -(434.4 - 453.4), offset = 434.4, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_flow(duration = 10, height = -(19.75 - 21.71), offset = 19.75, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Utilities.PRV prv(p0_par = 399999.9999999999)  annotation(
        Placement(visible = true, transformation(origin = {46, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(source.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletB) annotation(
        Line(points = {{-50, 62}, {2, 62}, {2, 34}}, color = {85, 0, 255}));
      connect(source_refrigerant.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletA) annotation(
        Line(points = {{42, -30}, {42, -28}, {16, -28}, {16, 24}}, color = {220, 38, 127}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletA, sink_refrigerant.inlet) annotation(
        Line(points = {{8, 24}, {8, 0}, {-22, 0}, {-22, -30}}, color = {220, 38, 127}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletB, steamAccumulator_pRulez.inlet) annotation(
        Line(points = {{22, 34}, {32, 34}, {32, 26}}, color = {85, 0, 255}));
      connect(Feed_flow.y, source.mdot0_var) annotation(
        Line(points = {{-78, 50}, {-50, 50}, {-50, 80}}, color = {0, 0, 127}));
      connect(WF_Ts.y, source_refrigerant.T0_var) annotation(
        Line(points = {{80, 30}, {54, 30}, {54, -30}}, color = {0, 0, 127}));
      connect(WF_flow.y, source_refrigerant.mdot0_var) annotation(
        Line(points = {{79, -30}, {60, -30}}, color = {0, 0, 127}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletB, prv.inlet) annotation(
        Line(points = {{22, 34}, {40, 34}, {40, 64}}, color = {85, 0, 255}));
      connect(prv.outlet, sink.inlet) annotation(
        Line(points = {{52, 64}, {98, 64}}, color = {85, 0, 255}));
      connect(Feed_p_table.y, source.p0_var) annotation(
        Line(points = {{-78, 90}, {-44, 90}, {-44, 74}}, color = {0, 0, 127}));
      annotation(
        uses(Modelica(version = "4.0.0")),
        Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));

    end Energies_Condenser_and_8barSA_TU;

    model Energies_Condenser_and_8barSA_TD

      NewHTHP.Boundaries.Sink_refrigerant sink_refrigerant(p0_par = 1713000) annotation(
        Placement(visible = true, transformation(origin = {-30, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source_refrigerant source_refrigerant(T0_par(displayUnit = "K") = 434.4, mdot0_par = 19.75, mdotFromInput = true, p0_par(displayUnit = "Pa") = 1713000, temperatureFromInput = true) annotation(
        Placement(visible = true, transformation(origin = {50, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      NewHTHP.Boundaries.Source source(T0_par(displayUnit = "K") = 416.75, mdot0_par = 2.33, mdotFromInput = true, p0_par(displayUnit = "bar") = 399999.9999999999) annotation(
        Placement(visible = true, transformation(origin = {-30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      NewHTHP.HeatExchangers.CtrFlow_IHX_Upwind_Refrigerant ctrFlow_IHX_Upwind_Refrigerant(A = 910, N_pA = 800, N_pB = 50, T_0 = 145 + 273.15) annotation(
        Placement(visible = true, transformation(origin = {12, 34}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
      NewHTHP.Boundaries.Sink sink(p0_par(displayUnit = "Pa") = 1) annotation(
        Placement(visible = true, transformation(origin = {106, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      NewHTHP.Volumes.SteamAccumulator_pRulez steamAccumulator_pRulez(V = 116,enableInlet = false, enableOutlet = true, initEnergy = false, k_flow_out = (2.33 - 1.21)/400000, p_0 = 799999.9999999999, phi_0 = 0.8)  annotation(
        Placement(visible = true, transformation(origin = {32, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      NewHTHP.Networks.Tee_pMatch tee_pMatch annotation(
        Placement(visible = true, transformation(origin = {74, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp Feed_flow(duration = 250, height = -(2.33 - 1.21), offset = 2.33, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {-90, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_Ts(duration = 250, height = -(434.4 - 433), offset = 434.4, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp WF_flow(duration = 250, height = -(19.75 - 10.28), offset = 19.75, startTime = 500) annotation(
        Placement(visible = true, transformation(origin = {90, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(source.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletB) annotation(
        Line(points = {{-30, 42}, {-30, 44}, {2, 44}, {2, 34}}, color = {85, 0, 255}));
      connect(source_refrigerant.outlet, ctrFlow_IHX_Upwind_Refrigerant.inletA) annotation(
        Line(points = {{42, -30}, {42, -28}, {16, -28}, {16, 24}}, color = {220, 38, 127}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletA, sink_refrigerant.inlet) annotation(
        Line(points = {{8, 24}, {8, 0}, {-22, 0}, {-22, -30}}, color = {220, 38, 127}));
      connect(tee_pMatch.outlet, sink.inlet) annotation(
        Line(points = {{80, 56}, {98, 56}, {98, 64}}, color = {85, 0, 255}));
      connect(source_refrigerant.mdot0_var, WF_flow.y) annotation(
        Line(points = {{60, -30}, {79, -30}}, color = {0, 0, 127}));
      connect(WF_Ts.y, source_refrigerant.T0_var) annotation(
        Line(points = {{80, 30}, {54, 30}, {54, -30}}, color = {0, 0, 127}));
      connect(Feed_flow.y, source.mdot0_var) annotation(
        Line(points = {{-78, 50}, {-30, 50}, {-30, 60}}, color = {0, 0, 127}));
      connect(ctrFlow_IHX_Upwind_Refrigerant.outletB, tee_pMatch.inlet1) annotation(
        Line(points = {{22, 34}, {70, 34}, {70, 56}}, color = {85, 0, 255}));
      connect(steamAccumulator_pRulez.outlet, tee_pMatch.inlet2) annotation(
        Line(points = {{42, 14}, {74, 14}, {74, 52}}, color = {85, 0, 255}));
      annotation(
        uses(Modelica(version = "4.0.0")),
        Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    
    end Energies_Condenser_and_8barSA_TD;
    annotation(
      Documentation(revisions = "<html>
  <p><img src=\"modelica:/ThermofluidStream/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</p>
  </html>", info = "<html>
  <p>
  This package contains several application examples that demonstarte the capabilities
  of the library and can act as a stating point for the user.
  </p>
  </html>"));
  end Examples;
  annotation(
    uses(Modelica(version = "4.0.0")));
end NewHTHP;
