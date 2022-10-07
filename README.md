This model was used in:
Jan F. Wiegner, Alexa Grimm, Lukas Weimann, and Matteo Gazzani (2022): Optimal Design and Operation of Solid Sorbent Direct Air Capture Processes at Varying Ambient Conditions, Industrial & Engineering Chemistry Research 61 (34), 12649-12667, DOI: 10.1021/acs.iecr.2c00681

# Dependencies
Obligatory:
-	YALMIP
- any MILP solver compatible with YALMIP (can be from the Matlab optimization toolbox)
- Matlab Statistics and Machine Learning Toolbox

Nice to have:
-	Gurobi or CLEX as external solvers


# Package content
Main: Loads and sets data and performs optimization.

## Data:

**AirInletCooling** contains data on the cooling effect of spraying water into the inlet air of the DAC unit

**DAC_Data_5_ConstOP** contains data on the process performance at constant operating parameters and different ambient conditions

**DAC_Data_5** contains data on the process performance at optimized operating parameters and different ambient conditions

**In Input_DAC** you can specify economic data and emission factors for the optimization:
  -	Demand is in t of CO2
  -	Module Investment cost is in EUR/module
  -	Lifetime is in years
  -	Maintenance Cost in % of annualized investment costs
  -	Carrier prices are in EUR/kWh
  -	Emission Factors are in t/kWh
  -	Eta_elth depicts the electric efficiency of the ohmic heating

**Settings** contains general data on the model assumptions:
  -	  Constant demand
      1: demand needs to be met as specified in demand
      0: demand needs to be met over the full period
  - D specifies the number of typical days. Input data is clustered by k-means respectively.
  - N specifies the number of total days
  - K specifies the number of time intervals per dayo
  - Size min is the minimal number of modules that need to be built
  - Size mas is the maximal number of modules that can be built
  - Optimization type (1: minimize cost, 2: minimize emissions, 3: minimize cost at emission constraint. You can set an emission constraints by specifying Settings.EmissionConstraint to a certain value)
  - ElThTradeOff (1: Ohmic heating allowed, 2: Ohmic heating turned off)
  - Nr_bp: number of breakpoint on DAC performance curve (leave it at 2)
  - ConstantOP (1: keep operating parameters constant, 0: allow for optimized operating parameters in each time slice)
  - RHadapt (0: no inlet water spraying, 1: allow inlet water spraying)
  - RHdiscrete: how many discrete point to consider for inlet water spraying

**Solver_Options**: sdpsetting from YALMIP adapted to the model

## Functions
**DAC_Standalone_6**: main function to perform the optimization

**DACOutput_4**: performs the fitting of the performance parameters based on the given RH and T vector (for variable operating parameters)

**DACOutput_4_ConsPar**: performs the fitting of the performance parameters based on the given RH and T vector (for constant operating parameters)

**PLR**: Performs a piece-wise linear regression

