#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())

module main_run

  using initialCondition,modifyFiles,computeSolution

  # Change parameters here ##################################
  nx = ny = 100; # lattice size

  rinit = 10.0; # initial skyrmion radius
  drinit = 0.0; # initial distance b/w skyrmion antiskyrmion

  # If you would like to input arguments from terminal, change 
  # desired inputs to parse(Float64,ARGS[1]), 
  # parse(Float64,ARGS[2]), etc. This parses the input argument
  j = 1.0;          # Exchange
  h = 0.01;         # External field
  dmi = 0.02;       # Dzyaloshinskii-Moriya
  pma = 0.0;        # Perpendicular magnetic anisotropy
  ed = 0.0;         # Dipole-dipole
  pbc = 1.0;        # Have to set pbc to false if using dipole-dipole

  tfree = 2.0;      # Time interval of RK
  nSteps = 10;      # Number of steps to take in RK solver 
  damping = 0.01    # LLG damping
  temp = 0.1        # Temperature
  ic = "skyrmion"

  # end of change parameter section
  #############################################################

  # Divide tfree and nSteps by 2 because in pulse noise we
  # compute half the time evolution at a time
  tRK = tfree/2; nRK = nSteps/2; 

  params = [j, h, dmi, pma, ed, pbc];
  evalParams = [tRK, nRK, damping, temp];

  # # Use the next three lines if you want to evaluate relaxation
  # # before pulse-noise dynamics.
  # # When running many iterations of pulse noise, it is better to
  # # run relaxation once, then import the result as a data file as
  # # done in the next chunk of code
  s0 = buildInitial(ic, rinit, drinit, nx, ny);
  runRelaxation!(s0, params, [1.0, 10, 1.0, 0.0])


  # If relaxation has already been completed and data has been saved,
  # comment out the previous runRelaxation!() and use the getDataH5()
  # command below
  # s0 = getDataH5(string("/data/spin_field_after_relaxation_T=0.0_H=",h,
  #   "_DMI=",dmi,"_DDI=",ed,"_PMA=",pma,"_.h5"))
  @time evaluateLL!(s0, params, evalParams)



end
