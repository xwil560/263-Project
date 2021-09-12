# ENGSCI 263 Computational Mechanics Group Project

Project 5 Group 15 - Thermal Recovery of Bitumen on the West Coast

We have used data from a pilot study in steam injection for bitumen extraction on the West Coast of the South Island to model and forecast how the temperature and pressure of a control volume would likely behave in a continued thermal recovery project proposed by Todd Energy. Within this folder we have:

- /data: A subdirectory containing the relevant data files (pressure & temperature levels as well as injection & extraction rates from the pilot study needed to fit the model to represent the system.

- data_vis.py: A script to display a time series of a pilot study we have been provided with, including the data used to fit the model. To see each plot, uncomment the data you wish to see i.e.pilot_plots(), exp2_plots() or exp3_plots() and run the file.

- benchmark.py: A script to carry out benchmarking of our method for pressure and temperature. Change the call in main to the benchmark you want to display and run the file to see the benchmarking plots.

- unit_test.py: A script containing multiple unit tests to check for correct outputs and error raises in our pressure and temperature ODE solvers. Run the file to test and if no assertion error appears in terminal, the unit tests have passed.

- LPM.py: This script contains all the function required to read data, solve odes, fit models, plot models, forecast future outcomes, and visualize uncertainty. To see each output, change the call in main to either:

- main.py: This script when run will call the correponding functions in LMP.py to display all our relevant plots, including the final model fit with misfit, benchmarking, forecast and uncertainty. The plots is calls are:

  - 'pilot_plots()': a plot for the time series of our provided data
  - 'plot_models()': a plot of the best fit temperature and pressure LPM models based on the past data provided.
  - 'temp_forecast()': a plot of the 300 day forecast that displays temperature best fit for past data as well as a future temperature prediction on the reservoir for 2 iterations of 4 proposed injection rates (i.e. 460 tonnes/day, 1000 tonnes/day, 250 tonnes/day & 0 tonnes/day).
  - 'uncertainty()': the plot above but with uncertainty illustrated.

  For these functions to work correctly, ensure the '/data' folder is within the same directory as the main.py, data_vis.py, benchmark.py, and unit_test.py files. The data folder should contain five .txt files containing the relevant data, titled:

  - 'tr_oil.txt' for production rate of oil in m^3/day
  - 'tr_p.txt' for pessure of well in kPa
  - 'tr_steam.txt' for injection rate of steam in tonnes/day
  - 'tr_T.txt' for temperature of well in degrees celcius 
  - 'tr_water.txt' for production rate of water in m^3/day

  Additionally, the data should be formatted with a single header, and two columns, the first containing the time(days), the second containing the data, seperated by a comma delimiter.
