# 263-Project

263 Group Project

We have used data from a pilot study in steam injection for bitumen extraction on the West Coast of the South Island to model and forecast how the temperature and pressure of a control volume would likely behave in a continued thermal recovery project proposed by Todd Energy. Within this folder we have:

- /data: A subdirectory containing the relevant data files (pressure & temperature levels as well as injection & extraction rates from the pilot study needed to fit the model to represent the system.

- data_vis.py: A script to display a time series of a pilot study we have been provided with, including the data used to fit the model. To see each plot, uncomment the data you wish to see i.e.pilot_plots(), exp2_plots() or exp3_plots() and run the file.

- benchmark.py: A script to carry out benchmarking of our method for pressure and temperature. Change the call in main to the benchmark you want to display and run the file to see the benchmarking plots.

- unit_test.py: A script containing multiple unit tests to check for correct outputs and error raises in our pressure and temperature ODE solvers. Run the file to test and if no assertion error appears in terminal, the unit tests have passed.

- main.py: This script contains all the function required to read data, solve odes, fit models, plot models, forecast future outcomes, and visualize uncertainty. To see each output, change the call in main to either:

- 'plot_models()': To see a plot of the best fit temperature and pressure LPM models based on the past data provided.
- 'temp_forecast()': To see a plot of the 300 day forecast that displays temperature best fit for past data as well as a future temperature prediction on the reservoir for 2 iterations of 4 proposed injection rates (i.e. 460 tonnes/day, 1000 tonnes/day, 250 tonnes/day & 0 tonnes/day).
- 'uncertainty()': To see the plot above but with uncertainty illustrated.

For these functions to work correctly, ensure the '/data' folder is within the same directory as the main.py, data_vis.py, benchmark.py, and unit_test.py files. The data folder should contain five .txt files containing the relevant data, titled:

- 'tr_oil.txt'
- 'tr_p.txt'
- 'tr_steam.txt'
- 'tr_T.txt'
- 'tr_water.txt'

Additionally, the data should be formated with a single header, and two columns, the first containing the data, the second containing the time (days), seperated by a comma deliminator.
