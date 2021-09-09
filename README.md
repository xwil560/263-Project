# 263-Project
263 Group Project 

We have used data from a pilot study in steam injection for bitumen extraction in the South Island to model and so forecast how the temperature and pressure of a control volume would likely behave in a new injection project proposed by Todd Energy. Within this folder we have:
     -data: a subdirectory containing the relevant data files (pressure, temperature, injection, and production rates)from the pilot study needed to fit the model to represent the system.

     -Data_vis.py: a script to display and data we have been provided with, including that used to fit the model, run to display this. To see each plot, in main uncomment the data you with to see, pilot_plots(), exp2_plots() or exp3_plots(). 

     -Benchmark.py: a script to carry out benchmarking for our method for pressure and temperature. change the call in main to the benchmark you want to display and run.

     -Unit_test.py: a script containing multiple unit test to check for correct outputs and error raises in our pressure and temperature solvers. Run to test.

     -main.py: This script contaitains all the function required to read data, solve odes, fit models, plot models, forecast future outcomes, and visualise uncertainty. To see each output, change the call in main to either : 
      - 'temp_forcast()' to see a plot of out final graph of temerature fit, and prediction with each different injection rate.
      - 'uncertainty()', to see the plot above but with uncertainty illustrated.
      - 'plot_models()'. to see a plot of out best fit temperature and pressure LMP models on the data provided.

    
    For these functions to work correctly, ensure the 'data' folder is within the same directory level as the main.py, data_vis.py, benchmark.py, and unit_test.py files. The data folder should contain five txt files containing the relevant data, titled:
        -tr_oil.txt
        -tr_p.txt
        -tr_steam.txt
        -tr_T.txt
        -tr_water.txt

    Additionally, the data should be formated with a single header, and two columns, the first containing the data, the second containing the time (days), seperated by a comma deliminator. 


