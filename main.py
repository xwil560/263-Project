from benchmark import *
from LPM import *
from data_vis import *

if __name__ == "__main__":
    # This is where all the plot functions should be called

    # Plotting the time series data for the pilot study
    pilot_plots()
    # Plot the models for pressure and temperature ODEs on the historical data
    plot_models()
    # Plot the pressure benchmark
    pressure_benchmark()
    # Plot the temperature benchmark
    temperature_benchmark()
    # Plot the temperature forecast for different injection rates
    temp_forecast()
    # Plot the forecast, factoring in uncertainty
    uncertainty()
