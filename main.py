from benchmark import *
from LPM import *

if __name__ == "__main__":
    # This is where all the plot functions should be called
    # Plot the model that we are using
    plot_models()
    # Plot the pressure benchmark
    pressure_benchmark()
    # Plot the temperature benchmark
    temperature_benchmark()
    # Plot the temperature forecast
    temp_forecast()
    # Plot uncertainty
    uncertainty()
