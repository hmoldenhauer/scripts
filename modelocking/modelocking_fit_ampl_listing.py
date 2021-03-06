# Script for fitting electron and hole spin rotation signals
# for mode-locking signals recorded by pump-probe experiments

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit
from matplotlib.ticker import FuncFormatter

# read file
def read_file(name):
    # the variables originate from the experiment
    mm, X_Calibrated, Y_Calibrated, ADC_1, ADC_2 = np.genfromtxt(name, unpack=True) 

    # change units to SI units
    m = mm * 1e-3
    amplitude = X_Calibrated * 1e-6

    # set peak position as zero
    #                           probably some work to do to find the correct peak
    # fist get index of peak
    max_ampl_index = np.argmax(amplitude)
    m_0 = m[max_ampl_index]
    # now subtract m_0 
    m = m - m_0

    # need to display times and not distances
    # -> calculate times by using the speed of light
    t =  m / const.c
    # as the light has to go back and forth so multiply by two
    t = 2 * t

    return t, amplitude, max_ampl_index

# cos with exponential decay
def f(t, ampl, omega, phase, tau):
    return ampl * np.cos(omega * t + phase) * np.exp(-t**2 / (2 * tau**2))

# sum function to fit electron and hole precision at
def f_eh(t, A_e, A_h, omega_e, omega_h, phi_e, phi_h, tau_e, tau_h, y):
    return f(t, A_e, omega_e, phi_e, tau_e) + f(t, A_h, omega_h, phi_h, tau_h) + y

# function to fit mode-locking at
def f_ml(t, ampl, omega, phase, tau, y):
    return f(t, ampl, omega, phase, tau) + y

#################################################
# Programm
#################################################

# insert all the files you want to work with
# for many files you might want to use:
#       ls *.dat >> file_list.txt
file_list = open("file_list.txt")
array_of_files = file_list.read().splitlines()

# file fit-data is written to
fit_file = open("ampl.txt", "w")

for data in array_of_files:

    # name of the data file
    data_name = data
    print()
    print(data_name)

    # create name for fit file
    fit_name = "fit/fit-" + data_name[data_name.index("/") + 1:data_name.index(".")] + ".pdf"
    print(fit_name)

    # do you want to fit?
    fit = True

    # do you want to plot starting parameter function of exponential sin after pulse?
    sum_ = False

    # do you want to fit the signal?
    eh = False

    # do you want do mode-locking?
    ml = True

    # do you want to plot only the interesting section of ml?
    print_ml = False

    # read file and get some parameters
    t, amplitude, max_index = read_file(data_name)

    # use these parameters to first find good starting parameters
    # for electron and hole frequencies, their relexation times
    # and their amplitude
    if not fit:
        x = np.linspace(0, 5e-10, 2000)
        ampl = 4.5e-4
        omega = 1e10
        phase = 0
        tau = 3e-10
        y = 2e-4
        if ml:
            # use these parameters to find good starting parameters
            # for hole mode-locking
            x_ml = np.linspace(-2e-9, 0, 2000)
            ampl_ml = 4e-6
            omega_ml = 1e9
            phase_ml = 0
            tau_ml = 4e-10
            y_ml = 2e-5

    # offsetparameter to erase left peak shoulder
    offset = 3

    # need for calculated starting values!!!
    # calculate minimum and maximum of modelocking
    ml_min = np.argmin(amplitude[:max_index - offset])
    ml_max = np.argmax(amplitude[:max_index - offset])

    ampl_ml_cal = amplitude[ml_max] + amplitude[ml_min]
    omega_ml_cal = 1 / abs(t[ml_max] - t[ml_min])

    print("initial guess:")
    print("ampl\t{:.3e}".format(ampl_ml_cal))
    print("omega\t{:.3e}".format(omega_ml_cal))

    # define starting parameters
    start_params = np.array([1e-6,  # electron amplitude
                            1e-6,   # hole amplitude
                            6e10,   # electron frequency
                            1.5e10, # hole frequency
                            1,      # electron phase
                            1,      # hole phase
                            4e-10,  # electron relaxation time
                            2e-10,  # hole relaxation time
                            0])     # y axis intersept

    # define starting parameters for mode-locking
    start_params_ml = np.array([ampl_ml_cal,  # hole amplitude
                               omega_ml_cal,  # hole frequency
                               0,             # hole phase
                               4e-10,         # hole relaxation time
                               2e-5])         # y axis intersept

    # fit function to data
    if fit:
        if eh:
            # fit signal
            params, cov = curve_fit(f_eh, t[max_index:], amplitude[max_index:], start_params)
            print()
            print("signal parameters")
            print("\te_params\th_params\tunit")
            print("ampl\t{:.3e}\t{:.3e}\tV".format(params[0], params[1]))
            print("freq\t{:.3e}\t{:.3e}\tHz".format(params[2], params[3]))
            print("phase\t{:.3e}\t{:.3e}\t".format(params[4], params[5]))
            print("relax\t{:.3e}\t{:.3e}\ts".format(params[6], params[7]))
        if ml:
            # fit mode-locking
            params_ml, cov_ml = curve_fit(f_ml, t[:max_index-offset], amplitude[:max_index-offset], start_params_ml)
            print()
            print("ml parameters")
            print("\th_params\tunit")
            print("ampl\t{:.3e}\tV".format(params_ml[0]))
            print("freq\t{:.3e}\tHz".format(params_ml[1]))
            print("phase\t{:.3e}\t".format(params_ml[2]))
            print("relax\t{:.3e}\ts".format(params_ml[3]))

            # write amplitude to ampl.txt
            fit_file.write("{:.3}\n".format(params_ml[0]))

    #################################################
    # plotting section
    #################################################

    # create figure and axes
    fig, ax = plt.subplots()

    # plot data
    if not print_ml:
        plt.plot(t, amplitude, "rx-", label=r"$\mathrm{Messwerte}$")
    else:
        plt.plot(t[:max_index-offset], amplitude[:max_index-offset], "rx-", label=r"$\mathrm{Messwerte}$")

    # plot functions to get starting parameters
    if not fit:
        if eh:
            plt.plot(x, f_ml(x, ampl, omega, phase, tau, y), "b", label=r"$\mathrm{Exponentialer}\ \sin$")
        if sum_:
            plt.plot(x, f_eh(x, *start_params), "b", label=r"$\mathrm{Summenfunktion}$")
        if ml:
            plt.plot(x_ml, f_ml(x_ml, ampl_ml, omega_ml, phase_ml, tau_ml, y_ml), "g", label=r"$\mathrm{Mode-Locking}$")

    # plot fit
    if fit:
        if eh:
            # plot signal after peak
            plt.plot(t[max_index:], f_eh(t[max_index:], *params), "b-", label=r"$\mathrm{Fit}$")
        if ml:
            # plot mode-locking fit
            plt.plot(t[:max_index], f_ml(t[:max_index], *params_ml), "g-", label=r"$\mathrm{Mode-Locking\ Fit}$")

    # label axis 
    plt.xlabel(r"$\mathrm{Pump-Probe\ Verzögerung} / \mathrm{ns}$")
    plt.ylabel(r"$\mathrm{Amplitude} / \mathrm{\mu V}$")

    # format axes to fit the data
    ax.xaxis.set_major_formatter(FuncFormatter(lambda t, pos: ('%.1f')%(t*1e9)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda t, pos: ('%.1f')%(t*1e6)))

    # plot legend
    plt.legend(loc="best")

    # some plot adjustments and saving the plot
    plt.tight_layout()
    plt.savefig(fit_name)

# close file
fit_file.close()
