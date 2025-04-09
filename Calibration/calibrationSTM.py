import numpy as np
import pandas as pd
import statsmodels.formula.api as sm
import sys

filename = "/mu2e/app/users/sisrael/adcPeaks.log"

'''
0.511 MeV is the beam flash
0.392 MeV is Zn113
0.662 MeV is the Cs137
0.898 MeV is the Y88
1.173 and 1.333 MeV are the Co60
'''
#energy_peaks_beam = np.array([0.511, 0.662,0.898,1.173,1.333])
energy_peaks_nobeam = np.array([0.392, 0.662,0.898,1.173,1.333])
run_num, detector = sys.argv[1:]

def chisquared(f_obs,f_exp,ndof, sigma):
    """
    Calculates the ratio between the chi squared value versus the number of
    degrees of freedom. Should be close to 1 for a good fit.

    Parameters
    ----------
    f_obs : Numpy array of floats or doubles
        The measured values from the fit means
    f_exp : Numpy array of floats or doubles
        The predicted values from inputting the energies into a linear equation from the regression parameters
    ndof : int
        Number of degrees of freedom, total measurements - number of fit parameters
    sigma : Numpy array of floats or doubles
        The error on the f_obs, in our case the error on the fit means, not the standard deviation

    Returns
    -------
    float
        chi**2/ndof
    """
    numerator = (f_obs - f_exp)**2
    denominator = sigma**2

    chi2 = np.sum(numerator/denominator)
    return chi2/ndof

def digi_regression(filename):
    data = pd.read_csv(filename, header=0)
    data[' par_name'].str.strip('')

    peaks = data[data[' par_name'] == ' Mean']
    adc_data = peaks[[' par_value', ' par_error']].sort_values(' par_value')
    adc_peaks = adc_data[' par_value']
    peak_error = adc_data[' par_error']

    data = {'adc' : adc_peaks, 'energy' : energy_peaks_nobeam, 'error' : peak_error}
    
    # Regression fit
    model = sm.ols(" adc~energy", data=data).fit()

    return model, adc_peaks, peaks, peak_error
# Linear function to take fitted parameters and plot regression line
def adc_pred(x, model):
    return x*model.params[1] + model.params[0]

def fit_result(model, adc_peaks, adc_predicted, peak_error):
    deg_freedom = len(adc_predicted) - len(model.params)
    chi2 = chisquared(adc_peaks, adc_predicted, deg_freedom, peak_error)
    gain = 1.0/model.params[1]
    cal = -model.params[0]/model.params[1]
    print(adc_predicted)
    print('Fit results: ', model.params)
    print('Gain: %3.3e' % (gain))
    print('Calibration: %3.3e' % (cal))
    print('R2 Value: ', model.rsquared)
    print('Chi2/dof: ', chi2)
    return gain, cal

def export_params(gain, cal, run_number):
    f = open('/mu2e/app/users/sisrael/STMCalibration/stmTest.txt', 'w+')
    f.write("TABLE STMEnergyPar " + str(run_number) + "\n")
    f.write("# E = p0 + p1*adc + p2*adc*adc\n")
    f.write("# ch,p0,p1,p2\n")
    f.write("HPGe," + str(cal) + "," + str(gain) + "," + str(0) + '\n')
    f.write("LaBr," + str(cal) + "," + str(gain) + "," + str(0) + '\n')
    f.close()
#%%
#####################################################################
model, adc_peaks, peaks, peak_error = digi_regression(filename)
adc_predicted = np.array([adc_pred(i, model) for i in energy_peaks_nobeam])

high_gain, high_cal = fit_result(model, adc_peaks, adc_predicted, peak_error)

'''
plt.errorbar(energy_peaks_nobeam, adc_peaks, yerr=peak_error, fmt='.', label='Data')
plt.plot(np.linspace(0,1.5,5000), np.linspace(0,1.5,5000)*model.params[1] + model.params[0], label='Fit', color='red')
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('ADC Number')
plt.grid()
plt.show()
'''
#%%
##############################################################################
export_params(high_gain, high_cal, run_number)
