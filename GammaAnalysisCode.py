# 12/02/2021
# Fada Guan
# modified from the example of gamma calculation for dose dicom files
# https://docs.pymedphys.com/lib/howto/gamma/from-dicom.html

import numpy as np
import csv
import matplotlib.pyplot as plt
import pymedphys
import pandas as pd

if __name__ == '__main__':
    print('1D gamma index calculation for input ASCII dose files')

#data1 = pd.read_excel("Data collection_update.xlsx", sheet_name="PDD 2x2")


data1 = np.genfromtxt('OCR5_10x10.csv',delimiter=',',skip_header=1)

axes_Measure = data1[:, 0] #1st column is x in mm
dose_Measure = data1[:, 1] #2nd column is dose in Gy/MU

axes_GBD = data1[:, 0]
dose_GBD = data1[:, 2]


gamma_options = {
    'dose_percent_threshold': 1,
    'distance_mm_threshold': 1,
    'lower_percent_dose_cutoff': 0,
    'interp_fraction': 10,  # Should be 10 or more for more accurate results
    'max_gamma': 2,
    'random_subset': None,
    'local_gamma': False,
    'ram_available': 2 ** 29  # 1/2 GB
}
#for global normalization, the maximum reference dose is used
#but in TG218, it said usually the prescribed dose or maximum dose in plan (evaluation) is used
gamma = pymedphys.gamma(
    axes_Measure, dose_Measure,
    axes_GBD, dose_GBD,
    **gamma_options)




####################################################################################################################
valid_gamma = gamma[~np.isnan(gamma)]
print('# of Measure points with valid gamma {0}'.format( len(valid_gamma)) )
num_bins = (
    gamma_options['interp_fraction'] * gamma_options['max_gamma'])
bins = np.linspace(0, gamma_options['max_gamma'], num_bins + 1)

flagProbDensity = True #FG
plt.hist(valid_gamma, bins, density=flagProbDensity) #y value is probability density in each bin
#plt.hist(valid_gamma, bins, density=False) #y value is counts in each bin
plt.xlim([0, gamma_options['max_gamma']])
plt.xlabel('gamma value') #FG
if flagProbDensity: #FG
    plt.ylabel('probability density')
else:
    plt.ylabel('counts')

pass_ratio = np.sum(valid_gamma <= 1) / len(valid_gamma)

if gamma_options['local_gamma']:
    gamma_norm_condition = 'Local gamma'
else:
    gamma_norm_condition = 'Global gamma'

plt.title(f"dose cut: {gamma_options['lower_percent_dose_cutoff']}%|{gamma_norm_condition} ({gamma_options['dose_percent_threshold']}%/{gamma_options['distance_mm_threshold']}mm)|Pass Rate(\u03B3<=1): {pass_ratio*100:.2f}%")
plt.vlines(x=[1], ymin=0, ymax=1, colors='purple', ls='--', lw=2,label='target')

#plt.title(f"Local Gamma ({gamma_options['dose_percent_threshold']}%/{gamma_options['distance_mm_threshold']}mm) | Percent Pass: {pass_ratio*100:.2f} %")
plt.savefig('1D_{0}_histogram.png'.format(gamma_norm_condition), dpi=300) #savefig must be before show()

####################################################################################################################
#plot gamma index all Measure points
max_Measure_dose = np.max(dose_Measure) #max Measure dose
max_GDB_dose = np.max(dose_GBD) #max GDB dose

lower_dose_cutoff = gamma_options['lower_percent_dose_cutoff'] / 100 * max_Measure_dose

#because the Measure dose data and GDB dose data may not have the same size and matched coordinates
#it is not convenient to directly get the dose difference unless we do interpolation
#therefore, we only plot Measure dose data, GDB dose data, and gamma index of each Measure dose point
#but note: in the gamma index calculation, the gamma function in PyMedPhys indeed interpolates the dose data
figure_2 = plt.figure(2, figsize=(8, 6), dpi=120, facecolor='w', edgecolor='k')
figure_2.suptitle('Measure and GDB dose 1D, and {0} index of Measure points'.format(gamma_norm_condition),fontsize=12)

ax_1 = figure_2.add_subplot(111)
ax_1.tick_params(which='minor',direction='in')  # , length=6, width=2, colors='b')
ax_1.tick_params(axis='x', bottom='on', top='on')  # show x ticks for top and bottom
ax_1.tick_params(labeltop='on')
ax_1.minorticks_on()  # if drawing speed is an issue, please use _off()
ax_1.set_xlabel('off-axis distance (mm)')

ax_2 = ax_1.twinx()
ax_2.minorticks_on()
ax_2.tick_params(which='minor',direction='in')

curve_0 = ax_1.plot(axes_Measure, dose_Measure,'x', label='Measure dose') #Measure dose 1D, 'k-' is black solid line
curve_1 = ax_1.plot(axes_GBD, dose_GBD,'bo', markersize=3, label='GBD dose') #GBD dose 1D, 'b--' is blue dashed line, 'bo' means blue o marker

curve_2 = ax_2.plot(axes_Measure, gamma,'r*', label=f"gamma ({gamma_options['dose_percent_threshold']}%/{gamma_options['distance_mm_threshold']}mm)") #gamma index of each Measure point
#'r*' is red star marker

ax_1.set_ylabel('profile (%)')
ax_2.set_ylabel('gamma index')

ax_1.set_ylim([0, max(max_Measure_dose,max_GDB_dose) * 1.1])
ax_2.set_ylim([0, gamma_options['max_gamma'] * 2.0 ])

curves = curve_0 + curve_1 + curve_2  # put curves together
labels = [l.get_label() for l in curves]  # find labels(in legend) for all curves
ax_1.legend(curves, labels, loc='best')  # show legend for first sub plot
ax_1.grid(True)

# save figure first and show it
figureName = '1D dose_Measure_GBD_{0} index.png'.format(gamma_norm_condition)
plt.savefig(figureName)  # save figure must be before plt.show()
plt.show()

#save gamma index in a .csv file
gammaDataFile = '1D {0} index.csv'.format(gamma_norm_condition)
aDataList = []
aDataList.append(axes_Measure)
aDataList.append(gamma)
aLabelList = []
aLabelList.append('x(mm)')
aLabelList.append('gamma index')

# save figure data to csv file
np.savetxt(gammaDataFile, np.transpose(aDataList), delimiter=',')

with open(gammaDataFile, newline='') as f:  # read in data
    r = csv.reader(f)
    data = [line for line in r]
with open(gammaDataFile, 'w', newline='') as f:  # write header strings and data
    w = csv.writer(f)
    w.writerow(aLabelList)
    w.writerows(data)




####################################################################################################################
