import numpy as np
import matplotlib.pyplot as plt
from uncertainties import unumpy
from uncertainties import ufloat
# ----------------------------------------------------------------------------------------------------------------------

# taking input parameters

m_data = "C:\\Users\\sree118\\AppData\\Local\\MASTER-8.0.0\\output\\KesslerTetherinertial_d.dia"  # MASTER output file
s_data = "C:\\Users\\sree118\\AppData\\Local\\MASTER-8.0.0\\output\\KesslerTetherinertial_d.dia.sigm" # Uncertanities
tether_volume = 3.122e-4   # volume available for tether storage
l_step_size = 10          # we iterate over tether length
l_min = 100
l_max = 2000
l_test_range = np.arange(l_min, l_max+ l_step_size, l_step_size)
tether_radius = np.sqrt(tether_volume/ (np.pi * l_test_range))        # generated to conserve volume
kc = 0.3                  # lethality coeefiecient -  the ratio of minimum lethal impactor diameter to radius
dt = 0.7                     # critical distance - ignores glancing blows
t = 2                      # simulated time, in months
t_y = t / 12               # simulated time, in years
alt = "400km"             # altitude of sat - only used for plotting, input to MASTER
# ----------------------------------------------------------------------------------------------------------------------

# defining functions


def a_eff(w, d, l , crit=dt):

    """
        effective area function - w for width of body, l for length, t for critical distance,
        and d for impactor diameter.
        This function accounts for macroscopic size of impactors and discounts glancing blows
    """

    a = l * ((crit * w) + d)
    return a


# ----------------------------------------------------------------------------------------------------------------------

# initialising data arrays

# input data arrays
diam = np.array([])  # impactor diameters, in m
flux = np.array([])  # differential impactor flux at this diameter, in impactors per square meter per year
sigma = np.array([])

# generated data arrays - used for plotting
cum_flux = np.array([])
rate_impacts_data = np.array([])
# ----------------------------------------------------------------------------------------------------------------------

# processing data input

for line in open(m_data):       # taking MASTER input
    li = line.strip()           # removing redundant white space
    if not li.startswith("#"):  # ignoring commented lines - simulation information

        li_array = np.array([float(x) for x in li.split()])  # takes data  for each particle size

        # taking relevant data
        diam = np.append(diam, li_array[0])
        flux = np.append(flux, li_array[-1])


for line in open(s_data):
    li = line.strip()
    if not li.startswith("#"):

        li_array = np.array([float(x) for x in li.split()])
        sigma = np.append(sigma, li_array[-1])

# ----------------------------------------------------------------------------------------------------------------------=

# initialising testing


for i, r in enumerate(tether_radius):  # testing radii of tether

    mask = diam > (kc * r)                         # defining a mask to only care about lethal impactors
    diam_lethal = diam[mask]
    flux_lethal = flux[mask]
    diam_sigma = diam_lethal[1:] - diam_lethal[:-1]
    flux_lethal_sigma = sigma[mask]

    diams = unumpy.uarray(diam_lethal[:-1].tolist(), diam_sigma.tolist())
    fluxs = unumpy.uarray(flux_lethal.tolist(), flux_lethal_sigma.tolist())

    cum_flux = np.append(cum_flux, np.sum(fluxs))

    areas = a_eff(r, diams, l_test_range[i])

    diff_rate_impacts = areas * fluxs[:-1]        # the differential rate of impacts, particles per year per meter

    rate_impacts = np.sum(diff_rate_impacts)  # impacts per year

    rate_impacts_data = np.append(rate_impacts_data, rate_impacts)

# ----------------------------------------------------------------------------------------------------------------------

# plotting results

plt.figure()
plt.suptitle(
    "Factors Relevant to Impact Survival given a Constant Volume and Changing Dimensions \n"
    "Altitude ={}, Mission Time ={} months, Volume = {} m^3".format(alt, t, tether_volume))
plt.subplot(2, 2, 1)
plt.title('Flux of Lethal Impactors')
plt.errorbar(l_test_range, unumpy.nominal_values(cum_flux), yerr=unumpy.std_devs(cum_flux), fmt="b.", ecolor='r')
plt.ylabel('Flux, phi, m^-2 year^-1')
plt.xlabel('Tether Length, l, m')

plt.subplot(2, 2, 2)
plt.title('Lethal Impacts per Mission')
plt.errorbar(l_test_range, unumpy.nominal_values(rate_impacts_data * t_y), yerr=unumpy.std_devs(rate_impacts_data * t_y), fmt="b.", ecolor='r')
plt.ylabel('Lethal Impacts ')
plt.xlabel('Tether Length, l, m')

plt.subplot(2, 2, 3)
plt.title('Probability of Surviving the Mission')
prob = unumpy.exp(-1 * rate_impacts_data * t_y)
plt.errorbar(l_test_range, unumpy.nominal_values(prob), yerr=unumpy.std_devs(prob), fmt="b.", ecolor='r')
plt.ylabel('Probability of surviving the mission, P, []')
plt.xlabel('Tether Length, l, m')

plt.subplot(2, 2, 4)
plt.title('Mean survival time')
tau = (1 / rate_impacts_data * t_y) * 12
plt.errorbar(l_test_range, unumpy.nominal_values(tau), yerr=unumpy.std_devs(tau), fmt="b.", ecolor='r')
plt.xlabel('Tether Length, l, m')
plt.ylabel('Mean survival Time, tau, months')

plt.show()
print(prob)
