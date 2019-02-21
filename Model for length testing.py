import numpy as np
import matplotlib.pyplot as plt
from uncertainties import unumpy
from uncertainties import ufloat
# ----------------------------------------------------------------------------------------------------------------------

# taking input parameters
m_data = "C:\\Users\\sree118\\AppData\\Local\\MASTER-8.0.0\\output\\KesslerTetherinertial_d.dia"  # MASTER output file
s_data = "C:\\Users\\sree118\\AppData\\Local\\MASTER-8.0.0\\output\\KesslerTetherinertial_d.dia.sigm" # Uncertanities
l_low = 100             # setting upper and lower bounds for the testing radius of the tether
l_high = 2000
l_step_size = 100
test_range = np.arange(l_low, l_high + l_step_size, l_step_size)
tether_radius = 0.3e-3     # in units of m
kc = 0.3                   # lethality coeefiecient -  the ratio of minimum lethal impactor diameter to radius
dt = 1                     # critical distance - ignores glancing blows
t = 8                      # simulated time, in months
t_y = t / 12               # simulated time, in years
alt = "400km"             # altitude of sat - only used for plotting, input to MASTER
# ----------------------------------------------------------------------------------------------------------------------

# defining functions


def a_eff(d, l, crit=dt, w=tether_radius):

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


for l in test_range:
    mask = diam > (kc * tether_radius)                         # defining a mask to only care about lethal impactors
    diam_lethal = diam[mask]
    flux_lethal = flux[mask]
    diam_sigma = diam_lethal[1:] - diam_lethal[:-1]
    flux_lethal_sigma = sigma[mask]

    diams = unumpy.uarray(diam_lethal[:-1].tolist(), diam_sigma.tolist())
    fluxs = unumpy.uarray(flux_lethal.tolist(), flux_lethal_sigma.tolist())

    cum_flux = np.append(cum_flux, np.sum(fluxs))

    areas = a_eff(diams, 1) * l

    diff_rate_impacts = areas * fluxs[:-1]        # the differential rate of impacts, particles per year per meter

    rate_impacts = np.sum(diff_rate_impacts)  # impacts per year

    rate_impacts_data = np.append(rate_impacts_data, rate_impacts)

# ----------------------------------------------------------------------------------------------------------------------

# plotting results

plt.figure()
plt.suptitle(
    "Factors Relevant to Impact Survival given a changing length and constant radius \n"
    "Altitude ={}, Mission Time ={} months, Radius = {} m".format(alt, t, tether_radius))
plt.subplot(2, 2, 1)
plt.title('Flux of Lethal Impactors')
plt.errorbar(test_range, unumpy.nominal_values(cum_flux), yerr=unumpy.std_devs(cum_flux), fmt="b.", ecolor='r')
plt.ylabel('Flux, phi, m^-2 year^-1')
plt.xlabel('Tether Length, l, m')

plt.subplot(2, 2, 2)
plt.title('Lethal Impacts per Mission')
plt.errorbar(test_range, unumpy.nominal_values(rate_impacts_data * t_y), yerr=unumpy.std_devs(rate_impacts_data * t_y), fmt="b.", ecolor='r')
plt.ylabel('Lethal Impacts ')
plt.xlabel('Tether Length, l, m')

plt.subplot(2, 2, 3)
plt.title('Probability of Surviving the Mission')
prob = unumpy.exp(-1 * rate_impacts_data * t_y)
plt.errorbar(test_range, unumpy.nominal_values(prob), yerr=unumpy.std_devs(prob), fmt="b.", ecolor='r')
plt.ylabel('Probability of surviving the mission, P, []')
plt.xlabel('Tether Length, l, m')

plt.subplot(2, 2, 4)
plt.title('Mean survival time')
tau = (1 / rate_impacts_data * t_y) * 12
plt.errorbar(test_range, unumpy.nominal_values(tau), yerr=unumpy.std_devs(tau), fmt="b.", ecolor='r')
plt.xlabel('Tether Length, l, m')
plt.ylabel('Mean survival Time, tau, months')

plt.show()
print(prob)
