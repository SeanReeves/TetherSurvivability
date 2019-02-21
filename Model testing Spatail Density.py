import numpy as np
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------------------------------------------------------

# taking input parameters

m_data = "C:\\Users\\sree118\\AppData\\Local\\MASTER-8.0.0\\output\\KesslerTetherSpatial_d.sdm"            # MASTER output file
d_low = 0.1e-3             # setting upper and lower bounds for the testing radius of the tether
d_high = 1.e-3
d_step_size = 1.e-4
test_range = np.arange(d_low, d_high, d_step_size)
tether_length = 1000        # in units of m
kc = 0.3                   # lethality coeefiecient -  the ratio of minimum lethal impactor diameter to radius
dt = 1                     # critical distance - ignores glancing blows
t = 8                      # simulated time, in months
t_y = t / 12               # simulated time, in years
alt = 400e3                  # altitude of sat
mu = 3.986e14              # standard gravitational parameter of Earth
er = 6371e3                # radius of Earth
# ----------------------------------------------------------------------------------------------------------------------

# defining functions


def a_eff(w, d, l=tether_length, crit=dt):

    """
        effective area function - w for width of body, l for length, t for critical distance,
        and d for impactor diameter.
        This function accounts for macroscopic size of impactors and discounts glancing blows
    """

    a = l * ((crit * w) + d)
    return a


def find_v(altitude=alt, mu=mu, earthrad= er):
    v = np.sqrt(mu/(altitude + earthrad))
    return v

x = find_v()
# ----------------------------------------------------------------------------------------------------------------------

# initialising data arrays

# input data arrays
diam = np.array([])  # impactor diameters, in m
flux = np.array([])  # differential impactor flux at this diameter, in impactors per square meter per year

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

# ----------------------------------------------------------------------------------------------------------------------=

# initialising testing


for r in test_range:  # testing radii of tether

    mask = diam > (kc * r)                         # defining a mask to only care about lethal impactors
    diam_lethal = diam[mask]
    flux_lethal = flux[mask] * find_v() / 1000
    diam_intervals_relevant = 1#diam_lethal[1:] - diam_lethal[:-1]

    cum_flux = np.append(cum_flux, np.sum(flux_lethal))

    areas = a_eff(r, diam_lethal)
    print(areas)
    diff_rate_impacts = areas * flux_lethal        # the differential rate of impacts, particles per year per meter

    rate_impacts = np.sum(diam_intervals_relevant * diff_rate_impacts[:-1])  # impacts per year

    rate_impacts_data = np.append(rate_impacts_data, rate_impacts)

# ----------------------------------------------------------------------------------------------------------------------

# plotting results

plt.figure()
plt.suptitle(
    "Factors relevant to survival impact as a Function of Tether Diameter \n"
    "Altitude ={}, Mission Time ={} months, Length = {} m".format( alt , t, tether_length))
plt.subplot(2, 2, 1)
plt.title('Flux of Lethal Impactors')
plt.scatter(test_range, cum_flux)
plt.ylabel('Flux, phi, m^-2 year^-1')
plt.xlabel('Tether Diameter, D, m')

plt.subplot(2, 2, 2)
plt.title('Lethal Impacts per Mission')
plt.scatter(test_range, rate_impacts_data * t_y)
plt.ylabel('Lethal Impacts ')
plt.xlabel('Tether Diameter, D, m')

plt.subplot(2, 2, 3)
plt.title('Probability of Surviving the Mission')
prob = np.exp(-1 * rate_impacts_data * t_y)
plt.scatter(test_range, prob)
plt.ylabel('Probability of surviving the mission, P, []')
plt.xlabel('Tether Diameter, D, m')

plt.subplot(2, 2, 4)
plt.title('Mean survival time')
tau = (1 / rate_impacts_data * t_y) * 12
plt.scatter(test_range, tau)
plt.xlabel('Tether Radius, r, m')
plt.ylabel('Mean survival Time, tau, months')

plt.show()
print(prob)
