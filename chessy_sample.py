


# =============================================================================
# Curve fitting tutorial
# =============================================================================
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.special import erf #importing the error function
plt.style.use('seaborn-poster')
from scipy.optimize import curve_fit #importing the curve fit function



#Importing the data

file_0V = pd.read_csv('0V_chessy_momentum m3.0_400 nm_FA 100_CA open_ext 10k_bias 0V_200 s_1_003 (2).csv' ) #importing a file having OV as biasing
file_1V = pd.read_csv('1V_chessy_momentum_m3.0_400_nm_FA_100_CA_open_ext_10k_bias_1V_200_s_005.csv') #importing a file having 1V as biasing
file_neg_1V = pd.read_csv('-1V_chessy_momentum m3.0_400 nm_FA 100_CA open_ext 10k_bias -1V_200 s_1_009 (2).csv') #importing a file having -1V as biasing
file_2V = pd.read_csv('2V_chessy_momentum m3.0_400 nm_FA 100_CA open_ext 10k_bias 2V_200 s_006.csv') #importing a file having 2V as biasing
file_neg_2V = pd.read_csv('-2V_chessy_momentum m3.0_400 nm_FA 100_CA open_ext 10k_bias -2V_200 s_1_012.csv') #importing a file having -2V as biasing


x= np.linspace(350,378,512 ) # the slice represents the time in ns






# =============================================================================
# ploting 0V biasing voltage
# =============================================================================




y_center_0V = (file_0V['Mean_center'])*np.pi*15*15*20  # the time represents the mean count of central part
y_full_0V = (file_0V['Mean_full'])*750*750  # the time represents the mean count of central part


plt.figure()  #  second plot
plt.xlabel('time of flight (ns)')
plt.ylabel('photoemission intensity (counts)')
plt.title('OV ')
plt.plot(x, y_full_0V, 'o', label='full')

plt.plot(x, y_center_0V,'o' ,color='gray', label='center')
plt.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

# define the function to fit the curve for full curve (left side of the graph which denotes the fermi edge)
x_data_full = x[(x>354) & (x<357.7) ]

y_data_full = y_full_0V[(x>354) & (x<357.7)]

def gauss_f(x,A, mu, sig): #gauss_f denotes the gaussian function having center at mu, sig denotes width, A denotes the amplitude
    return A*np.exp(-(x-mu)**2/sig**2)


#  Fit the gaussian function for full integrated curve
popt_left, pcov_left = curve_fit(gauss_f, x_data_full, y_data_full, p0=(175000,357,4))


# Plot fitted left curve
plt.plot(x_data_full, gauss_f(x_data_full, *popt_left), 'r-')


# Step : Define the function for curve fitting of central curve (right side of curve which defines the vacuum edge)
x_data = x[(x>358.4) & (x<360) ]
y_data = y_center_0V [(x>358.4) & (x<360)]


def func(x, a, b, c, d): # definng parameters
    return a * erf(b * (x -c)) + d #erf denotes the error function, a denotes the amplitude,


#  Fit the error function for center curve
popt_right, pcov_right = curve_fit(func, x_data, y_data, p0=(80000,0.05,360,66)) # 80000 is the maximum value  and 66 is the lower value of the curve fitting function in y cordinate. 360 is the maximum value in x cordinate.



# Plot fitted right curve
plt.plot(x_data, func(x_data, *popt_right), 'r-', label='fit')

#midpoint of vertical line for center curve
V_line = func(x_data, *popt_right)

# Calculate the center y value
center_y = (np.max(V_line) + np.min(V_line)) / 2

# Find the nearest points in y to the center_y
nearest_points = np.abs(V_line - center_y)
index_lower = np.argmin(nearest_points)
index_upper = index_lower + 1 if index_lower < len(V_line) - 1 else index_lower - 1

# Interpolate to find the x coordinate corresponding to center_y
x_lower = x_data[index_lower]
x_upper = x_data[index_upper]
y_lower = V_line[index_lower]
y_upper = V_line[index_upper]

slope = (y_upper - y_lower) / (x_upper - x_lower)
center_x = x_lower + (center_y - y_lower) / slope


 
plt.legend()
plt.show()

# =============================================================================
# ploting 1V biasing voltage
# =============================================================================

y_center_1V = (file_1V['Mean_center'])*np.pi*15*15*20 # the time represents the mean count of central part
y_full_1V = (file_1V['Mean_full'])*750*750  # the time represents the mean count of central part







plt.figure()  #  third plot
plt.xlabel('time of flight (ns)')
plt.ylabel('photoemission intensity (counts)')
plt.title('+1 V')

plt.plot(x, y_full_1V, 'o', label='full')

plt.plot(x, y_center_1V, 'o', color='gray' ,label= ' center')






#plt.ticklabel_format(style='sci', scilimits=(5,7))

# define the function to fit the curve for full curve (left side of the curve which denotes the fermi edge)
x_data_full = x[(x>356) & (x<360) ]

y_data_full = y_full_1V[(x>356) & (x<360)]

def gauss_f(x,A, mu, sig):
    return A*np.exp(-(x-mu)**2/sig**2)


#  Fit the error function for full integrated curve
popt_left, pcov_left = curve_fit(gauss_f, x_data_full, y_data_full, p0=(140000,360,4))


# Plot fitted right curve
plt.plot(x_data_full, gauss_f(x_data_full, *popt_left), 'r-')



# Step : Define the function for curve fitting of central curve (right side of curve which defines the vacuum edge)
x_data = x[(x>360.7) & (x<361.8) ]
y_data = y_center_1V [(x>360.7) & (x<361.8)]
x2 = np.linspace(min(x_data), max(x_data),1000)

def func(x, a, b, c, d):
    return a * erf(b * (x -c)) + d


#  Fit the error function for center curve
popt_right, pcov_right = curve_fit(func, x_data, y_data, p0=(60000,0.05,361,59))



# Plot fitted right curve
plt.plot(x_data, func(x_data, *popt_right), 'r-', label='fit')

#midpoint of vertical line for center curve
V_line = func(x_data, *popt_right)

# Calculate the center y value
center_y = (np.max(V_line) + np.min(V_line)) / 2

# Find the nearest points in y to the center_y
nearest_points = np.abs(V_line - center_y)
index_lower = np.argmin(nearest_points)
index_upper = index_lower + 1 if index_lower < len(V_line) - 1 else index_lower - 1

# Interpolate to find the x coordinate corresponding to center_y
x_lower = x_data[index_lower]
x_upper = x_data[index_upper]
y_lower = V_line[index_lower]
y_upper = V_line[index_upper]

slope = (y_upper - y_lower) / (x_upper - x_lower)
center_x = x_lower + (center_y - y_lower) / slope



#major and minor ticks on x axis
plt.minorticks_on()
plt.tick_params(axis='x', which='major',length=10)
plt.tick_params(axis='x', which='minor',length=5)

plt.xlim(352, 365)
plt.legend()



plt.show()






# =============================================================================
# ploting -1V biasing voltage
# =============================================================================

y_center_neg_1V = file_neg_1V['Mean_center']*np.pi*15*15*20  # the time represents the mean count of central part
y_full_neg_1V = file_neg_1V['Mean_full']*750*750  # the time represents the mean count of central part


plt.figure()  #  fourth plot
plt.xlabel('Time of flight (ns)')
plt.ylabel('Counts')

plt.plot(x, y_full_neg_1V, 'o', label='full')

plt.plot(x, y_center_neg_1V, 'o', color='gray', label= ' center')

plt.title('-1 V')
# define the function to fit the curve for full curve (left side of the curve which denotes the fermi edge)
x_data_full = x[(x>352) & (x<355.6) ]

y_data_full = y_full_neg_1V[(x>352) & (x<355.6)]

def gauss_f(x,A, mu, sig):
    return A*np.exp(-(x-mu)**2/sig**2)


#  Fit the error function for full integrated curve
popt_left, pcov_left = curve_fit(gauss_f, x_data_full, y_data_full, p0=(20000,355.5,4))


# Plot fitted right curve
plt.plot(x_data_full, gauss_f(x_data_full, *popt_left), 'r-')



# Step : Define the function for curve fitting of central curve (right side of curve which defines the vacuum edge)
x_data = x[(x>356.4) & (x<357.5) ]
y_data = y_center_neg_1V [(x>356.4) & (x<357.5)]
x2 = np.linspace(min(x_data), max(x_data),1000)

def func(x, a, b, c, d):
    return a * erf(b * (x -c)) + d


#  Fit the error function for center curve
popt_right, pcov_right = curve_fit(func, x_data, y_data, p0=(25000,0.05,356.5,59))



# Plot fitted right curve
plt.plot(x_data, func(x_data, *popt_right), 'r-', label='fit')

#midpoint of vertical line for center curve
#midpoint of vertical line for center curve
V_line = func(x_data, *popt_right)

# Calculate the center y value
center_y = (np.max(V_line) + np.min(V_line)) / 2

# Find the nearest points in y to the center_y
nearest_points = np.abs(V_line - center_y)
index_lower = np.argmin(nearest_points)
index_upper = index_lower + 1 if index_lower < len(V_line) - 1 else index_lower - 1

# Interpolate to find the x coordinate corresponding to center_y
x_lower = x_data[index_lower]
x_upper = x_data[index_upper]
y_lower = V_line[index_lower]
y_upper = V_line[index_upper]

slope = (y_upper - y_lower) / (x_upper - x_lower)
center_x = x_lower + (center_y - y_lower) / slope



#major and minor ticks on x axis
plt.minorticks_on()
plt.tick_params(axis='x', which='major',length=10)
plt.tick_params(axis='x', which='minor',length=5)

plt.xlim(352, 365)


plt.legend()
plt.show()

# =============================================================================
# #ploting 2V biasing voltage
# =============================================================================

y_center_2V = file_2V['Mean_center']*np.pi*15*15*20  # the time represents the mean count of central part
y_full_2V = file_2V['Mean_full']*750*750  # the time represents the mean count of central part


plt.figure()  #  fifth plot
plt.xlabel('Time of flight (ns)')
plt.ylabel('Counts')
plt.title('+2 V')
plt.plot(x, y_full_2V, 'o', label='full')

plt.plot(x, y_center_2V, 'o',color='gray', label= ' center')


# define the function to fit the curve for full curve (left side of the curve which denotes the fermi edge)
x_data_full = x[(x>358) & (x<362.3) ]

y_data_full = y_full_2V[(x>358) & (x<362.3)]

def gauss_f(x,A, mu, sig):
    return A*np.exp(-(x-mu)**2/sig**2)


#  Fit the error function for full integrated curve
popt_left, pcov_left = curve_fit(gauss_f, x_data_full, y_data_full, p0=(140000,362,3))


# Plot fitted right curve
plt.plot(x_data_full, gauss_f(x_data_full, *popt_left), 'r-')



# Step : Define the function for curve fitting of central curve (right side of curve which defines the vacuum edge)
x_data = x[(x>363.1) & (x<364.2) ]
y_data = y_center_2V [(x>363.1) & (x<364.2)]


def func(x, a, b, c, d): 
    return a * erf(b * (x -c)) + d #erf is the error function


#  Fit the error function for center curve
popt_right, pcov_right = curve_fit(func, x_data, y_data, p0=(40000,0.05,363,50))



# Plot fitted right curve
plt.plot(x_data, func(x_data, *popt_right), 'r-', label='fit')

#midpoint of vertical line for center curve
#midpoint of vertical line for center curve
V_line = func(x_data, *popt_right)

# Calculate the center y value
center_y = (np.max(V_line) + np.min(V_line)) / 2

# Find the nearest points in y to the center_y
nearest_points = np.abs(V_line - center_y)
index_lower = np.argmin(nearest_points)
index_upper = index_lower + 1 if index_lower < len(V_line) - 1 else index_lower - 1

# Interpolate to find the x coordinate corresponding to center_y
x_lower = x_data[index_lower]
x_upper = x_data[index_upper]
y_lower = V_line[index_lower]
y_upper = V_line[index_upper]

slope = (y_upper - y_lower) / (x_upper - x_lower)
center_x = x_lower + (center_y - y_lower) / slope



#major and minor ticks on x axis
plt.minorticks_on()
plt.tick_params(axis='x', which='major',length=10)
plt.tick_params(axis='x', which='minor',length=5)

plt.xlim(352, 365)

plt.legend()
plt.show()

# =============================================================================
# ploting -2V biasing voltage
# =============================================================================

y_center_neg_2V = file_neg_2V['Mean_center']*np.pi*15*15*20  # the time represents the mean count of central part
y_full_neg_2V = file_neg_2V['Mean_full']*750*750  # the time represents the mean count of central part




plt.figure()  #  sixth plot
plt.xlabel('Time of flight (ns)')
plt.ylabel('Counts')
plt.title('-2 V')
plt.plot(x, y_full_neg_2V, 'o', label='full')

plt.plot(x, y_center_neg_2V, 'o',color='gray', label= ' center')

# define the function to fit the curve for full curve (left side of the curve which denotes the fermi edge)
x_data_full = x[(x>352) & (x<352.6) ]

y_data_full = y_full_neg_2V[(x>352) & (x<352.6)]

def linear_f(x,M, C): #linear function which onset gives the x-cordinate
    return M*x +C


#  Fit the linear function for full integrated curve
popt_left, pcov_left = curve_fit(linear_f, x_data_full, y_data_full, p0=(1,352))


# Plot fitted left graph
plt.plot(x_data_full, linear_f(x_data_full, *popt_left), 'r-')



# Step : Define the function for curve fitting of right graph (right side of curve which defines the vacuum edge)
x_data = x[(x>354.3) & (x<355.4) ]
y_data = y_center_neg_2V [(x>354.3) & (x<355.4)]


def func(x, a, b, c, d): #func denotes the error function
    return a * erf(b * (x -c)) + d #erf is the error function



popt_right, pcov_right = curve_fit(func, x_data, y_data, p0=(100000,0.05,354,50))




# Plot fitted right curve
plt.plot(x_data, func(x_data, *popt_right), 'r-', label='fit')

#midpoint of vertical line for center curve
V_line = func(x_data, *popt_right)

# Calculate the center y value
center_y = (np.max(V_line) + np.min(V_line)) / 2

# Find the nearest points in y to the center_y
nearest_points = np.abs(V_line - center_y)
index_lower = np.argmin(nearest_points)
index_upper = index_lower + 1 if index_lower < len(V_line) - 1 else index_lower - 1


# Interpolate to find the x coordinate corresponding to center_y
x_lower = x_data[index_lower]
x_upper = x_data[index_upper]
y_lower = V_line[index_lower]
y_upper = V_line[index_upper]

slope = (y_upper - y_lower) / (x_upper - x_lower)
center_x = x_lower + (center_y - y_lower) / slope



#major and minor ticks on x axis
plt.minorticks_on()
plt.tick_params(axis='x', which='major',length=10)
plt.tick_params(axis='x', which='minor',length=5)



plt.legend()
plt.show()


# Plotting
