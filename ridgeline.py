"""
@author: christophernovitsky
"""
import pandas as pd
import matplotlib.pyplot as plt 
import myFunctions 
import scipy.optimize
import math as m 
import numpy as np 
import statsmodels.api as sm

#%% function to find the least cos(2*theta) for the data 

def sin2fit(y, x):
    # sin least-square fit function 
    # The elements of output parameter vector, s ( b in the function ) are:

    # s(1): sine wave amplitude (in units of y)
    # s(2): phase (phase is s(2)/(2*s(3)) in units of x)
    # s(3): offset (in units of y)

    yu = max(y)
    yl = min(y)
    yr = (yu-yl)                    # Range of y
    ym = np.mean(y)
    fit = lambda b,x: b[0]*(np.sin(2*m.pi*x/180 + 2*m.pi/b[1])) + b[2]
    fcn = lambda b: sum((fit(b,x) - y)**2)
    s =  scipy.optimize.fmin(fcn, [yr,  -1,  ym]);       # Minimise Least-Squares
    xp = np.linspace(min(x),max(x),1036);
    yp = fit(s,xp);
    
    return [xp,yp,s];

#%%
## Read in UTM GPS and travel-time data for circle shots 
filename1 = 'utm.xlsx'
utm = pd.read_excel(filename1, sheetname='Sheet1')

filename2 = 'new_arrivals_20m_1_11.xlsx'
num = pd.read_excel(filename2, sheetname='Sheet1') 

## Correct starting angles to true north and finds angles from 0 to 360 for 
## all gps data 
xAxis = myFunctions.yAxis(utm)

#%% Plot 
plt.figure(1)
id1 = [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5];
id2 = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1];
fig, axs = plt.subplots(nrows=6, ncols=2)
fig.subplots_adjust(bottom=0.06, hspace=0.5)
fig.suptitle('% Seismic Anisotropy Circle Shots', fontsize='large')

ang, ansp, anspError, angError = [], [], [], [] 

for i in range(1,12):

   y, x = [], []
   x = xAxis[:,i]
   y = num.loc[0:95,i]*100

   ax = axs[id1[i-1],id2[i-1]]
   ax.set_xlim(0, 360)
   

   ax.plot(x,y,'.')
   
   output = sin2fit(y, x) # output[0] - xfit output[1] - yfit output[2] - s 
   ax.plot(output[0],output[1],'k--')
   
   index_min  = np.argmin(output[1])
   if output[0][index_min] > 180: output[0][index_min] = output[0][index_min] - 180
   ang.append(output[0][index_min])
   
   ansp.append(abs(((max(output[1])-min(output[1])))/((max(output[1])+min(output[1])/2)))*100)
   
   ## Bootstrapping to find the error in angles and ansp
   angTemp, anspTemp = [], []
   for i in range(0,40):
     sub =np.random.permutation(95)[0:72] 
     outputTemp = sin2fit(y[sub], x[sub]) 
     indexMin  = np.argmin(outputTemp[1])
     if(outputTemp[0][indexMin]>180): outputTemp[0][indexMin] = outputTemp[0][indexMin]-180
     angTemp.append(outputTemp[0][indexMin])
     anspTemp.append(abs(((max(outputTemp[1])-min(outputTemp[1])))/((max(outputTemp[1])+min(outputTemp[1])/2)))*100)
   anspError.append(np.std(anspTemp)) 
   angError.append(np.std(angTemp))
   
ax = axs[5,1]    
ax.set_visible(False)     
plt.show()
  
 #%% Calculate the error in the depth of regolith 
 
depthRegolith = [12.944,11.974,6.775,13.454,9.162,9.496,9.322, 15.981,16.712,15.273,7.937,5.357,7.484,6.932,4.446]

filename3 = 'Brady_data_10m.txt'
data = np.loadtxt(filename3)
depthStd = [np.std(np.reshape(data[:,2],(11, 72))[i]) for i in range(0,11)]

#%% remove outlier
del ansp[2]
del depthRegolith[2]
del depthStd[2]
del anspError[2]

#%% plot regoth vs. Depth 
l = 9

plt.figure(3)
plt.suptitle('% Seiemic Anisotropy vs. Depth', fontsize='large')
plt.ylabel("% Seismic Anisotropy")
plt.xlabel("Depth of the Regolith (m)")

plt.errorbar(depthRegolith[0:l], ansp[0:l], xerr= depthStd[0:l], yerr= anspError[0:l], fmt='.k')

## Ordinary Least Squares
model = sm.OLS(ansp[0:l], depthRegolith[0:l]).fit()
print(model.summary())
plt.plot(depthRegolith[0:l], model.fittedvalues, '-.k')
plt.ylim(0, 14)
plt.xlim(6, 19)
plt.show()




print(model.fittedvalues)
#
