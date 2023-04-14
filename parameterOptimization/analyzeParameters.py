# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 23:17:49 2023

@author: laniq
"""
import comparsionPlot as compPlot
import copy

results = compPlot.np.load('output/allParameter_Results.npy', allow_pickle=True)

mask1 = 'transformingTestImage/Mask_0_withPosts.png'
mask2 = 'transformingTestImage/Mask_45_withPosts.png'

contractedFrame = 41
analytical_x = compPlot.np.load('comparisionPlots' + '/analytical_dispx_%i' % contractedFrame + '_withPosts.npy')
analytical_y = compPlot.np.load('comparisionPlots' + '/analytical_dispy_%i' % contractedFrame + '_withPosts.npy')

analytical_disp = {}
analytical_disp['x'] = analytical_x[:,:]
analytical_disp['y'] = analytical_y[:,:]

analytical_disp['x'][mask2==0] = compPlot.np.nan
analytical_disp['y'][mask2==0] = compPlot.np.nan

results_rsquared = []
allRSquaredValues = []

assert len(results) != 0, "List is empty."

for x in results:
    displacementField = x[1]
    rSquaredVal = compPlot.rSquared(displacementField[:,:,0], displacementField[:,:,1], analytical_disp['x'], analytical_disp['y'], mask2, contractedFrame)
    rSquaredAvg = (rSquaredVal[0]+rSquaredVal[1])/2
    data = copy.deepcopy([x[0], rSquaredAvg])
    allRSquaredValues.append(rSquaredAvg)
    results_rsquared.append(data)

max_index = allRSquaredValues.index(max(allRSquaredValues))
bestParameters = results_rsquared[max_index][0]
bestRSquared = results_rsquared[max_index][1]

print("parameters: ", bestParameters)
print("best R Squared: ", bestRSquared)
    

    







    