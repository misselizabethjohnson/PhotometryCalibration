#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 17:01:13 2019

@author: eljohnson
"""

import numpy as np
import matplotlib.pylab as plt
from matplotlib import gridspec
import pandas as pd
import math

dansFile = pd.read_csv('scolnic15_supercal_transformations.txt',delim_whitespace=True)
SDSSFile = pd.read_csv('Stripe82CrossMatch.csv')
DESFile = pd.read_csv('PS1CrossMatchedToDESSecondRegion.csv')


class Calibration: 
    
    def __init__(self):
        
        self.danFile = dansFile
        self.SDSSFile = SDSSFile
        self.DESFile = DESFile
        
        ## you can add more
        
    def TransformationPlot(self,PS1,otherSystem,xAxisColorOne,xAxisColorTwo,yAxisColor):
        
        if PS1 == 'PS1':
            x1 = '%sMeanPSFMag' % xAxisColorOne
            x2 = '%sMeanPSFMag' % xAxisColorTwo
            
            x1err = '%sMeanPSFMagErr' % xAxisColorOne
            x2err = '%sMeanPSFMagErr' % xAxisColorTwo
            
            danx1 = 'PS1-%s' % xAxisColorOne
            danx2 = 'PS1-%s' % xAxisColorTwo
            
            danxvals = self.danFile[danx1]-self.danFile[danx2]
            
            
        else:
            
            return "I'm confused, don't you want to compare to PS1?"
        
        if otherSystem == 'SDSS':
            
            y1 = 'psfMag_%s' % yAxisColor
            y1err = 'psfMagErr_%s' % yAxisColor
            
            danyvals = self.danFile['PS1-%s' % yAxisColor] - self.danFile['SDSS-%s' % yAxisColor]
            
            compareFile = self.SDSSFile
            
            ## curated cut
            #compareFile = compareFile.loc[compareFile['psfMag_r']-compareFile['psfMag_i'] > 0.2]
            
            
        if otherSystem == 'DES':
            
            if yAxisColor == 'g':
                y1 = 'G_PSF'
                y1err = 'G_PSFERR'
            if yAxisColor == 'r':
                y1 = 'R_PSF'
                y1err = 'R_PSFERR'
            if yAxisColor == 'i':
                y1 = 'I_PSF'
                y1err = 'I_PSFERR'
            if yAxisColor == 'z':
                y1 = 'Z_PSF'
                y1err = 'Z_PSFERR'
                
            danyvals = self.danFile['PS1-%s' % yAxisColor] - self.danFile['DES-%s' % yAxisColor]
                
            compareFile = self.DESFile
                
        ## you can add other systems with their naming conventions here
        ## compareFile = ...
        
        if yAxisColor == 'g':
            pltcolor = 'g'
        if yAxisColor == 'r':
            pltcolor = 'r'
        if yAxisColor == 'i':
            pltcolor = 'b'
        if yAxisColor == 'z':
            pltcolor = 'purple'
        
        ## now create the plot pieces: fitting function 
        ## put in a condition for error to be less than the mean error 
        
        xplotvals = []
        yplotvals = []
        
        for i,row in compareFile.iterrows():
            if (compareFile[x1][i] != -999 and compareFile[x2][i] != -999 
                and compareFile[x1err][i] < np.mean(compareFile.loc[compareFile[x1err] != -999, [x1err]]).item()
                and compareFile[x2err][i] < np.mean(compareFile.loc[compareFile[x2err] != -999, [x2err]]).item()
                and abs(compareFile[x1][i] - compareFile[y1][i]) < 1.5 ## a rough idea of a sigma cut 
                
                
                ):
                
                xplotvals.append(compareFile[x1][i] - compareFile[x2][i])
                yplotvals.append(compareFile[x1][i] - compareFile[y1][i])
        
        x = []
        y = []
        for i in range(0,len(xplotvals)):
            if math.isnan(xplotvals[i]) or math.isnan(yplotvals[i]):
                continue
            else:
                x.append(xplotvals[i])
                y.append(yplotvals[i])        

        from scipy.optimize import curve_fit

        def f(x, A, B): # this is your 'straight line' y=f(x)
            return A*x + B

        A,B = curve_fit(f, x, y)[0] # your data x, y to fit

        a = np.linspace(-2,2,100)
        
        Adan,Bdan = curve_fit(f,danxvals,danyvals)[0]
        
        ## make the plot 
        
        fig = plt.figure(figsize=(10,7))
        gs = gridspec.GridSpec(2,1,height_ratios=[5,2],hspace=0.2)
        
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
        
        ax0.plot(danxvals,danyvals,'.', color='k', label='Scolnic 2015')
        ax0.plot(xplotvals, yplotvals, '.', alpha=0.1, color=pltcolor, label='Acquried from Databases 2019')
        ax0.plot(a,f(a,A,B),color='limegreen',label='scipy curve_fit to database data')
        ax0.plot(a,f(a,Adan,Bdan),color='m',label='scipy curve_fit to Scolnic data')
        
        ax0.set_xlabel('%s %s - %s' % (PS1, x1, x2))
        ax0.set_ylabel('%s %s - %s %s' % (PS1,x1,otherSystem,y1))
        
        ax0.legend()
        plt.show()
        
        dataTable = [['Databases 2019','y = %sx + %s' % (np.round(A,3),np.round(B,3))],
                     ['Scolnic 2015','y = %sx + %s' % (np.round(Adan,3),np.round(Bdan,3))]]
        
        
        columnsTable = ('LOCATION', 'LINEAR SCIPY FIT')
        
        ax1.axis('off')
        tableFits = ax1.table(cellText=dataTable,cellLoc='center',colLabels=columnsTable,loc='center')
        tableFits.auto_set_font_size(True)
            
    def ColorColor(self,otherSystem,xAxisColorOne,xAxisColorTwo,yAxisColorOne,yAxisColorTwo):
        
        if otherSystem == 'SDSS':
            
            x1 = 'psfMag_%s' % xAxisColorOne
            x2 = 'psfMag_%s' % xAxisColorTwo
            
            x1err = 'psfMagErr_%s' % xAxisColorOne
            x2err = 'psfMagErr_%s' % xAxisColorTwo
            
            y1 = 'psfMag_%s' % yAxisColorOne
            y2 = 'psfMag_%s' % yAxisColorTwo    
        
            y1err = 'psfMagErr_%s' % yAxisColorOne
            y2err = 'psfMagErr_%s' % yAxisColorTwo
            
            compareFile = self.SDSSFile
            
        if otherSystem == 'DES':
            
            d = {'g' : 'G', 'r' : 'R', 'i' : 'I', 'z' : 'Z'}
            
            x1 = '%s_PSF' % d['%s' % xAxisColorOne]
            x2 = '%s_PSF' % d['%s' % xAxisColorTwo]
            
            x1err = '%s_PSFERR' % d['%s' % xAxisColorOne]
            x2err = '%s_PSFERR' % d['%s' % xAxisColorTwo]
            
            y1 = '%s_PSF' % d['%s' % yAxisColorOne]
            y2 = '%s_PSF' % d['%s' % yAxisColorTwo]
            
            y1err = '%s_PSFERR' % d['%s' % yAxisColorOne]
            y2err = '%s_PSFERR' % d['%s' % yAxisColorTwo]
                
            compareFile = self.DESFile
            
        ## you can add more systems here 
        
        xplotvals = []
        yplotvals = []
        
        for i in range(0,len(compareFile)):
            if (compareFile[x1][i] != -999 and compareFile[x2][i] != -999 
                and compareFile[y1][i] != -99 and compareFile[y2][i] != -99
                and compareFile[x1err][i] < np.mean(compareFile.loc[compareFile[x1err] != -999, [x1err]]).item()
                and compareFile[x2err][i] < np.mean(compareFile.loc[compareFile[x2err] != -999, [x2err]]).item()
                #and abs(compareFile[x1][i] - compareFile[y1][i]) < 1.5 ## a rough idea of a sigma cut 
                ):
                
                xplotvals.append(compareFile[x1][i] - compareFile[x2][i])
                yplotvals.append(compareFile[y1][i] - compareFile[y2][i])
            
        ## make the color-color plot
            
        plt.plot(xplotvals,yplotvals,'.',color='k',alpha=0.1)
        
        plt.xlabel('%s %s - %s' % (otherSystem, x1, x2))
        plt.ylabel('%s %s - %s' % (otherSystem, y1, y2))
        plt.tight_layout()
        
    def HRDiagram(self,otherSystem,xAxisColorOne,xAxisColorTwo,yAxis):
        
        if otherSystem == 'SDSS':
            
            x1 = 'psfMag_%s' % xAxisColorOne
            x2 = 'psfMag_%s' % xAxisColorTwo
            
            x1err = 'psfMagErr_%s' % xAxisColorOne
            x2err = 'psfMagErr_%s' % xAxisColorTwo
            
            y1 = 'psfMag_%s' % yAxis
        
            y1err = 'psfMagErr_%s' % yAxis
            
            compareFile = self.SDSSFile
        
        if otherSystem == 'DES':
            
            d = {'g' : 'G', 'r' : 'R', 'i' : 'I', 'z' : 'Z'}
            
            x1 = '%s_PSF' % d['%s' % xAxisColorOne]
            x2 = '%s_PSF' % d['%s' % xAxisColorTwo]
            
            x1err = '%s_PSFERR' % d['%s' % xAxisColorOne]
            x2err = '%s_PSFERR' % d['%s' % xAxisColorTwo]
            
            y1 = '%s_PSF' % d['%s' % yAxis]
            
            y1err = '%s_PSFERR' % d['%s' % yAxis]
                
            compareFile = self.DESFile

        ## you can add more systems here 
        
        xplotvals = []
        yplotvals = []
        
        for i in range(0,len(compareFile)):
            if (compareFile[x1][i] != -999 and compareFile[x2][i] != -999 
                and compareFile[y1][i] != -99 and compareFile[y2][i] != -99
                and compareFile[x1err][i] < np.mean(compareFile.loc[compareFile[x1err] != -999, [x1err]]).item()
                and compareFile[x2err][i] < np.mean(compareFile.loc[compareFile[x2err] != -999, [x2err]]).item()
                #and abs(compareFile[x1][i] - compareFile[y1][i]) < 1.5 ## a rough idea of a sigma cut 
                ):
                
                xplotvals.append(compareFile[x1][i] - compareFile[x2][i])
                yplotvals.append(compareFile[y1][i])



