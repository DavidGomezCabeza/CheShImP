
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tkinter.filedialog import askdirectory
import tkinter as tk
from ipywidgets import interact, widgets, fixed, interactive, HBox, Layout
import pickle
from sklearn.metrics import auc
import seaborn as sns
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import odeint
from scipy.optimize import differential_evolution

######################################################### Functions to analyse the CSI data for the BioMultiWell Project #########################################################


# Function to plot and visualise the CSI data and a single spectra. 
# Inputs: 
#   - ppms: Dictionary with chemical shift vectors
#   - CSIs: Dictionary with CSI arrays
#   - exp: Experiment selected for the visualisation
#   - j2: Column selection for spectra visualisation
#   - i2: Row selection for spectra visualisation
#   - scF: Scalling factor for the overall plot
#   - lw1: Linewith for the spectra plots of the CSI
#   - lw2: Linewith for the selected spectra visualisation
#   - mnn: Bottom Y limit for spectra visualisation
#   - mxx: Top Y limit for spectra visualisation
#   - ppmsel: Ppm selector

def plotCSIDat(ppms, CSIs, exp = 0, j2 = 0, i2 = 0, scF = 1, lw1 = 1, lw2 = 1, mnn = -0.1, mxx = 1.1, ppmsel = 171):
    sups = {}

    sizDt = np.shape(CSIs[str(exp-1)])
    mxdt = np.max(CSIs[str(exp-1)])
    mndt = np.min(CSIs[str(exp-1)])
    
    fig = plt.figure(figsize=(15*scF,8*scF))

    for i in range(sizDt[2]):
        for j in range(sizDt[1]):
            sups[str(i)+'-'+str(j)] = plt.subplot2grid((sizDt[2], sizDt[1]*3), (i, j), colspan=1, rowspan=1)
            
            sups[str(i)+'-'+str(j)].get_xaxis().set_visible(False)
            sups[str(i)+'-'+str(j)].get_yaxis().set_visible(False)
            sups[str(i)+'-'+str(j)].set_ylim(-0.1,1.1)
            sups[str(i)+'-'+str(j)].plot((CSIs[str(exp-1)][::,j,i] - mndt) / (mxdt - mndt),  color="#69dd52ff", linewidth = lw1)
            sups[str(i)+'-'+str(j)].spines[['right', 'top', 'left', 'bottom']].set_color('#e78829ff')
            if i == i2 and j == j2:
                sups[str(i)+'-'+str(j)].set_facecolor('#f07878ff')
            else:
                sups[str(i)+'-'+str(j)].set_facecolor('#474747ff')
            
    sups["Select"] = plt.subplot2grid((sizDt[2], sizDt[1]*3), (1, sizDt[1]+1), colspan=sizDt[1]*2, rowspan=sizDt[2]-2)
    sups["Select"].get_yaxis().set_visible(False)
    sups["Select"].spines[['right', 'top']].set_visible(False)
    sups["Select"].plot(ppms[str(exp-1)], (CSIs[str(exp-1)][::,j2,i2] - mndt) / (mxdt - mndt) ,  color="#69dd52ff", linewidth = lw2)
    # sups["Select"].plot([ppmsel,ppmsel], [-0.1,1.1], color = "Red", linewidth = 0.5)
    sups["Select"].set_xlim(sups["Select"].get_xlim()[::-1])
    sups["Select"].set_ylim(mnn,mxx)
    sups["Select"].set_ylabel("Intensity")
    sups["Select"].set_xlabel("$^{13}$C Chemical Shift (ppm)")
    sups["Select"].set_facecolor('#474747ff')
    
    plt.subplots_adjust(left=0.1, right=0.9, 
                    top=0.9, bottom=0.1, 
                    wspace=0, hspace=0)
    plt.show()
    
    

def inspecCSI(ppms, CSIs):    
    
    jsm = np.max([np.shape(CSIs[str(s)])[1] for s in range(len(CSIs))])   
    ism = np.max([np.shape(CSIs[str(s)])[2] for s in range(len(CSIs))])
        
    mm=interact(
        plotCSIDat,
        ppms = fixed(ppms),
        CSIs = fixed(CSIs),
        exp = widgets.Dropdown(options=[list(range(len(ppms)))[g]+1 for g in range(len(ppms))], value = 1, description = "Exp.:"), 
        j2 = widgets.IntSlider(min=0,max=jsm-1,step=1,value=0, description ='Column:', layout=Layout(width='1400px')), 
        i2 = widgets.IntSlider(min=0,max=ism-1,step=1,value=0, description ='Row:', layout=Layout(width='1400px')), 
        scF = widgets.FloatSlider(min=0.01,max=5,step=0.000001,value=1, readout = True, description ='Zoom', layout=Layout(width='1400px')), 
        lw1 = widgets.FloatSlider(min=0.01,max=5,step=0.000001,value=1, readout = True, description ='LW CSI', layout=Layout(width='1400px')), 
        lw2 = widgets.FloatSlider(min=0.01,max=5,step=0.000001,value=1, readout = True, description ='LW Spec.', layout=Layout(width='1400px')), 
        mnn = widgets.FloatSlider(min=-0.01,max=1.1,step=0.000001,value=-0.01, readout = True, description ='Z. Bot.', layout=Layout(width='1400px')), 
        mxx = widgets.FloatSlider(min=-0.01,max=1.1,step=0.000001,value=1.1, readout = True, description ='Z. Top.', layout=Layout(width='1400px'))
    )
    
    
    
    
##################################################################################################################################################################
    
# Funtions for the selection of integration regions of Pyruvate and Lactate
# Inputs: 
#   - ppms: Dictionary with chemical shift vectors
#   - CSIs: Dictionary with CSI arrays
#   - exp: Experiment selected for the visualisation
#   - j2: Column selection for spectra visualisation
#   - i2: Row selection for spectra visualisation
#   - scF: Scalling factor for the overall plot
#   - lw2: Linewith for the selected spectra visualisation 
#   - zmmY: Bottom and top Y limit for spectra visualisation
#   - zmmX: Min and Max ppm for the display (zoom) 
#   - intP: Vector containing the negative of the maximum and the negative of the minimum ppm for the Pyruvate integration region 
#   - intL: Vector containing the negative of the maximum and the negative of the minimum ppm for the Lactate integration region 


def SelectIntegReg(ppms, CSIs, exp = 1, j2 = 0, i2 = 0, scF = 1, lw2 = 1, zmmY = [0,0.1], zmmX = [160,200], intP=[-180, -160], intL=[-190, -181], intNM = 'IntegrPars8WCon'):

    mxdt = np.max(CSIs[str(exp-1)])
    mndt = np.min(CSIs[str(exp-1)])

    fig2, ax2 = plt.subplots(figsize=(15*scF,8*scF))
            
    ax2.spines[['right', 'top']].set_visible(False)
    ax2.plot(ppms[str(exp-1)], (CSIs[str(exp-1)][::,j2,i2] - mndt) / (mxdt - mndt) ,  color="#69dd52ff", linewidth = lw2)
    # sups["Select"].plot([ppmsel,ppmsel], [-0.1,1.1], color = "Red", linewidth = 0.5)

    ax2.set_ylim(zmmY[0],zmmY[1])
    ax2.set_xlim(zmmX[0],zmmX[1])
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylabel("Intensity")
    ax2.set_xlabel("$^{13}$C Chemical Shift (ppm)")
    ax2.set_facecolor('#474747ff')


    ax2.plot([abs(intP[1]), abs(intP[1])],[0, 1], color = '#e7be29ff', alpha = 0.4)
    ax2.plot([abs(intP[0]), abs(intP[0])],[0, 1], color = '#e7be29ff', alpha = 0.4)
    ax2.fill_between(ppms[str(exp-1)], 0, 1, where= (ppms[str(exp-1)] > abs(intP[1])) & (ppms[str(exp-1)] < abs(intP[0])),
                    facecolor='#e7be29ff', alpha=0.2)

    ax2.plot([abs(intL[1]), abs(intL[1])],[0, 1], color = '#e72984ff', alpha = 0.4)
    ax2.plot([abs(intL[0]), abs(intL[0])],[0, 1], color = '#e72984ff', alpha = 0.4)
    ax2.fill_between(ppms[str(exp-1)], 0, 1, where= (ppms[str(exp-1)] > abs(intL[1])) & (ppms[str(exp-1)] < abs(intL[0])),
                    facecolor='#e72984ff', alpha=0.2)

    if os.path.isdir('Processed') == False:
        os.mkdir('Processed')
        
    if os.path.isfile('Processed/'+intNM+'.p') == True:
        with open('Processed/'+intNM+'.p', 'rb') as fp:
            intPars = pickle.load(fp)
        intPars[str(exp-1)]= [intP, intL]
    else:
        intPars = {}
        intPars[str(exp-1)]= [intP, intL]

    with open('Processed/'+intNM+'.p', 'wb') as fp:
            pickle.dump(intPars, fp, protocol=pickle.HIGHEST_PROTOCOL)

    plt.show()
    
# SelectIntegReg(ppms, CSIs, exp = 1, j2 = 0, i2 = 0, scF = 1, lw2 = 1, zmmY = [0,0.1], zmmX = [160,200], intP=[-180, -170], intL=[-190, -181])


def dispInts(ppms, CSIs, intNM = 'IntegrPars8WCon'):
    phagan = True
    if os.path.isfile('Processed/'+intNM+'.p'):
        root = tk.Tk()
        root.attributes('-topmost',True)
        root.iconify()
        phagan = tk.messagebox.askyesno("WARNING", "You already selected the integration regions for this data. If you proceed the these values will be deleted and you will have to re-do it. Are you sure you want to proceed?")
        root.destroy()
            
    if phagan == True:
        intPars = {}
        for i in range(len(CSIs)):
            intPars[str(i)] = [[-180, -170], [-190, -181]]                
                
        jsm = np.max([np.shape(CSIs[str(s)])[1] for s in range(len(CSIs))])   
        ism = np.max([np.shape(CSIs[str(s)])[2] for s in range(len(CSIs))])

        minppm = np.min([np.min(ppms[str(s)]) for s in range(len(ppms))])
        maxppm = np.max([np.max(ppms[str(s)]) for s in range(len(ppms))])

        mm=interact(
            SelectIntegReg,
            ppms = fixed(ppms),
            CSIs = fixed(CSIs),
            intNM = fixed(intNM),
            exp = widgets.Dropdown(options=[list(range(len(ppms)))[g]+1 for g in range(len(ppms))], value = 1, description = "Exp.:"), 
            j2 = widgets.IntSlider(min=0,max=jsm-1,step=1,value=0, description ='Column:', layout=Layout(width='1400px')), 
            i2 = widgets.IntSlider(min=0,max=ism-1,step=1,value=0, description ='Row:', layout=Layout(width='1400px')), 
            scF = widgets.FloatSlider(min=0.01,max=5,step=0.000001,value=1, readout = True, description ='Zoom', layout=Layout(width='1400px')), 
            lw2 = widgets.FloatSlider(min=0.01,max=5,step=0.000001,value=1, readout = True, description ='LW Spec.', layout=Layout(width='1400px')), 
            zmmY = widgets.FloatRangeSlider(min=-0.1,max=1.1,step=0.0000001,value=[-0.1, 1.1], readout = False, layout=Layout(width='1400px'), description ='ZoomY'), 
            zmmX = widgets.FloatRangeSlider(min=minppm,max=maxppm,step=0.0000001,value=[minppm, maxppm], readout = False, layout=Layout(width='1400px'), description ='ZoomX'), 
            intP = widgets.FloatRangeSlider(min=-maxppm,max=-minppm,step=0.0000001,value=[-180, -170], readout = False, layout=Layout(width='1400px'), description ='Int1'), 
            intL = widgets.FloatRangeSlider(min=-maxppm,max=-minppm,step=0.0000001,value=[-190, -181], readout = False, layout=Layout(width='1400px'), description ='Int2')
        )
        
 
 ##################################################################################################################################################################


# Functions to extrat integrals from pyruvate and lactate and generate a colour map for each and the ratio. 

# def intMaps(ppms, CSIs_FC_NRMRW, scF2=0.5, exD = 1, cmp1 = 'YlOrBr', cmp2 = 'GnBu', cmp3 = 'YlGn', intNM = 'IntegrPars8WCon', savNM = '2_RatiosMap_8WellConditions', logg = False):

#     with open('Processed/'+intNM+'.p', 'rb') as fp:
#         intPars = pickle.load(fp)
    
#     pyrINT = {}
#     lacINT = {}

#     for ex in range(len(CSIs_FC_NRMRW)):
#         pyrIntR = np.logical_and(ppms[str(ex)] >= np.abs(intPars[str(ex)][0])[1],  ppms[str(ex)] <= np.abs(intPars[str(ex)][0])[0])
#         lacIntR = np.logical_and(ppms[str(ex)] >= np.abs(intPars[str(ex)][1])[1],  ppms[str(ex)] <= np.abs(intPars[str(ex)][1])[0])
        
#         sizDt = np.shape(CSIs_FC_NRMRW[str(ex)])
#         pyrMat = np.zeros([sizDt[1], sizDt[2]])
#         lacMat = np.zeros([sizDt[1], sizDt[2]])

#         for i in range(sizDt[2]):
#             for j in range(sizDt[1]):
#                 pyrMat[j,i] = auc(ppms[str(ex)][pyrIntR], CSIs_FC_NRMRW[str(ex)][pyrIntR, j,i])
#                 lacMat[j,i] = auc(ppms[str(ex)][lacIntR], CSIs_FC_NRMRW[str(ex)][lacIntR, j,i])
                
#         pyrINT[str(ex)] = pyrMat
#         lacINT[str(ex)] = lacMat
        
        
#         plt.ioff()
        
#         scF = 0.5
#         fig, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(sizDt[1]*3*scF, sizDt[2]*scF)) 
#         g1 = sns.heatmap(np.transpose(pyrINT[str(ex)]),ax=ax1, cmap=cmp1)
#         g1.get_xaxis().set_visible(False)
#         g1.get_yaxis().set_visible(False)
#         g1.set_title('Int1')

#         g2 = sns.heatmap(np.transpose(lacINT[str(ex)]),ax=ax2, cmap=cmp2)
#         g2.get_xaxis().set_visible(False)
#         g2.get_yaxis().set_visible(False)
#         g2.set_title('Int2')

#         if logg == False:
#             g3 = sns.heatmap(np.transpose(lacINT[str(ex)]/pyrINT[str(ex)]),ax=ax3, cmap=cmp3)
#             g3.set_title('Int2/Int1')
#         elif logg == True:
#             g3 = sns.heatmap(np.log(np.transpose(lacINT[str(ex)]/pyrINT[str(ex)])),ax=ax3, cmap=cmp3)
#             g3.set_title('log(Int2/Int1)')
#         g3.get_xaxis().set_visible(False)
#         g3.get_yaxis().set_visible(False)

#         plt.savefig('Processed/'+savNM+'_Exp'+str(ex)+'.svg')
#         plt.savefig('Processed/'+savNM+'_Exp'+str(ex)+'.png')
#         plt.close(fig)
        
        
#     # plt.ioff()
#     sizDt2 = np.shape(CSIs_FC_NRMRW[str(exD-1)])
#     fig2, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(sizDt2[1]*3*scF2, sizDt2[2]*scF2)) 
#     g1 = sns.heatmap(np.transpose(pyrINT[str(exD-1)]),ax=ax1, cmap=cmp1)
#     g1.get_xaxis().set_visible(False)
#     g1.get_yaxis().set_visible(False)
#     g1.set_title('Int1')

#     g2 = sns.heatmap(np.transpose(lacINT[str(exD-1)]),ax=ax2, cmap=cmp2)
#     g2.get_xaxis().set_visible(False)
#     g2.get_yaxis().set_visible(False)
#     g2.set_title('Int2')

#     if logg == False:
#         g3 = sns.heatmap(np.transpose(lacINT[str(ex)]/pyrINT[str(ex)]),ax=ax3, cmap=cmp3)
#         g3.set_title('Int2/Int1')
#     elif logg == True:
#         g3 = sns.heatmap(np.log(np.transpose(lacINT[str(ex)]/pyrINT[str(ex)])),ax=ax3, cmap=cmp3)
#         g3.set_title('log(Int2/Int1)')
#     g3.get_xaxis().set_visible(False)
#     g3.get_yaxis().set_visible(False)
    
#     plt.show()
    
def intMaps(ppms, CSIs_FC_NRMRW, scF2=0.5, exD = 1, cmp1 = 'YlOrBr', cmp2 = 'GnBu', cmp3 = 'YlGn', intNM = 'IntegrPars8WCon', savNM = '2_RatiosMap_8WellConditions', logg = False):

    with open('Processed/'+intNM+'.p', 'rb') as fp:
        intPars = pickle.load(fp)
    
    pyrINT = {}
    lacINT = {}
    
    

    # for ex in range(len(CSIs_FC_NRMRW)):
    ex = exD-1
    pyrIntR = np.logical_and(ppms[str(ex)] >= np.abs(intPars[str(ex)][0])[1],  ppms[str(ex)] <= np.abs(intPars[str(ex)][0])[0])
    lacIntR = np.logical_and(ppms[str(ex)] >= np.abs(intPars[str(ex)][1])[1],  ppms[str(ex)] <= np.abs(intPars[str(ex)][1])[0])
    
    sizDt = np.shape(CSIs_FC_NRMRW[str(ex)])
    pyrMat = np.zeros([sizDt[1], sizDt[2]])
    lacMat = np.zeros([sizDt[1], sizDt[2]])
    
    

    for i in range(sizDt[2]):
        for j in range(sizDt[1]):
            pyrMat[j,i] = auc(ppms[str(ex)][pyrIntR], CSIs_FC_NRMRW[str(ex)][pyrIntR, j,i])
            lacMat[j,i] = auc(ppms[str(ex)][lacIntR], CSIs_FC_NRMRW[str(ex)][lacIntR, j,i])
            
    if np.min(pyrMat) <=0:
        pyrMat = pyrMat+ (abs(np.min(pyrMat))+(abs(np.min(pyrMat))*0.1))
        lacMat = lacMat+ (abs(np.min(pyrMat))+(abs(np.min(pyrMat))*0.1))
    if np.min(lacMat) <=0:
        pyrMat = pyrMat+ (abs(np.min(lacMat))+(abs(np.min(lacMat))*0.1))
        lacMat = lacMat+ (abs(np.min(lacMat))+(abs(np.min(lacMat))*0.1))
    
    
    pyrINT[str(ex)] = pyrMat
    lacINT[str(ex)] = lacMat
    
    
    plt.ioff()
    
    scF = 0.5
    fig, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(sizDt[1]*3*scF, sizDt[2]*scF)) 
    g1 = sns.heatmap(np.transpose(pyrINT[str(ex)]),ax=ax1, cmap=cmp1)
    g1.get_xaxis().set_visible(False)
    g1.get_yaxis().set_visible(False)
    g1.set_title('Int1')

    g2 = sns.heatmap(np.transpose(lacINT[str(ex)]),ax=ax2, cmap=cmp2)
    g2.get_xaxis().set_visible(False)
    g2.get_yaxis().set_visible(False)
    g2.set_title('Int2')

    
    if logg == False:
        g3 = sns.heatmap(np.transpose(lacINT[str(ex)]/pyrINT[str(ex)]),ax=ax3, cmap=cmp3)
        g3.set_title('Int2/Int1')
    elif logg == True:
        g3 = sns.heatmap(np.log(np.transpose(lacINT[str(ex)]/pyrINT[str(ex)])),ax=ax3, cmap=cmp3)
        g3.set_title('log(Int2/Int1)')
    g3.get_xaxis().set_visible(False)
    g3.get_yaxis().set_visible(False)

    plt.savefig('Processed/'+savNM+'_Exp'+str(ex)+'.svg')
    plt.savefig('Processed/'+savNM+'_Exp'+str(ex)+'.png')
    plt.close(fig)
        
        
    # plt.ioff()
    sizDt2 = np.shape(CSIs_FC_NRMRW[str(exD-1)])
    fig2, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(sizDt2[1]*3*scF2, sizDt2[2]*scF2)) 
    g1 = sns.heatmap(np.transpose(pyrINT[str(exD-1)]),ax=ax1, cmap=cmp1)
    g1.get_xaxis().set_visible(False)
    g1.get_yaxis().set_visible(False)
    g1.set_title('Int1')

    g2 = sns.heatmap(np.transpose(lacINT[str(exD-1)]),ax=ax2, cmap=cmp2)
    g2.get_xaxis().set_visible(False)
    g2.get_yaxis().set_visible(False)
    g2.set_title('Int2')

    if logg == False:
        g3 = sns.heatmap(np.transpose(lacINT[str(ex)]/pyrINT[str(ex)]),ax=ax3, cmap=cmp3)
        g3.set_title('Int2/Int1')
    elif logg == True:
        g3 = sns.heatmap(np.log(np.transpose(lacINT[str(ex)]/pyrINT[str(ex)])),ax=ax3, cmap=cmp3)
        g3.set_title('log(Int2/Int1)')
    g3.get_xaxis().set_visible(False)
    g3.get_yaxis().set_visible(False)
    
    plt.show()
    
    
    # return(pyrINT, lacINT)


# intMaps(ppms, CSIs_FC_NRMRW, 0.5, 1, cmp1 = 'YlOrBr', cmp2 = 'GnBu', cmp3 = 'YlGn')

def pltRatsMaps(ppms, CSIs_FC_NRMRW, intNM = 'IntegrPars8WCon', savNM = '2_RatiosMap_8WellConditions', logg=False):
    lstCMPs = ['viridis', 'plasma', 'inferno', 'magma', 'cividis',
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
                'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
                'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper',
                'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
                'berlin', 'managua', 'vanimo',
                'twilight', 'twilight_shifted', 'hsv',
                'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
                'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b',
                'tab20c','flag', 'prism', 'ocean', 'gist_earth', 'terrain',
                'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap',
                'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet',
                'turbo', 'nipy_spectral', 'gist_ncar']

    mm=interact(
        intMaps,
        ppms = fixed(ppms),
        CSIs_FC_NRMRW = fixed(CSIs_FC_NRMRW),
        intNM = fixed(intNM),
        savNM = fixed(savNM),
        logg = fixed(logg),
        exD = widgets.Dropdown(options=[list(range(len(ppms)))[g]+1 for g in range(len(ppms))], value = 1, description = "Exp.:"), 
        scF2=widgets.FloatSlider(min=0.001,max=2,step=0.000001,value=0.5, readout = True, description ='Zoom', layout=Layout(width='1400px')), 
        cmp1 = widgets.Dropdown(options=lstCMPs, value = 'YlOrBr', description = "Col. Pyr"), 
        cmp2 = widgets.Dropdown(options=lstCMPs, value = 'GnBu', description = "Col. Lac"), 
        cmp3 = widgets.Dropdown(options=lstCMPs, value = 'YlGn', description = "Col. Rat")
    )
    
    
    
##################################################################################################################################################################

# Functions to plot (and save if desired) the voxel sum up final matrix 
# Inputs: 
#   - ppms: Dictionary with chemical shift vectors
#   - CSI_Sum_Norm: Dictionary with CSI arrays
#   - exp: Experiment selected for the visualisation
#   - scF: Scalling factor for the overall plot
#   - lw1: Linewith for the spectra visualisation 
#   - zmX: Min and Max ppm for the display (zoom) 
#   - sav: to save the plot or not in png and svg

def desigSumUpM(ppms, CSI_Sum_Norm, exp = 1, scF = 1, lw1 = 1, zmX = [-200, -150], zmY = [-0.1, 1.1], sav = False, pltNM = '2_SumUpMatrix_8WellConditions'):

    sizDt = np.shape(CSI_Sum_Norm[str(exp-1)])
    indsD = np.logical_and(ppms[str(exp-1)] >= -zmX[1], ppms[str(exp-1)] <= -zmX[0])
    
    print(sizDt)
    fig, ax = plt.subplots(sizDt[2],sizDt[1], figsize=(sizDt[1]*scF,sizDt[2]*scF))
    for i in range(sizDt[2]):
        for j in range(sizDt[1]):
            if sizDt[2] == 1:
                ax[j].get_xaxis().set_visible(False)
                ax[j].get_yaxis().set_visible(False)
                ax[j].set_ylim(zmY[0],zmY[1])
                
                ax[j].plot(ppms[str(exp-1)][indsD], CSI_Sum_Norm[str(exp-1)][::,j,i][indsD],  color="#69dd52ff", linewidth = lw1)
                ax[j].set_facecolor('#474747ff')
                ax[j].spines[['right', 'top', 'left', 'bottom']].set_color('#e78829ff')
                # ax[i,j].set_xlim(-zmX[1], -zmX[0])
                ax[j].set_xlim(ax[j].get_xlim()[::-1])
            elif sizDt[1] == 1:
                ax[i].get_xaxis().set_visible(False)
                ax[i].get_yaxis().set_visible(False)
                ax[i].set_ylim(zmY[0],zmY[1])
                
                ax[i].plot(ppms[str(exp-1)][indsD], CSI_Sum_Norm[str(exp-1)][::,j,i][indsD],  color="#69dd52ff", linewidth = lw1)
                ax[i].set_facecolor('#474747ff')
                ax[i].spines[['right', 'top', 'left', 'bottom']].set_color('#e78829ff')
                # ax[i,j].set_xlim(-zmX[1], -zmX[0])
                ax[i].set_xlim(ax[i].get_xlim()[::-1])
            else:
                ax[i,j].get_xaxis().set_visible(False)
                ax[i,j].get_yaxis().set_visible(False)
                ax[i,j].set_ylim(zmY[0],zmY[1])
                
                ax[i,j].plot(ppms[str(exp-1)][indsD], CSI_Sum_Norm[str(exp-1)][::,j,i][indsD],  color="#69dd52ff", linewidth = lw1)
                ax[i,j].set_facecolor('#474747ff')
                ax[i,j].spines[['right', 'top', 'left', 'bottom']].set_color('#e78829ff')
                # ax[i,j].set_xlim(-zmX[1], -zmX[0])
                ax[i,j].set_xlim(ax[i,j].get_xlim()[::-1])
            
    plt.subplots_adjust(left=0.1, right=0.9, 
                        top=0.9, bottom=0.1, 
                        wspace=0, hspace=0)

    if sav == True:
        plt.savefig('Processed/'+pltNM+'_Exp'+str(exp-1)+'.png')
        plt.savefig('Processed/'+pltNM+'_Exp'+str(exp-1)+'.svg')
        
    plt.show()
    
    
    
    
def inspecCSISumUp(ppms, CSI_Sum_Norm, pltNM = '2_SumUpMatrix_8WellConditions'):    
    
    minppm = np.min([np.min(ppms[str(s)]) for s in range(len(ppms))])
    maxppm = np.max([np.max(ppms[str(s)]) for s in range(len(ppms))])
        
    mm=interact(
        desigSumUpM,
        ppms = fixed(ppms),
        CSI_Sum_Norm = fixed(CSI_Sum_Norm),
        pltNM = fixed(pltNM),
        exp = widgets.Dropdown(options=[list(range(len(ppms)))[g]+1 for g in range(len(ppms))], value = 1, description = "Exp.:"), 
        scF = widgets.FloatSlider(min=0.01,max=5,step=0.000001,value=2, readout = True, description ='Zoom', layout=Layout(width='1400px')),
        lw1 = widgets.FloatSlider(min=0.01,max=5,step=0.000001,value=1, readout = True, description ='LW', layout=Layout(width='1400px')), 
        zmX = widgets.FloatRangeSlider(min=-maxppm ,max=-minppm,step=0.0000001,value=[-maxppm, -minppm], readout = False, layout=Layout(width='1400px'), description ='ZoomX'),
        zmY = widgets.FloatRangeSlider(min=-0.1 ,max=1.1,step=0.0000001,value=[-0.1, 1.1], readout = False, layout=Layout(width='1400px'), description ='ZoomY'),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON")
    )
    
    
    
##################################################################################################################################################################

# Functions to extract the integral values from previously selected regions, then do the ratio and display it

def getFinInts(ppms, CSI_Sum_Norm, intNM = "IntegrPars8WCon", datNM = '2_ProcDat8WCon', csvNM = '8WellConditions'):
    
    with open('Processed/'+intNM+'.p', 'rb') as fp:
        intPars = pickle.load(fp)
        
    pyrINT = {}
    lacINT = {}
    ratINT = {}

    for ex in range(len(CSI_Sum_Norm)):
        pyrIntR = np.logical_and(ppms[str(ex)] >= np.abs(intPars[str(ex)][0])[1],  ppms[str(ex)] <= np.abs(intPars[str(ex)][0])[0])
        lacIntR = np.logical_and(ppms[str(ex)] >= np.abs(intPars[str(ex)][1])[1],  ppms[str(ex)] <= np.abs(intPars[str(ex)][1])[0])
        
        sizDt = np.shape(CSI_Sum_Norm[str(ex)])
        pyrMat = np.zeros([sizDt[1], sizDt[2]])
        lacMat = np.zeros([sizDt[1], sizDt[2]])

        for i in range(sizDt[2]):
            for j in range(sizDt[1]):
                pyrMat[j,i] = auc(ppms[str(ex)][pyrIntR], CSI_Sum_Norm[str(ex)][pyrIntR, j,i])
                lacMat[j,i] = auc(ppms[str(ex)][lacIntR], CSI_Sum_Norm[str(ex)][lacIntR, j,i])
                
        pyrINT[str(ex)] = pyrMat
        lacINT[str(ex)] = lacMat
        ratINT[str(ex)] = lacMat/pyrMat    
                
    ratsd = [pyrINT, lacINT, ratINT]

    with open('Processed/'+datNM+'.p', 'wb') as fp:
        pickle.dump(ratsd, fp, protocol=pickle.HIGHEST_PROTOCOL)

    pdF = pd.DataFrame()
    pdFP = pd.DataFrame()
    pdFL = pd.DataFrame()
    
    
    for ex in range(len(ratINT)):
        if csvNM == '8WellConditions':
            pd1 = pd.DataFrame(np.transpose(ratINT[str(ex)]), columns=["exp "+str(ex+1), "exp "+str(ex+1)], index = ['Cell Glu+Gln / Glu+Gln Lys.','NADH Lys. / Cell','Cell NADH / Lys.','- C / + C'])
            pdF = pd.concat([pdF, pd1], axis = 1)
            
            pd1P = pd.DataFrame(np.transpose(pyrINT[str(ex)]), columns=["exp "+str(ex+1), "exp "+str(ex+1)], index = ['Cell Glu+Gln / Glu+Gln Lys.','NADH Lys. / Cell','Cell NADH / Lys.','- C / + C'])
            pdFP = pd.concat([pdFP, pd1P], axis = 1)
            
            pd1L = pd.DataFrame(np.transpose(lacINT[str(ex)]), columns=["exp "+str(ex+1), "exp "+str(ex+1)], index = ['Cell Glu+Gln / Glu+Gln Lys.','NADH Lys. / Cell','Cell NADH / Lys.','- C / + C'])
            pdFL = pd.concat([pdFL, pd1L], axis = 1)
            
        elif csvNM == '4WHepVsHeLa':
            pd1 = pd.DataFrame(np.transpose(ratINT[str(ex)]), columns=["exp "+str(ex+1)], index = ['HepG2','HeLa','HepG2','HeLa'])
            pdF = pd.concat([pdF, pd1], axis = 1)
            
            pd1P = pd.DataFrame(np.transpose(pyrINT[str(ex)]), columns=["exp "+str(ex+1)], index = ['HepG2','HeLa','HepG2','HeLa'])
            pdFP = pd.concat([pdFP, pd1P], axis = 1)
            
            pd1L = pd.DataFrame(np.transpose(lacINT[str(ex)]), columns=["exp "+str(ex+1)], index = ['HepG2','HeLa','HepG2','HeLa'])
            pdFL = pd.concat([pdFL, pd1L], axis = 1)
            
        elif csvNM == '4WDrugAZD':
            pd1 = pd.DataFrame(np.transpose(ratINT[str(ex)]), columns=["exp "+str(ex+1)], index = ['AZD 500 nM','AZD 5 nM','AZD 0.5 nM','AZD 0 nM'])
            pdF = pd.concat([pdF, pd1], axis = 1)
            
            pd1P = pd.DataFrame(np.transpose(pyrINT[str(ex)]), columns=["exp "+str(ex+1)], index = ['AZD 500 nM','AZD 5 nM','AZD 0.5 nM','AZD 0 nM'])
            pdFP = pd.concat([pdFP, pd1P], axis = 1)
            
            pd1L = pd.DataFrame(np.transpose(lacINT[str(ex)]), columns=["exp "+str(ex+1)], index = ['AZD 500 nM','AZD 5 nM','AZD 0.5 nM','AZD 0 nM'])
            pdFL = pd.concat([pdFL, pd1L], axis = 1)     
        elif csvNM == '4WDrugFX11':
            pd1 = pd.DataFrame(np.transpose(ratINT[str(ex)]), columns=["exp "+str(ex+1)], index = ['FX11 100 uM','FX11 40 uM','FX11 10 uM','FX11 0 uM'])
            pdF = pd.concat([pdF, pd1], axis = 1)
            
            pd1P = pd.DataFrame(np.transpose(pyrINT[str(ex)]), columns=["exp "+str(ex+1)], index = ['FX11 100 uM','FX11 40 uM','FX11 10 uM','FX11 0 uM'])
            pdFP = pd.concat([pdFP, pd1P], axis = 1)
            
            pd1L = pd.DataFrame(np.transpose(lacINT[str(ex)]), columns=["exp "+str(ex+1)], index = ['FX11 100 uM','FX11 40 uM','FX11 10 uM','FX11 0 uM'])
            pdFL = pd.concat([pdFL, pd1L], axis = 1)     
        
    pdF.to_csv(r'Processed\RatiosDataAll_'+csvNM+'.csv')
    pdFP.to_csv(r'Processed\PyrDataAll_'+csvNM+'.csv')
    pdFL.to_csv(r'Processed\LacDataAll_'+csvNM+'.csv')

    return(pyrINT, lacINT, ratINT)




def plotFinMatsRats(CSI_Sum_Norm, pyrINT, lacINT, ratINT, exD = 1, scF2 = 0.5, cmp1 = 'YlOrBr' , cmp2 = 'GnBu' , cmp3 = 'YlGn', figNM = '2_RatiosMapSumUp_8WellConditions'):


    for ex in range(len(CSI_Sum_Norm)):            
        sizDt = np.shape(CSI_Sum_Norm[str(ex)])
            
        plt.ioff()
            
        scF = 0.5
        fig, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(sizDt[1]*3*scF, sizDt[2]*scF)) 
        g1 = sns.heatmap(np.transpose(pyrINT[str(ex)]),ax=ax1, cmap=cmp1)
        g1.get_xaxis().set_visible(False)
        g1.get_yaxis().set_visible(False)
        g1.set_title('Pyr')

        g2 = sns.heatmap(np.transpose(lacINT[str(ex)]),ax=ax2, cmap=cmp2)
        g2.get_xaxis().set_visible(False)
        g2.get_yaxis().set_visible(False)
        g2.set_title('Lac')

        g3 = sns.heatmap(np.transpose(lacINT[str(ex)]/pyrINT[str(ex)]),ax=ax3, cmap=cmp3)
        g3.get_xaxis().set_visible(False)
        g3.get_yaxis().set_visible(False)
        g3.set_title('Lac/Pyr')

        plt.savefig('Processed/'+figNM+'_Exp'+str(ex)+'.svg')
        plt.savefig('Processed/'+figNM+'_Exp'+str(ex)+'.png')
        plt.close(fig)
        
        
    sizDt2 = np.shape(CSI_Sum_Norm[str(exD-1)])
    fig2, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(sizDt2[1]*3*scF2, sizDt2[2]*scF2)) 
    g1 = sns.heatmap(np.transpose(pyrINT[str(exD-1)]),ax=ax1, cmap=cmp1)
    g1.get_xaxis().set_visible(False)
    g1.get_yaxis().set_visible(False)
    g1.set_title('Pyr')

    g2 = sns.heatmap(np.transpose(lacINT[str(exD-1)]),ax=ax2, cmap=cmp2)
    g2.get_xaxis().set_visible(False)
    g2.get_yaxis().set_visible(False)
    g2.set_title('Lac')

    g3 = sns.heatmap(np.transpose(ratINT[str(exD-1)]),ax=ax3, cmap=cmp3)
    g3.get_xaxis().set_visible(False)
    g3.get_yaxis().set_visible(False)
    g3.set_title('Lac/Pyr')
    plt.show()


def pltRatsMapsSumUp(CSI_Sum_Norm, pyrINT, lacINT, ratINT, figNM):
    lstCMPs = ['viridis', 'plasma', 'inferno', 'magma', 'cividis',
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
                'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
                'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper',
                'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
                'berlin', 'managua', 'vanimo',
                'twilight', 'twilight_shifted', 'hsv',
                'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
                'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b',
                'tab20c','flag', 'prism', 'ocean', 'gist_earth', 'terrain',
                'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap',
                'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet',
                'turbo', 'nipy_spectral', 'gist_ncar']

    mm=interact(
        plotFinMatsRats,
        CSI_Sum_Norm = fixed(CSI_Sum_Norm),
        pyrINT = fixed(pyrINT), 
        lacINT = fixed(lacINT), 
        ratINT = fixed(ratINT),
        figNM = fixed(figNM),
        exD = widgets.Dropdown(options=[list(range(len(CSI_Sum_Norm)))[g]+1 for g in range(len(CSI_Sum_Norm))], value = 1, description = "Exp.:"), 
        scF2=widgets.FloatSlider(min=0.001,max=2,step=0.000001,value=0.5, readout = True, description ='Zoom', layout=Layout(width='1400px')), 
        cmp1 = widgets.Dropdown(options=lstCMPs, value = 'YlOrBr', description = "Col. Pyr"), 
        cmp2 = widgets.Dropdown(options=lstCMPs, value = 'GnBu', description = "Col. Lac"), 
        cmp3 = widgets.Dropdown(options=lstCMPs, value = 'YlGn', description = "Col. Rat")
    )
    
##################################################################################################################################################################













































##################################################################################################################################################################
# Functions to simulate a simple model for Pyr->Lac and estimate Kpl. 
def LDH(Y, t, p):
    
    P0 = Y[0]
    L0 = Y[1]
    
    k_f1 = p[0] 
    k_r1 = p[1]
    T1_P = p[2]
    T1_X = p[3] 
    
    
    d1 = k_r1 * L0 - ((1/T1_P) + k_f1) * P0
    d2 = k_f1 * P0 - ((1/T1_X) + k_r1) * L0
    
    # d1 = p[1] * Y[1] - ((1 / p[2]) + p[0]) * Y[0]
    # d2 = p[0] * Y[0] - ((1 / p[3]) + p[1]) * Y[1]
    
    dYdt = [d1, d2]
    
    return(dYdt)


def calc_method(theta, datsSin, timsSin):
    mxP = np.max(datsSin['Pyr'])

    Y0 = [datsSin['Pyr'][0], datsSin['Lac'][0]]

    result = odeint(LDH, Y0, timsSin, args=(theta,))

    pyr = [result[i][0] for i in range(len(result))]
    lac = [result[i][1] for i in range(len(result))]
    
    obP = np.sum( np.sqrt(((pyr/mxP) - (datsSin['Pyr']/mxP))**2) ) / len(pyr)
    obL = np.sum( np.sqrt(((lac/mxP) - (datsSin['Lac']/mxP))**2) ) / len(pyr)

    Obj = obP + obL
    
    return(Obj)


def estimParsModel(intNM2, datsAll, timsAll, mits = 100, bns = [(0,1), (0, 1e-5), (20, 80), (10, 70)], sts = 'randtobest1bin'):

    redo = True

    if os.path.isfile('Processed/'+intNM2+'.p'):
        root = tk.Tk()
        root.attributes('-topmost',True)
        root.iconify()
        redo = tk.messagebox.askyesno("WARNING", "You already estimated the model parameters for this set of experiments. If you proceed the these values will be deleted and you will be re-estimated. Are you sure you want to proceed?")
        root.destroy()
        
        
    if redo == True:
        resParFit = {}

        bounds = bns

        for i in range(len(datsAll)):
            resParFit[str(i)] = {}
            for j in range(3):
                result = differential_evolution(calc_method, bounds, args=(datsAll[str(i)][str(j)], timsAll[str(i)][str(j)]), maxiter=mits, strategy=sts, disp = False)
                resParFit[str(i)][str(j)] = result
                
        
        with open('Processed/'+intNM2+'.p', 'wb') as fp:
            pickle.dump(resParFit, fp, protocol=pickle.HIGHEST_PROTOCOL)

        kpls = np.zeros((len(resParFit)*3, 4))

        tmp = 0
        for i in range(len(datsAll)):
            for j in range(3):
                kpls[j+tmp,::] = resParFit[str(i)][str(j)].x
            tmp += 3
            
        clnm = []
        for i in range(len(datsAll)):
            for j in range(3):
                clnm.append('Exp '+str(i+1)+', Rep '+str(j+1))

        kps = pd.DataFrame(np.transpose(kpls), index=['k_PL','k_LP','T1_P','T1_L'], columns=clnm)
        
        kps.to_csv(r'Processed\\'+intNM2+'.csv')
        
    else:

        with open('Processed/'+intNM2+'.p', 'rb') as fp:
            resParFit = pickle.load(fp)
            
        kpls = np.zeros((len(resParFit)*3, 4))

        tmp = 0
        for i in range(len(datsAll)):
            for j in range(3):
                kpls[j+tmp,::] = resParFit[str(i)][str(j)].x
            tmp += 3
            
        clnm = []
        for i in range(len(datsAll)):
            for j in range(3):
                clnm.append('Exp '+str(i+1)+', Rep '+str(j+1))

        kps = pd.DataFrame(np.transpose(kpls), index=['k_PL','k_LP','T1_P','T1_L'], columns=clnm)
                    
    return(resParFit, kpls, kps)

def plotSims(ipre, datsAll, timsAll, resParFit, exps):

    i = ipre-1
    
    fig, ax = plt.subplots(3, 1, figsize=(7,12))

    for j in range(3):

        # Simulate Experiments
        Y0 = [datsAll[str(i)][str(j)]['Pyr'][0], datsAll[str(i)][str(j)]['Lac'][0]]
        tsim = np.arange(timsAll[str(i)][str(j)][0],  timsAll[str(i)][str(j)][-1::], 0.5)
            
        result = odeint(LDH, Y0, tsim, args=(resParFit[str(i)][str(j)].x,))
        pyr = [result[x][0] for x in range(len(result))]
        lac = [result[x][1] for x in range(len(result))]
        
        
        # PLOT
        ax[j].plot(tsim, (pyr- np.min(datsAll[str(i)][str(j)]['Pyr'])) / (datsAll[str(i)][str(j)]['Pyr'][0] - np.min(datsAll[str(i)][str(j)]['Pyr'])) *100, 
                color = 'blue', alpha=0.5, label = 'Simulation Pyr')
        
        ax[j].scatter(timsAll[str(i)][str(j)],  
                        (datsAll[str(i)][str(j)]['Pyr'] - np.min(datsAll[str(i)][str(j)]['Pyr'])) / (datsAll[str(i)][str(j)]['Pyr'][0] - np.min(datsAll[str(i)][str(j)]['Pyr'])) *100, 
                        s = 25, color='blue')

        ax[j].set_xlabel('time (s)')
        ax[j].set_ylabel('Norm. Pyr (A.U.)', color='blue')

        if j < 3:
            ax[j].get_xaxis().set_visible(False)
            
        if j == 0:
            ax[j].set_title(exps[i])

        ax[j].set_xlim(-10, 300)
        ax[j].legend()

        ax2 = ax[j].twinx()
        
        ax2.plot(tsim, (lac- np.min(datsAll[str(i)][str(j)]['Pyr'])) / (datsAll[str(i)][str(j)]['Pyr'][0] - np.min(datsAll[str(i)][str(j)]['Pyr'])) *100, 
                color = 'red', alpha=0.5, label = 'Simulation Lac')
        
        ax2.scatter(timsAll[str(i)][str(j)],  
                    (datsAll[str(i)][str(j)]['Lac'] - np.min(datsAll[str(i)][str(j)]['Pyr'])) / (datsAll[str(i)][str(j)]['Pyr'][0] - np.min(datsAll[str(i)][str(j)]['Pyr'])) *100, 
                    s = 25, color='red')
        
        mxlm = np.max([(np.max(datsAll[str(i)][str(j)]['Lac']) - np.min(datsAll[str(i)][str(j)]['Pyr'])) / (datsAll[str(i)][str(j)]['Pyr'][0] - np.min(datsAll[str(i)][str(j)]['Pyr'])) *105,
                      (np.max(lac) - np.min(datsAll[str(i)][str(j)]['Pyr'])) / (datsAll[str(i)][str(j)]['Pyr'][0] - np.min(datsAll[str(i)][str(j)]['Pyr'])) *105])
        
        ax2.set_ylim(-0.2, mxlm)
        ax2.set_ylabel('Norm. Lac (A.U.)', color='red') 
        

    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    plt.show()



def inspecSims(datsAll, timsAll, resParFit, exps):    
    
    mm=interact(
        plotSims,
        ipre  = widgets.IntSlider(min=1,max=len(datsAll),step=1,value=0, readout = True, description ='Exp.:', layout=Layout(width='1000px')),
        datsAll = fixed(datsAll), 
        timsAll = fixed(timsAll),
        resParFit = fixed(resParFit), 
        exps = fixed(exps)
    )
    

#####################################################################################################################


def plotCSIDat2(CSIs, savNM="Test", exp = 1, scF = [1,1], lw1 = 1, mnn = -0.1, mxx = 1.1, xzm = [0,1], sav = False):
    sups = {}

    sizDt = np.shape(CSIs[str(exp-1)])
    mxdt = np.max(CSIs[str(exp-1)])
    mndt = np.min(CSIs[str(exp-1)])
    
    fig = plt.figure(figsize=(scF[0],scF[1]))

    for i in range(sizDt[2]):
        for j in range(sizDt[1]):
            sups[str(i)+'-'+str(j)] = plt.subplot2grid((sizDt[2], sizDt[1]), (i, j), colspan=1, rowspan=1)
            
            sups[str(i)+'-'+str(j)].get_xaxis().set_visible(False)
            sups[str(i)+'-'+str(j)].get_yaxis().set_visible(False)
            sups[str(i)+'-'+str(j)].set_ylim(-0.1,1.1)
            sups[str(i)+'-'+str(j)].plot((CSIs[str(exp-1)][xzm[0]:xzm[1],j,i] - mndt) / (mxdt - mndt),  color="#69dd52ff", linewidth = lw1)
            sups[str(i)+'-'+str(j)].spines[['right', 'top', 'left', 'bottom']].set_color('#e78829ff')
            sups[str(i)+'-'+str(j)].set_ylim(mnn,mxx)

    
    plt.subplots_adjust(left=0.1, right=0.9, 
                    top=0.9, bottom=0.1, 
                    wspace=0, hspace=0)
    
    
    if sav == True:
        plt.savefig('Processed/'+savNM+'CSIGrid_Exp'+str(exp)+'.svg', transparent=True)
        plt.savefig('Processed/'+savNM+'CSIGrid_Exp'+str(exp)+'.png', transparent=True)
    plt.show()
    
    
    
def genCSI(CSIs, savNM):    
    
    jsm = np.max([np.shape(CSIs[str(s)])[1] for s in range(len(CSIs))])   
    ism = np.max([np.shape(CSIs[str(s)])[2] for s in range(len(CSIs))])
        
    mm=interact(
        plotCSIDat2,
        savNM = fixed(savNM),
        CSIs = fixed(CSIs),
        exp = widgets.Dropdown(options=[list(range(len(CSIs)))[g]+1 for g in range(len(CSIs))], value = 1, description = "Exp.:"), 
        scF = widgets.FloatRangeSlider(min=1,max=100,step=0.000001,value=[5,5], readout = True, description ='Im. Size', layout=Layout(width='1400px')), 
        lw1 = widgets.FloatSlider(min=0.01,max=5,step=0.000001,value=1, readout = True, description ='LW CSI', layout=Layout(width='1400px')), 
        mnn = widgets.FloatSlider(min=-0.01,max=1.1,step=0.000001,value=-0.01, readout = True, description ='Z. Bot.', layout=Layout(width='1400px')), 
        mxx = widgets.FloatSlider(min=-0.01,max=1.1,step=0.000001,value=1.1, readout = True, description ='Z. Top.', layout=Layout(width='1400px')),
        xzm = widgets.IntRangeSlider(min=0 ,max=np.shape(CSIs[str(0)])[0],step=1,value=[0, np.shape(CSIs[str(0)])[0]], readout = False, layout=Layout(width='1400px'), description ='ZoomX'),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON")
    )
    

#####################################################################################################################




































def sumVox1(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 2, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,0:4,::], 1)
    tmp2[::,1,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,4::,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,0:4],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,4:8],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,8:12],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,12::],2)
    
    return(tmp)

def sumVox2(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 2, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,0:3,::], 1)
    tmp2[::,1,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,3::,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,0:4],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,4:8],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,8:12],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,12::],2)
    
    return(tmp)

def sumVox3(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 2, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,0:3,::], 1)
    tmp2[::,1,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,3::,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,0:3],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,3:7],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,7:11],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,11::],2)
    
    return(tmp)

def sumVox4(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 2, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,0:4,::], 1)
    tmp2[::,1,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,4::,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,0:3],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,3:8],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,8:11],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,11::],2)
    
    return(tmp)


def sumVox5(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 1, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,::,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,0:4],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,4:8],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,8:12],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,12::],2)
    
    return(tmp)

def sumVox6(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 1, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,::,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,0:5],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,5:9],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,9:13],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,13::],2)
    
    return(tmp)

def sumVox7(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 1, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,::,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,0:5],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,5:9],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,9:12],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,12::],2)
    
    return(tmp)

def sumVox8(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 1, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,::,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,0:4],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,4:8],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,8:12],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,12::],2)
    
    return(tmp)

def sumVox9(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 2, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,0:5,::], 1)
    tmp2[::,1,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,5::,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,0:3],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,3:7],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,7:11],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,12::],2)
    
    return(tmp)

def sumVox10(siM, tmp, CSIs_FC_NRMRW, ex):
    tmp2 = np.zeros([siM[0], 1, siM[2]])
    tmp2[::,0,::] = np.sum(CSIs_FC_NRMRW[str(ex)][::,2:-1,::], 1)

    tmp[::,::,0] = np.sum(tmp2[::,::,1:5],2)
    tmp[::,::,1] = np.sum(tmp2[::,::,5:9],2)
    tmp[::,::,2] = np.sum(tmp2[::,::,9:12],2)
    tmp[::,::,3] = np.sum(tmp2[::,::,12::],2)
    
    return(tmp)