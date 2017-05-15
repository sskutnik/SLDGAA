import matplotlib 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import itertools
import math

sns.set(style="ticks",font_scale=1.5)

IMAGE_DEST='./images/'
GLOB_COL_WRAP = 2 # number of columns per line for small multiples
OUTPUT_FORMAT = 'pdf'
SHOW_PLOTS = False

def plot_scatterBox(df,xData,yData,title,fileName,plotAspect=1,colorVal=None):
    plt.figure(figsize=(6*plotAspect,6))
    if(colorVal):
        sns_plot = sns.boxplot(x=xData,y=yData,data=df,color=colorVal)
        sns.stripplot(x=xData,y=yData,size=9,data=df,color=colorVal,edgecolor='gray',linewidth=1)
    else:
        sns_plot = sns.boxplot(x=xData,y=yData,data=df)
        sns.stripplot(x=xData,y=yData,size=9,data=df,edgecolor='gray',linewidth=1)
    
    plt.title(title)
    fig = sns_plot.get_figure()
    process_plot(fileName)

def outline_markers(g):
    for ax in g.axes:
        for mk in ax.collections:
            mk.set_edgecolor('black')
            mk.set_linewidth(0.075)
            
def process_plot(fileName):

    sns.despine(top=True)
    plt.tight_layout()

    if(SHOW_PLOTS):
        plt.show()
    else:
        plt.savefig(IMAGE_DEST + fileName + '.' + OUTPUT_FORMAT)
        plt.close()        

def facetErrorBars(x, y, sigma, cDict, color=None, label=None, **kwargs):
    # Add error bars to a FacetGrid / FactorPlot and match the bar color to the hue
    data = kwargs.pop("data")
    hue = kwargs.pop("hue")

    ebColors = []
    for eval in data[hue].unique():
        ebColor = cDict[eval]
        xVals = data[x].loc[data[hue] == eval].values
        yVals = data[y].loc[data[hue] == eval].values
        yErrs = data[sigma].loc[data[hue] == eval].values    
        plt.errorbar(xVals, yVals, xerr=None, yerr=yErrs, ecolor=ebColor, color=None,**kwargs)
        
#================================================
# Singles data plotting
#================================================

# Fetch in data and do some cleanup
singles = pd.read_csv('singles.csv')
singles['C/E'] = singles['Calculated mass (ng)']/singles['Certified mass (ng)']
singles['UncertaintyProp'] = singles['C/E']*np.sqrt(singles['Certified uncertainty']**2 + singles['Calculated uncertainty']**2)

# Not using category type right now, but may in the future
#singles['Isotope'] = singles['Isotope'].astype('category')

singles['C/E-1'] = singles['C/E']-1.0
# Renormalize samples to sample 1-4 for each isotope 
# (Note that I should probably sort on category if the isotopes aren't sorted...)

# Rework numbering for this
sampleNums = []
for index, row in singles.iterrows():
    sampleNums.append(index%4+1)

singles['SampleNum'] = sampleNums
    
# Plot out small multiples of samples by isotope, with uncertainty
g = sns.FacetGrid(data=singles,col='Isotope',hue='Location',hue_kws=dict(marker=['s','^','o','v']),
                  hue_order=["PT-1", "PT-2 (Uncorrected)", "PT-2 (Corrected)"],
                  size=4,col_wrap=GLOB_COL_WRAP,legend_out=False)
g = (g.map(plt.errorbar,'SampleNum',"C/E",'UncertaintyProp',ls='',ms=9,capsize=6,elinewidth=1.5).add_legend())

g.fig.get_axes()[0].legend(loc='lower left')
(g.set_titles("{col_name}").set_xlabels("Sample"))

# Fix up x-axes to just show samples of interest
axes = g.axes
axes[0].set_xlim(0.1,4.9)

#outline_markers(g)
process_plot('singles')

isCalibrated = (singles["Location"] == "PT-1") | (singles["Location"] == "PT-2 (Corrected)")
SinglesCalibrated = singles[isCalibrated]

# Plot an overlay of a box plot + stripplot (individual measurements) to show spread in data by isotope
plot_scatterBox(SinglesCalibrated,"Isotope","C/E","Single-isotope measurements","singles_box",plotAspect=1.25)


#================================================
# Binary mixtures data plotting
#================================================
binaries = pd.read_csv('binaries.csv')
isCalibrated = (binaries["Location"] == "PT-2 (Corrected)")
isPT1 = (binaries["Location"] == "PT-1")

binaries['C/E'] = binaries['Calculated mass (ng)']/binaries['Certified mass (ng)']
binaries['UncertaintyProp'] = binaries['C/E']*np.sqrt(binaries['Certified mass uncert.']**2 + binaries['Calculated mass uncert.']**2)

binaries = binaries.sort_values(by=['Isotope','Sample ID'],ascending=True)

# Sum up fissile mass of each binary set
isFissile = (binaries['Isotope'] == "U-235") | (binaries['Isotope'] == "Pu-239") | (binaries['Isotope'] == 'U-233')
binSamples = binaries[isFissile].groupby(['Sample ID', 'Location'])['Certified mass (ng)'].sum().reset_index(name='Fissile mass (ng)')
#print(binSamples)

# Define isotope IDs to map out plots
isoDefs = dict([('Pu-239',0), ('U-233',1), ('U-235',2), ('U-238',3)])
isoSyms = {'Pu-239' : 's', 'U-233' : '^', 'U-235' : 'o', 'U-238' : 'v'} 
isoNumbers = []
fissileMasses = []

# Define a dictionary to consistently map colors / symbols to isotopes
binSyms = []
palette = itertools.cycle(sns.color_palette())
ebDict = { }    
for iso in binaries['Isotope'].unique():
    ebDict[iso] = next(palette)
    binSyms.append(isoSyms[iso])

# Find the fissile masses corresponding to each unique sample ID
for index, row in binaries.iterrows():
    sampleMass = binSamples.loc[binSamples['Sample ID'] == row['Sample ID']]['Fissile mass (ng)'].iloc[0]
    fissileMasses.append(sampleMass)
    isoNumbers.append(isoDefs[row['Isotope']])

# Add isotope ID number (for indexing / mapping plots) to table
binaries['IsoNum'] = pd.Series(isoNumbers,index=binaries.index)
# Add total fissile mass to table (for plot titles)
binaries['Fissile mass (ng)'] = pd.Series(fissileMasses,index=binaries.index)

#print(binaries.columns.values)
print(binaries[['Sample ID', 'Fissile mass (ng)']])
binariesCalib = binaries.loc[isCalibrated]

g = sns.factorplot(data=binariesCalib,kind='point',col='Sample ID',x='Isotope',y='C/E',hue='Isotope',
   legend=True,markers=binSyms,scale=1.0,size=4,col_wrap=GLOB_COL_WRAP)

g.map_dataframe(facetErrorBars,"IsoNum","C/E","UncertaintyProp",ebDict,data=binariesCalib, hue='Isotope',
       markeredgecolor='k',ls='',elinewidth=1.5,capsize=6,zorder=0)       

g.set_xlabels('Isotope').set_titles("{col_name}")

titles = binaries['Fissile mass (ng)'].unique()
#titles = binSamples
for ax, massTitle in zip(g.axes.flat, titles):
    ax.set_title("Sample {:s}\nFissile mass: {:.2f} ng".format(ax.get_title(),massTitle) )

#outline_markers(g)
process_plot('binaries_PT2')

#==========================
# PT-1 binaries data
#===========================

# WORKING
binariesPT1 = binaries.loc[isPT1]
pt1Syms = []
for iso in binariesPT1['Isotope'].unique():
    pt1Syms.append(isoSyms[iso])
    
g = sns.factorplot(data=binariesPT1,kind='point',col='Sample ID',x='Isotope',y='C/E',hue='Isotope',
   legend=True,markers=pt1Syms,scale=1.0,size=4,col_wrap=GLOB_COL_WRAP)

g.map_dataframe(facetErrorBars,"IsoNum","C/E","UncertaintyProp",ebDict,data=binariesPT1, hue='Isotope',
       markeredgecolor='k',ls='',elinewidth=1.5,capsize=6,zorder=0)
       
g.set_xlabels('Isotope').set_titles("{col_name}")
  
#titles = binaries['Fissile mass (ng)'].unique()
for ax, massTitle in zip(g.axes.flat, titles):
    ax.set_title("Sample {:s}\nFissile mass: {:.2f} ng".format(ax.get_title(),massTitle) )

#outline_markers(g)
process_plot('binaries_PT1')


#================================================
# IAEA swipe samples plotting
#================================================

df_swipes = pd.read_csv('IAEA_swipes.csv')

df_swipes['C/E'] = df_swipes['Calculated mass (ng)']/df_swipes['Certified mass (ng)']
df_swipes['UncertaintyProp'] = df_swipes['C/E']*np.sqrt(df_swipes['Certified mass uncertainty']**2 + df_swipes['Calculated mass uncertainty']**2)

df_swipes = df_swipes.sort_values(by=['Isotope','Sample'],ascending=True)

iaeaPalette = []
iaeaSyms = []
for iso in df_swipes["Isotope"].unique():
    iaeaPalette.append(ebDict[iso])
    iaeaSyms.append(isoSyms[iso])
    
# Sum up fissile mass of each binary set
isFissile = (df_swipes['Isotope'] == "U-235") | (df_swipes['Isotope'] == "Pu-239") | (df_swipes['Isotope'] == 'U-233')
iaeaSamples = df_swipes[isFissile].groupby('Sample')['Certified mass (ng)'].sum()
iaeaSamples.name = 'Fissile mass (ng)'

# Define isotope IDs to map out plots
isoDefs = dict([('U-235',0), ('U-238',1)])
isoNumbers = []
fissileMasses = []

for index, row in df_swipes.iterrows():
    fissileMasses.append(iaeaSamples[row['Sample']])
    isoNumbers.append(isoDefs[row['Isotope']])

# Add isotope ID number (for indexing / mapping plots) to table
df_swipes['IsoNum'] = pd.Series(isoNumbers,index=df_swipes.index)
# Add total fissile mass to table (for plot titles)
df_swipes['Fissile mass (ng)'] = pd.Series(fissileMasses,index=df_swipes.index)
    
g = sns.factorplot(data=df_swipes,kind='point',col='Sample',x='Isotope',y='C/E',
   palette=iaeaPalette,markers=iaeaSyms,hue='Isotope',legend=True,scale=1.0,size=4,col_wrap=GLOB_COL_WRAP)
g.map_dataframe(facetErrorBars,"IsoNum","C/E","UncertaintyProp",ebDict,data=df_swipes,hue='Isotope',
       markeredgecolor='k',ls='',elinewidth=1.5,capsize=6,zorder=0)
g.set_xlabels('Isotope').set_titles("{col_name}")
  
titles = df_swipes['Fissile mass (ng)'].unique()
for ax, massTitle in zip(g.axes.flat, titles):
    ax.set_title("Sample {:s}\nFissile mass: {:.2f} ng".format(ax.get_title(),massTitle) )

#outline_markers(g)
process_plot('IAEA_swipes')

#================================================
# Fission product recoveries plotting
#================================================

FPs = pd.read_csv('FP_recoveries.csv')

FPs['C/E'] = FPs['Comparator counts/ng']/FPs['Sample counts/ng']
isCalibrated = (FPs["Location"] == "PT-1") | (FPs["Location"] == "PT-2 (Corrected)")
FP_calibrated = FPs[isCalibrated]

# Plot an overlay of a box plot + stripplot (individual measurements) to show spread in data by isotope

isU235 = (FP_calibrated["Parent"] == "U-235")
isU233 = (FP_calibrated["Parent"] == "U-233")
isPu239 = (FP_calibrated["Parent"] == "Pu-239")

FP_U235 = FP_calibrated[isU235]
FP_U233 = FP_calibrated[isU233]
FP_Pu239 = FP_calibrated[isPu239]

palette = itertools.cycle(sns.color_palette())

plot_scatterBox(FP_Pu239,"Nuclide","C/E","Pu-239 fission product indicators","Pu239_FP",plotAspect=2.0,colorVal=next(palette))
plot_scatterBox(FP_U233,"Nuclide","C/E","U-233 fission product indicators","U233_FP",plotAspect=2.0,colorVal=next(palette))
plot_scatterBox(FP_U235,"Nuclide","C/E","U-235 fission product indicators","U235_FP",plotAspect=2.0,colorVal=next(palette))

dataSets = [FP_Pu239, FP_U233, FP_U235]

CE_cols = ["Nuclide","Parent","C/E mean","Uncertainty"]
dfCE = pd.DataFrame(columns=CE_cols)

for df in dataSets:
    measuredFPs = df["Nuclide"].unique()
    avgCE = []
    uncertCE = []
    tmpParent = []
    tmpFPs = []
    
    # Get the mean C/E and propagated uncertainty for each fission product
    for nuc in measuredFPs:
        nucLocs = (df["Nuclide"] == nuc)
        CE_mean = np.mean(df[nucLocs]["C/E"])
        # Sum the square of each uncertainty value for each measurement; 
        # then take the square root to get propagated uncertainty
        CE_propUnc = CE_mean*math.sqrt(sum(map(lambda x:x*x,df[nucLocs]["Uncertainty"])))

        tmpUnc = CE_propUnc
        dfTmp = pd.DataFrame([[nuc,df[nucLocs]["Parent"].iloc[0],CE_mean,tmpUnc]],columns=CE_cols)
        dfCE = dfCE.append(dfTmp,ignore_index=True)

# Define isotope IDs to map out plots
idx = 0
isoDefs = { }
for nuc in dfCE["Nuclide"].unique():
        isoDefs[nuc] = idx
        idx = idx + 1

isoNumbers = []

for index, row in dfCE.iterrows():
    isoNumbers.append(isoDefs[row['Nuclide']])

# Add isotope ID number (for indexing / mapping plots) to table
dfCE['IsoNum'] = pd.Series(isoNumbers,index=dfCE.index)

ceSyms = []
for iso in dfCE['Parent'].unique():
    ceSyms.append(isoSyms[iso])

sns.set(style="ticks",font_scale=2.0)
g = sns.factorplot(data=dfCE,kind='point',row='Parent',x='Nuclide',y='C/E mean',
   hue='Parent',legend=True,markers=ceSyms,scale=1.5,size=5,aspect=3,linestyles='None')
g.map_dataframe(facetErrorBars,"IsoNum","C/E mean","Uncertainty",ebDict,data=dfCE, hue='Parent',
       markeredgecolor='k',ls='',elinewidth=1.5,capsize=6,zorder=0)

#outline_markers(g)               
g.set_xlabels('Nuclide').set_titles("{row_name}")

process_plot('FP_errors')
