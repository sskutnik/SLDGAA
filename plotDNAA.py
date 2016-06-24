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
OUTPUT_FORMAT = 'png'
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

def process_plot(fileName):

	sns.despine(top=True)
	plt.tight_layout()

	if(SHOW_PLOTS):
		plt.show()
	else:
		plt.savefig(IMAGE_DEST + fileName + '.' + OUTPUT_FORMAT)
		plt.close()		
		
#================================================
# Singles data plotting
#================================================

# Fetch in data and do some cleanup
singles = pd.read_csv('singles.csv')
singles['C/E'] = singles['Calculated mass (ng)']/singles['Certified mass (ng)']
singles['UncertaintyProp'] = np.sqrt(singles['Certified uncertainty']**2 + singles['Calculated uncertainty']**2)

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
g = (g.map(plt.errorbar,'SampleNum',"C/E",'UncertaintyProp',ls='',ms=9).add_legend())
g.fig.get_axes()[0].legend(loc='lower left')
(g.set_titles("{col_name}").set_xlabels("Sample"))

# Fix up x-axes to just show samples of interest
axes = g.axes
axes[0].set_xlim(0.1,4.9)
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

binaries['C/E'] = binaries['Calculated mass (ng)']/binaries['Certified mass (ng)']
binaries['UncertaintyProp'] = np.sqrt(binaries['Certified mass uncert.']**2 + binaries['Calculated mass uncert.']**2)

binaries = binaries.sort_values(by=['Isotope','Sample ID'],ascending=True)

# Sum up fissile mass of each binary set
isFissile = (binaries['Isotope'] == "U-235") | (binaries['Isotope'] == "Pu-239") | (binaries['Isotope'] == 'U-233')
binSamples = binaries[isFissile].groupby('Sample ID')['Certified mass (ng)'].sum()
binSamples.name = 'Fissile mass (ng)'

# Define isotope IDs to map out plots
isoDefs = dict([('Pu-239',0), ('U-233',1), ('U-235',2), ('U-238',3)])
isoNumbers = []
fissileMasses = []

for index, row in binaries.iterrows():
	fissileMasses.append(binSamples[row['Sample ID']])
	isoNumbers.append(isoDefs[row['Isotope']])

# Add isotope ID number (for indexing / mapping plots) to table
binaries['IsoNum'] = pd.Series(isoNumbers,index=binaries.index)
# Add total fissile mass to table (for plot titles)
binaries['Fissile mass (ng)'] = pd.Series(fissileMasses,index=binaries.index)


# WORKING
binariesCalib = binaries.loc[isCalibrated]
g = sns.factorplot(data=binariesCalib,kind='point',col='Sample ID',x='Isotope',y='C/E',hue='Isotope',legend=True,markers=['s','^','o','v'],scale=1.0,size=4,col_wrap=GLOB_COL_WRAP)

myColors = sns.color_palette()
hexColors = []
for color in myColors:
	hexColors.append(matplotlib.colors.rgb2hex(color))

binariesPu = binariesCalib.loc[binaries['Isotope'] == "Pu-239"]
#print(binariesPu['C/E'])

# This is plotting the whole data set and not just the selection; why?
# IDEA: Iterate over hexColors and IsoNums; plot each isonum as a different *series*
g.map(plt.errorbar,"IsoNum","C/E",'UncertaintyProp',data=binariesPu,ls='',ms=0,zorder=0,elinewidth=2,color="gray")

g.set_xlabels('Isotope').set_titles("{col_name}")
  
titles = binaries['Fissile mass (ng)'].unique()
for ax, massTitle in zip(g.axes.flat, titles):
    ax.set_title("Sample {:s}\nFissile mass: {:.2f} ng".format(ax.get_title(),massTitle) )

process_plot('binaries')

#================================================
# IAEA swipe samples plotting
#================================================

df_swipes = pd.read_csv('IAEA_swipes.csv')

df_swipes['C/E'] = df_swipes['Calculated mass (ng)']/df_swipes['Certified mass (ng)']
df_swipes['UncertaintyProp'] = np.sqrt(df_swipes['Certified mass uncertainty']**2 + df_swipes['Calculated mass uncertainty']**2)

df_swipes = df_swipes.sort_values(by=['Isotope','Sample'],ascending=True)

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


# WORKING
g = sns.factorplot(data=df_swipes,kind='point',col='Sample',x='Isotope',y='C/E',hue='Isotope',legend=True,markers=['s','^','o','v'],scale=1.0,size=4,col_wrap=GLOB_COL_WRAP)

myColors = sns.color_palette()
hexColors = []
for color in myColors:
	hexColors.append(matplotlib.colors.rgb2hex(color))


g.map(plt.errorbar,"IsoNum","C/E",'UncertaintyProp',data=df_swipes,ls='',ms=0,zorder=0,elinewidth=2,color="gray")
g.set_xlabels('Isotope').set_titles("{col_name}")
  
titles = df_swipes['Fissile mass (ng)'].unique()
for ax, massTitle in zip(g.axes.flat, titles):
    ax.set_title("Sample {:s}\nFissile mass: {:.2f} ng".format(ax.get_title(),massTitle) )

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
palette = itertools.cycle(sns.color_palette())

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
		CE_propUnc = math.sqrt(sum(map(lambda x:x*x,df[nucLocs]["Uncertainty"])))

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
sns.set(style="ticks",font_scale=2.0)
g = sns.factorplot(data=dfCE,kind='point',row='Parent',x='Nuclide',y='C/E mean',hue='Parent',legend=True,markers=['s','^','o','v'],scale=1.5,size=5,aspect=3,linestyles='None')
g = (g.map(plt.errorbar,'IsoNum',"C/E mean",'Uncertainty',ls='',color='grey',ms=0,zorder=0,elinewidth=2.5))
g.set_xlabels('Nuclide').set_titles("{row_name}")
		
process_plot('FP_errors')