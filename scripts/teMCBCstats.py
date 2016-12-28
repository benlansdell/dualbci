import scipy.io
import numpy as np
import scipy.stats 
import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as mpl

a = scipy.io.loadmat('./worksheets/2016_06_10-resultsforpaper/teMCBCstats.mat')
img_out = './worksheets/2016_06_10-resultsforpaper/'

sns.set_style('white')

maxte = 0.0005
minte = -maxte

dfMCBC = pd.DataFrame({'MC':pd.Series(np.squeeze(a['teMC'])),\
	'BC':pd.Series(np.squeeze(a['teBC']))})
dfMCBC = dfMCBC[dfMCBC.MC < maxte]
dfMCBC = dfMCBC[dfMCBC.BC < maxte]

dfMCMC2 = pd.DataFrame({'MC':pd.Series(np.squeeze(a['teMCb'])),\
	'MC2':pd.Series(np.squeeze(a['teMC2']))})
dfMCMC2 = dfMCMC2[dfMCMC2.MC < maxte]
dfMCMC2 = dfMCMC2[dfMCMC2.MC2 < maxte]

dfMCDC = pd.DataFrame({'MC':pd.Series(np.squeeze(a['teMCc'])),\
	'DC':pd.Series(np.squeeze(a['teDC']))})
dfMCDC = dfMCDC[dfMCDC.MC < maxte]
dfMCDC = dfMCDC[dfMCDC.DC < maxte]

g1 = sns.jointplot("MC", "BC", data = dfMCBC, kind = "kde", color="b", shade = \
	True, xlim=(0,maxte), ylim=(0,maxte), n_levels = 10, z_order = 0)
sns.plt.show()
g1.savefig(img_out + 'seaborn_MCBC_scatter.eps')

g2 = sns.jointplot("MC", "MC2", data = dfMCMC2, kind = "kde", color="b",shade =\
	True, xlim=(0,maxte), ylim=(0,maxte), n_levels = 10, z_order = 0)
#sns.plt.show()
g2.savefig(img_out + 'seaborn_MCMC2_scatter.eps')

g3 = sns.jointplot("MC", "DC", data = dfMCDC, kind = "kde", color="b", shade = \
	True, xlim=(0,maxte), ylim=(0,maxte), n_levels = 10, z_order = 0)
#sns.plt.show()
g3.savefig(img_out + 'seaborn_MCDC_scatter.eps')

diffTEs = pd.DataFrame({'MCMC2':pd.Series(np.squeeze(a['teMCb']-a['teMC2'])),\
	'MCBC':pd.Series(np.squeeze(a['teMC']-a['teBC'])), 'MCDC':pd.Series(np.squeeze(a['teMCc'] - a['teDC']))})
diffTEs = diffTEs[abs(diffTEs.MCMC2) < maxte]
diffTEs = diffTEs[abs(diffTEs.MCBC) < maxte]
diffTEs = diffTEs[abs(diffTEs.MCDC) < maxte]

#diffTEs.MCMC2 = np.log(abs(diffTEs.MCMC2))
#diffTEs.MCBC = np.log(abs(diffTEs.MCBC))
#diffTEs.MCDC = np.log(abs(diffTEs.MCDC))

#Violin plots of changes
f, ax = mpl.subplots()
sns.boxplot(data = diffTEs);
sns.despine(offset = 10, trim = True);
sns.plt.show()

N = 3
ind = np.arange(N)
means = np.array(np.mean(diffTEs))
fig, ax = plt.subplots()
rects1 = ax.bar(ind, means, width, color='r', yerr=stds)