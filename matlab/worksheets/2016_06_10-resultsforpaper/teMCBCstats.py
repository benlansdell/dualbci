import scipy.io
import numpy as np
import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as mpl

a = scipy.io.loadmat('./worksheets/2016_06_10-resultsforpaper/teMCBCstats.mat')
img_out = './worksheets/2016_06_10-resultsforpaper/'

maxte = 0.0015
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

#g1 = sns.jointplot("MC", "BC", data = dfMCBC, kind = "kde", color="b", shade = False, xlim=(0,maxte), ylim=(0,maxte), n_levels = 10, z_order = 0)
g1 = sns.jointplot("MC", "BC", data = dfMCBC, color="b", xlim=(0,maxte), ylim=(0,maxte))
sns.plt.show()
g1.savefig(img_out + 'seaborn_MCBC_scatter.eps')

g2 = sns.jointplot("MC", "MC2", data = dfMCMC2, color="b", shade = False, xlim=(0,maxte), ylim=(0,maxte), n_levels = 10)
#sns.plt.show()
g2.savefig(img_out + 'seaborn_MCMC2_scatter.eps')

g3 = sns.jointplot("MC", "DC", data = dfMCDC, color="b", shade = False, xlim=(0,maxte), ylim=(0,maxte), n_levels = 10)
#sns.plt.show()
g3.savefig(img_out + 'seaborn_MCDC_scatter.eps')
