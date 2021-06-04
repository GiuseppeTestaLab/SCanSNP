#!/usr/bin/env python


import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import warnings
from scipy import stats
import numpy as np
import seaborn as sns



# def FittedMixturePlot(gm, outdir,comb, NoisedIDs):
# 	with warnings.catch_warnings():
# 		warnings.simplefilter('ignore')
		
# 		#Mixture Plot
# 		plt.clf()
# 		figure(num=None, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
		
# 		C1 = stats.norm(gm.means_[0][0], np.sqrt(gm.covariances_[0][0][0]))
# 		C2 = stats.norm(gm.means_[1][0], np.sqrt(gm.covariances_[1][0][0]))
		
		
# 		mc = gm.weights_
		
# 		x = np.linspace(min(NoisedIDs["logFC"]), max(NoisedIDs["logFC"]), 501)
		
		
# 		C1 = C1.pdf(x) * mc[0]
# 		C2 = C2.pdf(x) * mc[1]
		
		
# 		gm.means_[0][0] > gm.means_[1][0]
		
# 		sns.distplot(NoisedIDs["logFC"].to_numpy(), hist=False, label='Mixture')
# 		plt.plot(x, C1,'--', label='Fitted LowQuality' if gm.means_[0][0] < gm.means_[1][0] else 'Fitted_HighQuality' )
# 		plt.plot(x, C2,'--',  label='Fitted LowQuality' if gm.means_[0][0] > gm.means_[1][0] else 'Fitted_HighQuality' )
		
# 		plt.legend(prop={'size': 16}, title = 'Variants', fontsize=25)
# 		plt.title('logFC Density Plot ', fontsize=25)
# 		plt.xlabel("ID1-ID2_logFC", fontsize=25)
# 		plt.ylabel('Density', fontsize=25)
		
# 		plt.savefig(outdir+"/"+"_".join(list(comb))+ ".LowQ_MixtureFitting.svg")




def FittedMixturePlot(gm, outdir, NoisedIDs):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		
		#Mixture Plot
		plt.clf()
		figure(num=None, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
		
		C1 = stats.norm(gm.means_[0][0], np.sqrt(gm.covariances_[0][0][0]))
		C2 = stats.norm(gm.means_[1][0], np.sqrt(gm.covariances_[1][0][0]))
		
		
		mc = gm.weights_
		
		x = np.linspace(min(NoisedIDs["logFC"]), max(NoisedIDs["logFC"]), 501)
		
		
		C1 = C1.pdf(x) * mc[0]
		C2 = C2.pdf(x) * mc[1]
		
		
		gm.means_[0][0] > gm.means_[1][0]
		
		sns.distplot(NoisedIDs["logFC"].to_numpy(), hist=False, label='Mixture')
		plt.plot(x, C1,'--', label='Fitted LowQuality' if gm.means_[0][0] < gm.means_[1][0] else 'Fitted_HighQuality' )
		plt.plot(x, C2,'--',  label='Fitted LowQuality' if gm.means_[0][0] > gm.means_[1][0] else 'Fitted_HighQuality' )
		
		plt.legend(prop={'size': 16}, title = 'Variants', fontsize=25)
		plt.title('logFC Density Plot ', fontsize=25)
		plt.xlabel("ID1-ID2_logFC", fontsize=25)
		plt.ylabel('Density', fontsize=25)
		
		plt.savefig(outdir+"/MixtureFitting.svg")