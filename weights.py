# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 14:31:44 2021

@author: zahhr
"""

import pandas as pd, os, seaborn as sns, matplotlib.pyplot as plt, numpy as np
from datetime import date
src = "/home/kepecs/Documents/balbc_w2_pr2_cad_weights_zmd_may-july2021.xlsx"
dst = "/home/kepecs/Documents"
df = pd.read_excel(src, None)

#%%
#pr2w2
pr2w2 = df["pr2_w2"]
plt.figure()
sns.lineplot(x = "date", y = "weight (g)", hue = "animal name", data = pr2w2)
plt.xticks(np.unique(pr2w2["date"]), rotation = 90, fontsize = 7)
plt.tick_params('both', length=10, width=0.8, which='major')
# plt.savefig(os.path.join(dst, "pr2w2_weights_plot_{}.pdf".format(date.today().strftime("%Y%m%d"))), bbox_inches = "tight")
plt.title("PR2_W2 cohort\nwater deprived")
plt.savefig(os.path.join(dst, "pr2w2_weights_plot_{}.jpg".format(date.today().strftime("%Y%m%d"))), bbox_inches = "tight")

#%%
#cadw1
cadw1 = df["cadw1"]
plt.figure()
sns.lineplot(x = "date", y = "weight (g)", hue = "animal name", data = cadw1)
plt.xticks(np.unique(cadw1["date"]), rotation = 90, fontsize = 7)
plt.tick_params('both', length=10, width=0.8, which='major')
# plt.savefig(os.path.join(dst, "cadw1_weights_plot_{}.pdf".format(date.today().strftime("%Y%m%d"))), bbox_inches = "tight")
plt.title("CADW1 cohort")
plt.savefig(os.path.join(dst, "cadw1_weights_plot_{}.jpg".format(date.today().strftime("%Y%m%d"))), bbox_inches = "tight")
#%%
#cadw2
cadw2 = df["cadw2"]
plt.figure()
sns.lineplot(x = "date", y = "weight (g)", hue = "animal name", data = cadw2)
plt.xticks(np.unique(cadw2["date"]), rotation = 90, fontsize = 7)
plt.tick_params('both', length=10, width=0.8, which='major')
# plt.savefig(os.path.join(dst, "cadw2_weights_plot_{}.pdf".format(date.today().strftime("%Y%m%d"))), bbox_inches = "tight")
plt.title("CADW2 cohort")
plt.savefig(os.path.join(dst, "cadw2_weights_plot_{}.jpg".format(date.today().strftime("%Y%m%d"))), bbox_inches = "tight")
