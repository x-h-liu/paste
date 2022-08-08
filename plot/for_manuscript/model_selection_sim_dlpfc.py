import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

"""
Code for plotting model selection for simulated data
"""
plt.rcParams.update({'font.size': 22})


res_df = pd.DataFrame(columns=['Overlap','Delta', 'Estimation'])


res_df.loc[len(res_df)] = ['90%', '0.1', 0.9]
res_df.loc[len(res_df)] = ['90%', '1.0', 0.8]
res_df.loc[len(res_df)] = ['90%', '2.0', 0.8]
res_df.loc[len(res_df)] = ['90%', '3.0', 0.8]

res_df.loc[len(res_df)] = ['70%', '0.1', 0.7]
res_df.loc[len(res_df)] = ['70%', '1.0', 0.7]
res_df.loc[len(res_df)] = ['70%', '2.0', 0.6]
res_df.loc[len(res_df)] = ['70%', '3.0', 0.7]

res_df.loc[len(res_df)] = ['50%', '0.1', 0.5]
res_df.loc[len(res_df)] = ['50%', '1.0', 0.5]
res_df.loc[len(res_df)] = ['50%', '2.0', 0.8]
res_df.loc[len(res_df)] = ['50%', '3.0', 0.8]

res_df.loc[len(res_df)] = ['30%', '0.1', 0.3]
res_df.loc[len(res_df)] = ['30%', '1.0', 0.3]
res_df.loc[len(res_df)] = ['30%', '2.0', 0.9]
res_df.loc[len(res_df)] = ['30%', '3.0', 0.7]


g = sns.catplot(x="Delta", y="Estimation", col="Overlap",data=res_df,kind="bar", ci=None, aspect=1,legend=False,
                palette=sns.color_palette())
ax1, ax2, ax3, ax4 = g.axes[0]
ax1.axhline(0.9, color='r', ls='--')
ax2.axhline(0.7, color='r', ls='--')
ax3.axhline(0.5, color='r', ls='--')
ax4.axhline(0.3, color='r', ls='--')


plt.savefig('test.png')




"""
Code for plotting model selection for dlpfc data
"""
# plt.rcParams.update({'font.size': 22})

# res_df = pd.DataFrame(columns=['Sample','Pair','Subslice','Estimation'])

# res_df.loc[len(res_df)] = [1, 'AB', 'Horizontal', 0.5]
# res_df.loc[len(res_df)] = [1, 'AB', 'Vertical', 0.5]
# res_df.loc[len(res_df)] = [1, 'BC', 'Horizontal', 0.3]
# res_df.loc[len(res_df)] = [1, 'BC', 'Vertical', 0.5]
# res_df.loc[len(res_df)] = [1, 'CD', 'Horizontal', 0.6]
# res_df.loc[len(res_df)] = [1, 'CD', 'Vertical', 0.5]

# res_df.loc[len(res_df)] = [2, 'AB', 'Horizontal', 0.5]
# res_df.loc[len(res_df)] = [2, 'AB', 'Vertical', 1.0]
# res_df.loc[len(res_df)] = [2, 'BC', 'Horizontal', 0.5]
# res_df.loc[len(res_df)] = [2, 'BC', 'Vertical', 0.5]
# res_df.loc[len(res_df)] = [2, 'CD', 'Horizontal', 0.7]
# res_df.loc[len(res_df)] = [2, 'CD', 'Vertical', 0.5]

# res_df.loc[len(res_df)] = [3, 'AB', 'Horizontal', 0.8]
# res_df.loc[len(res_df)] = [3, 'AB', 'Vertical', 0.8]
# res_df.loc[len(res_df)] = [3, 'BC', 'Horizontal', 0.8]
# res_df.loc[len(res_df)] = [3, 'BC', 'Vertical', 0.6]
# res_df.loc[len(res_df)] = [3, 'CD', 'Horizontal', 0.8]
# res_df.loc[len(res_df)] = [3, 'CD', 'Vertical', 0.8]


# g = sns.catplot(x="Pair", y="Estimation", hue='Subslice', col="Sample", data=res_df,kind="bar", ci=None, aspect=1,legend=True,
#                 palette=sns.color_palette())
# ax1, ax2, ax3 = g.axes[0]
# ax1.axhline(0.7, color='r', ls='--')
# ax2.axhline(0.7, color='r', ls='--')
# ax3.axhline(0.7, color='r', ls='--')


# plt.savefig('test.png')


