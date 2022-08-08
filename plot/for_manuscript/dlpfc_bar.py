import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt




res_df = pd.DataFrame(columns=['Sample','Pair','Kind', 'ARI'])


res_df.loc[len(res_df)] = [1, 'AB', 'PASTE2', 0.5982726987358716]
res_df.loc[len(res_df)] = [1, 'BC', 'PASTE2', 0.06780848340640959]
res_df.loc[len(res_df)] = [1, 'CD', 'PASTE2', 0.5650396516478512]
res_df.loc[len(res_df)] = [2, 'AB', 'PASTE2', 0.428116044914222]
res_df.loc[len(res_df)] = [2, 'BC', 'PASTE2', 0.6637812938623632]
res_df.loc[len(res_df)] = [2, 'CD', 'PASTE2', 0.36385628822237787]
res_df.loc[len(res_df)] = [3, 'AB', 'PASTE2', 0.5587703813466455]
res_df.loc[len(res_df)] = [3, 'BC', 'PASTE2', 0.4724182915603778]
res_df.loc[len(res_df)] = [3, 'CD', 'PASTE2', 0.39816624188317135]


g = sns.catplot(x="Pair", y="ARI", hue='Kind', col="Sample",data=res_df,kind="bar", ci=None, aspect=1,legend=False,
                palette=sns.color_palette("Set2")[:1]+sns.color_palette("Set2")[2:3]+sns.color_palette("Set2")[4:6],hue_order=['PASTE2'])
g.set_axis_labels("Pair", "ARI").set(ylim=(0, 1))

plt.show()