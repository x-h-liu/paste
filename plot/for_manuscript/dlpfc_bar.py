import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt




plt.rcParams.update({'font.size': 22})
# plt.figure(figsize=(20,5))

res_df = pd.DataFrame(columns=['Sample','Pair','Kind', 'ARI'])


res_df.loc[len(res_df)] = [1, 'AB', 'PASTE2', 0.5982726987358716]
res_df.loc[len(res_df)] = [1, 'BC', 'PASTE2', 0.06780848340640959]
res_df.loc[len(res_df)] = [1, 'CD', 'PASTE2', 0.5650396516478512]
res_df.loc[len(res_df)] = [2, 'AB', 'PASTE2', 0.6479499641985342]
res_df.loc[len(res_df)] = [2, 'BC', 'PASTE2', 0.6637812938623632]
res_df.loc[len(res_df)] = [2, 'CD', 'PASTE2', 0.36385628822237787]
res_df.loc[len(res_df)] = [3, 'AB', 'PASTE2', 0.5587703813466455]
res_df.loc[len(res_df)] = [3, 'BC', 'PASTE2', 0.4724182915603778]
res_df.loc[len(res_df)] = [3, 'CD', 'PASTE2', 0.39816624188317135]

res_df.loc[len(res_df)] = [1, 'AB', 'PASTE', 0.1807979710315713]
res_df.loc[len(res_df)] = [1, 'BC', 'PASTE', 0.04244505600543241]
res_df.loc[len(res_df)] = [1, 'CD', 'PASTE', 0.18117888497485246]
res_df.loc[len(res_df)] = [2, 'AB', 'PASTE', 0.428116044914222]
res_df.loc[len(res_df)] = [2, 'BC', 'PASTE', 0.6496472030858088]
res_df.loc[len(res_df)] = [2, 'CD', 'PASTE', 0.12121316738977594]
res_df.loc[len(res_df)] = [3, 'AB', 'PASTE', 0.1989594553422047]
res_df.loc[len(res_df)] = [3, 'BC', 'PASTE', 0.20166888798045457]
res_df.loc[len(res_df)] = [3, 'CD', 'PASTE', 0.1924570862190806]

res_df.loc[len(res_df)] = [1, 'AB', 'Pamona', 0.1062529968268354]
res_df.loc[len(res_df)] = [1, 'BC', 'Pamona', 0.09997692561863605]
res_df.loc[len(res_df)] = [1, 'CD', 'Pamona', 0.23667879508822984]
res_df.loc[len(res_df)] = [2, 'AB', 'Pamona', 0.17103816665343557]
res_df.loc[len(res_df)] = [2, 'BC', 'Pamona', 0.1499653462148793]
res_df.loc[len(res_df)] = [2, 'CD', 'Pamona', 0.002926864537740061]
res_df.loc[len(res_df)] = [3, 'AB', 'Pamona', 0.15037046488108843]
res_df.loc[len(res_df)] = [3, 'BC', 'Pamona', 0.12804688846617474]
res_df.loc[len(res_df)] = [3, 'CD', 'Pamona', 0.18540193182333908]

res_df.loc[len(res_df)] = [1, 'AB', 'Tangram', 0.004754446594825131]
res_df.loc[len(res_df)] = [1, 'BC', 'Tangram', 0.06370320700808636]
res_df.loc[len(res_df)] = [1, 'CD', 'Tangram', 0.04665607074322447]
res_df.loc[len(res_df)] = [2, 'AB', 'Tangram', 0.10419877335390436]
res_df.loc[len(res_df)] = [2, 'BC', 'Tangram', 0.14582357052085756]
res_df.loc[len(res_df)] = [2, 'CD', 'Tangram', 0.0029407457628681756]
res_df.loc[len(res_df)] = [3, 'AB', 'Tangram', 0.056725607904686295]
res_df.loc[len(res_df)] = [3, 'BC', 'Tangram', 0.04791785266554903]
res_df.loc[len(res_df)] = [3, 'CD', 'Tangram', 0.05813060253275463]




# g = sns.catplot(x="Pair", y="ARI", hue='Kind', col="Sample",data=res_df,kind="bar", ci=None, aspect=1,legend=False,
#                 palette=sns.color_palette("Set2")[:1]+sns.color_palette("Set2")[2:3]+sns.color_palette("Set2")[4:6],hue_order=['PASTE2', 'PASTE'])
g = sns.catplot(x="Pair", y="ARI", hue='Kind', col="Sample",data=res_df,kind="bar", ci=None, aspect=1,legend=False,
                palette=sns.color_palette(),hue_order=['PASTE2', 'PASTE', 'Pamona', 'Tangram'])
g.set_axis_labels("Pair", "ARI").set(ylim=(0, 0.7))
g.add_legend(title='',labels=['PASTE2', 'PASTE', 'Pamona', 'Tangram'],bbox_to_anchor=(0.5, -0.025),ncol=4)#(1.0, 0.5))(0.52,-0.1)

plt.savefig("test.png")