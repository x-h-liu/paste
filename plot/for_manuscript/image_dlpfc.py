import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



# """
# Top down
# """
# plt.rcParams.update({'font.size': 22})


# res_df = pd.DataFrame(columns=['Sample','Pair','Method', 'ARI'])



# res_df.loc[len(res_df)] = [1, 'AB', 'Expression only', 0.5982726987358716]
# res_df.loc[len(res_df)] = [1, 'BC', 'Expression only', 0.06780848340640959]
# res_df.loc[len(res_df)] = [1, 'CD', 'Expression only', 0.5650396516478512]
# res_df.loc[len(res_df)] = [2, 'AB', 'Expression only', 0.6479499641985342]
# res_df.loc[len(res_df)] = [2, 'BC', 'Expression only', 0.6637812938623632]
# res_df.loc[len(res_df)] = [2, 'CD', 'Expression only', 0.36385628822237787]
# res_df.loc[len(res_df)] = [3, 'AB', 'Expression only', 0.5587703813466455]
# res_df.loc[len(res_df)] = [3, 'BC', 'Expression only', 0.4724182915603778]
# res_df.loc[len(res_df)] = [3, 'CD', 'Expression only', 0.39816624188317135]

# res_df.loc[len(res_df)] = [1, 'AB', 'Expression + Image', 0.596964448718805]
# res_df.loc[len(res_df)] = [1, 'BC', 'Expression + Image', 0.2618641600945121]
# res_df.loc[len(res_df)] = [1, 'CD', 'Expression + Image', 0.48693969912449453]
# res_df.loc[len(res_df)] = [2, 'AB', 'Expression + Image', 0.6423396191946658]
# res_df.loc[len(res_df)] = [2, 'BC', 'Expression + Image', 0.6884785330750198]
# res_df.loc[len(res_df)] = [2, 'CD', 'Expression + Image', 0.3690141159270531]
# res_df.loc[len(res_df)] = [3, 'AB', 'Expression + Image', 0.4959730461101442]
# res_df.loc[len(res_df)] = [3, 'BC', 'Expression + Image', 0.483613865942996]
# res_df.loc[len(res_df)] = [3, 'CD', 'Expression + Image', 0.4080224400539087]



# g = sns.catplot(x="Pair", y="ARI", hue='Method', col="Sample",data=res_df,kind="bar", ci=None, aspect=1,legend=True,
#                 palette=sns.color_palette())
# g.set_axis_labels("Pair", "ARI").set(ylim=(0, 1))


# plt.savefig('test.png')



# """
# Left right
# """
# plt.rcParams.update({'font.size': 22})


# res_df = pd.DataFrame(columns=['Sample','Pair','Method', 'ARI'])



# res_df.loc[len(res_df)] = [1, 'AB', 'Expression only', 0.5763368119922595]
# res_df.loc[len(res_df)] = [1, 'BC', 'Expression only', 0.4312785685555477]
# res_df.loc[len(res_df)] = [1, 'CD', 'Expression only', 0.5379887967330693]
# res_df.loc[len(res_df)] = [2, 'AB', 'Expression only', 0.6914330208160833]
# res_df.loc[len(res_df)] = [2, 'BC', 'Expression only', 0.6422029895361098]
# res_df.loc[len(res_df)] = [2, 'CD', 'Expression only', 0.687653609937276]
# res_df.loc[len(res_df)] = [3, 'AB', 'Expression only', 0.5189486914217226]
# res_df.loc[len(res_df)] = [3, 'BC', 'Expression only', 0.4623086489723203]
# res_df.loc[len(res_df)] = [3, 'CD', 'Expression only', 0.3426678283641906]

# res_df.loc[len(res_df)] = [1, 'AB', 'Expression + Image', 0.6146856922928862]
# res_df.loc[len(res_df)] = [1, 'BC', 'Expression + Image', 0.2904057264861441]
# res_df.loc[len(res_df)] = [1, 'CD', 'Expression + Image', 0.5749537011524328]
# res_df.loc[len(res_df)] = [2, 'AB', 'Expression + Image', 0.6418699373588566]
# res_df.loc[len(res_df)] = [2, 'BC', 'Expression + Image', 0.5814365098235177]
# res_df.loc[len(res_df)] = [2, 'CD', 'Expression + Image', 0.686696877729989]
# res_df.loc[len(res_df)] = [3, 'AB', 'Expression + Image', 0.49229579222215664]
# res_df.loc[len(res_df)] = [3, 'BC', 'Expression + Image', 0.4584654001109203]
# res_df.loc[len(res_df)] = [3, 'CD', 'Expression + Image', 0.4632599209342587]



# g = sns.catplot(x="Pair", y="ARI", hue='Method', col="Sample",data=res_df,kind="bar", ci=None, aspect=1,legend=True,
#                 palette=sns.color_palette())
# g.set_axis_labels("Pair", "ARI").set(ylim=(0, 1))


# plt.savefig('test.png')




"""
Sample 3
"""
plt.rcParams.update({'font.size': 22})


res_df = pd.DataFrame(columns=['Subslice','Pair','Method', 'ARI'])



# res_df.loc[len(res_df)] = ['Horizontal', 'AB', 'Expression only', 0.5587703813466455]
# res_df.loc[len(res_df)] = ['Horizontal', 'BC', 'Expression only', 0.4724182915603778]
# res_df.loc[len(res_df)] = ['Horizontal', 'CD', 'Expression only', 0.39816624188317135]
# res_df.loc[len(res_df)] = ['Horizontal', 'AB', 'Expression + Image', 0.4959730461101442]
# res_df.loc[len(res_df)] = ['Horizontal', 'BC', 'Expression + Image', 0.483613865942996]
# res_df.loc[len(res_df)] = ['Horizontal', 'CD', 'Expression + Image', 0.4080224400539087]


res_df.loc[len(res_df)] = ['Vertical', 'AB', 'Expression only', 0.5189486914217226]
res_df.loc[len(res_df)] = ['Vertical', 'BC', 'Expression only', 0.4623086489723203]
res_df.loc[len(res_df)] = ['Vertical', 'CD', 'Expression only', 0.3426678283641906]
res_df.loc[len(res_df)] = ['Vertical', 'AB', 'Expression + Image', 0.49229579222215664]
res_df.loc[len(res_df)] = ['Vertical', 'BC', 'Expression + Image', 0.4584654001109203]
res_df.loc[len(res_df)] = ['Vertical', 'CD', 'Expression + Image', 0.4632599209342587]



g = sns.catplot(x="Pair", y="ARI", hue='Method', col="Subslice",data=res_df,kind="bar", ci=None, aspect=1,legend=False,
                palette=sns.color_palette())
g.set_axis_labels("Pair", "Label Transfer ARI").set(ylim=(0, 0.8))


plt.savefig('test.png')



