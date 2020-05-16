import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go

movments = pd.read_csv(r"/Users/daniellemiller/Downloads/internal_movements_count_verbosed.csv")



isr_cities = movments['country_self'].unique()
world = list(set([c for c  in movments['country_first_change'] if c not in isr_cities]))

# the first layer
layer_1 = {}
for i, w in enumerate(world):
    layer_1[w] = i

layer_2 = {}
for i, w in enumerate(isr_cities):
    layer_2[w] = i + len(world)

layer_3 = {}
for i, w in enumerate(isr_cities):
    layer_3[w] = i + len(world) + len(isr_cities)

movments['source'] = movments['country_first_change'].apply(lambda x: layer_1[x] if x in world else layer_2[x])
movments['target'] = movments.apply(lambda row: layer_3[row['country_self']] if
row['country_first_change'] in isr_cities else layer_2[row['country_self']], axis=1)

movments.rename(columns={'size':'value'}, inplace=True)
pal = sns.color_palette('Set2', len(world) + len(isr_cities)).as_hex()

all_names = list(isr_cities) + world
label_2_color = {}
for i, name in enumerate(all_names):
    label_2_color[name] = pal[i]

label_list = list(layer_1.keys()) + list(layer_2.keys()) + list(layer_3.keys())


fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = list(layer_1.keys()) + list(layer_2.keys()) + list(layer_3.keys()),
      color = [label_2_color[n] for n in label_list]
    ),
    link = dict(
      source = movments['source'], # indices correspond to labels, eg A1, A2, A2, B1, ...
      target = movments['target'],
      value = movments['value'],
      #color = [new_pal[i] for i in g['source']]
  ))])

#fig.write_image(r"/Users/daniellemiller/Google Drive/covid19/paper/COVID19_phylodynamics/figures/transmission_patterns.png", scale=10, width=1000)
fig.show()