import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math


values = np.loadtxt('bottom_values_NONE.txt')

values = values[(values != 0) & (values != 1)]
log_values = np.log2(values)

plt.figure(figsize=(10, 6))
plt.hist(log_values, bins=100, edgecolor='black')
plt.title('Log-Scaled Distribution of Bottom Values for NONE')
plt.xlabel('Log2(Bottom Value)')
plt.ylabel('Frequency')
plt.savefig('bottom_values_NONE.png')
plt.show()


values = np.loadtxt('bottom_values_SUM1.txt')

values = values[(values != 0) & (values != 1)]
log_values = np.log2(values)


plt.figure(figsize=(10, 6))
plt.hist(log_values, bins=100, edgecolor='black')
plt.title('Log-Scaled Distribution of Bottom Values for SUM1')
plt.xlabel('Log2(Bottom Value)')
plt.ylabel('Frequency')
plt.savefig('bottom_values_SUM1.png')
plt.show()


values = np.loadtxt('bottom_values_EXP.txt')

values = values[(values != 0) & (values != 1)]
log_values = np.log2(values)


plt.figure(figsize=(10, 6))
plt.hist(log_values, bins=100, edgecolor='black')
plt.title('Log-Scaled Distribution of Bottom Values for EXP')
plt.xlabel('Log2(Bottom Value)')
plt.ylabel('Frequency')
plt.savefig('bottom_values_EXP.png')
plt.show()


values = np.loadtxt('bottom_values_MAX.txt')

values = values[(values != 0) & (values != 1)]
log_values = np.log2(values)


plt.figure(figsize=(10, 6))
plt.hist(log_values, bins=100, edgecolor='black')
plt.title('Log-Scaled Distribution of Bottom Values for MAX')
plt.xlabel('Log2(Bottom Value)')
plt.ylabel('Frequency')
plt.savefig('bottom_values_MAX.png')
plt.show()






import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math


df = pd.read_csv('bottom_values_by_layer_NONE.csv')

df = df[(df['bottom_value'] != 0) & (df['bottom_value'] != 1)]

# 计算log2(bottom_value)
df['log2_bottom_value'] = np.log2(df['bottom_value'])

layers = sorted(df['layer'].unique())
print(f"Layers found: {layers}")

sns.set(style="whitegrid")

num_layers = len(layers)
cols = 2  
rows = math.ceil(num_layers / cols)

fig, axes = plt.subplots(rows, cols, figsize=(12, rows * 4))
axes = axes.flatten()

for idx, layer in enumerate(layers):
    df_layer = df[df['layer'] == layer]
    log_values = df_layer['log2_bottom_value']

    ax = axes[idx]
    sns.histplot(log_values, bins=100, edgecolor='black', ax=ax)
    ax.set_title(f'Layer {layer}')
    ax.set_xlabel('Log2(Bottom Value)')
    ax.set_ylabel('Frequency')
    # ax.set_xlim(df['log2_bottom_value'].min(), df['log2_bottom_value'].max())


for idx in range(len(layers), len(axes)):
    fig.delaxes(axes[idx])

fig.suptitle('Log-Scaled Value Distribution By Layer for MAX', fontsize=16)

plt.tight_layout(rect=[0, 0, 1, 0.96]) 
plt.savefig('bottom_values_by_layer_log2_none.png')
plt.show()

df = pd.read_csv('bottom_values_by_layer_MAX.csv')

df = df[(df['bottom_value'] != 0) & (df['bottom_value'] != 1)]

df['log2_bottom_value'] = np.log2(df['bottom_value'])

layers = sorted(df['layer'].unique())
print(f"Layers found: {layers}")

sns.set(style="whitegrid")

num_layers = len(layers)
cols = 2  
rows = math.ceil(num_layers / cols)

fig, axes = plt.subplots(rows, cols, figsize=(12, rows * 4))
axes = axes.flatten()

for idx, layer in enumerate(layers):
    df_layer = df[df['layer'] == layer]
    log_values = df_layer['log2_bottom_value']

    ax = axes[idx]
    sns.histplot(log_values, bins=100, edgecolor='black', ax=ax)
    ax.set_title(f'Layer {layer}')
    ax.set_xlabel('Log2(Bottom Value)')
    ax.set_ylabel('Frequency')
    # ax.set_xlim(df['log2_bottom_value'].min(), df['log2_bottom_value'].max())

for idx in range(len(layers), len(axes)):
    fig.delaxes(axes[idx])

fig.suptitle('Log-Scaled Value Distribution By Layer for MAX', fontsize=16)

plt.tight_layout(rect=[0, 0, 1, 0.96]) 
plt.savefig('bottom_values_by_layer_log2_max.png')
plt.show()



df = pd.read_csv('bottom_values_by_layer_SUM1.csv')

df = df[(df['bottom_value'] != 0) & (df['bottom_value'] != 1)]

df['log2_bottom_value'] = np.log2(df['bottom_value'])

layers = sorted(df['layer'].unique())
print(f"Layers found: {layers}")

sns.set(style="whitegrid")


num_layers = len(layers)
cols = 2  
rows = math.ceil(num_layers / cols)

fig, axes = plt.subplots(rows, cols, figsize=(12, rows * 4))
axes = axes.flatten()

for idx, layer in enumerate(layers):
    df_layer = df[df['layer'] == layer]
    log_values = df_layer['log2_bottom_value']

    ax = axes[idx]
    sns.histplot(log_values, bins=100, edgecolor='black', ax=ax)
    ax.set_title(f'Layer {layer}')
    ax.set_xlabel('Log2(Bottom Value)')
    ax.set_ylabel('Frequency')
    # ax.set_xlim(df['log2_bottom_value'].min(), df['log2_bottom_value'].max())

for idx in range(len(layers), len(axes)):
    fig.delaxes(axes[idx])

fig.suptitle('Log-Scaled Value Distribution By Layer for SUM1', fontsize=16)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('bottom_values_by_layer_log2_sum1.png')
plt.show()





df = pd.read_csv('bottom_values_by_layer_EXP.csv')

df = df[(df['bottom_value'] != 0) & (df['bottom_value'] != 1)]

df['log2_bottom_value'] = np.log2(df['bottom_value'])

layers = sorted(df['layer'].unique())
print(f"Layers found: {layers}")

sns.set(style="whitegrid")

num_layers = len(layers)
cols = 2  
rows = math.ceil(num_layers / cols)

fig, axes = plt.subplots(rows, cols, figsize=(12, rows * 4))
axes = axes.flatten()

for idx, layer in enumerate(layers):
    df_layer = df[df['layer'] == layer]
    log_values = df_layer['log2_bottom_value']

    ax = axes[idx]
    sns.histplot(log_values, bins=100, edgecolor='black', ax=ax)
    ax.set_title(f'Layer {layer}')
    ax.set_xlabel('Log2(Bottom Value)')
    ax.set_ylabel('Frequency')
    # ax.set_xlim(df['log2_bottom_value'].min(), df['log2_bottom_value'].max())

for idx in range(len(layers), len(axes)):
    fig.delaxes(axes[idx])

fig.suptitle('Log-Scaled Value Distribution By Layer for EXP', fontsize=16)

plt.tight_layout(rect=[0, 0, 1, 0.96]) 
plt.savefig('bottom_values_by_layer_log2_exp.png')
plt.show()


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns




values = np.loadtxt('addition_values_NONE.txt')

values = values[(np.isfinite(values))&(values != 0)]

values = values[values != 1]
log_values = np.log2(values)

plt.figure(figsize=(10, 6))
sns.histplot(log_values, bins=100, kde=True)
plt.title('Log2 Distribution of Multiplication Values for EXP')
plt.xlabel('Log2(Product Value)')
plt.ylabel('Frequency')
plt.savefig('log2_addition_values_distribution_NONE.png')
plt.show()



values = np.loadtxt('addition_values_EXP.txt')

values = values[(np.isfinite(values))&(values != 0)]

values = values[values != 1]
log_values = np.log2(values)

plt.figure(figsize=(10, 6))
sns.histplot(log_values, bins=100, kde=True)
plt.title('Log2 Distribution of Multiplication Values for EXP')
plt.xlabel('Log2(Product Value)')
plt.ylabel('Frequency')
plt.savefig('log2_addition_values_distribution_EXP.png')
plt.show()



values = np.loadtxt('addition_values_SUM1.txt')

values = values[(np.isfinite(values))&(values != 0)]

values = values[values != 1]
log_values = np.log2(values)

plt.figure(figsize=(10, 6))
sns.histplot(log_values, bins=100, kde=True)
plt.title('Log2 Distribution of Multiplication Values for SUM1')
plt.xlabel('Log2(Product Value)')
plt.ylabel('Frequency')
plt.savefig('log2_addition_values_distribution_SUM1.png')
plt.show()


values = np.loadtxt('addition_values_MAX.txt')

values = values[(np.isfinite(values))&(values != 0)]

values = values[values != 1]
log_values = np.log2(values)

plt.figure(figsize=(10, 6))
sns.histplot(log_values, bins=100, kde=True)
plt.title('Log2 Distribution of Multiplication Values for MAX')
plt.xlabel('Log2(Product Value)')
plt.ylabel('Frequency')
plt.savefig('log2_addition_values_distribution_MAX.png')
plt.show()