import matplotlib.pyplot as plt
from pandas.plotting import scatter_matrix
import numpy as np
import pandas as pd 

data= pd.read_csv('prob.csv')

print(data.describe())
print(data)
#data.plot(subplots=True, layout=(2, 2))

scatter_matrix(data, alpha=0.2, figsize=(6, 6), diagonal="kde")

plt.show()