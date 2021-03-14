import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


df = pd.read_csv("samples.csv")

if(True):
    walker_id = df.walker.drop_duplicates().array


    fig, ax = plt.subplots(3,1,constrained_layout=True)    
    for walker in walker_id:
        df_aux = df.loc[df["walker"]==walker]
        size = len(df_aux["th13"].index)
        x = np.linspace(0,size,size)
        ax[0].plot(x, df_aux["th13"], alpha=0.3, color = 'black')
        ax[1].plot(x, df_aux["dm31"], alpha=0.3, color = 'black')
        ax[2].plot(x, df_aux["real"], alpha=0.3, color = 'black')
        #ax[3].plot(x, df_aux["imag"], alpha=0.3, color = 'black')
        #ax[4].plot(x, df_aux["pull1"], alpha=0.3, color = 'black')
        #ax[5].plot(x, df_aux["pull2"], alpha=0.3, color = 'black')
        #ax[6].plot(x, df_aux["pull3"], alpha=0.3, color = 'black')
        #ax[7].plot(x, df_aux["pull4"], alpha=0.3, color = 'black')
    ax[0].set_xlabel("th13")
    ax[1].set_xlabel("dm31")
    ax[2].set_xlabel("real")
    #ax[3].set_xlabel("imag")
    plt.savefig("SBL-ReSet-convergence-500walks.png", dpi=300)
    plt.show()


if(True):
    df = df[['th13','dm31','real']]
    sns.set_theme(style="dark")
    g = sns.PairGrid(df, diag_sharey=False)
    g.map_upper(sns.scatterplot, s=5, alpha = 0.1, color = 'red')
    g.map_lower(sns.kdeplot, fill=True, levels=[0.01, 0.1 ,0.64, 0.95, 0.99, 1], cmap = 'Reds')
    g.map_diag(sns.histplot, kde=True, color = 'red')
    plt.savefig("SBL-ReSet-500walks.png", dpi=300)
    plt.show()
