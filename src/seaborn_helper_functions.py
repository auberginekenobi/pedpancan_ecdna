import seaborn as sns
import matplotlib.pyplot as plt

# Some helper functions
def savefig(fig, basename):
    pngName = basename + ".png"
    svgName = basename + ".svg"
    fig.savefig(pngName,format='png')
    fig.savefig(svgName,format='svg')
def set_plot_defaults(linewidth=0.75,fontsize=10):
    sns.set(rc={'svg.fonttype':'none',
                'pdf.use14corefonts':True,
                'figure.figsize':(7,2),
                'font.size': fontsize,'axes.labelsize': fontsize,'axes.titlesize': fontsize,'xtick.labelsize': fontsize,
                'ytick.labelsize': fontsize,'legend.fontsize': fontsize,
                'font.family':'sans-serif',
                'font.sans-serif':'Arial',
                'axes.linewidth':linewidth,
                'ytick.major.width':linewidth,
                })
    plt.tight_layout()
    sns.set_style("white")
    sns.despine()
    return