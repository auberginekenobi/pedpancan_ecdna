import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

# Some helper functions
def savefig(fig, basename,dpi=300):
    pngName = basename + ".png"
    svgName = basename + ".svg"
    fig.savefig(pngName,format='png', dpi=dpi, bbox_inches='tight')
    fig.savefig(svgName,format='svg')
def set_plot_defaults(linewidth=0.75,fontsize=7,figsize=(7,2)):
    sns.set_style("white")
    mpl.rcParams.update({
        'svg.fonttype':'none',
        'pdf.use14corefonts':True,
        'figure.figsize':figsize,
        'font.size': fontsize,'axes.labelsize': fontsize,'axes.titlesize': fontsize,'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,'legend.fontsize': fontsize,
        'font.family':'sans-serif',
        'font.sans-serif':'Arial',
        'axes.linewidth':linewidth,'ytick.major.width':linewidth,'ytick.minor.width':linewidth,'xtick.major.width':linewidth,
        'xtick.minor.width':linewidth,
        'xtick.major.size':2,
        'patch.force_edgecolor':False,
        'ytick.left':True,
    })
    return