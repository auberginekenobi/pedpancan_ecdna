import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

# Some helper functions
def savefig(fig, basename,dpi=300):
    pngName = basename + ".png"
    svgName = basename + ".svg"
    fig.savefig(pngName,format='png', dpi=dpi, bbox_inches='tight')
    fig.savefig(svgName,format='svg')
def set_plot_defaults(linewidth=0.75,fontsize=7):
    sns.set_style("white")
    mpl.rcParams.update({
        'svg.fonttype':'none',
        'pdf.use14corefonts':True,
        'figure.figsize':(7,2),
        'font.size': fontsize,'axes.labelsize': fontsize,'axes.titlesize': fontsize,'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,'legend.fontsize': fontsize,
        'font.family':'sans-serif',
        'font.sans-serif':'Arial',
        'axes.linewidth':linewidth,'ytick.major.width':linewidth,'xtick.major.width':linewidth,
        'patch.force_edgecolor':False,
        'ytick.left':True,
    })
    return