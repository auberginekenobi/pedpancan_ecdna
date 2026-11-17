import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import yaml

COLORS_YML = Path(__file__).resolve().parents[1] / 'config' / 'colors.yml'

def get_palette(name, shade=None, path=COLORS_YML):
    '''
    Dict (label -> hex) for palette `name` in config/colors.yml, usable as a seaborn palette=.
    shade: for palettes with dark/light entries (tumor_type, tumor_subtype), which one to return;
           None returns the {'dark': ..., 'light': ...} dicts as stored (missing shades are None).
    '''
    with open(path) as f:
        palettes = yaml.safe_load(f)
    if name not in palettes:
        raise KeyError(f"no palette '{name}' in {path}")
    pal = palettes[name]
    if shade is not None:
        pal = {k: v[shade] for k, v in pal.items()}
    return pal

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