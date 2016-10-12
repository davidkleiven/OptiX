ONE_COLUMN = False
if ( ONE_COLUMN ):
    figsize = [6.0, 4.5]
else:
    figsize = [3.5, 2.25]
params = {'backend': 'ps',
       'text.latex.preamble': ['\usepackage{gensymb}\IfFileExists{siunitx.sty}{\usepackage{siunitx}}{\providecommand{\SI}[2]{#1\mathrm{#2}}}'],
       'axes.labelsize': 8, # fontsize for x and y labels (was 10)
       'axes.titlesize': 8,
       'font.size': 8, # was 10
       'legend.fontsize': 8, # was 10
       'xtick.labelsize': 8,
       'ytick.labelsize': 8,
       'text.usetex': True,
       'font.family': 'serif',
       'axes.linewidth':0.1,
       'figure.figsize': figsize
}
