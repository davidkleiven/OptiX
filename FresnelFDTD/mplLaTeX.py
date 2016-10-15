ONE_COLUMN = False
LATEX_PREAMBLE = "\usepackage{gensymb}"
LATEX_PREAMBLE += "\IfFileExists{siunitx.sty}{\usepackage{siunitx}}{ \providecommand{\SI}[2]{{#1} \; {#2}} }"
LATEX_PREAMBLE += "\providecommand{\micro}{\ensuremath{\mathrm{\mu}}}"
LATEX_PREAMBLE += "\providecommand{\milli}{\ensuremath{\mathrm{m}}}"
LATEX_PREAMBLE += "\providecommand{\meter}{\ensuremath{\mathrm{m}}}"

if ( ONE_COLUMN ):
    figsize = [6.0, 4.5]
else:
    figsize = [3.5, 2.25]
params = {'backend': 'ps',
       'text.latex.preamble': [LATEX_PREAMBLE],
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
