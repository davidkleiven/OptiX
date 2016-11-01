# *-* coding: utf-8 *-*
#from __future__ import unicode_literals
ONE_COLUMN = False
LATEX_PREAMBLE = u"\usepackage{gensymb}"
LATEX_PREAMBLE += u"\IfFileExists{siunitx.sty}{\usepackage{siunitx}}{}"

# Provide an SI command and macros if the siunitx-package is not found
LATEX_PREAMBLE += u"\providecommand{\SI}[2]{{#1} \; {#2}}"
LATEX_PREAMBLE += u"\providecommand{\micro}{\ensuremath{\mathrm{\mu}}}"
LATEX_PREAMBLE += u"\providecommand{\milli}{\ensuremath{\mathrm{m}}}"
LATEX_PREAMBLE += u"\providecommand{\meter}{\ensuremath{\mathrm{m}}}"

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
