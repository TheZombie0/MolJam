import sys

import matplotlib

plt = None
sns = None
Rectangle = None
mpatches = None
GridSpec = None
PdfPages = None
inset_axes = None


def ensure_plotting_imports(backend="Agg"):
    global plt, sns, Rectangle, mpatches, GridSpec, PdfPages, inset_axes

    if plt is not None and sns is not None:
        return

    if backend and "matplotlib.pyplot" not in sys.modules:
        try:
            matplotlib.use(backend)
        except Exception:
            pass

    import matplotlib.pyplot as _plt
    from matplotlib.patches import Rectangle as _Rectangle
    import matplotlib.patches as _mpatches
    from matplotlib.gridspec import GridSpec as _GridSpec
    from matplotlib.backends.backend_pdf import PdfPages as _PdfPages
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes as _inset_axes
    import seaborn as _sns

    plt = _plt
    Rectangle = _Rectangle
    mpatches = _mpatches
    GridSpec = _GridSpec
    PdfPages = _PdfPages
    inset_axes = _inset_axes
    sns = _sns

