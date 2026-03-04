try:
    from pyecharts import options as opts
    from pyecharts.charts import Sankey as EchartsSankey
    from pyecharts.render import make_snapshot

    try:
        from snapshot_selenium import snapshot

        SNAPSHOT_AVAILABLE = True
    except ImportError:
        SNAPSHOT_AVAILABLE = False

    PYECHARTS_AVAILABLE = True
except ImportError:
    opts = None
    EchartsSankey = None
    make_snapshot = None
    snapshot = None
    SNAPSHOT_AVAILABLE = False
    PYECHARTS_AVAILABLE = False

