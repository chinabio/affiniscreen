"""Launch-time guard against production-timestep regressions (v2.5.62)."""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
from __future__ import annotations
MAX_PROD_DT_PS = 0.001
class UnsafeTimestepError(RuntimeError): pass
def _check(label,value):
    try: v=float(value)
    except (TypeError,ValueError): return
    if v>MAX_PROD_DT_PS+1e-12:
        raise UnsafeTimestepError(
            f"Production timestep {label}={v} ps exceeds GPU-stable ceiling "
            f"{MAX_PROD_DT_PS} ps. ABFE restraint mid-band blows up at 2 fs. "
            f"Set dt=0.001 and double the matching nstlim.")
def assert_prod_dt_safe(cfg=None):
    fields=("dt_ps","prod_dt_ps")
    if cfg is not None:
        for f in fields:
            if hasattr(cfg,f): _check(f,getattr(cfg,f))
        md=getattr(cfg,"md",None)
        if md is not None:
            for f in fields:
                if hasattr(md,f): _check("md."+f,getattr(md,f))
        return
    import dataclasses
    from . import config as _c
    for name in dir(_c):
        o=getattr(_c,name)
        if dataclasses.is_dataclass(o):
            for fld in dataclasses.fields(o):
                if fld.name in fields and fld.default is not dataclasses.MISSING:
                    _check(f"{name}.{fld.name}",fld.default)
if __name__=="__main__":
    assert_prod_dt_safe(); print(f"[dt-guard] OK: all production dt defaults <= {MAX_PROD_DT_PS} ps")
