from ._hiclassy import (
    __version__,
    HiClass,
    CosmoSevereError,
    CosmoComputationError
)

# Alias kept so that code written for CLASS's `classy.Class` keeps working
# with `from hiclassy import Class`.
Class = HiClass

__all__ = ["__version__", "HiClass", "Class", "CosmoSevereError", "CosmoComputationError"]
