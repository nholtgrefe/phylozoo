"""
Inference module.

This module provides classes and functions for phylogenetic network inference.
"""

from .inference import NetworkInferrer, infer_network_from_msa

__all__ = [
    "NetworkInferrer",
    "infer_network_from_msa",
]

