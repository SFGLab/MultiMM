from enum import Enum


class InitialStructureType(str, Enum):
    RW = "rw"
    CONFINED_RW = "confined_rw"
    KNOT = "knot"
    SELF_AVOIDING_RW = "self_avoiding_rw"
    CIRCLE = "circle"
    HELIX = "helix"
    SPIRAL = "spiral"
    SPHERE = "sphere"
    HILBERT = "hilbert"
