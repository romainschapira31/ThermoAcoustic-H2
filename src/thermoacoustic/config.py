"""Loading and validation of YAML configuration files."""

from pathlib import Path
import yaml

REQUIRED_KEYS = {
    "fuel": ["name", "mechanism", "oxidizer", "equivalence_ratio", "T_fresh", "pressure"],
    "geometry": ["tube_length", "flame_position", "speed_of_sound"],
    "flame_model": ["characteristic_distance", "coupling_factor"],
    "sweep": ["n_min", "n_max", "n_points"],
}


def load_config(path):
    """Read a YAML configuration file and check that all required entries exist."""
    with open(Path(path)) as f:
        config = yaml.safe_load(f)
    for section, keys in REQUIRED_KEYS.items():
        if section not in config:
            raise KeyError(f"Missing section '{section}' in {path}")
        for key in keys:
            if key not in config[section]:
                raise KeyError(f"Missing key '{section}.{key}' in {path}")
    return config


def check_same_setup(cfg_a, cfg_b):
    """Raise if two configs differ in anything other than the fuel itself,
    so that a fuel comparison is always like-for-like."""
    for section in ("geometry", "flame_model", "sweep"):
        if cfg_a[section] != cfg_b[section]:
            raise ValueError(f"Configs differ in '{section}': comparison would be invalid")
    for key in ("T_fresh", "pressure"):
        if cfg_a["fuel"][key] != cfg_b["fuel"][key]:
            raise ValueError(f"Configs differ in fuel.{key}: comparison would be invalid")