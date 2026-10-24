import pytest
from pathlib import Path

from thermoacoustic.config import load_config, check_same_setup

CONFIG_DIR = Path(__file__).resolve().parent.parent / "configs"


def test_load_real_configs():
    """Both shipped YAML files load and contain all required entries."""
    for name in ("ch4_rijke.yaml", "h2_rijke.yaml"):
        cfg = load_config(CONFIG_DIR / name)
        assert cfg["fuel"]["name"] in ("CH4", "H2")
        assert cfg["geometry"]["tube_length"] > 0


def test_missing_key_raises(tmp_path):
    """A config with a missing entry must fail loudly, not later and obscurely."""
    bad = tmp_path / "bad.yaml"
    bad.write_text("fuel:\n  name: CH4\n")
    with pytest.raises(KeyError):
        load_config(bad)


def test_check_same_setup_accepts_shipped_configs():
    cfg_ch4 = load_config(CONFIG_DIR / "ch4_rijke.yaml")
    cfg_h2 = load_config(CONFIG_DIR / "h2_rijke.yaml")
    check_same_setup(cfg_ch4, cfg_h2)   # must not raise


def test_check_same_setup_rejects_different_geometry():
    cfg_ch4 = load_config(CONFIG_DIR / "ch4_rijke.yaml")
    cfg_h2 = load_config(CONFIG_DIR / "h2_rijke.yaml")
    cfg_h2["geometry"]["tube_length"] = 2.0
    with pytest.raises(ValueError):
        check_same_setup(cfg_ch4, cfg_h2)
        
def test_check_same_setup_rejects_different_equivalence_ratio():
       """A fuel comparison must be made at the same equivalence ratio."""
       a = load_config(CONFIG_DIR / "ch4_rijke.yaml")
       b = load_config(CONFIG_DIR / "h2_rijke.yaml")
       b["fuel"]["equivalence_ratio"] = 0.9
       with pytest.raises(ValueError):
           check_same_setup(a, b)