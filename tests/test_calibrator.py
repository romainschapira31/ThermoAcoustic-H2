from pathlib import Path

from thermoacoustic.config import load_config
from thermoacoustic.flame import FlameCalibrator

CONFIG_DIR = Path(__file__).resolve().parent.parent / "configs"


def test_calibrator_from_config_uses_the_configuration():
    cfg = load_config(CONFIG_DIR / "h2_rijke.yaml")
    calibrator = FlameCalibrator.from_config(cfg)
    assert calibrator.fuel == cfg["fuel"]["name"]
    assert calibrator.mechanism == cfg["fuel"]["mechanism"]
    assert calibrator.T_in == cfg["fuel"]["T_fresh"]
    assert calibrator.P_in == cfg["fuel"]["pressure"]
    assert calibrator.characteristic_distance == cfg["flame_model"]["characteristic_distance"]