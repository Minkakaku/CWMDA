from pathlib import Path

from omics.general.scripts.config_utils import load_config
from omics.general.scripts.qc_common import write_qc_report


def test_load_and_qc(tmp_path):
    repo_root = Path(__file__).resolve().parents[3]
    common_config = repo_root / "omics" / "general" / "config" / "common_config.yaml"
    genome_config = repo_root / "omics" / "genomics" / "config" / "genome_config.yaml"
    config = load_config(str(common_config), str(genome_config))
    assert "logging" in config

    output_dir = tmp_path / "qc"
    write_qc_report({"sample": "demo"}, output_dir, "stats.tsv", "report.html")
    assert (output_dir / "stats.tsv").exists()
    assert (output_dir / "report.html").exists()
