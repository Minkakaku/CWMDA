"""Genome pipeline entrypoint (example skeleton)."""
from __future__ import annotations

from pathlib import Path

from omics.general.scripts.config_utils import load_config, resolve_path
from omics.general.scripts.file_utils import ensure_dir, validate_file
from omics.general.scripts.log_utils import setup_logger
from omics.general.scripts.qc_common import write_qc_report


def main() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    common_config = repo_root / "omics" / "general" / "config" / "common_config.yaml"
    genome_config = repo_root / "omics" / "genomics" / "config" / "genome_config.yaml"
    config = load_config(str(common_config), str(genome_config))
    output_dir = Path(config["genome"]["output_dir"])
    ensure_dir(output_dir)

    logger = setup_logger("genome", config["logging"]["log_dir"], config["logging"]["level"])
    logger.info("开始基因组流程示例")

    input_bam = resolve_path(config["genome"]["input_bam"])
    reference_fasta = resolve_path(config["genome"]["reference_fasta"])
    validate_file(input_bam, "输入 BAM")
    validate_file(reference_fasta, "参考基因组")

    qc_stats = {
        "sample": config["genome"]["sample_id"],
        "qc_tool": config["qc"]["tool"],
        "input_bam": str(input_bam),
        "reference": str(reference_fasta),
    }

    write_qc_report(
        qc_stats,
        output_dir,
        config["qc"]["tsv_name"],
        config["qc"]["html_name"],
    )

    logger.info("质控报告已生成: %s", output_dir)


if __name__ == "__main__":
    main()
