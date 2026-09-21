"""General project visualisations."""

from __future__ import annotations

from pathlib import Path


def save_experimental_design(output_path: str | Path) -> None:
    """Save a compact schematic of the analysis design."""
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(8, 5))
    axis.axis("off")
    text = (
        "DMSO: perturbation + non-targeting\n"
        "Dasatinib: perturbation + non-targeting\n\n"
        "RNA + ATAC\n↓\nBackground correction\n↓\n"
        "Cross-modal scoring\n↓\nCandidate ranking"
    )
    axis.text(0.5, 0.5, text, ha="center", va="center", fontsize=13)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)

