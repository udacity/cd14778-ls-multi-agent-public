"""
vina_helpers_bio.py — Tiny wrappers for BLAST+, Meeko, AutoDock Vina, and SDF generation.

All usage is documented inline in each function docstring. Keep this file in your project
and import what you need. No README necessary.

Requires command-line tools in PATH:
- BLAST+: blastp, blastdbcmd      https://blast.ncbi.nlm.nih.gov/
- Meeko: mk_prepare_ligand.py, mk_prepare_receptor.py   https://github.com/forlilab/Meeko
- AutoDock Vina: vina             http://vina.scripps.edu/
- (Optional) Open Babel: obabel   https://openbabel.org/
- (Optional) RDKit (Python)       https://www.rdkit.org/

Tip: set BLASTDB env to your database dir, e.g. `export BLASTDB=/home/ubuntu/blastdb`
"""

from __future__ import annotations
import csv
import os
import re
import shutil
import subprocess
import tempfile
import textwrap
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

__all__ = [
    "BLAST_OUTFMT",
    "blastp_from_seq",
    "meeko_ligand",
    "meeko_receptor",
    "smiles_to_sdf_3d",
    "vina_score_only",
    "vina_dock",
]


def _run(cmd: Sequence[str]) -> str:
    """Run a subprocess and return stdout (raises on non-zero).
    """
    return subprocess.run(list(cmd), check=True, capture_output=True, text=True).stdout


# ----------------------------- BLAST+ ----------------------------------------

BLAST_OUTFMT = "6 sacc pident length evalue bitscore sscinames"


def blastp_from_seq(
    seq: str,
    db: str = "pdbaa",
    *,
    max_targets: int = 10,
    threads: int = 2,
    blastdb: str | None = None,
) -> List[Dict[str, str]]:
    """Run BLASTP against a local DB and return tabular hits as dicts.

    Args:
        seq: Protein sequence (FASTA body).
        db: BLAST protein DB name (e.g., "pdbaa", "swissprot").
        max_targets: Max target sequences.
        threads: Number of CPU threads.
        blastdb: Optional BLASTDB directory to export for this call.

    Returns:
        List of dictionaries with keys: sacc, pident, length, evalue, bitscore, sscinames.

    Notes:
        - Ensure BLAST+ is installed and DBs are in $BLASTDB; taxonomy fields require taxdb.*.
        - CPU-bound; set `threads` as appropriate.

    Example:
        >>> hits = blastp_from_seq("MTEYKLVVVGAGGVGKSALTIQLIQNHF...",
        ...                         db="pdbaa", max_targets=5, threads=4,
        ...                         blastdb="/home/ubuntu/blastdb")
        >>> [(h["sacc"], h["pident"]) for h in hits][:3]
        [('sp|P01112|', '95.1'), ('sp|Q61PZ0|', '93.7'), ...]
    """
    if blastdb:
        os.environ["BLASTDB"] = blastdb
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
        fh.write(textwrap.dedent(f">q\n{seq}\n"))
        q = fh.name
    out = _run(
        [
            "blastp",
            "-query",
            q,
            "-db",
            db,
            "-max_target_seqs",
            str(max_targets),
            "-num_threads",
            str(threads),
            "-outfmt",
            BLAST_OUTFMT,
        ]
    )
    reader = csv.DictReader(
        out.splitlines(),
        fieldnames=["sacc", "pident", "length", "evalue", "bitscore", "sscinames"],
        delimiter="\t",
    )
    return list(reader)


# ----------------------------- Meeko -----------------------------------------


def meeko_ligand(in_sdf: str | os.PathLike, out_pdbqt: str | os.PathLike) -> None:
    """Convert a small-molecule SDF to PDBQT using Meeko.
    Args:
        in_sdf: Input .sdf path.
        out_pdbqt: Output .pdbqt path.
    Example:
        >>> meeko_ligand("lig.sdf", "lig.pdbqt")
    """
    _run(["mk_prepare_ligand.py", "-i", str(in_sdf), "-o", str(out_pdbqt)])


def meeko_receptor(in_pdb: str | os.PathLike, out_pdbqt: str | os.PathLike) -> None:
    """Convert a receptor PDB to PDBQT using Meeko.
    Args:
        in_pdb: Input .pdb path.
        out_pdbqt: Output .pdbqt path.
    Example:
        >>> meeko_receptor("protein.pdb", "receptor.pdbqt")
    """
    _run(["mk_prepare_receptor.py", "-i", str(in_pdb), "-o", str(out_pdbqt)])


# ----------------------------- SDF generation --------------------------------


def smiles_to_sdf_3d(smiles: str, out_sdf: str | os.PathLike) -> None:
    """Build a 3D SDF from SMILES via RDKit, falling back to Open Babel if RDKit missing.
    Args:
        smiles: SMILES string.
        out_sdf: Output .sdf path.
    Notes:
        - RDKit path (preferred): Embed with ETKDG + UFF optimize.
        - Fallback requires `obabel` in PATH: `obabel -:'<smiles>' -O out.sdf --gen3d`
    Example:
        >>> smiles_to_sdf_3d("c1ccccc1", "benzene.sdf")
    """
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem  # type: ignore
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)
        w = Chem.SDWriter(str(out_sdf))
        w.write(mol)
        w.close()
    except Exception:
        if not shutil.which("obabel"):
            raise RuntimeError("Need RDKit or Open Babel (obabel) to build 3D SDF.")
        _run(["obabel", f"-:{smiles}", "-O", str(out_sdf), "--gen3d"])


# ----------------------------- AutoDock Vina ---------------------------------


def vina_score_only(
    receptor_pdbqt: str | os.PathLike,
    ligand_pdbqt: str | os.PathLike,
    config_txt: str | os.PathLike | None = None,
) -> Tuple[float, str]:
    """Run Vina in score-only mode and return (best_affinity_kcal_per_mol, raw_stdout).
    Args:
        receptor_pdbqt: Receptor .pdbqt path.
        ligand_pdbqt: Ligand .pdbqt path.
        config_txt: Optional Vina config file (box, etc.).
    Notes:
        - Affinity is in kcal/mol (more negative is better).
        - If parsing fails, best affinity is NaN but full stdout is returned.
    Example:
        >>> best, log = vina_score_only("receptor.pdbqt", "ligand.pdbqt", "config.txt")
        >>> print(best)  # -7.4 (example)
    """
    cmd = ["vina", "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt), "--score_only"]
    if config_txt:
        cmd += ["--config", str(config_txt)]
    out = _run(cmd)
    m = re.search(r"Affinity:\s*(-?\d+(?:\.\d+)?)", out)
    best = float(m.group(1)) if m else float("nan")
    return best, out


def vina_dock(
    receptor_pdbqt: str | os.PathLike,
    ligand_pdbqt: str | os.PathLike,
    out_pdbqt: str | os.PathLike,
    *,
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    size: Tuple[float, float, float] = (20.0, 20.0, 20.0),
    num_modes: int = 9,
    exhaustiveness: int = 8,
) -> str:
    """Dock with Vina and write poses to `out_pdbqt`. Returns Vina stdout.
    Args:
        receptor_pdbqt: Receptor .pdbqt.
        ligand_pdbqt: Ligand .pdbqt.
        out_pdbqt: Output poses .pdbqt.
        center: (cx, cy, cz) Å box center.
        size: (sx, sy, sz) Å box size.
        num_modes: Number of poses to output.
        exhaustiveness: Search thoroughness (time ↑ with value).
    Example:
        >>> log = vina_dock("receptor.pdbqt","ligand.pdbqt","poses.pdbqt",
        ...                 center=(10,12,5), size=(22,22,22), num_modes=10)
    """
    cx, cy, cz = map(str, center)
    sx, sy, sz = map(str, size)
    cmd = [
        "vina",
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        cx,
        "--center_y",
        cy,
        "--center_z",
        cz,
        "--size_x",
        sx,
        "--size_y",
        sy,
        "--size_z",
        sz,
        "--num_modes",
        str(num_modes),
        "--exhaustiveness",
        str(exhaustiveness),
        "--out",
        str(out_pdbqt),
    ]
    return _run(cmd)

