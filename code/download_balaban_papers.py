# -*- coding: utf-8 -*-
"""
Download Balaban Papers and Dimock Pedagogical Series

This script downloads the foundational papers for the rigorous
Balaban multi-scale renormalization group analysis.
"""

import os
import sys
import urllib.request
import ssl
from pathlib import Path

# Create output directories
BASE_DIR = Path(__file__).parent.parent
PAPERS_DIR = BASE_DIR / "papers"
BALABAN_DIR = PAPERS_DIR / "balaban"
DIMOCK_DIR = PAPERS_DIR / "dimock"

BALABAN_DIR.mkdir(parents=True, exist_ok=True)
DIMOCK_DIR.mkdir(parents=True, exist_ok=True)

# Dimock papers on arXiv (freely accessible) - CORRECT arXiv IDs
DIMOCK_PAPERS = [
    {
        "arxiv_id": "1108.1335",
        "title": "The Renormalization Group According to Balaban I - Small Fields",
        "year": 2013,
        "filename": "dimock_balaban_I_small_fields.pdf"
    },
    {
        "arxiv_id": "1212.5562",
        "title": "The Renormalization Group According to Balaban II - Large Fields",
        "year": 2013,
        "filename": "dimock_balaban_II_large_fields.pdf"
    },
    {
        "arxiv_id": "1304.0705",
        "title": "The Renormalization Group According to Balaban III - Convergence",
        "year": 2013,
        "filename": "dimock_balaban_III_convergence.pdf"
    },
]

# Balaban papers - Springer links (require institutional access)
BALABAN_PAPERS = [
    {
        "doi": "10.1007/BF01206179",
        "title": "Propagators and renormalization transformations I",
        "journal": "Commun. Math. Phys. 95",
        "year": 1984,
        "pages": "17-40",
        "filename": "balaban_1984_cmp95_propagators_I.pdf"
    },
    {
        "doi": "10.1007/BF01211907",
        "title": "Propagators and renormalization transformations II",
        "journal": "Commun. Math. Phys. 96",
        "year": 1984,
        "pages": "223-250",
        "filename": "balaban_1984_cmp96_propagators_II.pdf"
    },
    {
        "doi": "10.1007/BF01206178",
        "title": "Averaging operations for lattice gauge theories",
        "journal": "Commun. Math. Phys. 98",
        "year": 1985,
        "pages": "17-51",
        "filename": "balaban_1985_cmp98_averaging.pdf"
    },
    {
        "doi": "10.1007/BF01212893",
        "title": "Ultraviolet stability of 3D lattice pure gauge field theories",
        "journal": "Commun. Math. Phys. 102",
        "year": 1985,
        "pages": "255-275",
        "filename": "balaban_1985_cmp102_uv_stability_3d.pdf"
    },
    {
        "doi": "10.1007/BF01215324",
        "title": "RG approach to lattice gauge theories I: Cluster expansion",
        "journal": "Commun. Math. Phys. 109",
        "year": 1987,
        "pages": "249-301",
        "filename": "balaban_1987_cmp109_rg_cluster.pdf"
    },
    {
        "doi": "10.1007/BF01205037",
        "title": "Spaces of regular gauge field configurations",
        "journal": "Commun. Math. Phys. 116",
        "year": 1988,
        "pages": "1-22",
        "filename": "balaban_1988_cmp116_regular_spaces.pdf"
    },
    {
        "doi": "10.1007/BF01207382",
        "title": "Large field renormalization I: The basic step",
        "journal": "Commun. Math. Phys. 122",
        "year": 1989,
        "pages": "175-202",
        "filename": "balaban_1989_cmp122_large_field_I.pdf"
    },
    {
        "doi": "10.1007/BF01217969",
        "title": "Large field renormalization II: Localization, exponentiation",
        "journal": "Commun. Math. Phys. 122",
        "year": 1989,
        "pages": "355-392",
        "filename": "balaban_1989_cmp122_large_field_II.pdf"
    },
]


def download_arxiv_pdf(arxiv_id: str, output_path: Path) -> bool:
    """Download PDF from arXiv."""
    url = f"https://arxiv.org/pdf/{arxiv_id}.pdf"

    try:
        # Create SSL context that doesn't verify (for some institutional networks)
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE

        print(f"  Downloading {arxiv_id}...")

        # Use context with urllib
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)')

        with urllib.request.urlopen(req, context=ctx, timeout=60) as response:
            with open(output_path, 'wb') as f:
                f.write(response.read())

        print(f"  -> Saved to {output_path.name}")
        return True
    except Exception as e:
        print(f"  ERROR: {e}")
        return False


def create_bibtex_entry(paper: dict, paper_type: str) -> str:
    """Create BibTeX entry for a paper."""
    if paper_type == "arxiv":
        return f"""@article{{dimock:{paper['arxiv_id']},
    author = {{Dimock, Jonathan}},
    title = {{{paper['title']}}},
    year = {{{paper['year']}}},
    eprint = {{{paper['arxiv_id']}}},
    archivePrefix = {{arXiv}},
    primaryClass = {{math-ph}}
}}
"""
    else:  # springer
        return f"""@article{{balaban:{paper['doi'].split('/')[-1]},
    author = {{Balaban, Tadeusz}},
    title = {{{paper['title']}}},
    journal = {{{paper['journal']}}},
    year = {{{paper['year']}}},
    pages = {{{paper['pages']}}},
    doi = {{{paper['doi']}}}
}}
"""


def main():
    print("=" * 60)
    print("BALABAN PAPERS DOWNLOAD SCRIPT")
    print("=" * 60)

    # Download Dimock papers from arXiv (freely accessible)
    print("\n[1/2] Downloading Dimock papers from arXiv...")
    print("-" * 40)

    dimock_success = 0
    for paper in DIMOCK_PAPERS:
        output_path = DIMOCK_DIR / paper["filename"]
        if output_path.exists():
            print(f"  {paper['arxiv_id']}: Already exists, skipping")
            dimock_success += 1
        elif download_arxiv_pdf(paper["arxiv_id"], output_path):
            dimock_success += 1

    print(f"\nDimock: {dimock_success}/{len(DIMOCK_PAPERS)} papers downloaded")

    # Create reference for Balaban papers (require institutional access)
    print("\n[2/2] Creating Balaban paper references...")
    print("-" * 40)
    print("NOTE: Balaban papers are on Springer and require institutional access.")
    print("      Creating reference file with DOI links...")

    refs_file = BALABAN_DIR / "PAPER_REFERENCES.md"
    with open(refs_file, "w") as f:
        f.write("# Balaban Papers - Reference Guide\n\n")
        f.write("These papers require institutional access via Springer.\n\n")
        f.write("## Download Links\n\n")

        for paper in BALABAN_PAPERS:
            f.write(f"### {paper['title']}\n")
            f.write(f"- **Journal**: {paper['journal']} ({paper['year']})\n")
            f.write(f"- **Pages**: {paper['pages']}\n")
            f.write(f"- **DOI**: https://doi.org/{paper['doi']}\n")
            f.write(f"- **Direct Link**: https://link.springer.com/article/{paper['doi']}\n")
            f.write(f"- **Save as**: `{paper['filename']}`\n\n")

    print(f"  Reference file created: {refs_file}")

    # Create combined BibTeX file
    bibtex_file = PAPERS_DIR / "yang_mills_references.bib"
    with open(bibtex_file, "w") as f:
        f.write("% Yang-Mills Mass Gap Proof - Bibliography\n")
        f.write("% Balaban and Dimock papers\n\n")

        for paper in DIMOCK_PAPERS:
            f.write(create_bibtex_entry(paper, "arxiv"))

        for paper in BALABAN_PAPERS:
            f.write(create_bibtex_entry(paper, "springer"))

    print(f"  BibTeX file created: {bibtex_file}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\nDimock papers (arXiv): {dimock_success}/{len(DIMOCK_PAPERS)}")
    print(f"Balaban references:    {len(BALABAN_PAPERS)} DOI links created")
    print(f"\nOutput directories:")
    print(f"  - Dimock: {DIMOCK_DIR}")
    print(f"  - Balaban: {BALABAN_DIR}")
    print(f"  - BibTeX: {bibtex_file}")

    print("\n" + "=" * 60)
    print("NEXT STEPS")
    print("=" * 60)
    print("""
1. If you have institutional access to Springer:
   - Use the DOI links in PAPER_REFERENCES.md to download Balaban papers
   - Save them with the specified filenames

2. Alternative sources:
   - Check your university library's electronic resources
   - Request via interlibrary loan
   - Check Project Euclid for some older papers

3. The Dimock papers provide a more accessible exposition of
   Balaban's methods and may be sufficient for understanding
   the key ideas.
""")


if __name__ == "__main__":
    main()
