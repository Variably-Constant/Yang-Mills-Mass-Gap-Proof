# -*- coding: utf-8 -*-
"""Extract text from Dimock PDFs for analysis."""

import sys
from pathlib import Path
from PyPDF2 import PdfReader

DIMOCK_DIR = Path(r"E:\POL\Millennium Prize\papers\dimock")
OUTPUT_DIR = Path(r"E:\POL\Millennium Prize\papers\dimock_text")
OUTPUT_DIR.mkdir(exist_ok=True)

def extract_pdf(pdf_path: Path) -> str:
    """Extract text from a PDF file."""
    reader = PdfReader(pdf_path)
    text = []
    for i, page in enumerate(reader.pages):
        page_text = page.extract_text()
        if page_text:
            text.append(f"\n--- PAGE {i+1} ---\n")
            text.append(page_text)
    return "\n".join(text)

def main():
    print("Extracting text from Dimock PDFs...")

    for pdf_file in DIMOCK_DIR.glob("*.pdf"):
        print(f"  Processing: {pdf_file.name}")
        try:
            text = extract_pdf(pdf_file)
            output_file = OUTPUT_DIR / f"{pdf_file.stem}.txt"
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(text)
            print(f"    -> {output_file.name} ({len(text)} chars)")
        except Exception as e:
            print(f"    ERROR: {e}")

    print("\nDone! Text files saved to:", OUTPUT_DIR)

if __name__ == "__main__":
    main()
