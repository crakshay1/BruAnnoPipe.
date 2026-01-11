# BruAnno_Pipe

**BruAnno_Pipe** is an automated bioinformatics pipeline designed for the comprehensive annotation of genomic sequences. It integrates local command-line tools with web-based prediction and validation services to identify transposable elements, predict gene structures, and validate protein sequences.

---

## 🔍 Features

* **Interactive GUI**: Easy-to-use selection for species identifiers and input files via `easygui`.
* **Transposable Element (TE) Analysis**:
    * Local alignment against `TrepDB` using local BLASTn.
    * Automated web-submission to **Censor** (Giri Institute) for soft-masking and TE detection.
* **Gene & Protein Prediction**:
    * **Augustus**: Local prediction executed on both raw and masked sequences.
    * **FGENESH**: Automated web-based prediction (optimized for plant genomes) executed on both raw and masked sequences.
* **Automated Validation**:
    * **BLASTn**: Validation of mRNA against the TSA database.
    * **BLASTp**: Protein validation against NR and SwissProt databases.
    * **BLASTx**: Final check of the full sequence against SwissProt to catch missed coding regions.
* **Anonymity & Access**: Integrates with the **Tor** network to manage web requests.

---

## 📖 Prerequisites

### System Requirements
* **Linux/Unix** environment (recommended).
* **Mozilla Firefox** and **Geckodriver** (for Selenium automation).
* **Tor Service** (installed and configurable via `sudo service tor`).

### External Tools
The following must be installed and available in your system `$PATH`:
* **NCBI BLAST+** (makeblastdb, blastn, blastp, blastx).
* **Augustus**.

### Python Dependencies
```bash
pip install -r requirements.txt
```

---

## 📁 Setup & Directory Structure

To function correctly, the script requires a specific environment setup and will generate organized results as follows:

### 1. Pre-run Requirements
* **TE Database**: You must have a directory named `TE_DB` in the parent folder of the script.
* **Reference File**: This folder must contain the file `trep_db_renamed.fasta` for local BLAST searches.

### 2. Output Organization
The script automatically creates a `Results/` directory. Within that folder, a sub-directory is created for each run using the format `[Species]_[Filename]`.



### 3. Folder Contents
Inside each result folder, you will find:
* **`TE/`**: Contains the local BLAST output against TrepDB and the `result.pdf` masking report from Censor.
* **`pred_by_augustus/`**: Includes the full Augustus text output, as well as sub-folders for predicted proteins and mRNA in `.fasta` format.
* **`pred_by_fgenesh/`**: (If selected) Contains the parsed PDF results and corresponding `.fasta` predictions.
* **`blastp_results/`**: Organized by database (`nr/` and `swiss/`), containing the alignment reports for predicted proteins.
* **`blastn_results/`**: Contains the validation reports for predicted mRNA against the TSA database.
* **`blastx_result/`**: Contains a final validation check of the original sequence against protein databases.
