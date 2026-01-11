import os
import time
import re
import easygui as eg
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.firefox.options import Options
from selenium.common.exceptions import TimeoutException
from pdfreader import PDFDocument, SimplePDFViewer

# =========================
#   INTERACTIVE WINDOW
# =========================
# Interactive window for species choice 
species_str = [
    f"{scientific} : {short}"
    for short, scientific in [
        ('human', 'Homo sapiens'),
        ('fly', 'Drosophila melanogaster'),
        ('arabidopsis', 'Arabidopsis thaliana'),
        ('brugia', 'Brugia malayi'),
        ('aedes', 'Aedes aegypti'),
        ('tribolium', 'Tribolium castaneum'),
        ('schistosoma', 'Schistosoma mansoni'),
        ('tetrahymena', 'Tetrahymena thermophila'),
        ('galdieria', 'Galdieria sulphuraria'),
        ('maize', 'Zea mays'),
        ('toxoplasma', 'Toxoplasma gondii'),
        ('caenorhabditis', 'Caenorhabditis elegans'),
        ('(elegans)', 'Caenorhabditis elegans'),
        ('aspergillus_fumigatus', 'Aspergillus fumigatus'),
        ('aspergillus_nidulans', 'Aspergillus nidulans'),
        ('(anidulans)', 'Aspergillus nidulans'),
        ('aspergillus_oryzae', 'Aspergillus oryzae'),
        ('aspergillus_terreus', 'Aspergillus terreus'),
        ('botrytis_cinerea', 'Botrytis cinerea'),
        ('candida_albicans', 'Candida albicans'),
        ('candida_guilliermondii', 'Candida guilliermondii'),
        ('candida_tropicalis', 'Candida tropicalis'),
        ('chaetomium_globosum', 'Chaetomium globosum'),
        ('coccidioides_immitis', 'Coccidioides immitis'),
        ('coprinus', 'Coprinus cinereus'),
        ('coprinus_cinereus', 'Coprinus cinereus'),
        ('coyote_tobacco', 'Nicotiana attenuata'),
        ('cryptococcus_neoformans_gattii', 'Cryptococcus neoformans gattii'),
        ('cryptococcus_neoformans_neoformans_B', 'Cryptococcus neoformans neoformans'),
        ('cryptococcus_neoformans_neoformans_JEC21', 'Cryptococcus neoformans neoformans'),
        ('(cryptococcus)', 'Cryptococcus neoformans'),
        ('debaryomyces_hansenii', 'Debaryomyces hansenii'),
        ('encephalitozoon_cuniculi_GB', 'Encephalitozoon cuniculi'),
        ('eremothecium_gossypii', 'Eremothecium gossypii'),
        ('fusarium_graminearum', 'Fusarium graminearum'),
        ('(fusarium)', 'Fusarium graminearum'),
        ('histoplasma_capsulatum', 'Histoplasma capsulatum'),
        ('(histoplasma)', 'Histoplasma capsulatum'),
        ('kluyveromyces_lactis', 'Kluyveromyces lactis'),
        ('laccaria_bicolor', 'Laccaria bicolor'),
        ('lamprey', 'Petromyzon marinus'),
        ('leishmania_tarentolae', 'Leishmania tarentolae'),
        ('lodderomyces_elongisporus', 'Lodderomyces elongisporus'),
        ('magnaporthe_grisea', 'Magnaporthe grisea'),
        ('neurospora_crassa', 'Neurospora crassa'),
        ('(neurospora)', 'Neurospora crassa'),
        ('phanerochaete_chrysosporium', 'Phanerochaete chrysosporium'),
        ('(pchrysosporium)', 'Phanerochaete chrysosporium'),
        ('pichia_stipitis', 'Pichia stipitis'),
        ('rhizopus_oryzae', 'Rhizopus oryzae'),
        ('saccharomyces_cerevisiae_S288C', 'Saccharomyces cerevisiae'),
        ('saccharomyces_cerevisiae_rm11-1a_1', 'Saccharomyces cerevisiae'),
        ('(saccharomyces)', 'Saccharomyces cerevisiae'),
        ('schizosaccharomyces_pombe', 'Schizosaccharomyces pombe'),
        ('thermoanaerobacter_tengcongensis', 'Thermoanaerobacter tengcongensis'),
        ('trichinella', 'Trichinella spiralis'),
        ('ustilago_maydis', 'Ustilago maydis'),
        ('(ustilago)', 'Ustilago maydis'),
        ('yarrowia_lipolytica', 'Yarrowia lipolytica'),
        ('nasonia', 'Nasonia vitripennis'),
        ('tomato', 'Solanum lycopersicum'),
        ('chlamydomonas', 'Chlamydomonas reinhardtii'),
        ('amphimedon', 'Amphimedon queenslandica'),
        ('pneumocystis', 'Pneumocystis jirovecii'),
        ('wheat', 'Triticum aestivum'),
        ('chicken', 'Gallus gallus'),
        ('zebrafish', 'Danio rerio'),
        ('E_coli_K12', 'Escherichia coli'),
        ('s_aureus', 'Staphylococcus aureus'),
        ('volvox', 'Volvox carteri'),
    ]
]

msg ="Choose an identifier"
title = "Identifier | Species"
chosen_species = eg.choicebox(msg, title, species_str)
c_s = chosen_species.split(":")[1].replace(" ","")

# Interactive window for file selection
sequence = eg.fileopenbox(multiple=False)
file_name = os.path.basename(sequence).split(".")[0]
if not os.path.isdir(f"{c_s}_{file_name}"):
    os.system(f"mkdir {c_s}_{file_name}") # Creating a folder for the given sequence
os.chdir(f"{c_s}_{file_name}")

# =========================
#   TRANSPOSABLE ELEMENTS
#   DETECTION AND SOFT MASKING 
# ========================
if not os.path.isdir(f"TE"):
    os.system(f"mkdir TE")
os.chdir("../TE_DB")

# Using Blast with TrepDB to detect TE
os.system(f"makeblastdb -in trep_db_renamed.fasta -dbtype nucl -parse_seqids")
os.system(f"blastn -db trep_db_renamed.fasta -query {sequence} -outfmt 7 -out ../{c_s}_{file_name}/TE/Blast_TE_output7") 
os.system(f"blastn -db trep_db_renamed.fasta -query {sequence} -out ../{c_s}_{file_name}/TE/Blast_TE_output") 
print("Blast with TrepDB done.")

# Using Censor to detect and softmask TE
os.system("sudo service tor stop")
os.system("sudo service tor start") # OH HELL NAH I'M NOT PAYING
opts = Options()
opts.headless = True
opts.set_preference("print.print_headerleft", "")
opts.set_preference("print.print_headerright", "")
opts.set_preference("print.print_footerleft", "")
opts.set_preference("print.print_footerright", "")

output_pdf = os.path.join(os.path.abspath(f"../{c_s}_{file_name}"), "result.pdf")
opts.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/pdf")
opts.set_preference("print.always_print_silent", True)
opts.set_preference("print.print_to_file", True)
opts.set_preference("print.print_to_filename", output_pdf)
opts.set_preference("print_printer", "Mozilla Save to PDF")
opts.set_preference("print.printer_Mozilla_Save_to_PDF.print_to_file", True)
opts.set_preference("print.printer_Mozilla_Save_to_PDF.print_to_filename", output_pdf)
opts.set_preference("pdfjs.disabled", True) 

# Basically I don't want to pay so I change my proxy with Tor (Kinda illegal btw)
opts.set_preference("network.proxy.type", 1)
opts.set_preference("network.proxy.socks", "127.0.0.1")
opts.set_preference("network.proxy.socks_port", 9050)
opts.set_preference("network.proxy.socks_remote_dns", True)

driver = webdriver.Firefox(options=opts)
driver.get("https://www.girinst.org/censor/")
driver.set_page_load_timeout(30)

# Choosing the taxon
name = chosen_species.split(":")[0][:-1]
name1 = name
if name == "Triticum aestivum":
    name = "Triticum"
driver.find_element(By.CSS_SELECTOR,"select[name='taxon']").click()
options = driver.find_elements(By.CSS_SELECTOR, "option")
for option in options:
    if name in option.text:
        option.click()
driver.find_element(By.CSS_SELECTOR, "div[id='content']").click()

# Okay time to submit the file
file_input = driver.find_element(By.CSS_SELECTOR, "input[type='file']")
file_input.send_keys(os.path.abspath(sequence))
driver.find_element(By.CSS_SELECTOR, "input[type='submit']").click()

# Waiting for the new page to load
time.sleep(5)
wait = WebDriverWait(driver, 40)
link = wait.until(
    EC.presence_of_element_located((By.XPATH, "//a[@href]"))
)
link.click()

# Saving the results
with open (f'../{c_s}_{file_name}/{file_name}_masked.fasta', 'w') as sq:
    sq.write(driver.find_element(By.CSS_SELECTOR, "pre").text) # OMG ? The masked sequence is here
driver.execute_script("window.print();")
os.chdir(f"../{c_s}_{file_name}")
time.sleep(7)
os.system(f"mv result.pdf TE")

# Perfect we detected the TE !
print("Softmasking + detection done.")


# =========================
#   GENE PREDICTION
# =========================
def augustus(masked: bool):
    # Using Augustus to predict genes
    folder = "augustus_pred"
    if masked :
        os.system(f"augustus --species={c_s} --codingseq=on --exonnames=on {file_name}_masked.fasta > augustus_pred_masked/augustus_result.txt")
        folder = "augustus_pred_masked"
        os.system('grep -v "^#" augustus_pred_masked/augustus_result.txt > augustus_pred_masked/augustus_result.gff3')
    else:
        os.system(f"augustus --species={c_s} --codingseq=on --exonnames=on {sequence} > augustus_pred/augustus_result.txt")
        os.system('grep -v "^#" augustus_pred/augustus_result.txt > augustus_pred/augustus_result.gff3')

    # Extracting every predicted proteins, predicted CDS and predicted genes
    with open(f"{folder}/augustus_result.txt", 'r', encoding='utf-8') as res:
        with open(f"{folder}/predicted_genes_by_augustus.txt", 'w', encoding="utf-8") as pred :
            prots = ""
            cds = ""
            gene_found, protein_found, CDS_found = False, False, False
            for line in res:
                if line.startswith("# start gene"):
                    gene_found = True
                if line.startswith("# end gene"):
                    pred.write(line)
                    gene_found = False
                if gene_found:
                    if line.startswith("# coding sequence"):
                        CDS_found = True
                    if line.startswith("# protein sequence"):
                        CDS_found = False
                        protein_found = True
                    if line.startswith("# Evidence "):
                        protein_found = False
                    pred.write(line)
                if CDS_found:
                    line = line.replace("coding sequence = ", "").replace("# ", "")
                    cds += line
                if protein_found:
                    line = line.replace("protein sequence = ", "").replace("# ", "")
                    prots += line

    # We'll use Blastp to check the protein predictions
    if not os.path.isdir(f"{folder}/predicted_proteins"):
        os.system(f"mkdir {folder}/predicted_proteins")

    # Creating .fasta files for proteins
    seqs = re.findall(r"\[(.*?)\]", prots, flags=re.DOTALL)
    gene_number = 0
    for seq in seqs : 
        gene_number += 1
        with open(f"{folder}/predicted_proteins/Protein{gene_number}.fasta", "w", encoding='utf-8') as pro :
            pro.write(f">Protein{gene_number}\n")
            pro.write(seq)

    # We'll use Blastn to check the mRNA predictions
    if not os.path.isdir(f"{folder}/predicted_mRNA"):
        os.system(f"mkdir {folder}/predicted_mRNA")

    # Creating .fasta files for mRNA
    seqs1 = re.findall(r"\[(.*?)\]", cds, flags=re.DOTALL)
    gene_number = 0
    for seq in seqs1 : 
        gene_number += 1
        with open(f"{folder}/predicted_mRNA/mRNA{gene_number}.fasta", "w", encoding='utf-8') as mrna :
            mrna.write(f">mRNA{gene_number}\n")
            mrna.write(seq)

def fgenesh(masked: bool):
    driver.get("http://www.softberry.com/berry.phtml?topic=fgenesh&group=programs&subgroup=gfind")
    folder = "FGENESH"
    if masked:
        folder = "FGENESH_masked"

    # Submitting the fasta file
    file_input2 = driver.find_element(By.CSS_SELECTOR, "input[type='file']")
    if masked :
        file_input2.send_keys(os.path.abspath(f"{file_name}_masked.fasta"))
    else :
        file_input2.send_keys(os.path.abspath(f"{sequence}"))

    try:
        # Choosing the organism
        driver.find_element(By.XPATH, "//span[contains(text(), 'Organism')]").click()
        driver.find_element(By.CSS_SELECTOR, "input[tabindex='5']").send_keys(name)
        driver.find_element(By.CSS_SELECTOR, "li.active-result").click() # First result

        # Getting the .pdf file of the results
        driver.find_element(By.CSS_SELECTOR, "input[type='submit']").click()
        wait.until(
            EC.presence_of_element_located((By.XPATH, "//a[contains(text(), 'Show picture of predicted genes in PDF file')]"))
        ).click()
        download_dir = "/home/crakshay/Téléchargements"
        files = [os.path.join(download_dir, f) for f in os.listdir(download_dir)]
        latest_pdf = max(files, key=os.path.getmtime)
        os.system(f"mv {latest_pdf} {folder}")

        # Parsing the .pdf file of the results
        with open(f"{folder}/getfile.pdf", "rb") as f:
            viewer = SimplePDFViewer(f)
            doc = PDFDocument(f)
            whole_file = ""
            nb_pages = len(list(doc.pages()))+1
            for page in range(1, nb_pages): 
                viewer.navigate(page)
                viewer.render()
                whole_file += "".join(viewer.canvas.strings) + "\n"

        # Creating .fasta files for proteins and mRNA
        seqs = re.findall(r'>(FGENESH[^\n>]*)([\s\S]*?)(?=>FGENESH|$)', whole_file, flags=re.DOTALL)
        prot_number = 0
        mRNA_number = 0
        for seq in seqs : 
            if not 'mRNA' in seq[0]:
                with open(f"{folder}/predicted_proteins/FGENESH_prot{prot_number}.fasta", "w", encoding='utf-8') as pro :
                    prot_number += 1
                    comma = seq[0].index(",")
                    try:
                        oof = seq[0].index("+", comma)
                    except:
                        oof = seq[0].index("-", comma)
                    se = seq[0][:oof+1] + "\n" + seq[0][oof+1:]
                    pro.write(f">{se+seq[1]}")
            else:
                with open(f"{folder}/predicted_mRNA/FGENESH_mRNA{mRNA_number}.fasta", "w", encoding='utf-8') as mRNA :
                    mRNA_number += 1
                    comma = seq[0].index(",")
                    try:
                        oof = seq[0].index("+", comma)
                    except:
                        oof = seq[0].index("-", comma)
                    se = seq[0][:oof+1] + "\n" + seq[0][oof+1:]
                    mRNA.write(f">{se+seq[1]}")
    except:
        pass

# Using Augustus to predict genes and proteins
if not os.path.isdir("augustus_pred"): # Creating a folder for the predictions
    os.system("mkdir augustus_pred")
augustus(False)
if not os.path.isdir("augustus_pred_masked"): 
    os.system("mkdir augustus_pred_masked")
augustus(True)

# Using FGENESH to predict genes and proteins
use_fgenesh = eg.boolbox("Do you want to use FGENESH? (Plant Genome only !)", choices=("YES", "NO"))
if use_fgenesh:
    # Non softmasked
    if not os.path.isdir("FGENESH"):
        os.system("mkdir FGENESH")
    if not os.path.isdir("FGENESH/predicted_proteins"):
        os.system("mkdir FGENESH/predicted_proteins")
    if not os.path.isdir("FGENESH/predicted_mRNA"):
        os.system("mkdir FGENESH/predicted_mRNA")
    fgenesh(False)

    # Softmasked
    if not os.path.isdir("FGENESH_masked"):
        os.system("mkdir FGENESH_masked")
    if not os.path.isdir("FGENESH_masked/predicted_proteins"):
        os.system("mkdir FGENESH_masked/predicted_proteins")
    if not os.path.isdir("FGENESH_masked/predicted_mRNA"):
        os.system("mkdir FGENESH_masked/predicted_mRNA")
    fgenesh(True)
else :
    print("No prediction made by FGENESH")

print("Proteins and genes prediction done. Check the folders.")

# Creating folders for the specific databases
if not os.path.isdir("blastp_results"): # Creating a folder for the predictions
    os.system("mkdir blastp_results")
if not os.path.isdir("blastp_results/nr"):
    os.system("mkdir blastp_results/nr")
if not os.path.isdir("blastp_results/swiss"):
    os.system("mkdir blastp_results/swiss")

if not os.path.isdir("blastn_results"): # Creating a folder for the mRNA predictions
    os.system("mkdir blastn_results")
if not os.path.isdir("blastx_result"): 
    # Creating a folder for overall protein sequences validation
    os.system("mkdir blastx_result")

# =========================
#   PROTEIN PREDICTIONS
#   mRNA PREDICTIONS
#   VALIDATION
# =========================
def blast(dr, fi, type, directory, BLOSUM, DATABASE):
    try:
        dr.get(f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blast{type}&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome")
    except TimeoutException:
        dr.refresh()
    file_input1 = dr.find_element(By.CSS_SELECTOR, "input[type='file']")
    if type in "px":
        if type == "p":
            file_input1.send_keys(os.path.abspath(f"{directory}/predicted_proteins/{fi}"))
        if type == "x":
            file_input1.send_keys(os.path.abspath(fi))
        dr.find_element(By.CSS_SELECTOR,"select[name='DATABASE']").click()
        dr.find_element(By.CSS_SELECTOR, f"option[value='{DATABASE}']").click() # nr_cluster_seq / swissprot
        dr.find_element(By.CSS_SELECTOR, "button[class='usa-accordion-button']").click()
        dr.find_element(By.CSS_SELECTOR,"select[name='MATRIX_NAME']").click()
        dr.find_element(By.CSS_SELECTOR, f"option[value='{BLOSUM}']").click() #BLOSUM62
        dr.find_element(By.CSS_SELECTOR, "input[class='blastbutton']").click()
    else :
        # TSA = archive of computationally assembled transcript sequences from primary data such as ESTs and NGS.
        file_input1.send_keys(os.path.abspath(f"{directory}/predicted_mRNA/{fi}"))
        dr.find_element(By.CSS_SELECTOR,"select[name='DATABASE']").click()
        dr.find_element(By.CSS_SELECTOR, f"option[value='{DATABASE}']").click() # TSA
        dr.find_element(By.CSS_SELECTOR, "input[id='qorganism']").send_keys(name1)
        first_li = wait.until(
            EC.presence_of_element_located(
                (By.CSS_SELECTOR, "ul.ui-ncbiautocomplete-options li") # First result
            )
        )
        first_li.click()
        dr.find_element(By.CSS_SELECTOR, "input[class='blastbutton']").click()

def download(dr, file, destination, suffix):
    wait.until(
        EC.presence_of_element_located((By.XPATH, "//h1[contains(text(), 'results for')]"))
    )
    try :
        dr.find_element(By.XPATH, "//div[contains(text(), 'No significant similarity found')]")
    except:
        time.sleep(3)
        rid = dr.find_element(By.CSS_SELECTOR, "a[title='BLAST search request ID']").text
        dr.set_page_load_timeout(5)
        try:
            dr.get(f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&RID={rid}&FORMAT_TYPE=SAM_SQ&FORMAT_OBJECT=Alignment& \
               DESCRIPTIONS=100&ALIGNMENTS=100&CMD=Get&DOWNLOAD_TEMPL=Results_All&ADV_VIEW=on")
        except:
            download_dir = "/home/crakshay/Téléchargements"
            files = [os.path.join(download_dir, f) for f in os.listdir(download_dir)]
            latest_file = max(files, key=os.path.getmtime)
            os.system(f"mv {latest_file} {os.path.basename(file).split(".")[0]}{suffix}.txt")
            os.system(f"mv {os.path.basename(file).split(".")[0]}{suffix}.txt {destination}")
            time.sleep(3)
            dr.set_page_load_timeout(30)

def validate_mRNA(file, directory, suffix):
    # Blasting with TSA database
    blast(driver, file, "n", directory, '', 'tsa_nt')
    # Getting the results
    download(driver, file, "blastn_results", suffix) 

def validate_proteins(file, directory, suffix):
    # Blasting with NR database
        # Blosum62
    blast(driver, file, "p", directory, 'BLOSUM62', 'nr_cluster_seq')
    driver.execute_script("window.open('');")  
    driver.switch_to.window(driver.window_handles[1])   

    # Blasting with Swissprot database
        # Blosum62
    blast(driver, file, "p", directory, 'BLOSUM62', 'swissprot')
    driver.execute_script("window.open('');")  

    # Getting the results
    driver.switch_to.window(driver.window_handles[0])
    download(driver, file, f"blastp_results/nr", suffix) # Blosum60
    driver.switch_to.window(driver.window_handles[1])
    download(driver, file, f"blastp_results/swiss", suffix) # Blosum60
    driver.switch_to.window(driver.window_handles[0]) 

# Please Blast I need this... My protein is kinda homeless

for file in os.listdir("augustus_pred/predicted_mRNA"):
    validate_mRNA(file, "augustus_pred", "_ap")
for file in os.listdir("augustus_pred_masked/predicted_mRNA"):
    validate_mRNA(file, "augustus_pred", "_apm")

for file in os.listdir("augustus_pred/predicted_proteins"):
    validate_proteins(file, "augustus_pred", "_ap")
for file in os.listdir("augustus_pred_masked/predicted_proteins"):
    validate_proteins(file, "augustus_pred", "_apm")

if use_fgenesh:
    for file1 in os.listdir("FGENESH/predicted_mRNA"):
        validate_mRNA(file1, "FGENESH", "_fg")
    for file1 in os.listdir("FGENESH_masked/predicted_mRNA"):
        validate_mRNA(file1, "FGENESH_masked", "_fgm")

    for file1 in os.listdir("FGENESH/predicted_proteins"):
        validate_proteins(file1, "FGENESH", "_fg")
    for file1 in os.listdir("FGENESH_masked/predicted_proteins"):
        validate_proteins(file1, "FGENESH_masked", "_fgm")

print("BlastP and BlastN done. Check the folders.")

# Let's check if some proteins are missing...
blast(driver, sequence, "x", "", 'BLOSUM62', 'swissprot')
download(driver,sequence,"blastx_result", "_") 
print("BlastX done.")

# ARI ARI ARI ARIVEDERCCI
driver.quit()
os.system("sudo service tor stop")

# Little correction
if use_fgenesh:
    if not os.path.isdir("pred_by_fgenesh"):
        os.system("mkdir pred_by_fgenesh")
    os.system("mv FGENESH pred_by_fgenesh")
    os.system("mv FGENESH_masked pred_by_fgenesh")

if not os.path.isdir("pred_by_augustus"):
    os.system("mkdir pred_by_augustus")
os.system("mv augustus_pred pred_by_augustus")
os.system("mv augustus_pred_masked pred_by_augustus")

os.chdir("..")
if not os.path.isdir("Results"):
    os.system("mkdir Results")
os.system(f"mv {c_s}_{file_name} Results")
