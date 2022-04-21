from . import miRNA_normalisers as miR_norm
from . import check_input_files as chk_files

xl_file = miR_norm.get_filename()
sheet_name = input("Enter Sheet Name with Your Data: ")
normaliser_file = miR_norm.get_filename("normaliser")

normalisers = miR_norm.get_candidate_normalisers(normaliser_file)
if chk_files.all_checks(miR_norm.read_xl(xl_file), normaliser_file):
    print("All file checks passed: Proceeding to analysis")
else:
    print("Correct files where indicated, and re-run")

