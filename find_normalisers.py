from . import miRNA_normalisers as miR_norm
from . import check_input_files as chk_files

xl_file = miR_norm.get_filename()
sheet_name = input("Enter Sheet Name with Your Data: ")
normaliser_file = miR_norm.get_filename("normaliser")

normalisers = miR_norm.get_candidate_normalisers(normaliser_file)

if miR_norm.check_xl_file(xl_file, sheet_name, normaliser_file):
    pass
else:
    break

