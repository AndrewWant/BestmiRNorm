# Utils file for generally useful functions, etc.

import os
from datetime import datetime

def generate_output_filename(filename_stub,
                             location=""
                             ):
    """
    Creates an output filename for saving plots, allowing systematic saving of plots/images from analysis.
    :param filename_stub:
    :param location:
    :return:
    """
    if isinstance(location, str):
        if (location == "") or (not os.path.isdir(location)):
            location = os.getcwd()
    else:
        print("Invalid argument: {} supplied for location. Only strings allowed.".format(location))    
        
    today_string = datetime.today().strftime("%Y-%m-%d_%H_%M")
    return os.sep.join((location, "_".join((today_string, filename_stub))))
  
  
