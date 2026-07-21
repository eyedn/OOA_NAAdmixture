###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           log_msg.py
###############################################################################
# format timestamped workflow log messages.


##### set up ##################################################################
from datetime import datetime


##### main function ###########################################################
'''
print one workflow message with a sortable wall-clock timestamp. Output is
flushed immediately so Slurm logs preserve stage progress during long jobs.
'''
def log_msg(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{timestamp} | {message}", flush=True)
