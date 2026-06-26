###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           log_msg.py
###############################################################################


from datetime import datetime


# Print a timestamped workflow message for shell-script logs.
def log_msg(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{timestamp} | {message}", flush=True)
