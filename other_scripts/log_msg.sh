#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           log_msg.sh
###############################################################################

# shared logging helper for timestamped pipeline messages.


##### main function ###########################################################
# print one pipeline message with a sortable wall-clock timestamp.
log_msg() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') | $1"
}
