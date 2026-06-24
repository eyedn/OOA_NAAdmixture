###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           simple_table.py
###############################################################################


# pattern: Functional Core


class SimpleTable:
    def __init__(self, rows):
        self.rows = rows
        self.columns = list(rows[0].keys()) if rows else []

    def to_dict(self, orient="dict"):
        if orient != "records":
            raise ValueError("SimpleTable only supports orient='records'")
        return self.rows
