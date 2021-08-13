"""
Module with the definitions of the processing errors
"""

class Registration_errors():
    """Class used to define the error that can happen during the registration
    of the reference channel.
    
    Errors and corresponding values:
    - registration file has no counts: -1
    - registration file with the counts is missing: -2
    - registration file with the counts cannot be loaded: -3
    - fish channel doesn't have counts: -4
    - fish file cannot be loaded: -5
    - registration is below extraction resolution: -6
    
    """
  
    def __init__(self):

        # The registration file has no counts
        self.missing_counts_reg_channel = -1

        # The registration file with the counts is missing
        self.missing_file_reg_channel = -2

        # The file with the counts cannot be loaded
        self.cannot_load_file_reg_channel = -3

        # Fish channel doesn't have counts
        self.missing_counts_fish_channel = -4

        # Fish channel cannot be loaded
        self.cannot_load_file_fish_channel = -5

        # Registration is below extraction resolution
        self.registration_below_extraction_resolution = -6
