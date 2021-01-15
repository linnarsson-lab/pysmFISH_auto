class Registration_errors():
    """
    Class used to define the error that characterize the
    possible errors in the registration

    """

    def __init__(self):

        # The registration file has no counts
        self.missing_counts_reg_channel = -1

        # The registration file with the counts is missing
        self.missing_counts_file_reg_channel = -2

        # The file with the counts cannot be loaded
        self.cannot_load_file_reg_channel = -3

        # Fish channel doesn't have counts
        self.missing_counts_fish_channel = -4

        # Fish channel doesn't have counts
        self.cannot_load_file_fish_channel = -5



