class ReadLinesError(Exception):

    def __init__(self, filename, msg=None):
        if msg is None:
            msg = "File '{}' doesn't exist".format(filename)
        self.msg = msg
        self.filename = filename


class ReadConfigError(Exception):

    def __init__(self, msg=None):
        if msg is None:
            msg = "Error in config file"
        self.msg = msg


class LabelPowerError(ReadConfigError):

    def __init__(self, msg=None):
        if msg is None:
            msg = "ERROR in config-file!Label power of nitrogen labels is insufficient." \
                  " At least one of 13C-spectra and one of 13C-labels should be used."
        super().__init__(msg=msg)