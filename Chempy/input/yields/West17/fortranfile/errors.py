"""
Contains all error classes
"""

############################################################
# read errors
############################################################

class RecordSizeError(Exception):
    """
    Exception raised when assert_eor finds EOR == False.

    This means that the data read does not match the record length.
    """

    def __init__(self, reclen, pos):
        self.reclen = reclen
        self.pos = pos

    def __str__(self):
        """Return error message."""
        return ("\n" +
                "Not at end of data record in\n" +
                "Poistion = {:d}\n" +
                "RecLen   = {:d}"
                ).format(
            self.pos,
            self.reclen)

class ReadError(Exception):
    """
    Exception raised when data could not be read correctly.
    """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        """
        Return error message.
        """
        return (
            "\n" +
            "Error reading record from file: {}".format(
                self.message))



class DataSizeError(Exception):
    """
    Exception raised when assert_eor finds EOR == False.

    This means that the data read does not match the record length.
    """

    def __init__(self, filename, reclen, pos):
        self.filename = filename
        self.reclen = reclen
        self.pos = pos

    def __str__(self):
        """Return error message."""
        return ("Not at end of data record in file {}.\n"+\
               "Poistion = {:d}\n"+\
               "RecLen   = {:d}")\
               .format(self.filename,
                       self.pos,
                       self.reclen)

class DataFileError(Exception):
    """
    Exception raised when assert_eof finds EOF == False.

    This means that the data file size does not match the file length.
    """

    def __init__(self, filename, filesize, pos):
        self.filename = filename
        self.filesize = filesize
        self.pos = pos

    def __str__(self):
        """Return error message."""
        return ("Not at end of data file in file {}.\n"+\
               "Poistion = {:d}\n"+\
               "Filesize   = {:d}")\
               .format(self.filename,
                       self.pos,
                       self.filesize)

class RecLenError(Exception):
    """
    Exception raised when a record was not read correctly.

    This means that the record header does not match the
    record trailer (both should contain the length of the
    data area).
    """

    def __init__(self, filename, message):
        self.filename = filename # name of file which caused an error
        self.message = message

    def __str__(self):
        """Return error message."""
        return "Error reading record from file {}.\n{}"\
               .format(self.filename,self.message)

class FileReadError(Exception):
    """
    Exception raised when data could not be read correctly.
    """

    def __init__(self, filename, message):
        self.filename = filename # name of file which caused an error
        self.message = message

    def __str__(self):
        """
        Return error message.
        """
        return "Error reading record from file {}: {}"\
               .format(self.filename,self.message)


class StringReadError(Exception):
    """
    Exception raised string length is not specified.
    """

    def __init__(self, message):
        self.message = message

    def __str__(self):
        """
        Return error message.
        """
        return self.message

############################################################
# Write Errors
############################################################


class RecordBeginningError(Exception):
    """
    Exception raised when assert_bor finds BOR == False.

    This means that the data buffer is not empty.
    """

    def __init__(self, pos):
        self.pos = pos

    def __str__(self):
        """Return error message."""
        return ("\n" +
                "Not at beginning of data record:\n" +
                "Poistion = {:d}"
                ).format(
            self.pos)


class WriteError(Exception):
    """
    Exception raised when data could not be written correctly.
    """

    def __init__(self, filename, message):
        self.filename = filename
        self.message = message

    def __str__(self):
        """
        Return error message.
        """
        return "Error writing record to file {}: {}"\
               .format(self.filename,
                       self.message)
