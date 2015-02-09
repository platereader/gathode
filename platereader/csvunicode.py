"""
This module implements the classes to read and write csv supporting different encodings, both python 2 and 3.
"""

import sys
import csv
if sys.version < '3': 
    import codecs
    import cStringIO

class UTF8Recoder(object):
    """
    Iterator that reads an encoded stream and reencodes the input to UTF-8
    """
    def __init__(self, f, encoding):
        self.reader = codecs.getreader(encoding)(f)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.reader).encode("utf-8")

    def next(self):
        return self.__next__()

class CsvUnicodeReader(object):
    """
    A CSV reader which will iterate over lines in the CSV file "f",
    which is encoded in the given encoding.
    """

    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
        f = UTF8Recoder(f, encoding)
        self.reader = csv.reader(f, dialect=dialect, **kwds)

    def __iter__(self):
        return self

    def __next__(self):
        row = next(self.reader)
        return [unicode(s, "utf-8") for s in row]

    def next(self):
        return self.__next__()

class CsvFileUnicodeReader(object):
    """
    A context manager interface for a csv reader reading from file, handling python 2 vs 3 differences.
    """

    def __init__(self,filename,dialect=csv.excel,encoding="utf-8",**kwds):
        self.filename = filename
        self.dialect = dialect
        self.encoding = encoding
        self.kwds = kwds

    def __enter__(self):
        if sys.version < '3': 
            self.f = open(self.filename, 'rb')
            self.reader = CsvUnicodeReader(self.f,dialect=self.dialect,
                                           encoding=self.encoding,**self.kwds)
        else:
            self.f = open(self.filename, 'r', newline='', encoding=self.encoding)
            self.reader = csv.reader(self.f,dialect=self.dialect,**self.kwds)
        return self.reader

    def __exit__(self,type,value,traceback):
        self.f.close()

    def __iter__(self):
        return self.reader

    def __next__(self):
        return next(self.reader)

    def next(self):
        return self.__next__()


class CsvUnicodeWriter:
    """
    A CSV writer which will write rows to CSV file "f",
    which is encoded in the given encoding.
    """

    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
        # Redirect output to a queue
        self.queue = cStringIO.StringIO()
        self.writer = csv.writer(self.queue, dialect=dialect, **kwds)
        self.stream = f
        self.encoder = codecs.getincrementalencoder(encoding)()

    def writerow(self, row):
        encodedrow=[]
        for s in row:
            # do not encode numerical types, as we may loose precision
            # FIXME this list is not exhaustive
            if type(s) == str or type(s) == unicode:
                encodedrow.append(s.encode("utf-8"))
            else:
                encodedrow.append(s)
        self.writer.writerow(encodedrow)
        # Fetch UTF-8 output from the queue ...
        data = self.queue.getvalue()
        data = data.decode("utf-8")
        # ... and reencode it into the target encoding
        data = self.encoder.encode(data)
        # write to the target stream
        self.stream.write(data)
        # empty queue
        self.queue.truncate(0)

    def writerows(self, rows):
        for row in rows:
            self.writerow(row)

class CsvFileUnicodeWriter(object):
    """
    A context manager interface for a csv writer writing to a file, handling python 2 vs 3 differences.
    """

    def __init__(self,filename,dialect=csv.excel,encoding="utf-8",**kwds):
        self.filename = filename
        self.dialect = dialect
        self.encoding = encoding
        self.kwds = kwds

    def __enter__(self):
        if sys.version < '3': 
            self.f = open(self.filename, 'wb')
            self.writer = CsvUnicodeWriter(self.f,dialect=self.dialect,
                                           encoding=self.encoding,**self.kwds)
        else:
            self.f = open(self.filename, 'w', newline='', encoding=self.encoding)
            self.writer = csv.writer(self.f,dialect=self.dialect,**self.kwds)
        return self.writer

    def __exit__(self,type,value,traceback):
        self.f.close()

    def writerow(self, row):
        self.writer.writerow(row)

    def writerows(self, rows):
        self.writer.writerows(rows)
