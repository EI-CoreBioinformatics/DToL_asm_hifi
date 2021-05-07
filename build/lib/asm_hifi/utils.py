import os
import sys
import errno
from pathlib import Path

def createDirIfNotExist(path_as_string):
	try:
		os.makedirs(path_as_string)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise
	return path_as_string

def assertFilenameValid(path_as_string,err_msg=None):
	try:
		return str(Path(path_as_string).resolve(strict=True))
	except FileNotFoundError:
		if err_msg == None:
			print("File Not Found Error: Not found " + path_as_string)
		else:
			print(err_msg + ": " + path_as_string)
		sys.exit(1)
