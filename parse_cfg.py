#!/home/observer/miniconda2/bin/python

import sys, os
from collections import namedtuple

def ctp(val):
    try:
      ans=int(val)
      return ans
    except ValueError:
      try:
        ans=float(val)
        return ans
      except ValueError:
        if val.lower()=="false":
          return False
        elif val.lower()=="true":
          return True
        else:
          return val


def parse_cfg(cfg_file, return_dict = False, comment_identifier = ("#"), delim = None, ret_comments=False, ret_studs=False):
  '''
  Reads any text configuration file into a namedtuple (or dict).

  The first column becomes the key and the second/rest column(s) become(s) the value. If there are more than one values for a key, all the values stay in a string as a single value, unless 'delim' is specified in which case the values will be split into an array using 'delim' as the delimeter.
  (Note: Any occurence of the specified 'delim' in the first column (i.e. keys) will be ignored)

  You may specify a comment_identifier string, which by default is '#'. All lines starting with the 'comment_identifier' are ignored unless ret_comments is true when all the comments are appended to a list with key name 'Comments_'. Empty comments will be ignored. Comments towards the end of a line (e.g. :  x  5  #sample comment) will be appended to comments list if ret_comments is True.

  Any lines with only one column will be ignored unless ret_studs is True when they are appended to a list with key name 'Studs_'.

  Only alphanumeric characters are supported as valid first characters of keys if returning a named_tuple (default). Any string will work as key if return_dict is True.
  '''
  if not os.path.exists(cfg_file):
    raise IOError("File {0} does not exist".format(cfg_file))
  
  if not os.access(cfg_file, os.R_OK):
    raise IOError("Do not have read permissions on {0}".format(cfg_file))
  
  conf_dict = {}
  if ret_comments:
    conf_dict["Comments_"] = []
  if ret_studs:
    conf_dict["Studs_"] = []

  with open(cfg_file) as c:
    while True:
      line = c.readline()
      if line=="":
        break
      elif line.strip() == "":
        continue
      elif len(line.split()) < 2:
        if ret_studs:
          conf_dict["Studs_"].append(line.strip())
        continue
      elif line.strip().startswith(comment_identifier):
        if len(line.strip().strip(delim)) > 0 and ret_comments:
          conf_dict["Comments_"].append(line.strip())
        continue
      else:
        key = line.split()[0].strip()
        val = line.strip()[len(key):].strip().split(comment_identifier)[0].strip()
        val = ctp(val)

        if comment_identifier in line and ret_comments:
          conf_dict["Comments_"].append(line.split(comment_identifier)[1].strip())

        if delim is not None and isinstance(val, str):
          x = val.split(delim)
          if len(x) > 1:
            x = [ctp(i) for i in x]
            val = x
        conf_dict[key] = val
  
  if return_dict:
    return conf_dict.copy()

  tmp = namedtuple("CONF", conf_dict.keys())
  parsed_namedtuple = tmp(*conf_dict.values())

  return parsed_namedtuple


def parse_str(cfg_str, return_dict = False, comment_identifier = ("#"), delim = None, ret_comments=False, ret_studs=False, newline_char = '\n', rogue_chars = '\x00'):
  '''
  Reads any config string into a namedtuple (or dict). newline_char may be specified to something else instead of '\\n'

  The first column becomes the key and the second/rest column(s) become(s) the value. If there are more than one values for a key, all the values stay in a string as a single value, unless 'delim' is specified in which case the values will be split into an array using 'delim' as the delimeter.
  (Note: Any occurence of the specified 'delim' in the first column (i.e. keys) will be ignored)

  You may specify a comment_identifier string, which by default is '#'. All lines starting with the 'comment_identifier' are ignored unless ret_comments is true when all the comments are appended to a list with key name 'Comments_'. Empty comments will be ignored. Comments towards the end (e.g. :  x  5  #sample comment) will be appended to the comment list if ret_comments is True.

  Any lines with only one column will be ignored unless ret_studs is True when they are appended to a list with key name 'Studs_'.

  Only alphanumeric characters are supported as valid first characters of keys if returning a named_tuple (default). Any string will work as key if return_dict is True.
  '''
  
  conf_dict = {}
  if ret_comments:
    conf_dict["Comments_"] = []
  if ret_studs:
    conf_dict["Studs_"] = []
  
  for rc in rogue_chars:
    cfg_str = cfg_str.replace(rc, '')

  if len(cfg_str) == 0:
    raise ValueError("Empty string!")

  c = cfg_str.split(newline_char)
  ll = -1
  while True:
      ll += 1
      line = c[ll]
      
      if ll == (len(c)-1):
        break
      elif line.strip() == "":
        continue
      elif len(line.split()) < 2:
        if ret_studs:
          conf_dict["Studs_"].append(line.strip())
        continue
      elif line.strip().startswith(comment_identifier):
        if len(line.strip().strip(delim)) > 0 and ret_comments:
          conf_dict["Comments_"].append(line.strip())
        continue
      else:
        key = line.split()[0].strip()
        val = line.strip()[len(key):].strip().split(comment_identifier)[0].strip()
        val = ctp(val)

        if comment_identifier in line:
          conf_dict["Comments_"].append(line.split(comment_identifier)[1].strip())

        if delim is not None and isinstance(val, str):
          x = val.split(delim)
          if len(x) > 1:
            x = [ctp(i) for i in x]
            val = x
        conf_dict[key] = val

  if return_dict:
    return conf_dict.copy()

  tmp = namedtuple("CONF", conf_dict.keys())
  parsed_namedtuple = tmp(*conf_dict.values())

  return parsed_namedtuple

