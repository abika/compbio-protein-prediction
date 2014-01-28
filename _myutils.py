#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Aug  9 16:24:00 2012

@author: Alexander Bikadorov

Tested with python2!
"""

import os
import shutil
import re
import glob
import string
import random
import logging
import zipfile
import hashlib
import datetime
import importlib
import subprocess

def split_path(path):
    """Return a list of componenents of string 'path'"""
    folders=[]
    while True:
        path, folder=os.path.split(path)
        if not folder:
            if path:
                folders.append(path)
            break        
        folders.append(folder)

    folders.reverse()
    return folders

def _rename(file_path_str):
    dir_name, fname = os.path.split(file_path_str)
    bname, ext = os.path.splitext(fname)
    bname_pre = re.split('_[0-9]+$', bname)[0]
    num_files = len(files_in_dir(dir_name, bname_pre+'_[0-9]+'+ext))
    while True:
        file_path_str = os.path.join(dir_name, bname_pre+'_'+str(num_files)+ext)
        if not os.path.exists(file_path_str):
            break
        num_files += 1
    return file_path_str

def move(src_str, dest_str, rename=True, into_folder=True):
    """Move an existing file or directory from 'src_str' to 'dest_str'.
       If 'rename' is true and the destination exists already the target will be renamed.
       'into_folder' is used for that to determine if source should be moved into a destination directory or source
       should be renamed instead.
    """
    if not os.path.exists(src_str):
        logging.warning('does not exist (can not move): '+src_str)
        return False
    if into_folder and os.path.isdir(dest_str):
        if os.path.samefile(src_str, dest_str):
            logging.warning('can not move directory into itself: '+src_str)
            return False
        dest_str = os.path.join(dest_str, os.path.basename(src_str))
    if os.path.exists(dest_str):
        if rename:
            logging.info('does already exist (renaming): '+dest_str)
            dest_str = _rename(dest_str)
        else:
            logging.warning('does already exist (not moving):'+dest_str)
            return False
    logging.info('moving '+src_str+' to '+dest_str)
    shutil.move(src_str, dest_str)
    return True

def write_file(file_path_str, str_, append=False, rename=True, verbose=True):
    """Write string to file.
       If 'append' is false an existing file will not be overwritten. Instead an index is added to the filename.
       If 'append' is true the string is appended to file if it already exists, the file is created otherwise.
       Returns the file path the file was written to.
    """
    dir_name, fname = os.path.split(file_path_str)
    if not append and os.path.exists(file_path_str):
        if rename:
            logging.info('file already exists (renaming): '+file_path_str)
            file_path_str = _rename(file_path_str)
        else:
            logging.warning('file already exists (not overwriting): '+file_path_str)
            return None
    if dir_name and not os.path.isdir(dir_name):
        logging.warning('directory does not exist (can not create file)): '+file_path_str)
        return None
    file_ = open(file_path_str, 'a') if append else open(file_path_str, 'w')
    file_.write(str_)
    file_.close()
    if verbose:
        logging.info("wrote "+str(len(str_))+" bytes to file: "+file_path_str)
    return file_path_str

def read_file_lines(file_):
    """Read content of file
       file_: either path to file or a file (like) object.
       Return list of lines read.
    """
    if type(file_) is str:
        try:
            file_ = open(file_, 'r')
        except IOError:
            logging.warning('file does not exist (can not read): ' + file_)
            return []
        cont_str = file_.read()
        file_.close()
    else:
        cont_str = file_.read()
    return [url_str.rstrip() for url_str in cont_str.splitlines()]

def files_in_dir(dir_, regex='*.*'):
    """Return a list of all files in 'dir_' matching 'regex' with absolute path"""
    abs_path = os.path.abspath(dir_)
    if not os.path.isdir(abs_path):
        logging.warning('does not exist/is not a directory: '+abs_path)
    return glob.glob(os.path.join(abs_path, regex))

def open_archive(arch_path_str, verbose=False):
    """Open an existing archive for reading.
       Return a ZipFile object.
    """
    if not os.path.exists(arch_path_str):
        if verbose:
            logging.warning('archive does not exist: '+os.path.abspath(arch_path_str))
        return None
    try:
        # bug in python zip lib: can't open empty archive
        if os.path.getsize(arch_path_str) <= 22:
            return None
        archive = zipfile.ZipFile(arch_path_str,'r')
    except IOError:
        logging.warning('could not open archive '+os.path.abspath(arch_path_str))
        return None
    else:
        return archive

def get_from_archive(archive, filename):
    """Read file from archive.
       archive: ZipFile object or path to zip archive. 
       Return tuple (zipInfo, data) or (None, None) if something went wrong.
    """
    if type(archive) is str:
        if not os.path.exists(archive):
            return (None, None)
        archive = zipfile.ZipFile(archive, 'r')
    try:
        info_entry = archive.getinfo(filename)
    except (KeyError, AttributeError): # not in archive, archive does not exist
        return (None, None) 
    data_str = archive.read(filename)
    return (info_entry, data_str)

def add_to_archive(archive, file_str, file_is_path=True, name=None):
    """Add file or string to archive. 
       archive: ZipFile or path to ZipFile. If archive does not exist it is created.
       file_str: path to file or string to be written to file.  If 'file_str' is not a path 
       'file_is_path' must be set to "False".
       name: will be used to set the name of the file in archive. It can be a 
       string or ZipInfo instance. If not specified and 'file_' is a string to be written ValueError is raised.
    """
    if not file_is_path and name is None:
        raise ValueError('No filename given (can not add to archive)')        
    archive = zipfile.ZipFile(archive, 'a', zipfile.ZIP_DEFLATED) if type(archive) is str else archive    
    if type(name) is zipfile.ZipInfo:
        arch_filename_str = name.filename
    elif type(name) is str:
        arch_filename_str = name
    else:
        arch_filename_str = os.path.basename(file_str)    
    # writing a file that already exists does NOT OVERWRITE!
    # the file is archived two times and overwrites itself when extracting
    if arch_filename_str in archive.namelist():
        logging.warning('file already in archive (skipping): '+arch_filename_str)
        return False    
    if file_is_path:
        archive.write(file_str, arcname=arch_filename_str)
    else:
        archive.writestr(name, file_str)        
    return True

def create_dirs(dir_list):
    """Create all directories listed as path string in 'dir_list'."""
    for dir_ in dir_list:
        if not os.path.exists(dir_):
            logging.info('INIT: create directory: '+dir_)
            os.makedirs(dir_)

def datetime_from_zipinfo(z_info):
    """Return date/time string from a zipinfo file object."""
    #(year, month, day, h, m, s) = z_info.date_time if z_info else (2000, 1, 1, 0, 0, 0)
    #date_time = datetime.datetime(year, month, day, h, m, s)    
    return datetime.datetime(*(z_info.date_time)).strftime('%Y-%m-%d %H:%M:%S')

def md5sum(str_):
    """Calculate the md5sum in hex of a string"""
    sanitized_str = ''.join([c for c in str_ if ord(c) < 128])             
    md5sum_str = hashlib.md5(sanitized_str).hexdigest()
    return md5sum_str

def group_it(l, n):
    """Yield successive 'n'-sized chunks from sequence 'l'.
       The last chunk contains len('l') modulo 'n' elements.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]

def load_object(path_str):
    """Load an object given its path as string"""
    mod_str, _, obj_str = path_str.partition('.')
    try:
        mod = importlib.import_module(mod_str)
    except ImportError as e:
        raise ImportError("Error loading module '%s': %s" % (path_str, e))
    try:
        obj = getattr(mod, obj_str)
    except AttributeError:
        raise NameError("Error loading object '%s' doesn't define any object named '%s'" % (obj_str, mod_str))
    return obj

def remove_dups(some_list, comp_item_index=None):
    """Remove duplicated items in list, preserves order.
       If 'comp_item' is given 'some_list' is treated as list of sequences, the duplicates
       will be found by comparing the items at 'comp_item_index' position in the sequences.
    """
    seen = set()
    seen_add = seen.add
    if comp_item_index is None:
        filt_list = [x for x in some_list if x not in seen and not seen_add(x)]
    else:
        filt_list = [t for t in some_list if t[comp_item_index] not in seen and not seen_add(t[comp_item_index])]
    print('removed '+str(len(some_list)-len(filt_list))+' duplicates')
    return filt_list

def is_executable(exe_name_str):
    """Return true if executable 'exe_name_str' can be executed, else false."""
    devnull = open(os.devnull, 'w')
    try:
        popen = subprocess.Popen(exe_name_str, stdout=devnull, stderr=devnull)
    except OSError as e:
        logging.warning('"'+exe_name_str+'" is not executable, OSError was: '+str(e))
        return False
    finally:
        devnull.close()
    popen.kill()
    return True


# ############
# OLD STUFF...
# ############

def removeFile(file_str):
    if not os.path.exists(file_str):
        print('warning, file does not exists (can not remove): '+file_str)
        return False
    os.remove(file_str)
    return True

""" DEPRECATED
def group(iterable, size):
    sourceiter = iter(iterable)
    while True:
        batchiter = islice(sourceiter, size)
        yield chain([batchiter.next()], batchiter)
"""

def chunkIt(seq, num):
    """Return a list containing the elements in seq divided into num evenly sized chunks."""
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

def rand_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))
