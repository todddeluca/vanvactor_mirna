
import os
import subprocess

'''
Common functions for downloading and uncompressing files.
'''

def call(cmd):
    '''
    Run a shell command string using check_call.  Also print the command before
    running it.
    '''
    print cmd
    return subprocess.check_call(cmd, shell=True)


def makedirs(path, mode=0775):
    '''
    Idempotent function for making a directory.  If path is not a directory,
    make path and all of its non-existent parent directories.  Either way,
    return path.
    '''
    if not os.path.exists(path):
        os.makedirs(path, mode)
    return path


def download(url, dest, mode=0775):
    '''
    dest: the filename to download the url to.  Will create any parent
    directories of dest that do not already exist.
    '''
    makedirs(os.path.dirname(dest), mode)
    call('curl -o {} {}'.format(dest, url))
    return dest


def gunzip(path):
    if path.endswith('.gz'):
        call('gunzip {}'.format(path))
    return path.rstrip('.gz')


def unzip(path):
    path_dir = os.path.dirname(path)
    if path.endswith('.zip'):
        call('unzip -d {} {}'.format(path_dir, path))
    return path.rstrip('.zip')


def download_and_unzip(url, dest_dir):
    '''
    download a url to a file in dest_dir named with the basename of the url.
    If the file is named *.zip, unzip the file in the directory.  Finally,
    return the filename with any '.zip' suffix returned.
    '''
    dest = os.path.join(dest_dir, os.path.basename(url))
    download(url, dest)
    return unzip(dest)


