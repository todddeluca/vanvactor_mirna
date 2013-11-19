

import os
import subprocess

import sparql


def download_sparql_ntriples(qry, endpoint):
    ''' 
    Download triples from a sparql endpoint and return them as N-triples.

    :param qry: a SPARQL Construct query.
    :type qry: str.
    :param endpoint: a SPARQL endpoint URL.
    :type endpoint: str.
    '''
    print 'Downloading triples from endpoint: {}'.format(endpoint)
    print 'Query:', qry
    sp = sparql.Sparql(endpoint)
    try:
        # Use text/plain for N-triples
        nt = sp.query(qry, accept='text/plain')
    except Exception as e:
        print dir(e)
        print e.response.text
        print e.response.headers
        raise
    return nt


def count_and_write_ntriples(nt, filename):
    '''
    :param nt: N-triples, one per line.
    :type nt: str
    '''
    print 'Writing triples to file:', filename
    count = 0
    with open(filename, 'w') as fh:
        for line in nt.splitlines(True):
            if line.strip():
                count += 1
                fh.write(line)

    print 'Wrote {} triples.'.format(count)
    return count


def download_sparql_construct_rdf(qry, filename, endpoint):
    '''
    Run qry, a SPARQL query, on endpoint, a SPARQL endpoint, and save
    the constructed triples to a file as N-Triples.  Return the count of lines
    (triples) returned from the query.

    qry: a SPARQL Construct query.
    filename: where to save the constructed triples (as N-Triples).
    endpoint: a sparql endpoint url
    '''
    nt = download_sparql_ntriples(qry, endpoint)
    count = count_and_write_ntriples(nt, filename)
    return count


def map_flatten_set_list(func, *args_lists):
    '''
    Run a function, one that returns a sequence, on each element of items.
    Flatten the sequences into a single list, and return a list of the
    unique items.
    Example: 
        >>> map_flatten_set_list(lambda x: (x, x+1), [1, 2, 4])
        [1, 2, 3, 4, 5]

    '''
    result_lists = map(func, *args_lists)
    return list(set(i for list in results_lists for i in lst))


def map_flatten_set_list_old(func, items):
    '''
    Run a function, one that returns a sequence, on each element of items.
    Flatten the sequences into a single list, and return a list of the
    unique items.
    Example: 
        >>> map_flatten_set_list(lambda x: (x, x+1), [1, 2, 4])
        [1, 2, 3, 4, 5]
    '''
    return list(set(i for lst in [func(i) for i in items] for i in lst))


def example(message):

    print 'Example CLI function.'
    print 'message:', message


def call(cmd):
    '''
    Run a shell command string using check_call.  Also print the command before
    running it.
    '''
    print cmd
    return subprocess.check_call(cmd, shell=True)


def makefiledirs(path, mode=0775):
    '''
    Idempotent function for making any non-existent directories for the
    file path `path`.  Returns path.
    '''
    makedirs(os.path.dirname(path), mode=mode)
    return path


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


def download_and_unzip(url, dest_dir):
    '''
    download a url to a file in dest_dir named with the basename of the url.
    If the file is named *.zip, unzip the file in the directory.  Finally,
    return the filename with any '.zip' suffix returned.
    '''
    dest = os.path.join(dest_dir, os.path.basename(url))
    download(url, dest)
    if dest.endswith('.zip'):
        call('unzip -d {} {}'.format(dest_dir, dest))
    return dest.rstrip('.zip')


