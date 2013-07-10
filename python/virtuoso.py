
from __future__ import print_function
import os
import shutil
import subprocess
import time


def call(cmd):
    '''
    Run a shell command string using check_call.  Also print the command before
    running it.
    '''
    print(cmd)
    return subprocess.check_call(cmd, shell=True)


##########
# VIRTUOSO


class Virtuoso:
    '''
    Class implementing shared functionality between virtuoso 6 and virtuoso 7.
    '''

    def __init__(self, password, load_dir=None):
        '''
        password: the password for the virtuoso dba user
        load_dir: a directory assigned to 'DirsAllowed' in 'virtuoso.ini', used
          to copy files to when bulk loading data into virtuoso.
        '''
        self.user = 'dba'
        self.password = password
        self.host = 'localhost'
        self.port = '8890' # default virtuoso sparql port
        self.load_dir = load_dir

    def sparql_endpoint(self):
        '''
        Endpoint for SPARQL query requests.
        '''
        return 'http://{host}:{port}/sparql'.format(host=self.host,
                                                    port=self.port)

    def sparul_endpoint(self):
        '''
        Endpoint for SPARQL Update requests.
        http://www.w3.org/TR/2013/REC-sparql11-protocol-20130321/#update-operation
        '''
        return self.sparql_endpoint()

    def load_graph(self, graph, filename):
        '''
        Copy files to the virtuoso load dir, bulk load them into the named
        graph, and delete them from the load dir.

        On an EC2 m3.xlarge 64-bit EBS with moderate network performance, 15GB RAM,
        4 vCPUs and 13 ECUs, I loaded ~500k triples/minute load (i.e. ~8k triples/sec)

        The page http://virtuoso.openlinksw.com/dataspace/doc/dav/wiki/Main/VOSArticleLUBMBenchmark
        mentions loading about 1 gigatriple (8000 universities) at about 29-36
        Kt/s.

        filename: a file path, like '/path/to/foo.rdf' or './bar.rdf'.
        graph: a named graph URI, like 'http://example.com/mygraph'.  The
        triples in the file will be loaded into this graph.
        '''
        assert os.path.exists(self.load_dir)

        dest_fn = os.path.join(self.load_dir, os.path.basename(filename))
        global_fn = os.path.join(self.load_dir, 'global.graph')

        # clear out the load dir
        print('load_graph: clearing load dir')
        call('rm -rf {} {}'.format(dest_fn, global_fn))

        # copy file to the load dir
        print('load_graph: copying {} to {}'.format(filename, self.load_dir))
        shutil.copy(filename, dest_fn)

        # create global.graph file
        print('load_graph: creating global.graph file')
        with open(global_fn, 'w') as fh:
            fh.write('{}\n'.format(graph))

        # register files to be loaded
        print('load_graph: registering {} for loading'.format(dest_fn))
        cmd = "ld_dir('{dir}', '{file}', '{graph}')".format(
            dir=self.load_dir,
            file=os.path.basename(filename),
            graph=graph)
        self.run_isql_cmd(cmd)

        # load registered files
        print('load_graph: loading all files.')
        print('Use "isql-vt localhost:1111 dba <password> \'EXEC=select * from DB.DBA.load_list\'" to track progress.')
        self.run_isql_cmd('rdf_loader_run()')

        # clean up the files
        os.remove(dest_fn)
        os.remove(global_fn)

        # Remove file load records for the graph (so it can be reloaded)
        cmd = "delete from DB.DBA.load_list where ll_graph = '{graph}'".format(
            graph=graph)
        self.run_isql_cmd(cmd)

    def run_isql_cmd(self, cmd, password=None):
        '''
        Run a virtuoso isql command from the command line, by setting EXEC={cmd}.
        See http://docs.openlinksw.com/virtuoso/isql.html for details.

        cmd: a string used as a value for EXEC, an isql function to run.  For
        example, 'status()' or 'ld_dir("/home/ubuntu/tmp", "my.rdf", "my:graph:uri")'
        '''
        if password is None:
            password = self.password

        cmd = [self.isql_exe(), 'localhost:1111', self.user, password,
               'EXEC={cmd}'.format(cmd=cmd)]
        print(cmd)
        subprocess.check_call(cmd)

    def clear_virtuoso_db(self):
        '''
        Stop the virtuoso server, delete all data in the database (!!!) and
        restart the virtuoso server.
        '''
        self.stop()
        call('sudo rm -rf {}'.format(self.db_file()))
        self.start()

        # since deleting the database also deletes the dba password
        # reset the dba password using the default isql password
        # http://docs.openlinksw.com/virtuoso/newadminui.html
        default_password = 'dba'
        cmd = 'set password "{old}" "{new}"'.format(old=default_password, 
                                                    new=self.password)
        self.run_isql_cmd(cmd, password=default_password)

    def grant_select_to_sparql_user(self):
        '''
        After installing virtuoso, one must grant privileges to user SPARQL
        in order to run SPARQL queries via the web interface (or http
        endpoint?)

        How to fix the sparql query failure:
        http://answers.semanticweb.com/questions/10360/virtuoso-federated-query
        After fixing the quotes (they must be plain double quotes), I ran the
        queries to fix the SPARQL user permissions:

            isql-vt localhost:1111 dba <password>
            grant select on "DB.DBA.SPARQL_SINV_2" to "SPARQL";
            grant execute on "DB.DBA.SPARQL_SINV_IMP" to "SPARQL";
        '''
        self.run_isql_cmd('grant select on "DB.DBA.SPARQL_SINV_2" to "SPARQL"')
        self.run_isql_cmd('grant execute on "DB.DBA.SPARQL_SINV_IMP" to "SPARQL"')

    def db_file(self):
        raise NotImplementedError

    def isql_exe(self):
        raise NotImplementedError

    def stop(self):
        raise NotImplementedError

    def start(self):
        raise NotImplementedError


class Virtuoso7(Virtuoso):

    def db_file(self):
        return '/usr/local/virtuoso-opensource-7.0/var/lib/virtuoso/db/virtuoso.db'

    def isql_exe(self):
        return '/usr/local/virtuoso-opensource-7.0/bin/isql'

    def stop(self):
        call('sudo service virtuoso-opensource-7.0 stop')
        # Without a brief pause, if start() immediately follows stop(), start()
        # will fail.
        time.sleep(2)

    def start(self):
        call('sudo service virtuoso-opensource-7.0 start')


class Virtuoso6(Virtuoso):

    def db_file(self):
        return '/var/lib/virtuoso-opensource-6.1/db/virtuoso.db'

    def isql_exe(self):
        return '/usr/bin/isql-vt'

    def stop(self):
        call('sudo service virtuoso-opensource-6.1 stop')
        # Without a brief pause, if start() immediately follows stop(), start()
        # will fail.
        time.sleep(2)

    def start(self):
        call('sudo service virtuoso-opensource-6.1 start')


