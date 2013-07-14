

from __future__ import print_function
import subprocess

import secrets


def call(cmd):
    '''
    Run a shell command string using check_call.  Also print the command before
    running it.
    '''
    print(cmd)
    return subprocess.check_call(cmd, shell=True)


class Stardog:
    '''
    '''

    def __init__(self, admin_password, database, bin_dir):
        '''
        admin_password: the password for the admin user
        database: The database this Stardog instance is associated with.
        bin_dir: location of 'stardog' and 'stardog-admin' commands, e.g.
        /usr/local/stardog-1.2.3/bin.
        '''
        self.admin_user = 'admin'
        self.admin_password = admin_password
        self.host = 'localhost'
        self.port = '5822' # default sparql port
        self.sparql_user = 'anonymous'
        self.sparql_password = 'anonymous'
        self.db = database

    def sparql_user_password(self):
        '''
        Return a username and password for basic http authentication to the
        stardog sparql endpoint.
        '''
        return (self.sparql_user, self.sparql_password)

    def sparql_endpoint(self):
        '''
        Endpoint for SPARQL query requests.  Typically requires http
        authentication.
        '''
        return 'http://{host}:{port}/{db}/query'.format(
            host=self.host, port=self.port, db=self.db)
        # The Stardog SPARQL query endpoint is http://<host>[:<port]/<db>/query.

    def sparul_endpoint(self):
        '''
        Endpoint for SPARQL Update requests.
        http://www.w3.org/TR/2013/REC-sparql11-protocol-20130321/#update-operation
        '''
        msg = 'Stardog, as of version 1.2.3 does not support SPARQL Update.'
        raise NotImplementedError(msg)
        # return self.sparql_endpoint()

    def drop_database(self):
        cmd = 'time stardog-admin db drop --username {} --passwd {}'.format(
            secrets.admin_user, secrets.admin_password)
        cmd += ' ' + self.db
        call(cmd)

    def load_graph(self, graph, filename, mediatype):
        '''
        '''
        raise NotImplementedError()


