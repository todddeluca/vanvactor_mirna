

import json

import requests


class Sparql:
    '''
    Like SPARQLWrapper, instances of this class issue requests to SPARQL
    endpoints and return the results.
    '''
    def __init__(self, endpoint):
        self.endpoint = endpoint

    def query(self, qry, accept=None):
        '''
        Send an HTTP GET request to the SPARQL endpoint with the given SPARQL
        Query string. Return the response data.  JSON data will be
        parsed/decoded before being returned.

        qry: a SPARQL query (SELECT, ASK, CONSTRUCT, ...) string.
        accept: A media type for an accept header. The default, None, means
        that no accept header will be sent to the endpoint.
        For SELECT queries, the media types are typically one of:

            JSON: application/sparql-results+json,
            XML: application/sparql-results+xml,
            CSV: text/csv; charset=utf-8,
            TSV: text/tab-separated-values; charset=utf-8,

        For CONSTRUCT queries (which return RDF data), the media types are
        typically one of:

            RDF/XML: application/rdf+xml
            N-Triples: text/plain
            Turtle: application/x-turtle
            N3: text/rdf+n3
            TriX: application/trix
            TriG: application/x-trig
            N-Quads: text/x-nquads or application/n-quads
            JSON-LD: application/ld+json
        '''
        url = self.endpoint
        params = {'query': qry}
        if accept:
            headers = {'accept': accept}
        else:
            headers = {}
        r = requests.get(url, params=params, headers=headers)
        # , auth=('anonymous', 'anonymous'))
        r.raise_for_status()
        # Parse JSON for json response content types
        resp_ct = str(r.headers.get('content-type'))
        if (resp_ct.startswith('application/sparql-results+json') or
            resp_ct.startswith('application/ld+json')):
            return r.json()
        else:
            return r.text

    def ask(self, qry):
        '''
        Send a SPARQL ASK query to the endpoint and return the answer as a
        boolean value.

        qry: a SPARQL ASK string.
        '''
        # Sesame does not seem to like outfmt 'json' or 'tsv' for responding
        # to ASK queries.  If no 'accept' header is given, Sesame responds
        # with a plain text 'true' or 'false', which is json-compatible.
        return json.loads(self.query(qry))

    def update(self, qry):
        '''
        POST a SPARQL UPDATE request to the endpoint.  Raise an exception if
        the request fails.

        qry: a SPARQL UPDATE string.
        '''
        url = self.endpoint
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        data = {'update': qry}
        # Do I need to add "auth=('anonymous', 'anonymous')" to the post() call?
        r = requests.post(url, headers=headers, data=data)
        r.raise_for_status()
        return r.status_code

    def drop_graph(self, graph, silent=False):
        '''
        Drop the named graph from the graph store.
        http://www.w3.org/TR/2013/REC-sparql11-update-20130321/#drop

        graph: a URI like 'http://example.com/mygraph'
        silent: If True, the server will always return success, even if the
        named graph does not exist.
        '''
        qry = 'DROP {silent} GRAPH <{graph}>'.format(
            silent='SILENT' if silent else '',
            graph=graph)
        return self.update(qry)

    def exists_graph(self, graph):
        '''
        Return True iff the graph exists in the graph store.  In some
        implementations, e.g. OpenRDF Sesame 2.6.10, this returns true iff the
        graph exists AND is non-empty.

        http://answers.semanticweb.com/questions/1745/ask-whether-a-graph-exists-or-not
        The empty group pattern should match any graph, including the empty
        graph.  However in Sesame, it always returns true, even if the graph
        does not exist in the repository (at least for the OWLIM-Lite SAIL
        and Sesame In-Memory Store I tested.)
        Therefore, the best we can do is test if the graph is empty or not.
        '''
        qry = 'ASK {{ GRAPH <{graph}> {{ }} }}'.format(graph=graph)
        return self.ask(qry)


