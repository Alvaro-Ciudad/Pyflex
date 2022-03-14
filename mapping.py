#!/usr/bin/env python

import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

def mapping_to_pdb(uniprot_id):
    """This function mapps the databases identifiers in oder to obtain PDB ID from Uniprot ID."""

    params = {
    'from': 'ID',
    'to': 'PDB_ID',
    'format': 'tab',
    'query': ''
    }

    params['query'] = uniprot_id
    
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    print(response.decode('utf-8'))

id = mapping_to_pdb("P40926")
print(id)