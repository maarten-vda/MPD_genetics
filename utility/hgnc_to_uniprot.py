import httplib2 as http
import json
import argparse

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

def fetch_hgnc_data(hgnc_id):
    headers = {
        'Accept': 'application/json',
    }

    uri = 'https://rest.genenames.org'
    path = f'/fetch/hgnc_id/{hgnc_id}'  # Using f-string to insert the hgnc_id

    target = urlparse(uri + path)
    method = 'GET'
    body = ''

    h = http.Http()

    response, content = h.request(
        target.geturl(),
        method,
        body,
        headers)

    if response['status'] == '200':
        # Assume that content is a JSON reply
        data = json.loads(content)
        print(data['response']['docs'][0]['uniprot_ids'][0])
    else:
        print('Error detected: ' + response['status'])

if __name__ == '__main__':
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Fetch data from HGNC by HGNC ID")
    parser.add_argument('-i', '--hgnc_id', required=True, type=str, help="HGNC ID to fetch data for")
    
    # Parse arguments
    args = parser.parse_args()

    # Call the fetch function with the inputted hgnc_id
    fetch_hgnc_data(args.hgnc_id)
