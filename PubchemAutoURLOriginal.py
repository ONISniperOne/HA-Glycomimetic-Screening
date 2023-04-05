import requests, pandas as pd, pubchempy as pcp, time, json
from rdkit import Chem
from tqdm import tqdm

# Input search structure
search_query = 'CCO'

# Build the URL with a placeholder for the structure search query
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/smiles/{search_query}/JSON?MaxRecords=50000"
print("Structure search query URL:", url)

# Make the API request
response = requests.get(url)

# Get the Waiting.ListKey element from the JSON response
list_key = json.loads(response.content).get('Waiting', {}).get('ListKey')
print("ListKey:", list_key or "No compounds found matching the structure query.")

if list_key:
    # Build the URL to retrieve the search results
    url2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{list_key}/cids/JSON"

    # Wait for the search to complete
    while 'Waiting' in json.loads(requests.get(url2).content):
        print("Waiting for 10 seconds...")
        time.sleep(10)

    # Retrieve properties for each compound
    cid_elements = json.loads(requests.get(url2).content).get('IdentifierList', {}).get('CID', [])
    properties = [{'CID': cid, 'InChI': pcp.Compound.from_cid(cid).inchi}
                  for cid in tqdm(cid_elements, desc="Retrieving properties") if not time.sleep(0.2) and not ValueError]

    # Convert to pandas dataframe and filter out unwanted compounds
    df = pd.DataFrame(properties)
    df = df[~df['InChI'].str.contains('/i|\*|\.\d+|,')]

    # Convert InChI to SMILES and save the filtered dataframe to a new CSV file
    df['SMILES'] = df['InChI'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromInchi(x)))
    df.to_csv('Library_filtered.csv', index=False)

