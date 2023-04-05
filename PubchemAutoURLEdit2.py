import requests, pandas as pd, pubchempy as pcp, time, json
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm

# Input search structure
search_query = 'CO(H)'

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
    properties = []
    for cid in tqdm(cid_elements, desc="Retrieving properties"):
        try:
            # Try to retrieve the InChI for the current CID
            inchi = pcp.Compound.from_cid(cid).inchi
            # If successful, add the CID and InChI to the properties list
            properties.append({'CID': cid, 'InChI': inchi})
            # Wait for 0.2 seconds before continuing (5 requests per second)
            time.sleep(0.1)
        except ValueError:
            # If there's an error, skip the current CID and move on to the next one
            pass

    # Create a DataFrame from the properties list
    df = pd.DataFrame(properties)

    # Filter out unwanted compounds based on their InChI
    df = df[~df['InChI'].str.contains('/i|\*|\.\d+|,')]

    # Convert InChI to SMILES and save the filtered dataframe to a new CSV file
    df['SMILES'] = df['InChI'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromInchi(x)))
    df.to_csv(f'{search_query}_Library_filtered.csv', index=False)
    
    # Read the CSV file and filter compounds with MW > 500
    df_filtered = pd.read_csv(f'{search_query}_Library_filtered.csv')
    df_filtered['MW'] = df_filtered['SMILES'].apply(lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x)))
    df_filtered = df_filtered[df_filtered['MW'] <= 500]
    df_filtered.to_csv(f'{search_query}_Library_filtered_MW500.csv', index=False)

