import requests
import os
from stat import S_ISDIR
import traceback #  for error handling'
import re
from Bio.SeqUtils import seq1
import pandas as pd
import paramiko
from Bio import Align
import json


def get_pdb_ids(gene_symbol, entrez_gene_id, synonyms):
    api_url = f"http://mygene.info/v3/query?q={gene_symbol}&species=human&fields=entrezgene,ensembl,pdb"
    response = requests.get(api_url)
    # check response
    if response.status_code != 200:
        return response.status_code, None
    response = response.json()
    for hit in response['hits']:
        if ( str(hit['_id']) == str(int(entrez_gene_id)) ) or ( synonyms and 'ensembl' in hit and hit['ensembl']['gene'] in synonyms ): # (sometimes entrez_gene_id might not match but ensembl gene does
            if 'pdb' in hit:
                return 200, hit['pdb'] # status_code, pdb_ids
            else: # no structures found 
                continue
    return 200, None


def dwl_pdb_file(pdb_id):
    base_url = 'https://files.rcsb.org/download/'
    pdb_url = base_url + pdb_id + '.pdb'
    response = requests.get(pdb_url)
    # Check if the request was successful
    if response.status_code == 200:
        # Save the PDB file to pdb_files folder
        with open(f'../pdb_files/{pdb_id}.pdb', 'wb') as f:
            f.write(response.content)
            return True
    else:
        print(f"Failed to download PDB file {pdb_id}. Status code:", response.status_code)
        return False

def get_resolution(pdb_id):
    try:
        with open(f'../pdb_files/{pdb_id}.pdb', 'r') as f:
            pdb_text = f.read()

            # find line REMARK   2 RESOLUTION.
            resolution = None
            for line in pdb_text.split('\n'):
                if line.startswith('REMARK   2 RESOLUTION.'):
                    resolution = line.split("REMARK   2 RESOLUTION.")[1].strip()
                    break
        if "ANGSTROM" in resolution:
            resolution = resolution.replace(" ANGSTROMS.", "").strip()
        else:
            # print(f"Resolution not found for {pdb_id}")
            resolution = None
    except Exception as e:
        print(f"Failed to get resolution for {pdb_id} with error {e}")
        resolution = None
    return resolution


def get_dbref_data(pdb_id):
    with open(f'../pdb_files/{pdb_id}.pdb', 'r') as f:
        pdb_text = f.read()
        data = []
        for line in pdb_text.split('\n'):
            # contents of a DBREF line:  https://www.wwpdb.org/documentation/file-format-content/format33/sect3.html
            # only reading relevant contents    
            if line.startswith('DBREF') and line[26:32].strip()=="UNP": # if database is UNP 
                data.append({
                    'chain': line[12], # Chain identifier
                    'uniprot': line[33:41].strip(), # uniprot accession code
                    'start': int(line[14:18]), # Initial sequence number of the PDB sequence segment
                    'end': int(line[20:24]), # Ending sequence number of the PDB sequence segment
                    'uniprot_gene': line[42:54].strip() # uniprot gene name 
                })  
    return data


def calculate_sbna_and_download(client, sftp_client, pdb_id, chains=None):

    # pdb_chains is a dictionary with keys as pdb_id and values as chains
    # if chains is empty, all chains will be processed

    def download_directory(sftp, remote_dir, local_dir):
        for item in sftp.listdir_attr(remote_dir):
            remote_path = remote_dir + '/' + item.filename
            local_path = os.path.join(local_dir, item.filename)
            if S_ISDIR(item.st_mode):
                os.makedirs(local_path, exist_ok=True)
                download_directory(sftp, remote_path, local_path)
            else:
                sftp.get(remote_path, local_path)

    REMOTE_DIR = 'ym65_scratch/yliy0004/NetworkAnalysis'
    LOCAL_DIR = '../sbna_results'

    # get available chains, if empty then all chains will be processed
    if not chains:
        try:
            chains = get_dbref_data(pdb_id)
            chains = [chain['chain'] for chain in chains]
        except Exception as e:
            print(f"Failed to get chains for {pdb_id} with error {e}")
            return

    sftp_client.mkdir(f'{REMOTE_DIR}/{pdb_id}')

    for chain in chains: 
        local_dir = f"{LOCAL_DIR}/{pdb_id}/{chain}"
        if os.path.exists(f"{local_dir}/{pdb_id}_monomer/FinalSum"): # skip if analysis already done and downloaded
            continue
        try:
            print(f"\nProcessing {pdb_id} with chain {chain}...")

            # remove the directory if already exists and create new one
            client.exec_command(f"rm -r {REMOTE_DIR}/{pdb_id}/{pdb_id}_monomer")
            client.exec_command(f"rm -r {REMOTE_DIR}/{pdb_id}/{pdb_id}_multimer") 
            
            # upload pdb file to the remote directory
            sftp_client.put(f"../pdb_files/{pdb_id}.pdb", f'{REMOTE_DIR}/{pdb_id}/{pdb_id}.pdb') 
            
            # copy the sbna.sh file to the subdirectory
            client.exec_command(f"cp {REMOTE_DIR}/sbna.sh {REMOTE_DIR}/{pdb_id}/sbna.sh") 

            # run the sbna.sh script
            _,stdout,stderr=client.exec_command(f'cd {REMOTE_DIR}/{pdb_id}; sh sbna.sh {pdb_id}.pdb "Chain {chain}"') 
            print(stdout.read().decode())
            
            # check if Final_sum file exists
            print(stderr.read().decode('utf-8'))
            _,stdout,stderr=client.exec_command(f"ls {REMOTE_DIR}/{pdb_id}/{pdb_id}_monomer/FinalSum")
            # print error message if no output

            # download the sbna results
            if not stdout.read().decode():
                print(f"Final_sum file not found for {pdb_id}")
                # client.exec_command(f"rm -rf {base_path}/{pdb_id}") 
                continue
            
            # download the pdb_id directory
            os.makedirs(local_dir, exist_ok=True) # create the directory if it doesn't exist
            download_directory(sftp_client, f"{REMOTE_DIR}/{pdb_id}", local_dir)
            print(f"SBNA results for {pdb_id} succesfully downloaded")
        
        except Exception:
            # print full error message
            print(f"Failed to process {pdb_id} with error:")
            print(traceback.print_exc())

            # remove directory
            client.exec_command(f"rm -r {REMOTE_DIR}/{pdb_id}/{pdb_id}_monomer")
            client.exec_command(f"rm -r {REMOTE_DIR}/{pdb_id}/{pdb_id}_multimer")
            continue
    
    client.exec_command(f"rm -rf {REMOTE_DIR}/{pdb_id}") # remove the pdb_id remote directory
        
    return
    
def get_uniprot_id(pdb_id, chain):
    # read pdb file 
    with open(f'../pdb_files/{pdb_id}.pdb', 'r') as f:
        pdb_text = f.read()
        for line in pdb_text.split('\n'):
            vals = line.split()
            if line.startswith('DBREF') and chain in vals[2] and 'UNP' in vals:
                return vals[6] # uniprot id is the 6th value


def get_uniprot_sequence(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.ok:
        # Extract the sequence from the FASTA format
        lines = response.text.split('\n')
        sequence = ''.join(lines[1:])  # Skip the header line
        return sequence
    else:
        print("Failed to fetch sequence")
        return None


def get_gene_name(uniprot_id):
    response = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta")
    response = response.text
    match = re.search(r'GN=([^\s]+)', response)
    if match:
        return match.group(1)
    else:
        return None
    
    
def read_scoresEnergetics(pdb_id, chain, monomer=True, z_score=False):
    z = "Z" if z_score else ""
    n_mer = "monomer" if monomer else "multimer"
    file_path = f"../sbna_results/{pdb_id}/{chain}/{pdb_id}_{n_mer}/{pdb_id}_{n_mer}_nowaters_ScoresEnergetics{z}"
    with open(file_path) as f:
        lines = f.readlines()
        # first line is header
        header = lines[0].strip().split('\t')
        header.insert(0, 'res_num_chain') # add residue number to header as first element
        data = []
        for line in lines[1:]:
            data.append(line.strip().split('\t'))
        data = pd.DataFrame(data, columns=header)

    # convert 3 letter code to 1 letter code and get residue number
    data['num'] = data['res_num_chain'].apply(lambda x: ''.join(filter(str.isnumeric, x))).astype(int)
    # letters before the first number (last letter is usually the chain identifier)
    data['res'] = data['res_num_chain'].apply(lambda x: ''.join(filter(str.isalpha, x[:-1])))
    data['res_code'] = data['res'].apply(seq1)
    return data

def read_finalsum(pdb_id, chain):
    file_path = f'../sbna_results/{pdb_id}/{chain}/{pdb_id}_monomer/FinalSum'
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            # convert final_sum to dataframe
            final_sum = f.readlines()
            final_sum = [i.split() for i in final_sum]
            rows = []
            for row in final_sum:
                # get letters only from row[0]
                aa_3 = ''.join(filter(str.isalpha, row[0]))
                if row[1] == 'NA' or len(aa_3) != 3: 
                    # skip if residue number is NA or residue code is not 3 letters
                    continue
                
                aa_1 = seq1(aa_3) # convert to 1 letter code
                res_num = ''.join(filter(str.isnumeric, row[0]))
                score = float(row[1])
                # print(f"{aa_3} {aa_1} {res_num} {score}")
                rows.append([aa_3, aa_1, int(res_num), score])
            final_sum = pd.DataFrame(rows, columns=['res', 'res_code', 'num', 'network_score'])
        return final_sum
    else:
        print(f"FinalSum file not found for {pdb_id}")
        return None
    

def read_finalsum_decomp(pdb_id, chain):
    file_path = f"../sbna_results/{pdb_id}/{chain}/{pdb_id}_monomer/FinalSum_decomp"
    if os.path.exists(file_path):
        with open(file_path) as f:
            lines = f.readlines()
            # first line is header
            header = lines[0].strip().split('\t')
            data = []
            for line in lines[1:]:
                data.append(line.strip().split('\t'))
            data = pd.DataFrame(data, columns=header)
        # convert 3 letter code to 1 letter code and get residue number
        data['num'] = data['Acid'].apply(lambda x: ''.join(filter(str.isnumeric, x))).astype(int)
        data['res'] = data['Acid'].apply(lambda x: ''.join(filter(str.isalpha, x)))

        # convert SecondOrderIntermodularDegree_AVERAGE	NodeEdgeBetweennessSTRIDE_sidechain_MAX	LigandMULTIMERCENTROIDSC_MIN to float, if nan then set to 
        # for col in data.columns[2:]:
        data['res_code'] = data['res'].apply(seq1)

        # a quick check if this is correct is res_code shouldnt be empty
        data = data[data['res_code'] != '']
        
        # convert SecondOrderIntermodularDegree_AVERAGE	NodeEdgeBetweennessSTRIDE_sidechain_MAX	LigandMULTIMERCENTROIDSC_MIN to float
        for col in ['SecondOrderIntermodularDegree_AVERAGE', 'NodeEdgeBetweennessSTRIDE_sidechain_MAX', 'LigandMULTIMERCENTROIDSC_MIN']:
            data[col] = data[col].apply(lambda x: x if x != 'NA' else None).astype(float)
            
        return data
    else:
        print(f"FinalSum_decomp file not found for {pdb_id}")
        return None
    

def read_all_scores_and_average(pdb_id, chain):
    # read final scores and decomps and average across chains
    final_sum = read_finalsum(pdb_id, chain)
    final_sum_decomp = read_finalsum_decomp(pdb_id, chain)

    # merge the two dataframes
    final_sum = final_sum.merge(final_sum_decomp, on='num', suffixes=('', '_decomp'))
    return final_sum

def map_sequ_sbna_pdb(sbna, pdb):
    # maps the sbna sequence to the pdb sequence
    sbna_to_pdb_mapping = []
    # i starts at where the first letter is found (not -) in pdb_aligned
    i = 0
    for (s, p) in zip(sbna, pdb):
        if s == p:
            sbna_to_pdb_mapping.append(i)
            i += 1
        elif s == '-':
            i += 1
        else:
            sbna_to_pdb_mapping.append("?")
    assert(len(sbna_to_pdb_mapping) == len("".join(aa for aa in sbna if aa != '-')))  # verify that the length is correct
    return sbna_to_pdb_mapping


def convert_auth_to_pdb(pdb_id):
    # get the author chain id to pdb chain id mapping
    response = requests.get(f'''
    https://data.rcsb.org/graphql?query={{entry(entry_id:"{pdb_id}")
    {{polymer_entities {{
            polymer_entity_instances {{
                rcsb_polymer_entity_instance_container_identifiers {{
                auth_asym_id
                asym_id
                entry_id
                entity_id
                }}
            }}
            }}
        }}
    }}
    ''')
    # check response
    if response.status_code != 200:
        print("Failed to fetch data from RCSB: convert auth_asym_id to asym_id")
        return None
    response = response.json()

    auth_pdb_map = {}
    for polymer in response['data']['entry']['polymer_entities']:
        # represent one "macromolecule" section in PDB website
        for chain_instance in polymer['polymer_entity_instances']:
            # represents one chain
            # map author chain id to rcsb pdb chain id
            auth_pdb_map[chain_instance['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id']] = \
                chain_instance['rcsb_polymer_entity_instance_container_identifiers']['asym_id']
    return auth_pdb_map

def get_alignment_regions(uniprot_id, pdb_id, chain, auth_pdb_map):
    """ returns the aligned sequences like e.g. 
    [{'query_begin': 696, 'query_end': 772, 'target_begin': 2, 'target_end': 78}, {'query_begin': 773, 'query_end': 1022, 'target_begin': 82, 'target_end': 331}]
    and the PDB sequence"""

    def tmp(unp):
        return requests.get(f"""https://1d-coordinates.rcsb.org/graphql?query={{
        alignment(
            from:UNIPROT,
            to:PDB_INSTANCE,
            queryId:"{unp.split('-')[0]}",
        ){{
            query_sequence
            target_alignment {{
            target_id
            target_sequence
            coverage{{
                query_coverage
                query_length
                target_coverage
                target_length
            }}
            aligned_regions {{
                query_begin
                query_end
                target_begin
                target_end
            }}
            }}
        }}
        }}""")

    # maps UniProt - PDBs (one to many)
    response = tmp(uniprot_id)
    if response.status_code != 200 or response.json()['data']['alignment']['target_alignment'] is None:
        # print("Failed to fetch data from RCSB: get alignment regions")

        # try to get uniprot id another way
        response_tmp = requests.get(f"""https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}""")
        response_tmp = response_tmp.json()
        found = None
        for uniprot_id, data in response_tmp[pdb_id.lower()]['UniProt'].items():
            for mapping in data["mappings"]:
                if mapping["struct_asym_id"] == chain:
                    response = tmp(uniprot_id)
                    found = True
                    break
            if found:
                break

    response = response.json()
    for alignment in response['data']['alignment']['target_alignment']:
        # search for correct pdb and chain
        if alignment['target_id'] == f"{pdb_id}.{auth_pdb_map[chain]}": 
            target_sequence = alignment['target_sequence']
            aligned_regions = alignment['aligned_regions']
            return aligned_regions, target_sequence
    return None, None
    

def map_sequ_sbna_pdb(sbna, pdb):
    # maps the sbna sequence to the pdb sequence 
    sbna_to_pdb_mapping = []
    # i starts at where the first letter is found (not -) in pdb_aligned
    i = 1
    for (s, p) in zip(sbna, pdb):
        if s == p:
            sbna_to_pdb_mapping.append(i)
            i += 1
        elif s == '-':
            i += 1
        else:
            sbna_to_pdb_mapping.append("?")

    assert(len(sbna_to_pdb_mapping) == len("".join(aa for aa in sbna if aa != '-')))  # verify that the length is correct
    return sbna_to_pdb_mapping  

def align_finalsum_with_uniprot(final_sum, pdb_id, chain):
    """
    final_sum is the output from SBNA, where each row is a residue in the chain (although not always correct)
    pdb_id is the pdb id
    chain is the chain id

    This function aligns the SBNA sequence with the PDB sequence, and then maps the PDB sequence to the UniProt sequence
    """
    data = get_dbref_data(pdb_id) # obtain the DBREF lines details from the pdb file
    ranges = [range(i['start'], i['end']+1) for i in data if i['chain'] == chain] # get the start and end sequence number in the chain
    r = [item for sublist in ranges for item in sublist] # flatten the list. r is the list of all sequence numbers in the chain
    final_sum = final_sum[final_sum['num'].isin(r)].reset_index(drop=True)

    # get uniprot sequence (use a cache to avoid multiple requests for the same uniprot id)
    uniprot_id = get_uniprot_id(pdb_id, chain)
    with open("../cache/uniprot_seqs.json", "r") as f:
        uniprot_seqs = json.load(f) # saved uniprot sequences
    if uniprot_id in uniprot_seqs:
        uniprot_seq = uniprot_seqs[uniprot_id]
    else:
        uniprot_seq = get_uniprot_sequence(uniprot_id)
        uniprot_seqs[uniprot_id] = uniprot_seq
        with open("../cache/uniprot_seqs.json", "w") as f:
            json.dump(uniprot_seqs, f) # save uniprot sequences
    
    # convert author chain id to pdb chain id
    auth_pdb_map = convert_auth_to_pdb(pdb_id) # convert author chain id to pdb chain id
    aligned_regions, target_sequence = get_alignment_regions(uniprot_id, pdb_id, chain, auth_pdb_map)

    # sometimes SBNA has multiple residues for the same number (due to other chains/polymer entities)
    # where there are multiple residues for the same number, make final_sum use the n'th value
    # method is not perfect, but works for most cases
    if any(final_sum.groupby('num')['res_code'].count() > 1):
        print(f"Multiple residues for the same number detected for {pdb_id}{chain}")
        try:
            # non-duplicated rows
            non_dups = final_sum.drop_duplicates(subset='num', keep=False)

            # see how many chains there are for each number
            dbrefs = get_dbref_data(pdb_id)

            # duplicated rows
            dups = final_sum[final_sum.duplicated('num', keep=False)].reset_index(drop=True)

            # unique chains in dbrefs (dbrefs is a dict)
            unique_chains = []
            for chain_entity in dbrefs:
                if chain_entity['chain'] not in unique_chains:
                    unique_chains.append(chain_entity['chain'])
            i = unique_chains.index(chain)

            # get the i'th value for duplicated rows
            # Group by 'num' and create a counter for each group
            dups['counter'] = dups.groupby('num').cumcount()
            # Filter the DataFrame to get the ith occurrence (e.g., n=3)
            nth_items = dups[dups['counter'] == i]
            nth_items = nth_items.drop(columns=['counter'])

            # combine non-duplicated and duplicated rows
            final_sum = pd.concat([non_dups, nth_items])

            # sort by num
            final_sum = final_sum.sort_values(by='num')
        except:
            print("Error in handling multiple residues for the same number, using all values.")
            pass

    sbna_seq = ''.join(final_sum['res_code'].values)

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(sbna_seq, target_sequence)

    sbna_aligned, pdb_aligned  = alignments[0][0], alignments[0][1]

    mapped = map_sequ_sbna_pdb(sbna_aligned, pdb_aligned)

    # align pdb sequence to uniprot seqeunce using aligned regions
    pdb_to_uniprot_mapping = {}
    for region in aligned_regions:
        pdb_begin = region['target_begin']
        pdb_end = region['target_end']
        uniprot_begin = region['query_begin']
        uniprot_end = region['query_end']

        # map numbers from pdb to uniprot
        for i in range(pdb_begin, pdb_end+1):
            pdb_to_uniprot_mapping[i] = i - pdb_begin + uniprot_begin

    # map sbna sequence to uniprot sequence
    sbna_to_uniprot_mapping = []
    for i in mapped:
        if i != "?":
            try:
                sbna_to_uniprot_mapping.append(pdb_to_uniprot_mapping[i])
            except:
                sbna_to_uniprot_mapping.append("?")
        else:
            sbna_to_uniprot_mapping.append("?")
    final_sum['uniprot_num'] = sbna_to_uniprot_mapping
    final_sum['uniprot_res'] = [uniprot_seq[i-1] if i != "?" else "?" for i in sbna_to_uniprot_mapping]
        
    return final_sum