"""
# Exam 04

## Important!
For each question, provide:
- a callable function that returns the requested result
- a second function that tests the behaviour of the first function (that answers the question), by running the first function on a simple testcase, and using an assertion or raising an error if the test fails. (N.B.: you need to decide what a simple, reasonable testcase is! I provide a hint for each question)
"""


"""
## Q1: git repository

Write commands below that:
- create a git repository
- add a python file and a README.md file, to be staged before commit
- commit those changes (with a commit message that mention your name and what changes you made!)
- show the commit log


WRITE COMMANDS HERE:
----
# Create a git repository
git init

# Add a python file and a README.md file to the staging area
git add py_exam4.py README.md

# Commit the changes with a descriptive message
git commit -m "Add initial files by Aarti: Python file and README"

# Show the commit log
git log

----

"""

"""
## Q2: center of mass in a volume

You are tasked with writing a function to find the center of mass of a provided 3d volume.

You are provided with code to load a volume, and display projections of a volume,
now use scipy's center of mass function (or a function of your own!) to find the center of mass of the provided volume.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass

def load_volume_data(filename: str = "./1CS4.npz", arrayname: str = "map_1CS4"):
    """
    Load the 3D volume from a .npz file.

    Parameters:
    - filename: The path to the .npz file.
    - arrayname: The key in the .npz file for the 3D volume data.

    Returns:
    - volume: The 3D volume as a numpy array.
    """
    volume = np.load(filename)[arrayname]
    return volume

def plot_volume_projection(volume: np.array):
    """
    Display the center slice and the x-axis projection of the volume.

    Parameters:
    - volume: The 3D volume as a numpy array.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    
    axes[0].imshow(volume[volume.shape[0] // 2], cmap='gray')
    axes[0].set_title("Center Slice")
    
    axes[1].imshow(volume.sum(axis=0), cmap='gray')
    axes[1].set_title("X-axis Projection")
    plt.show()

def compute_volume_COM(volume: np.array):
    """
    Compute the center of mass of the volume.

    Parameters:
    - volume: The 3D volume as a numpy array.

    Returns:
    - com: The coordinates of the center of mass.
    """
    return center_of_mass(volume)

def test_compute_volume_COM():
    """
    Test function to load a volume, compute its center of mass, and visualize it.
    """
    volume = load_volume_data("./1CS4.npz", "map_1CS4")
    com = compute_volume_COM(volume)
    print(f"Center of Mass: {com}")
    plot_volume_projection(volume)

test_compute_volume_COM()


def test_compute_volume_COM_with_controlled_data():
    
    volume = np.ones((3, 3, 3))
    
    volume[0, 0, 0] = 2
    volume[0, 0, 1] = 2
    volume[0, 1, 0] = 2
    volume[1, 0, 0] = 2
    
    expected_com = center_of_mass(volume)
    
    computed_com = compute_volume_COM(volume)
    
    assert np.allclose(computed_com, expected_com), f"Computed COM {computed_com} does not match expected COM {expected_com}"
    
    return "Test passed. Computed center of mass matches the expected value."

test_compute_volume_COM_with_controlled_data()


    


"""
## Q3: finding motifs, distribution

You are provided with a list of motifs, and a clustal file containing aligned protein sequences.

Write a function that, using a pandas dataframe, store the number of times each motif appear in a given sequence, and return that dataframe.

**Only write a test for this function, not the next one**

Write a separate function that computes the average number of occurrences of the motif in the provided file.
"""
import pandas as pd
from Bio import AlignIO

def get_motif_counts_df(aligned_filename: str, motifs: list):
    """
    Reads a clustal file containing aligned protein sequences,
    counts the occurrences of each motif in each sequence, and
    stores these counts in a pandas DataFrame.

    Parameters:
    - aligned_filename: Path to the clustal file.
    - motifs: List of motifs to search for in the sequences.

    Returns:
    - A pandas DataFrame with one row per sequence and one column per motif,
      containing the count of each motif in each sequence.
    """

    alignment = AlignIO.read(aligned_filename, "clustal")
    counts = {motif: [] for motif in motifs}
    
   
    for record in alignment:
        sequence = str(record.seq)
        for motif in motifs:
           
            counts[motif].append(sequence.count(motif))
    
   
    df = pd.DataFrame(counts, index=[record.id for record in alignment])
    return df

def test_get_counts_df():
    """
    Test function for get_motif_counts_df.
    """
   
    aligned_filename = "./hAPP.clustal"
    motifs = ["LL", "EE", "W"]
    
   
    df = get_motif_counts_df(aligned_filename, motifs)
    print(df)


def get_median(df):
    """
    Computes the mean number of occurrences of each motif across all sequences.

    Parameters:
    - df: DataFrame with counts of motifs in sequences.

    Returns:
    - A pandas Series with the mean occurrences of each motif.
    """
    return df.mean()


test_get_counts_df()



class MockSeqRecord:
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq

def test_get_motif_counts_df_corrected():
    
    mock_records = [
        MockSeqRecord("Sequence1", "AACCTGEEELLKKLLWWKK"),
        MockSeqRecord("Sequence2", "AAGGTEEEKKLLMMWWKK"),
        MockSeqRecord("Sequence3", "ACCCTGEEEKKLLWWMMKK"),
    ]

    motifs = ["LL", "EE", "W"]

    expected_counts = {
        "LL": [2, 1, 1],
        "EE": [1, 1, 1],  
        "W": [2, 2, 2]
    }

    counts = {motif: [] for motif in motifs}
    for record in mock_records:
        sequence = record.seq
        for motif in motifs:
            counts[motif].append(sequence.count(motif))
    df = pd.DataFrame(counts, index=[record.id for record in mock_records])

    for motif, count_list in expected_counts.items():
        assert all(df[motif] == count_list), f"Counts for motif {motif} are incorrect: expected {count_list}, got {df[motif].tolist()}"

    print("All motif counts are correct.")
    return df

test_result_df_corrected = test_get_motif_counts_df_corrected()
test_result_df_corrected



"""
## Q4 EMDB data, pandas

You are provided with a SQLite database, `emdb.db`, containing a table `entries` with the emdb id, title, and date of deposition of various EMDB entries.

Write a function that takes the database and a search term as an input and:
- select the data from the database
- discard the entries whose title do not match the search term
- create and return a dataframe with the content from the database, sorted by descending order of deposition (newest first, oldest last)
"""

import pandas as pd
import sqlite3

def find_matching_entries(database, search_term: str):
    """
    Searches for entries in the EMDB database that match a given search term in their title,
    sorts them by deposition date in descending order, and returns them in a pandas DataFrame.
    
    Parameters:
    - database: The path to the SQLite database file or a sqlite3.Connection object.
    - search_term: The search term to match in the titles of database entries.
    
    Returns:
    - A pandas DataFrame containing the filtered and sorted entries.
    """
    if isinstance(database, str):
        conn = sqlite3.connect(database)
    elif isinstance(database, sqlite3.Connection):
        conn = database
    else:
        raise ValueError("database must be a path (str) or a sqlite3.Connection object")
    
    query = "SELECT id, title, date FROM entries WHERE title LIKE ? ORDER BY date DESC"
    df = pd.read_sql_query(query, conn, params=('%' + search_term + '%',))
    
    if isinstance(database, str):
        conn.close()
    
    return df

def test_find_matching_entries():
    """
    Tests the find_matching_entries function against controlled scenarios.
    """
    database_path = "./emdb.db" 

    # Test with a search term expected to return no results
    df_no_results = find_matching_entries(database_path, "ZZZZXYZ")
    assert df_no_results.empty, "Expected no results for 'ZZZZXYZ', but some were found."

    df_all_results = find_matching_entries(database_path, "")
    assert not df_all_results.empty, "Expected some results for an empty search term, but none were found."


    conn = sqlite3.connect(':memory:')
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE entries (id TEXT, title TEXT, date DATE)")
    test_data = [
        ('1', 'First Entry Title - Virus', '2024-03-23'),
        ('2', 'Second Entry Title - Life', '2024-03-22'),
        ('3', 'Third Entry Title - Bacteria', '2024-03-21'),
    ]
    cursor.executemany("INSERT INTO entries VALUES (?, ?, ?)", test_data)
    conn.commit()

   
    df_in_memory = find_matching_entries(conn, "Virus")
    expected_ids_in_memory = ['1']
    assert set(df_in_memory['id'].tolist()) == set(expected_ids_in_memory), "The in-memory database test did not return the expected entries."

    conn.close()

    print("All tests passed successfully.")

test_find_matching_entries()

if __name__ == '__main__':
    # Git repository commands are written as comments and are to be executed in a shell.
    """
    git init
    git add py_exam4.py README.md
    git commit -m "Add initial files by Aarti: Python file and README"
    git log
    """
    
    # Execute tests for Q2, Q3, and Q4
    test_compute_volume_COM()
    test_get_counts_df()
    test_find_matching_entries()
