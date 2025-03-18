import ROOT as R
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--dataset", default="0")

args = parser.parse_args()
DATA  = f"datos/loop_{args.dataset}.feather"

''' 
Require events to have exactly one hit in each board
And to have CALcode values that are close to the average CALcode
'''
def filterEvents(event):
    if not len(event) == 4: return False
    if not len(event['board']) == len(np.unique(event['board'])): return False
    return True

# Function to write Feather file data to ROOT
def feather_to_root(feather_file, root_file):
    # Read the Feather file into a pandas DataFrame
    df = pd.read_feather(feather_file)

    # Open the ROOT file for writing
    root_file = R.TFile(root_file, "RECREATE")

    # Create a TTree object
    tree = R.TTree("data_tree", "Tree storing event data")

    # Get the unique events (evt)
    events = df['evt'].unique()
    print(events)

    # Create ROOT branches for each variable
    bcid = R.std.vector('int')()
    l1a_counter = R.std.vector('int')()
    ea = R.std.vector('int')()
    board = R.std.vector('int')()
    row = R.std.vector('int')()
    col = R.std.vector('int')()
    toa = R.std.vector('int')()
    tot = R.std.vector('int')()
    cal = R.std.vector('int')()

    # Add branches to the tree
    #tree.Branch('evt', R.std.vector('int')()) # evt will be a single int for each entry
    tree.Branch('bcid', bcid)
    tree.Branch('l1a_counter', l1a_counter)
    tree.Branch('ea', ea)
    tree.Branch('board', board)
    tree.Branch('row', row)
    tree.Branch('col', col)
    tree.Branch('toa', toa)
    tree.Branch('tot', tot)
    tree.Branch('cal', cal)

    # Iterate over each unique event 'evt'
    for evt in tqdm(events, total=max(events)):
        # Filter the dataframe for this event
        event_data = df[df['evt'] == evt]
        if not len(event_data) == 4: continue
        if not len(event_data['board']) == len(np.unique(event_data['board'])): continue
        #print(event_data)
        
        # Clear the vectors for the current event
        bcid.clear()
        l1a_counter.clear()
        ea.clear()
        board.clear()
        row.clear()
        col.clear()
        toa.clear()
        tot.clear()
        cal.clear()

        # Fill the vectors with the values for this event
        for index, row_data in event_data.iterrows():
            bcid.push_back(int(row_data['bcid']))
            l1a_counter.push_back(int(row_data['l1a_counter']))
            ea.push_back(int(row_data['ea']))
            board.push_back(int(row_data['board']))
            row.push_back(int(row_data['row']))
            col.push_back(int(row_data['col']))
            toa.push_back(int(row_data['toa']))
            tot.push_back(int(row_data['tot']))
            cal.push_back(int(row_data['cal']))

        # Set the evt value for this entry
        #tree.GetBranch('evt').SetAddress(R.std.vector('int')([int(evt)]))
        
        # Fill the tree with the current event
        tree.Fill()

    # Write the TTree to the ROOT file
    root_file.Write()

    # Close the ROOT file
    root_file.Close()

if __name__=='__main__':
    ## Read dataset
    #df = pd.read_feather(DATA)
    #print(df)

    ## Filter events
    #print("Filtering events...")
    #start = time.time()
    #df = df.groupby('evt').filter(filterEvents)
    #end = time.time()
    #print(f"Elapsed time: {end-start:.3f} s")

    feather_to_root(DATA, f"rootFiles/loop_{args.dataset}.root")
