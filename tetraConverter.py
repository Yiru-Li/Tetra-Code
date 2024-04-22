import numpy as np
import pandas as pd

def tetraConverter(input_list, intype='auto', outtype=['auto']):
    """
    Convert between MNI coordinates on the scalp/brain and Tetra Code. Also recognizes and supports conversion of CPC coordinates, Beam-style X/Y percentages, and 10-10 EEG electrode locations.
    Tries to automatically detect the input type and returns the corresponding scalp MNI coordinates (if input is Tetra Code) or Tetra Code (if input is anything else). Each row in the input would have a corresponding output
    
    Parameters:
        - input_list (2D numpy array or list of strings): inputs you want to convert
        - intype (string): input type, options are 'auto' (default), 'scalp', 'brain', 'tetra', 'CPC', 'XY', or 'EEG'
        - outtype (list of strings): output type, options are 'auto' (default), 'all', or a combination of 'scalp', 'brain', 'tetra', 'CPC', 'XY', and/or 'EEG'
            
    Returns:
        out: converted list (for Tetra Code or EEG location) or numpy array (for all other output types). If multiple or 'all' outtypes are specified, returns a dictionary of outputs
    """
    assert intype in ['auto', 'scalp', 'brain', 'tetra', 'CPC', 'XY', 'EEG'], "intype must be 'auto' (default), 'scalp', 'brain', 'tetra', 'CPC', 'XY', or 'EEG'"
    assert isinstance(outtype, list) and (set(outtype) <= set(['auto', 'scalp', 'brain', 'tetra', 'CPC', 'XY', 'EEG', 'all'])), "outtype must be a list consisting of 'auto' (default), 'all', or a combination of 'scalp', 'brain', 'tetra', 'CPC', 'XY', and/or 'EEG'"
    if isinstance(input_list, str):
        input_list = [input_list]
    if (not all(isinstance(item, str) for item in input_list)) and (not (type(input_list) == np.ndarray and input_list.ndim == 2)):
        try:
            input_list = np.array(input_list, ndmin=2)
        except:
            raise Exception('input_list must be a list of strings or convertable to a 2D array')
    if intype == 'auto':
        if type(input_list) == np.ndarray: # numeric input has to be CPC, XY, or MNI
            if input_list.shape[1] == 2: # assume CPC or XY input
                if any(input_list[0]>=1.1084) & any(input_list[1]>1):
                    intype = 'XY'
                    print('Assuming input_list is Beam-style X/Y percentages')
                else:
                    intype = 'CPC'
                    print('Assuming input_list is CPC')
            elif not input_list.shape[1] == 3: # otherwise assume MNI input; if not:
                raise Exception('Numeric input_list must be convertable to a Nx3 array (MNI) or Nx2 array (CPC or Beam-style X/Y percentages')
        else:
            print('Assuming input is Tetra Code or EEG')
            intype = 'tetra' # EEG and tetra are parsed the same way
            
    ref = pd.read_excel('tet_code2MNI_lookup_extended.xlsx', sheet_name='Reference')
    
    istetra = [0]*len(input_list)
    if intype == 'auto':
        assert input_list.shape[1] == 3, 'MNI input_list must be convertable to a Nx3 array'
        matchID = findMinDisID(np.vstack((np.array(ref[['Scalp X', 'Scalp Y', 'Scalp Z']]), 
                                          np.array(ref[['Brain X', 'Brain Y', 'Brain Z']]))), 
                               input_list) % len(ref)
    elif intype == 'brain':
        assert input_list.shape[1] == 3, 'MNI input_list must be convertable to a Nx3 array'
        matchID = findMinDisID(np.array(ref[['Brain X', 'Brain Y', 'Brain Z']]), input_list)
    elif intype == 'scalp':
        assert input_list.shape[1] == 3, 'MNI input_list must be in convertable to Nx3 array'
        matchID = findMinDisID(np.array(ref[['Scalp X', 'Scalp Y', 'Scalp Z']]), input_list)
    elif intype == 'CPC':
        assert input_list.shape[1] == 2, 'CPC input_list must be convertable to a Nx2 array'
        assert (input_list>=0).all() & (input_list<1.1084).all(), 'CPC input out of range [0, 1.1084)'
        matchID = findMinDisID(np.array(ref[['Pnz', 'Pal']]), input_list)
    elif intype == 'XY':
        assert input_list.shape[1] == 2, 'Beam-style X/Y input_list must be convertable to a Nx2 array'
        assert (input_list>=0).all() & (input_list<152.0659).all(), 'Beam-style X/Y percentage input out of range [0, 152.0659)'
        matchID = findMinDisID(np.array(ref[['X %', 'Y %']]), input_list)
    elif intype in ('EEG', 'tetra'):
        matchID = [None]*len(input_list)
        for row in range(len(input_list)):
            try:
                matchID[row], istetra[row] = np.where(input_list[row]==ref[['EEG', 'Tetra Code']])
                matchID[row] = matchID[row][0]
                istetra[row] = istetra[row][0]
            except:
                raise Exception('Match not found for "'+input_list[row]+'", please check spelling and capitalization')
    if 'auto' in outtype:
        if any(istetra):
            outtype = ['scalp']
        else:
            outtype = ['tetra']
    
    out = {'scalp': np.array(ref[['Scalp X', 'Scalp Y', 'Scalp Z']].iloc[matchID]),
           'brain': np.array(ref[['Brain X', 'Brain Y', 'Brain Z']].iloc[matchID]),
           'tetra': np.array(ref[['Tetra Code']].iloc[matchID]).tolist(),
           'CPC': np.array(ref[['Pnz', 'Pal']].iloc[matchID]),
           'XY': np.array(ref[['X %', 'Y %']].iloc[matchID]),
           'EEG': np.array(ref[['EEG']].iloc[matchID]).tolist()}
    if 'all' not in outtype:
        if len(outtype)==1:
            out = out[outtype[0]]
        else:
            out = {key: out[key] for key in outtype}
    return out

def findMinDisID(ref, targets):
    assert type(ref) == np.ndarray and type(targets) == np.ndarray, 'inputs must be arrays'
    assert ref.shape[1] == targets.shape[1], 'input arrays must have the same width'
    ref = ref.astype(float)
    targets = targets.astype(float)
    matchID = [None]*len(targets)
    for row in range(len(targets)):
        matchID[row] = np.argmin(np.sum(np.square(ref-targets[row, :]), axis=1))
    return np.array(matchID)