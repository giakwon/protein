import math, sys, os
import numpy as np
import time
from tqdm import tqdm
# from latticeproteins.interactions import miyazawa_jernigan
try:
    import cPickle as pickle
except ImportError:
    import pickle

miyazawa_jernigan = {
'CC':-5.44, 'CM':-5.05, 'CF':-5.63, 'CI':-5.03, 'CL':-5.03, 'CV':-4.46, 'CW':-4.76, 'CY':-3.89, 'CA':-3.38, 'CG':-3.16, 'CT':-2.88, 'CS':-2.86, 'CQ':-2.73, 'CN':-2.59, 'CE':-2.08, 'CD':-2.66, 'CH':-3.63, 'CR':-2.70, 'CK':-1.54, 'CP':-2.92,
'MC':-5.05, 'MM':-6.06, 'MF':-6.68, 'MI':-6.33, 'ML':-6.01, 'MV':-5.52, 'MW':-6.37, 'MY':-4.92, 'MA':-3.99, 'MG':-3.75, 'MT':-3.73, 'MS':-3.55, 'MQ':-3.17, 'MN':-3.50, 'ME':-3.19, 'MD':-2.90, 'MH':-3.31, 'MR':-3.49, 'MK':-3.11, 'MP':-4.11,
'FC':-5.63, 'FM':-6.68, 'FF':-6.85, 'FI':-6.39, 'FL':-6.26, 'FV':-5.75, 'FW':-6.02, 'FY':-4.95, 'FA':-4.36, 'FG':-3.72, 'FT':-3.76, 'FS':-3.56, 'FQ':-3.30, 'FN':-3.55, 'FE':-3.51, 'FD':-3.31, 'FH':-4.61, 'FR':-3.54, 'FK':-2.83, 'FP':-3.73,
'IC':-5.03, 'IM':-6.33, 'IF':-6.39, 'II':-6.22, 'IL':-6.17, 'IV':-5.58, 'IW':-5.64, 'IY':-4.63, 'IA':-4.41, 'IG':-3.65, 'IT':-3.74, 'IS':-3.43, 'IQ':-3.22, 'IN':-2.99, 'IE':-3.23, 'ID':-2.91, 'IH':-3.76, 'IR':-3.33, 'IK':-2.70, 'IP':-3.47,
'LC':-5.03, 'LM':-6.01, 'LF':-6.26, 'LI':-6.17, 'LL':-5.79, 'LV':-5.38, 'LW':-5.50, 'LY':-4.26, 'LA':-3.96, 'LG':-3.43, 'LT':-3.43, 'LS':-3.16, 'LQ':-3.09, 'LN':-2.99, 'LE':-2.91, 'LD':-2.59, 'LH':-3.84, 'LR':-3.15, 'LK':-2.63, 'LP':-3.06,
'VC':-4.46, 'VM':-5.52, 'VF':-5.75, 'VI':-5.58, 'VL':-5.38, 'VV':-4.94, 'VW':-5.05, 'VY':-4.05, 'VA':-3.62, 'VG':-3.06, 'VT':-2.95, 'VS':-2.79, 'VQ':-2.67, 'VN':-2.36, 'VE':-2.56, 'VD':-2.25, 'VH':-3.38, 'VR':-2.78, 'VK':-1.95, 'VP':-2.96,
'WC':-4.76, 'WM':-6.37, 'WF':-6.02, 'WI':-5.64, 'WL':-5.50, 'WV':-5.05, 'WW':-5.42, 'WY':-4.44, 'WA':-3.93, 'WG':-3.37, 'WT':-3.31, 'WS':-2.95, 'WQ':-3.16, 'WN':-3.11, 'WE':-2.94, 'WD':-2.91, 'WH':-4.02, 'WR':-3.56, 'WK':-2.49, 'WP':-3.66,
'YC':-3.89, 'YM':-4.92, 'YF':-4.95, 'YI':-4.63, 'YL':-4.26, 'YV':-4.05, 'YW':-4.44, 'YY':-3.55, 'YA':-2.85, 'YG':-2.50, 'YT':-2.48, 'YS':-2.30, 'YQ':-2.53, 'YN':-2.47, 'YE':-2.42, 'YD':-2.25, 'YH':-3.33, 'YR':-2.75, 'YK':-2.01, 'YP':-2.80,
'AC':-3.38, 'AM':-3.99, 'AF':-4.36, 'AI':-4.41, 'AL':-3.96, 'AV':-3.62, 'AW':-3.93, 'AY':-2.85, 'AA':-2.51, 'AG':-2.15, 'AT':-2.15, 'AS':-1.89, 'AQ':-1.70, 'AN':-1.44, 'AE':-1.51, 'AD':-1.57, 'AH':-2.09, 'AR':-1.50, 'AK':-1.10, 'AP':-1.81,
'GC':-3.16, 'GM':-3.75, 'GF':-3.72, 'GI':-3.65, 'GL':-3.43, 'GV':-3.06, 'GW':-3.37, 'GY':-2.50, 'GA':-2.15, 'GG':-2.17, 'GT':-2.03, 'GS':-1.70, 'GQ':-1.54, 'GN':-1.56, 'GE':-1.22, 'GD':-1.62, 'GH':-1.94, 'GR':-1.68, 'GK':-0.84, 'GP':-1.72,
'TC':-2.88, 'TM':-3.73, 'TF':-3.76, 'TI':-3.74, 'TL':-3.43, 'TV':-2.95, 'TW':-3.31, 'TY':-2.48, 'TA':-2.15, 'TG':-2.03, 'TT':-1.72, 'TS':-1.59, 'TQ':-1.59, 'TN':-1.51, 'TE':-1.45, 'TD':-1.66, 'TH':-2.31, 'TR':-1.97, 'TK':-1.02, 'TP':-1.66,
'SC':-2.86, 'SM':-3.55, 'SF':-3.56, 'SI':-3.43, 'SL':-3.16, 'SV':-2.79, 'SW':-2.95, 'SY':-2.30, 'SA':-1.89, 'SG':-1.70, 'ST':-1.59, 'SS':-1.48, 'SQ':-1.37, 'SN':-1.31, 'SE':-1.48, 'SD':-1.46, 'SH':-1.94, 'SR':-1.22, 'SK':-0.83, 'SP':-1.35,
'QC':-2.73, 'QM':-3.17, 'QF':-3.30, 'QI':-3.22, 'QL':-3.09, 'QV':-2.67, 'QW':-3.16, 'QY':-2.53, 'QA':-1.70, 'QG':-1.54, 'QT':-1.59, 'QS':-1.37, 'QQ':-0.89, 'QN':-1.36, 'QE':-1.33, 'QD':-1.26, 'QH':-1.85, 'QR':-1.85, 'QK':-1.02, 'QP':-1.73,
'NC':-2.59, 'NM':-3.50, 'NF':-3.55, 'NI':-2.99, 'NL':-2.99, 'NV':-2.36, 'NW':-3.11, 'NY':-2.47, 'NA':-1.44, 'NG':-1.56, 'NT':-1.51, 'NS':-1.31, 'NQ':-1.36, 'NN':-1.59, 'NE':-1.43, 'ND':-1.33, 'NH':-2.01, 'NR':-1.41, 'NK':-0.91, 'NP':-1.43,
'EC':-2.08, 'EM':-3.19, 'EF':-3.51, 'EI':-3.23, 'EL':-2.91, 'EV':-2.56, 'EW':-2.94, 'EY':-2.42, 'EA':-1.51, 'EG':-1.22, 'ET':-1.45, 'ES':-1.48, 'EQ':-1.33, 'EN':-1.43, 'EE':-1.18, 'ED':-1.23, 'EH':-2.27, 'ER':-2.07, 'EK':-1.60, 'EP':-1.40,
'DC':-2.66, 'DM':-2.90, 'DF':-3.31, 'DI':-2.91, 'DL':-2.59, 'DV':-2.25, 'DW':-2.91, 'DY':-2.25, 'DA':-1.57, 'DG':-1.62, 'DT':-1.66, 'DS':-1.46, 'DQ':-1.26, 'DN':-1.33, 'DE':-1.23, 'DD':-0.96, 'DH':-2.14, 'DR':-1.98, 'DK':-1.32, 'DP':-1.19,
'HC':-3.63, 'HM':-3.31, 'HF':-4.61, 'HI':-3.76, 'HL':-3.84, 'HV':-3.38, 'HW':-4.02, 'HY':-3.33, 'HA':-2.09, 'HG':-1.94, 'HT':-2.31, 'HS':-1.94, 'HQ':-1.85, 'HN':-2.01, 'HE':-2.27, 'HD':-2.14, 'HH':-2.78, 'HR':-2.12, 'HK':-1.09, 'HP':-2.17,
'RC':-2.70, 'RM':-3.49, 'RF':-3.54, 'RI':-3.33, 'RL':-3.15, 'RV':-2.78, 'RW':-3.56, 'RY':-2.75, 'RA':-1.50, 'RG':-1.68, 'RT':-1.97, 'RS':-1.22, 'RQ':-1.85, 'RN':-1.41, 'RE':-2.07, 'RD':-1.98, 'RH':-2.12, 'RR':-1.39, 'RK':-0.06, 'RP':-1.85,
'KC':-1.54, 'KM':-3.11, 'KF':-2.83, 'KI':-2.70, 'KL':-2.63, 'KV':-1.95, 'KW':-2.49, 'KY':-2.01, 'KA':-1.10, 'KG':-0.84, 'KT':-1.02, 'KS':-0.83, 'KQ':-1.02, 'KN':-0.91, 'KE':-1.60, 'KD':-1.32, 'KH':-1.09, 'KR':-0.06, 'KK':0.13, 'KP':-0.67,
'PC':-2.92, 'PM':-4.11, 'PF':-3.73, 'PI':-3.47, 'PL':-3.06, 'PV':-2.96, 'PW':-3.66, 'PY':-2.80, 'PA':-1.81, 'PG':-1.72, 'PT':-1.66, 'PS':-1.35, 'PQ':-1.73, 'PN':-1.43, 'PE':-1.40, 'PD':-1.19, 'PH':-2.17, 'PR':-1.85, 'PK':-0.67, 'PP':-1.18}



class ConformationsError(Exception):
    """Error finding or storing a conformation."""
#----------------------------------------------------------------------
# '_loop_in_C' is a Boolean switch specifying whether we try to speed
# up the energy calculations by doing the looping with the C-extension
# 'contactlooper'.  'True' means we try to do this.
# _loop_in_C = True
# if _loop_in_C:
    # from latticeproteins.contactlooper import ContactLooper

class PickleProtocolError(Exception):
    """Error is pickle version is too old. """

PROTOCOL = pickle.HIGHEST_PROTOCOL
if PROTOCOL < 2:
    raise PickleProtocolError("Version of pickle is outdated for this package.")

def fold_energy(sequence, conformation, interactions=miyazawa_jernigan):
    """Calculate the energy of the sequence with the given conformation.

    Parameters
    ----------
    sequence : str
        Amino acid sequence to fold.
    conformation : str
        Conformation according to latticemodel's conformations format (e.g. 'UDLLDRU')

    Returns
    ------
    energy : float
        energy of the conformation (sum of all contact energies)
    """
    contacts = lattice_contacts(sequence, conformation)
    energy = sum([interactions[c] for c in contacts])
    return energy

def lattice_contacts(sequence, conformation):
    """Find all contacts in conformation.

    Parameters
    ----------
    sequence : str
        Amino acid sequence to fold.
    conformation : str
        Conformation according to latticemodel's conformations format (e.g. 'UDLLDRU')

    Returns
    -------
    contacts : list
        list of contact pairs
    """
    sites = list(sequence)
    length = len(sites)
    try:
        moves = list(conformation)
    except TypeError:
        raise Exception("""Protein conformation is None; is there a native state? """)
    # build a coordinate system, note that odd rotation of intuitive coordinates
    # since we are working in numpy array grid.
    coordinates = {"U": [-1,0], "D":[1,0], "L":[0,-1], "R":[0,1]}
    grid = np.zeros((2*length+1, 2*length+1), dtype=str)
    x = y = length # initial position on the grid is at the center of the 2d array
    grid[x,y] = sites[0]
    # move on grid, populate with amino acid at that site, and store all contacting neighbors.
    contacts = []
    for i in range(length-1):
        step = coordinates[moves[i]]
        x += step[0]
        y += step[1]
        grid[x,y] = sites[i+1]
        neighbors = [sites[i+1] + grid[x+c[0], y+c[1]] for c in coordinates.values()]
        contacts += [n for n in neighbors if len(n) == 2]
    # subtract the contacts that have bonds between them.
    for i in range(1,length):
        # Try forward
        try:
            contacts.remove(sequence[i-1:i+1])
        # Else try reverse
        except ValueError:
            contacts.remove(sequence[i] + sequence[i-1])
    return contacts


class PickleProtocolError(Exception):
    """Error is pickle version is too old. """

PROTOCOL = pickle.HIGHEST_PROTOCOL
if PROTOCOL < 2:
    raise PickleProtocolError("Version of pickle is outdated for this package.")


class Conformations(object):
    """Creates a database of conformations for a protein of specified length.

    The created 'Conformations' object 'c' stores the contact
        lists and the number of conformations with these contact sets
        for all self-avoiding walks of length 'length'.  It can then
        be used to compute the free energy of a protein folding to
        the lowest energy conformation.

    Parameters
    ----------
    length :
        is an integer specifying the length of the protein
        for which we are computing the contacts.  It must be >= 2.
    database_dir :
        specifies the name of the database directory storing
        existing conformations.  If the conformation instance
        already exists in this database we return the
        existing data, and if it doesn't we store it in the database.
    interaction_energies :
        specifies the interaction energies between
        residues. By default, this is interactions.miyazawa_jernigan.


    Attributes
    ----------
    _numconformations : dict
        A dictionary mapping the number of contact sets to the number of conformations
        with that contact set.
    _contactsets : list of lists
        'self._contactsets' is a list of contact sets.  'self._contactsets[i]'
        is the contact set for contact i.  It is a list of numbers.
        'x = self._contactsets[i]' describes the residues in contact
        in contact 'i'.  If this contact is between residues 'ires'
        and 'jres', then 'x = self._length * ires + jres' where
        0 <= ires, jres < 'self._length', and ires < jres + 1 contact sets
    _contactsetdegeneracy : list
        'self._contactsetdegeneracy' is a list of integers giving the
        degeneracy of the contact sets (the number of different conformations
        with this contact set).  'self_contactsetdegeneracy[i]' is the
        degeneracy of the contact set 'self._contactsets[i]'
    _contactsetconformation : list
        'self._contactsetconformation' is a list of the conformations
        associated with each contact set.  If contact set 'self._contactsets[i]'
        is degenerate ('self._contactsetdegeneracy[i]' > 1), the value
        'self._contactsetconformation[i]' is 'None'.  Otherwise, it
        is the string representing the conformation that gives rise
        to contact set 'self._contactsets[i]'.  The conformations are given
        such that 'self._contactsetconformation[i][j]'
        gives the conformation of bond 'j' (0 <= j < 'self._length' - 1)
        as 'U' (Up), 'R' (Right), 'D' (Down), or 'L' (Left).
        We require the first bond to be Up, and the first
        non-Up bond to be Right.
    _numcontactsets : dict
        'self._numcontactsets[i]' holds the number of different contact
        sets with 'i' contacts.
    """
    def __init__(self, length, database_dir="database/", interaction_energies=miyazawa_jernigan):
        self._interaction_energies = interaction_energies
        if not os.path.isdir(database_dir):
            os.makedirs(database_dir)
            #raise IOError("Cannot find database directory of %s." % database_dir)
        object_properties = ['_length', '_numconformations', '_contactsets', '_contactsetdegeneracy', '_contactsetconformation', '_numcontactsets']
        foundone = didntfindone = False
        for prop in object_properties:
            f = "%s/%d%s.pickle" % (database_dir, length, prop)
            if os.path.isfile(f):
                try:
                    val = pickle.load(open("%s/%d_length.pickle" % (database_dir, length), 'rb'))
                    foundone = True
                except ValueError:
                    foundone = False
            else:
                didntfindone = True
        if foundone and not didntfindone:
            # return existing values
            self._length = pickle.load(open("%s/%d_length.pickle" % (database_dir, length), 'rb'))
            if self._length != length:
                raise ValueError("Length mismatch.")
            self._numconformations = pickle.load(open("%s/%d_numconformations.pickle" % (database_dir, length), 'rb'))
            self._contactsets = pickle.load(open("%s/%d_contactsets.pickle" % (database_dir, length), 'rb'))
            self._contactsetdegeneracy = pickle.load(open("%s/%d_contactsetdegeneracy.pickle" % (database_dir, length), 'rb'))
            self._contactsetconformation = pickle.load(open("%s/%d_contactsetconformation.pickle" % (database_dir, length), 'rb'))
            self._numcontactsets = pickle.load(open("%s/%d_numcontactsets.pickle" % (database_dir, length), 'rb'))
            # sort the contact set information so that the conformations
            # with the most contacts come first
            n = len(self._contactsets)
            decorated_list = [(len(self._contactsets[i]), self._contactsets[i], self._contactsetdegeneracy[i], self._contactsetconformation[i]) for i in range(n)]
            decorated_list.sort()
            decorated_list.reverse()
            self._contactsets = [decorated_list[i][1] for i in range(n)]
            self._contactsetdegeneracy = [decorated_list[i][2] for i in range(n)]
            self._contactsetconformation = [decorated_list[i][3] for i in range(n)]
            self._foldedsequences = {}
            return
        elif foundone and didntfindone:
            raise ValueError("Found some but not all conformations for length %d in %s." % (length, database_dir))
        # If we have made it here, we are generating new information for this conformation
        # A listing of all class variables
        # 'self._length' is the length of the protein (self-avoiding walk)
        self._length = length
        # there are 'self._numconformations[i]' conformations with 'i' contacts
        self._numconformations = {}
        # 'self._contactsets' is a list of contact sets.  'self._contactsets[i]'
        # is the contact set for contact i.  It is a list of numbers.
        # 'x = self._contactsets[i]' describes the residues in contact
        # in contact 'i'.  If this contact is between residues 'ires'
        # and 'jres', then 'x = self._length * ires + jres' where
        # 0 <= ires, jres < 'self._length', and ires < jres + 1
        self._contactsets = []
        # 'self._contactsetdegeneracy' is a list of integers giving the
        # degeneracy of the contact sets (the number of different conformations
        # with this contact set).  'self_contactsetdegeneracy[i]' is the
        # degeneracy of the contact set 'self._contactsets[i]'
        self._contactsetdegeneracy = []
        # 'self._contactsetconformation' is a list of the conformations
        # associated with each contact set.  If contact set 'self._contactsets[i]'
        # is degenerate ('self._contactsetdegeneracy[i]' > 1), the value
        # 'self._contactsetconformation[i]' is 'None'.  Otherwise, it
        # is the string representing the conformation that gives rise
        # to contact set 'self._contactsets[i]'.  The conformations are given
        # such that 'self._contactsetconformation[i][j]'
        # gives the conformation of bond 'j' (0 <= j < 'self._length' - 1)
        # as 'U' (Up), 'R' (Right), 'D' (Down), or 'L' (Left).
        # We require the first bond to be Up, and the first
        # non-Up bond to be Right.
        self._contactsetconformation = []
        # 'self._numcontactsets[i]' holds the number of different contact
        # sets with 'i' contacts.
        self._numcontactsets = {}
        # 'self._foldedsequences = {}' holds recently folded sequences.
        # The keys are the arguments to 'FoldSequence', and the items
        # are how many times the saved sequence was accessed and the
        # information about the folded sequence.  Sequences in the keys
        # are strings.
        self._foldedsequences = {}
        #
        # Now do some error checking on input variables
        if not (isinstance(self._length, int) and self._length >= 2):
            raise ConformationsError("Invalid 'length' of %r." % self._length)
        #
        # The initial generation of conformations uses the variable
        # 'contactsets' which is a dictionary as described below.  Once
        # we have constructed this dictionary, we use it to create
        # 'self._contactsets', 'self._contactsetdegeneracy', and
        # 'self._contactsetconformation'.  The format of 'contactsets':
        # The keys are string representing the contact sets.  These
        # strings are the numbers that will appear in 'self._contactsets',
        # separated by spaces.
        # The values of 'contactsets' provides information
        # about the conformations encoding the contact sets.  If there is
        # a single conformation encoding a contact set 'cs', then
        # 'contactsets[cs]' is a string.  If the contact set stored as
        # 'contactsets[cs]' is coded for by multiple conformations,
        # then 'self._contactsets[cs]' is an integer > 1 that represents
        # the number of conformations coding for this contact set.
        contactsets = {}
        # Now being looping over all possible conformations
        # The initial conformation is all 'U'
        dx = {'U':0, 'R':1, 'D':0, 'L':-1}
        dy = {'U':1, 'R':0, 'D':-1, 'L':0}
        next = {'U':'R', 'R':'D', 'D':'L', 'L':'U'}
        opposite = {'U':'D', 'D':'U', 'R':'L', 'L':'R'}
        n = self._length - 2 # index of last bond in 'conformation'
        conformation = ['U' for i in range(n + 1)]
        first_R = n # index of the first 'R' in the conformation
        ncount = 0
        while True: # keep finding new conformations
            # See if the current conformation has overlap
            # If no overlap, store the contact set
            x = y = j = 0
            res_positions = {(x, y):j} # keyed by coords, items are residue numbers
            res_coords = [(x, y)] # 'res_coords[j]' is coords of residue 'j'
            for c in conformation:
                x += dx[c]
                y += dy[c]
                coord = (x, y)
                if coord in res_positions: # overlap
                    # increment at the step that gave the problem
                    for k in range(j + 1, n + 1):
                        conformation[k] = 'U'
                    conformation[j] = next[conformation[j]]
                    while conformation[j] == 'U':
                        j -= 1
                        conformation[j] = next[conformation[j]]
                    if j == first_R and conformation[j] not in ['R', 'U']:
                        first_R -= 1
                        conformation[first_R] = 'R'
                        for k in range(j, n + 1):
                            conformation[k] = 'U'
                    break
                j += 1
                res_positions[coord] = j
                res_coords.append(coord)
            else: # loop finishes normally, this is a valid conformation
                # generate the contact set
                cs = []
                numcontacts = 0
                for j in range(len(res_coords)):
                    (x, y) = res_coords[j]
                    partners_list = []
                    for c in ['U', 'R', 'D', 'L']:
                        try:
                            k = res_positions[(x + dx[c], y + dy[c])]
                            if k > j + 1:
                                partners_list.append(k)
                        except KeyError:
                            pass
                    partners_list.sort()
                    jtimeslength = j * self._length
                    for k in partners_list:
                        cs.append("%d " % (jtimeslength + k))
                        numcontacts += 1
                cs = ''.join(cs)
                # add conformation to count
                try:
                    self._numconformations[numcontacts] += 1
                except KeyError:
                    self._numconformations[numcontacts] = 1
                # store contact set
                try:
                    contactsets[cs] += 1
                except KeyError:
                    contactsets[cs] = ''.join(conformation)
                    try:
                        self._numcontactsets[numcontacts] += 1
                    except KeyError:
                        self._numcontactsets[numcontacts] = 1
                except TypeError:
                    contactsets[cs] = 2
                # generate the next conformation
                i = n
                conformation[i] = next[conformation[i]]
                while conformation[i] == 'U':
                    i -= 1
                    conformation[i] = next[conformation[i]]
                # make sure first non-'U' is 'R'
                if i == first_R and conformation[i] not in ['R', 'U']:
                    first_R -= 1
                    conformation[first_R] = 'R'
                    for j in range(i, n + 1):
                        conformation[j] = 'U'
            # see if we are done
            if first_R == 0:
                break
        #
        # Now use 'contactsets' to generate 'self._contactsets',
        # 'self._contactsetdegeneracy', and 'self._contactsetconformation'
        for (cs, n_or_conf) in contactsets.items():
            # convert 'cs' to the appropriate list
            clist = [int(x) for x in cs.split()]
            self._contactsets.append(clist)
            if isinstance(n_or_conf, str):
                self._contactsetdegeneracy.append(1)
                self._contactsetconformation.append(n_or_conf)
            else:
                self._contactsetdegeneracy.append(n_or_conf)
                self._contactsetconformation.append(None)
        #
        del contactsets
        #
        # to make the energy evaluations quicker, sort so that
        # contact sets with more conformations are first:
        decorated_list = [(len(self._contactsets[i]), self._contactsets[i], self._contactsetdegeneracy[i], self._contactsetconformation[i]) for i in range(len(self._contactsets))]
        del self._contactsets, self._contactsetdegeneracy, self._contactsetconformation
        decorated_list.sort()
        decorated_list.reverse()
        self._contactsets = [decorated_list[i][1] for i in range(len(decorated_list))]
        self._contactsetdegeneracy = [decorated_list[i][2] for i in range(len(decorated_list))]
        self._contactsetconformation = [decorated_list[i][3] for i in range(len(decorated_list))]
        #
        # store the conformations data in the database
        pickle.dump(self._length, open("%s/%d_length.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._numconformations, open("%s/%d_numconformations.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._contactsets, open("%s/%d_contactsets.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._contactsetdegeneracy, open("%s/%d_contactsetdegeneracy.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._contactsetconformation, open("%s/%d_contactsetconformation.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._numcontactsets, open("%s/%d_numcontactsets.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)

    def length(self):
        """Returns the length of the protein these conformations are for."""
        return self._length

    def k_lowest_confs(self, seq, temp, k):
        """Get the `k` lowest conformations in the sequence's conformational ensemble.
        """
        length = len(seq)
        # Calculate the kth lowest conformations
        ncontacts = self.max_contacts()
        confs = []
        nconfs, n = 0, ncontacts
        while nconfs < k and n > 1:
            confs += self.unique_conformations(n)
            nconfs = len(confs)
            n -= 1
        # Loop through conformations and find lowest energy.
        confs = np.array(confs)
        energies = np.empty(len(confs), dtype=float)
        for i, conf in enumerate(confs):
            energies[i] = fold_energy(seq, conf, interactions=self._interaction_energies)
        sorted_e = np.argsort(energies)
        states = confs[sorted_e[0:k]]
        return states

    def fold_sequence(self, seq, temp):
        """Folds a protein sequence; calculates native energy and partition sum.

        Parameters
        ----------
        seq : string
            is the sequence of the protein to be folded as one-letter amino
            acid codes.  It should be a string or list of length 'c.Length()'.
        temp :
            is the temperature at which the protein is to be folded.  It
            must be a number > 0.  It represents a reduced temperature, scaled
            so that a value of 1 represents 273 K.

        Returns
        -------
        minE : float
            The energy of the lowest energy conformation.
        conf : str
            Lowest energy conformation.
        partitionsum : float
            Total partition function sum.
        numcontacts : int
            Number of contacts in the native conformation.
        folds : bool
            True if lattice protein has a single lowest energy.
        """
        # do some error checking on the input variables
        if len(seq) != self._length:
            raise ConformationsError("Invalid 'seq' length: is %r but should be %r." % (len(seq), self._length))
        try:
            temp = float(temp)
            if temp <= 0.0:
                raise ConformationsError("'temp' is <= 0: %r." % temp)
        except KeyError:
            raise ConformationsError("Invalid 'temp' of %r." % temp)

        # create 'res_interactions'.  'res_interactions[j]' holds the energy
        # for the interaction 'j' as specified in 'self._contactsets[i][j]'
        res_interactions = [0.0 for i in range(self._length**2)]
        for ires in range(self._length):
            itimeslength = ires * self._length
            for jres in range(ires + 1, self._length):
                respair = "%s%s" % (seq[ires], seq[jres])
                res_interactions[itimeslength + jres] = self._interaction_energies[respair]

        # Use C function to loop through contacts
        (minE, ibest, partitionsum, folds) = ContactLooper(res_interactions, self._contactsets, self._contactsetdegeneracy, float(temp))

        # Set folds to boolean
        if folds == 1:
            folds = True
        else:
            folds = False

        # Get the native conformation
        conf = self._contactsetconformation[ibest]
        return minE, conf, partitionsum, folds

    def unique_conformations(self, numcontacts):
        """Gets all unique conformations with specified number of contacts.

        Parameters
        ----------
        numcontacts : int
            Number of contacts to include in unique conformations list.

        Returns
        -------
        clist : list
            is of all unique conformations
            with exactly 'numcontacts' contacts.  A conformation
            is "unique" if it is the only conformation that gives
            rise to its particular contact set.  If there are
            no unique conformations with 'numcontacts' contacts,
            'clist' is an empty list.  Conformations are specified
            as strings of 'U', 'R', 'L', and 'D' as described in
            'FoldSequence'."""
        if not (isinstance(numcontacts, int) and numcontacts >= 0):
            raise ConformationsError("Invalid 'numcontacts' of %r." % numcontacts)
        clist = []
        for i in range(len(self._contactsets)):
            if self._contactsetdegeneracy[i] == 1:
                if len(self._contactsets[i]) == numcontacts:
                    clist.append(self._contactsetconformation[i])
        return clist

    def max_contacts(self):
        """Gets the most contacts of any conformation.

        Returns
        -------
        n : int
            is returned as the number of contacts for the conformation
            with the most contacts."""
        n = 0
        for cs in self._contactsets:
            if len(cs) > n:
                n = len(cs)
        return n

    def num_conformations(self, contacts = None):
        """Returns the number of conformations.

        If 'contacts' has its default value of 'None', returns the total
            number of conformations (self-avoiding walks).
        If 'contacts' has an integer value, returns the number of conformations
            with 'contacts' contacts.  If there are no walks with this number
            of contacts, returns 0.
        """
        if contacts == None:
            n = 0
            for x in self._numconformations.values():
                n += x
        else:
            try:
                n = self._numconformations[contacts]
            except KeyError:
                if isinstance(contacts, int) and contacts >= 0:
                    n = 0
                else:
                    raise ConformationsError("Invalid 'contacts' of %r." % contacts)
        return n

    def num_contact_sets(self, contacts = None):
        """Returns the number of unique contact sets.

        If 'contacts' has its default value of 'None', returns the total
            number of unique contact sets (defined as the list of all
            contacts of non-adjacent residues).
        If 'contacts' has an integer value, returns the number of unique
            contact sets with 'contacts' contacts.  If there are no
            contact sets with this number of contacts, returns 0."""
        if contacts == None:
            return len(self._contactsets)
        else:
            try:
                return self._numcontactsets[contacts]
            except KeyError:
                if isinstance(contacts, int) and contacts >= 0:
                    n = 0
                else:
                    raise Conformationserror("Invalid 'contacts' of %r." % contacts)
                


def main():

    start_time = time.time()  
    c = Conformations(length=8, database_dir='database')
    # seq = "HPHPPH"
    seq = "HPPHHPHP"
    # seq = "HPHPPHHPHHPHHPPHPH"
    temp = 1.0  
    k = 100

    # print(c)

    # lowest_confs = c.k_lowest_confs(seq, temp, k)
    # print(lowest_confs)

    for i in tqdm(range(10)):
        print(i, c.unique_conformations(i))
        
main()