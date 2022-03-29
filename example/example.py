#  example.py
#  
#  Copyright 2017 Chichau Miau <zmiao@ebi.ac.uk>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import sys,os
import argparse
import RNA_normalizer

from operator import attrgetter

__doc__="Calculation of the RMSD and the Deformation Index"

def isPDBfile(path):
	if not os.path.isfile(path):
		if os.path.isdir(path):
			msg = "{0} is a directory (requires file)".format(path)
		else:
			msg = "{0} does not exist".format(path)
		raise argparse.ArgumentTypeError(msg)
	elif path[-4:] != ".pdb" and path[-4:] != ".ent":
		msg = "{0} is not a PDB file".format(path)
		raise argparse.ArgumentTypeError(msg)
	return path

def isdir(path):
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file (requires directory)".format(path)
        else:
            msg = "{0} does not exist".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """Retrieves the arguments of the program."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-n', dest='native_pdb', type=isPDBfile, required=True,
    					help = "PDB file of the native RNA structure")
    parser.add_argument('-d', dest='path_data', type=isdir, default=os.getcwd(),
    					help = "Dataset repository that contains PDB files of structure that should be compared\nto the native structure")
    parser.add_argument('-e', dest='exp_pdb', type=isPDBfile,
    					help = "PDB file of the experimental RNA structure (incompatible with argument -d)")
    parser.add_argument('-o', dest='out_file', type=str, default="RMSD_DI",
    					help = "Output file name")
    return parser.parse_args()



RESIDUES_LIST = "data/residues.list"
ATOMS_LIST = "data/atoms.list"

def CleanFormat(f):
	os.system( "mac2unix -q %s" %f )
	os.system( "dos2unix -q %s" %f )

def normalize_structure(struct, out_file = None, index_file=None, extract_file = None):
	pdb_normalizer = RNA_normalizer.PDBNormalizer( RESIDUES_LIST, ATOMS_LIST )
	ok = pdb_normalizer.parse( struct, out_file )
	if not ok:
		sys.stderr.write("ERROR: structure not normalized!\n")
	else:
		sys.stderr.write("INFO: Normalization succeded!\n")
	if not extract_file is None:
		coords=open(index_file).read()
		extract.extract_PDB(SOLUTION_NORMAL,coords, extract_file)
		sys.stderr.write("INFO:	structure extracted\n")

# PVALUE set according to Hajdin et al., RNA (7) 16, 2010, either "+" or "-"
def calc_RMSD(native_file, prediction_file, native_index = None, prediction_index = None, PVALUE = "-"):
	res_struct = RNA_normalizer.PDBStruct()
	res_struct.load( native_file, native_index )
	res_raw_seq = res_struct.raw_sequence()
	
	sol_struct = RNA_normalizer.PDBStruct()
	sol_struct.load( prediction_file, prediction_index )
	sol_raw_seq = sol_struct.raw_sequence()
	
	if( sol_raw_seq != res_raw_seq ):
		sys.stderr.write("ERROR Result sequence != Solution sequence!\n")
		sys.stderr.write("DATA Solution sequence --> '%s'\n" %sol_raw_seq )
		sys.stderr.write("DATA Result sequence   --> '%s'\n" %res_raw_seq )
		return(-1)
	# computes the RMSD
	comparer = RNA_normalizer.PDBComparer()
	rmsd = comparer.rmsd( sol_struct, res_struct )
	sys.stderr.write("INFO Partial RMSD --> %f\n" %rmsd )
	pvalue = comparer.pvalue( rmsd, len(sol_raw_seq), PVALUE )
	sys.stderr.write("INFO Partial P-Value --> %e\n" %pvalue )
	return(rmsd, pvalue)

def InteractionNetworkFidelity(native_file, prediction_file, native_index = None, prediction_index = None):
	res_struct = RNA_normalizer.PDBStruct()
	res_struct.load( native_file, native_index )
	res_raw_seq = res_struct.raw_sequence()
	
	sol_struct = RNA_normalizer.PDBStruct()
	sol_struct.load( prediction_file, prediction_index )
	sol_raw_seq = sol_struct.raw_sequence()
	
	if( sol_raw_seq != res_raw_seq ):
		sys.stderr.write("ERROR Result sequence != Solution sequence!\n")
		sys.stderr.write("DATA Solution sequence --> '%s'\n" %sol_raw_seq )
		sys.stderr.write("DATA Result sequence   --> '%s'\n" %res_raw_seq )
		return(-1)
	# computes the RMSD
	comparer = RNA_normalizer.PDBComparer()
	rmsd = comparer.rmsd( sol_struct, res_struct )
	INF_ALL = comparer.INF( sol_struct, res_struct, type="ALL" )
	DI_ALL = rmsd / INF_ALL
	#INF_WC = comparer.INF( sol_struct, res_struct, type="PAIR_2D" )
	#INF_NWC = comparer.INF( sol_struct, res_struct, type="PAIR_3D" )
	#INF_STACK = comparer.INF( sol_struct, res_struct, type="STACK" )
	return (rmsd,DI_ALL) # (INF_ALL, INF_WC, INF_NWC,INF_STACK)
	

if __name__ == '__main__':
	# Get arguments
	args = get_arguments()

	# Normalize PDB format, correct residue names and atom names. 
	#normalize_structure(native_pdb,native_pdb[:-4]+'_normalized.pdb')
	
	# calculate RMSD for RNA structures
	# require biopython
	#print(calc_RMSD(native_pdb, exp_pdb))

	# calculate InteractionNetworkFidelity and Deformation Index for RNA structures
	# need to have MA-annotate in the directory or set in mcannotate.py

	if args.exp_pdb:
		rmsd, DI_ALL = InteractionNetworkFidelity(args.native_pdb, args.exp_pdb)
		print("RMSD: {0}\tDeformation Index: {1}".format(rmsd, DI_ALL))
	else:
		out = open(args.out_file, "w")
		for file in os.listdir(args.path_data):
			if file[-4:] == ".pdb" or file[-4:] == ".ent":
				rmsd, DI_ALL = InteractionNetworkFidelity(args.native_pdb, args.path_data+file)
				out.write("{0}\t\t{1}\t{2}".format(file, rmsd, DI_ALL))
			else:
				sys.stderr.write("ERROR: File {0} is not a PDB file !\n".format(file))
		out.close()

