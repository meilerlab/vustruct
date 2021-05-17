"""
Given command line arguments in args, both ddg_launch.py and ddg_run.py must
be able to load a source structure, and prepare a cross reference dictionary
to the Rosetta ready structure, if not already available.  Place all that in one module for code sharing.
"""

import os
import datetime
import gzip
import lzma
import logging
import warnings

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio import BiopythonWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from psb_shared import psb_config
from psb_shared.ddg_repo import DDG_repo
from psb_shared import ddg_clean
from psb_shared.ddg_monomer import DDG_monomer
from psb_shared.psb_progress import PsbStatusManager

from lib import PDBMapSwiss
from lib import PDBMapModbase2020
from lib import PDB36Parser
# from lib import PDBMapModbase2016
# from lib import PDBMapModbase2013

LOGGER = logging.getLogger(__name__)



class LoadStructureError(Exception):
    """Exception raised in case of failure to load a structure file expected in the filesystem"""

    def __init__(self,message):
        self.message = message
        super().__init__(self.message)




def ddg_load_structure(args, config_dict):
    structure_info_dict = {}
    # We need to load the structure, clean it, and save the xref

    def _load_structure_file(id, coord_filename):
        tryParser36 = False
        with gzip.open(coord_filename, 'rt') if coord_filename.split('.')[-1] == "gz" else \
          lzma.open(coord_filename, 'rt') if coord_filename.split('.')[-1] == "xz" else \
          open(coord_filename,'r') as fin:
            warnings.filterwarnings('ignore', category=PDBConstructionWarning)
            try:
                structure = PDBParser().get_structure(id, fin)
            except ValueError:
                tryParser36 = True # We will make a last ditch effort to read this because of alpha in int columns
                structure = None
        if tryParser36: # Try hybrid36 format on modbase
           LOGGER.critical("ValueError with traditional parser - how trying hybrid36 parser on %s"%coord_filename)
           with gzip.open(coord_filename, 'rt') if coord_filename.split('.')[-1] == "gz" else \
             lzma.open(coord_filename, 'rt') if coord_filename.split('.')[-1] == "xz" else \
             open(coord_filename,'r') as fin:
                 warnings.filterwarnings('ignore', category=PDBConstructionWarning)
                 structure = PDB36Parser().get_structure(id, fin)
        return structure

    def mine_structure_info_from_mmCIF(mmCIF_dict):
        method = mmCIF_dict['_exptl.method'][0]
        pdb_id = method, mmcif_dict['_entry.id']
        deposition_date = None

        raw_deposition_date = mmCIF_dict.get('_pdbx_database_status.recvd_initial_deposition_date', None)
        if raw_deposition_date:  # Split into YYYY-MM-DD
            raw_YYYY_MM_DD = raw_deposition_date[0].split('-')
            if len(raw_YYYY_MM_DD) == 3:
                deposition_date = datetime.datetime(
                    int(raw_YYYY_MM_DD[0]), int(raw_YYYY_MM_DD[1]), int(raw_YYYY_MM_DD[2]))

            if not deposition_date:
                LOGGER.warning("Deposition date %s apparently not YYYY-MM-DD", raw_deposition_date)

        resolution = None  # Not applicable is correct for NMR structures
        resolution_keys = ['_refine.ls_d_res_high', '_reflns.d_resolution_high', '_refine_hist.d_res_high']
        if method == 'X-RAY DIFFRACTION':
            method = 'X-RAY'
        elif method == 'ELECTRON MICROSCOPY':
            method = 'EM'
            resolution_keys = ['_em_3d_reconstruction.resolution']
        elif method.find('NMR') > -1:  # resolution is not applicable for nmr
            resolution_keys = []
        elif method.find('SOLUTION SCATTERING') > -1:
            # SOLUTION SCATTERING does not have 'resolution' notion
            # But the keys are often in the .cif file - and we'll try to fish something out.
            # Generally, it will flow through like NMR
            method = 'SCAT'
        else:
            LOGGER.critical('%s Experimental method for %s is unknown' % (method, pdb_id))

        for resolution_key in resolution_keys:
            if resolution_key in mmCIF_dict:
                # Congratulations.  You found the old REMARK 2 RESOLUTION entry.... or so we hope.
                if len(mmCIF_dict[resolution_key]) > 0 and mmCIF_dict[resolution_key][0].strip() != '.':
                    resolution = float(mmCIF_dict[resolution_key][0])
                break
        if not resolution and method.find('X-RAY') > 0:
            LOGGER.warning("pdb %s is method=%s with no resolution entry" % (pdb_id, method))

        structure_info_dict = {'method': method, 'deposition_date': deposition_date}
        if resolution:
            structure_info_dict['resolution'] = resolution
        return structure_info_dict

    structure = None
    if args.pdb:
        # Load mmCIF from pdb repository
        structure_id = args.pdb.lower()
        original_structure_filename = os.path.join(
            config_dict['pdb_dir'],
            'structures',
            'divided',
            'mmCIF',
            structure_id[1:3],
            "%s.cif.gz" % structure_id)
        original_structure_type = 'mmCIF'
    elif args.swiss:
        PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'], config_dict['swiss_summary']);
        swiss_info = PDBMapSwiss.get_info(args.swiss)
        assert swiss_info, "%s was not found in the swiss INDEX JSON file" % args.swiss
        original_structure_type = 'pdb'
        original_structure_filename = PDBMapSwiss.get_coord_file(args.swiss)
        structure_id = swiss_info['modelid']
        structure_info_dict['pdb_template'] = swiss_info['template']
        # It is OK to run ddG on any chain in a multimeric swiss model
        # until Andrew Waterhouse gives more details
        # The swiss template info has source PDB chain in the last position
        # pdb_template_chain = structure_info_dict['pdb_template'][-1]
        # assert pdb_template_chain == args.chain, \
        #     "With swiss models, chain on command line %s must match modelled chain %s" % (
        #         args.chain, pdb_template_chain)
        remark3_metrics = PDBMapSwiss.load_REMARK3_metrics(structure_id)

        if not DDG_monomer.evaluate_swiss(structure_id, remark3_metrics):
            raise LoadStructureError(
                "Terminating as %s of insufficient quality for ddg monomer" % args.swiss)

        structure_info_dict['method'] = remark3_metrics['mthd']
        structure_info_dict['template_identity'] = float(remark3_metrics['sid'])
    elif args.modbase:  # This is a ENSP.... modbase file.  Only supporting modbase 2020 for time being
        structure_id = args.modbase.lower()
        original_structure_type = 'pdb'
        modbase20 = PDBMapModbase2020(config_dict)
        original_structure_filename = modbase20.get_coord_file(args.modbase)
        if os.path.exists(original_structure_filename):
            structure = _load_structure_file(args.modbase, original_structure_filename)
        else:
            raise LoadStructureError(
                  "Modbase model id %s not found in modbase2020 directory" % \
                  args.modbase)
        # modbase16 = PDBMapModbase2016(config_dict)
        # original_structure_filename = modbase16.get_coord_file(args.modbase)
        # if os.path.exists(original_structure_filename):
        #    structure = _load_structure_file(args.modbase, original_structure_filename)
        # else:
        #     modbase13 = PDBMapModbase2013(config_dict)
        #     original_structure_filename = modbase13.get_coord_file(args.modbase)
        #     if os.path.exists(original_structure_filename):
        #         structure = _load_structure_file(args.modbase, original_structure_filename)
        #     else:
        #         raise LoadStructureError(
        #             "Modbase model id %s not found in modbase2013/2016 directories" % \
        #            args.modbase)

        if not DDG_monomer.evaluate_modbase(original_structure_filename):
            raise LoadStructureError(
                "Terminating as %s of insufficient quality for ddg monomer" % args.modbase)
    else:
        assert args.usermodel is not None
        structure_id = args.usermodel.lower()
        original_structure_type = 'pdb'
        original_structure_filename = args.usermodel
        if os.path.exists(original_structure_filename):
            structure = _load_structure_file(args.modbase, original_structure_filename)
        else:
            raise LoadStructureError(
                  "UserModel %s not found" % original_structure_filename)

    if not os.path.exists(original_structure_filename):
        exit_str = "%s: No structure file %s" % (structure_id, original_structure_filename)
        raise LoadStructureError(exit_str)

    # Open a file handle to our structure
    # LOGGER.debug("Loading orignal structure %s", original_structure_filename)
    LOGGER.info("Loading structure id=%s from %s" % (structure_id, original_structure_filename))
    structure_fin = None
    if original_structure_filename.endswith(".gz"):
        structure_fin = gzip.open(original_structure_filename, 'rt')
    elif original_structure_filename.endswith(".xz"):
        structure_fin = lzma.open(original_structure_filename, 'rt')
    else:
        structure_fin = open(original_structure_filename, 'rt')

    mmcif_dict = None
    if original_structure_type == 'mmCIF':
        mmCIF_parser = MMCIFParser()  # QUIET=True)
        structure = mmCIF_parser.get_structure(structure_id, structure_fin)
        structure_fin.close()

        # DDG_Cleaner uses the mmcif information to remove non-species tags in PDB files
        mmcif_dict = mmCIF_parser._mmcif_dict  # << This is a little dangerous but the get_structure code above looks quite clean

        # Gather resolution information, etc, from the
        structure_info_dict = mine_structure_info_from_mmCIF(mmcif_dict)
    elif structure is None:
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure(structure_id, structure_fin)
        # If we ever return to deposited PDBs for ddg calculations (we now use mmcif), then this has to
        # be worked out.  For now, the only .pdb files we are processing are models that will not have
        # spurious residues that need to be cleaned out
        mmcif_dict = None

    return structure_id,structure,structure_info_dict,mmcif_dict
