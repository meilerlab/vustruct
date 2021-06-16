#!/usr/bin/env python
from configparser import SafeConfigParser

def ddg_commandline_parser(commandline_parser: SafeConfigParser) -> SafeConfigParser:
    """
    Initialize the ddg command line with arguments shared by ddg_monomer, ddg_cartesian, and ddg_relax programs.
    @param command_line_parser:
    @return:
    """

    commandline_parser.add_argument(
        "--ddg_config",
        help="ddG Configuration File specifying binary programs, parameters, and rep location",
        required=False, metavar="FILE")

    # Input parameters
    group = commandline_parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pdb', type=str,
                       help="4 character PDB ID with optional .chain suffix")
    group.add_argument('--biounit', type=str,
                       help="4 character PDB ID with optional .chain suffix")
    group.add_argument('--modbase', type=str,
                       help="Modbase 13 or Modbase 16 model ID with optional .chain suffix")
    group.add_argument('--swiss', type=str,
                       help="Swissmodel ID with optional .chain suffix")
    group.add_argument('--usermodel', type=str, metavar='FILE',
                       help="Filename of a user model.  Requires explicit transcript specifications")
    # cmdline_parser.add_argument("entity",type=str,
    #                   help="Gene ID, UniProt AC, PDB ID, or PDB filename")

    commandline_parser.add_argument("--chain", type=str,
                                help="Limit the ddG processing to one particular chain")

    commandline_parser.add_argument(
                "-o", "--UDNoutdir", type=str,
                help="Optional directory to echo output and results back to UDN pipeline")
    commandline_parser.add_argument(
                "--UDNuniquekey", type=str, required=False,
                help="Optional gene/refseq/mutation/structure/chain/flavor unique identifer for this pipeline ddG run")
    commandline_parser.add_argument(
                "--variant", type=str,
                help="Amino Acid Variant in modified HGVS format.  PDB insertion codes can be added to residue #")
    commandline_parser.add_argument(
                "--species_filter", type=str, required=False, default=None,
                help="Optionally retain only residues with DBREF of, example 'HUMAN'")

    commandline_parser.add_argument("--label", type=str, default='',
                                help="Optional analysis label (overrides entity inference)")

    return commandline_parser
# cmdline_parser.add_argument("--no-timestamp", "-nt", action="store_true", default=False,
#                            help="Disables output directory timestamping")
# cmdline_parser.add_argument("--overwrite", action="store_true", default=False,
#                             help="Overwrite previous results. Otherwise, exit with error")
