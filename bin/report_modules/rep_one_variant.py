#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# Filename       : rep_one_variant.py
# Authors        : Chris Moth and R. Michael Sivley
# project Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2023 Januar 27
# Description    : The function report_one_variant_one_isoform, excerpted from psb_rep.py,
#                : is called by psb_rep.py to create the webpage (and support files)
#                : for a single variant.  It is complex, in and of itself, and
#                : with the 2023 introduction of the web interface, has been factored
#                : into its own module
# =============================================================================#

import os
import sys
import pprint
from typing import Dict, List
import time
import logging
import numpy as np
import json
from jinja2 import Environment, FileSystemLoader
from lib import PDBMapProtein
from lib import PDBMapGlobals
from lib import PDBMapComplex
from lib import PDBMapTranscriptUniprot
from lib import PDBMapAlphaMissense
from logging.handlers import RotatingFileHandler

from .pfam_graphic import PfamDomainGraphics

from .load_one_variant import CalculationResultsLoader

LOGGER = logging.getLogger(__name__)



def _initialize_local_logging_in_variant_directory(
        case_root_dir: str,
        variant_directory: str,
        local_log_level: str) -> RotatingFileHandler:

    log_filename = os.path.join(CalculationResultsLoader.case_root_dir_variant_directory_segment(
        case_root_dir,
        variant_directory),
            "psb_rep.log")
    LOGGER.info("additionally logging to file %s" % log_filename)
    _need_roll = os.path.isfile(log_filename)
    local_fh = RotatingFileHandler(log_filename, backupCount=7)
    formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',
                                  datefmt="%H:%M:%S")
    local_fh.setFormatter(formatter)
    if local_log_level.lower() == "debug":
        local_fh.setLevel(logging.DEBUG)
    elif local_log_level.lower() == "info":
        local_fh.setLevel(logging.INFO)
    else:
        local_fh.setLevel(logging.WARN)

    LOGGER.addHandler(local_fh)
    if _need_roll:
        local_fh.doRollover()

    return local_fh  # Be sure to _end_local_logging(local_fh) on this later


def _end_local_logging(local_fh: RotatingFileHandler):
    local_fh.flush()
    local_fh.close()
    LOGGER.removeHandler(local_fh)
    local_fh = None



def report_one_variant_one_isoform(project: str,
                                   case_root_dir: str,
                                   variant_directory_segment: str,
                                   parent_report_row: Dict,
                                   local_log_level: str,
                                   vustruct_pipeline_launched: bool,
                                   config_pathprox_dict: dict
                                   ) -> ( Dict, Dict, List[str] ):
    """
    Generate extensive structure analysis report .html for a single variant of a single transcript by
    reading all PathPRox/ddG/etc calculations into dataframes, transforming the data a bit, and then
    writing out into .html via jinja2 templates.

    Variant_isoform summary data are returned to the caller in an extensive dictionary, which the caller
    uses to fill in tables on the main case-wide page.

    @case_root_dir: the parent directory for all the variant-specific directory, including the variant directory we
             are building a report for in this function.

    @param variant_directory: Typically a string like "Gene_transcriptId_aaVariant" that is the root subdirectory
            of the case directory that houses all calculation results.  The string is also a prefix to
             key files that record calculation status

    @param parent_report_row: In the event that no calculations were launched for the variant, then
           the unp from the parent row (on the main page report) populates variables so a minimal report can be
           completed

    @param local_log_level: If defined as "debug", "info", "warm", "error", or "critical", then an _Additional_
             log file is create in the variant_directory, with the entries specific to ONLY the variant for which
             this function is building a report.  This typically is in addition to the main psb_rep.py level
             log file which includes all logging for all variants.

    @param vustruct_pipeline_launched: IF it is very early, and the pipeline has NOT been launched,
            THEN the variant report will include
            only top-level information about the transcript context, and the tabke of workstatus will not be
            consoluted to mine PathProx, ddG, etc results.

    @return: 1 of 2) An variant_isoform_summary dictionary that can help the caller assemble the main report page.
             or {} if no caluclation results are available for any reason
             2 of 2) A list of the files that must be added to the final website tarball, by the caller.
    """

    # Start logging - after which logging entries are mirrored to variant_directory/psb_rep.log
    variant_directory_fullpath = CalculationResultsLoader.case_root_dir_variant_directory_segment(
        case_root_dir,
        variant_directory_segment)

    if not os.path.exists(variant_directory_fullpath):
        LOGGER.warning("The variant directory has not even been created via psb_plan yet,")
        return {}, {}, []

    local_logger_fh = _initialize_local_logging_in_variant_directory(
        case_root_dir,
        variant_directory_fullpath,
        local_log_level)

    variant = parent_report_row['mutation']
    gene = parent_report_row['gene']
    uniprot_id = None
    unp_transcript = None
    alphamissense_score = None

    if 'unp' in parent_report_row and parent_report_row['unp']:
        uniprot_id = parent_report_row['unp']
        unp_transcript = PDBMapTranscriptUniprot(uniprot_id)
        if PDBMapProtein.isCanonicalByUniparc(uniprot_id):
            alphamissense_score = PDBMapAlphaMissense.score_from_uniprot_id(uniprot_id,variant)
            LOGGER.info("Alphamissense for %s %s = %s" % (uniprot_id, variant, alphamissense_score))
    if alphamissense_score is None: # Try to use ENST transcripts to get alphafold
        # If the transcript has NO ensembl transcript IDs, then do not attempt this method
        if parent_report_row['transcript']:
            enst_transcript_ids = parent_report_row['transcript'].split(';')
            for transcript_id in enst_transcript_ids:
                alphamissense_score = PDBMapAlphaMissense.score_from_ENSEMBL_isoform_id(transcript_id, variant)
                LOGGER.info("Alphamissense for %s %s = %s" % (transcript_id, variant, alphamissense_score))
                if alphamissense_score is not None:
                    break
    
    if alphamissense_score is not None:            
        alphamissense_score = np.round(alphamissense_score,2)

       



    website_filelist = []

    # The rate4site score for the variant is shown in a table of transcripts near top, if rate4site data is found
    rate4site_dict = {}
    if 'transcript' in parent_report_row and parent_report_row['transcript']:
        ensembl_transcript_ids =  parent_report_row['transcript'].split(';')
        for ensembl_transcript_id in ensembl_transcript_ids:
            rate4site_norm_rates_filename = os.path.join(PDBMapGlobals.config['rate4site_dir'],
                                                     "%s_norm_rates.txt" % ensembl_transcript_id.split('.')[0])
            if os.path.exists(rate4site_norm_rates_filename):
                rate4site_df, alpha_parameter, average, std = PDBMapComplex._load_one_rate4site_file(
                    rate4site_norm_rates_filename)
                variant_pos = int(parent_report_row['mutation'][1:-1])
                # Now grab just th dataframe row of all the rate4site entries for this position
                rate4site_dict[ensembl_transcript_id] = rate4site_df[rate4site_df['pos'] == variant_pos].iloc[0].to_dict()
            else:
                LOGGER.warning("No rate4site results were found for %s: %s",
                               parent_report_row['unp'], ensembl_transcript_id)



    cosmis_dict = {}
    if uniprot_id:
        # _load_one_cosmis_set below will return empty if non-canonical
        # which is a very good thing
        cosmis_df = PDBMapComplex._load_one_cosmis_set(uniprot_id)
        if cosmis_df.empty:
            LOGGER.warning("No cosmis results were found for %s", uniprot_id)
        else:
            variant_pos = int(parent_report_row['mutation'][1:-1])
            # Now grab just th dataframe row of all the rate4site entries for this position
            cosmis_uniprot_pos_df = cosmis_df[cosmis_df['uniprot_pos'] == variant_pos]
            if cosmis_uniprot_pos_df.empty:
                LOGGER.warning("Cosmis results were found for %s, not not for position %d", uniprot_id, variant_pos)
            else:
                cosmis_dict[uniprot_id] = cosmis_uniprot_pos_df.iloc[0].to_dict()

    # If the psb_launch.py application was run, THEN we should have
    # workstatus csv tables and we can load these to make a much more complete website.
    # Though now, we need to think through the progress of jobs and how to display that
    # progress
    calculation_results_loader = {}
    if vustruct_pipeline_launched:
        calculation_results_loader = CalculationResultsLoader(
            case_root_dir,
            variant_directory_segment, 
            parent_report_row,
            config_pathprox_dict)
        calculation_results_loader.load_dataframes()
        additional_website_files = calculation_results_loader.load_structure_graphics_dicts()
        website_filelist.extend(additional_website_files)

    # Various data max/min/etc computations from the gathered data, which are great to return to the caller
    # and also useful to compose the report.
    variant_isoform_summary = {}
    variant_isoform_details = {} # Carries specific per-structure calculations to pass to caller

    # Load the template psb_report.html that is one directory above where this file is sourced from
    jinja2_environment = Environment(loader=FileSystemLoader(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
    variant_isoform_template = jinja2_environment.get_template("psb_report.html")
    structure_iframe_template = jinja2_environment.get_template("html/structureIframeTemplate.html")

    calculation_results_loader_as_dict = {}
    if calculation_results_loader:  # True only if pipeline got as far as 'launch'
        calculation_results_loader_as_dict = {key: value for key, value in vars(calculation_results_loader).items() \
                                              if not (
                key.startswith('_'))}  # -> don't exclude callable, you lose properties callable(value)))},

        for attribute in ['gene', 'unp', 'variant']: # No longer needed , 'project' and 'unp_transcript'] populated below
            calculation_results_loader_as_dict[attribute] = getattr(calculation_results_loader, attribute)
         

        calculation_results_loader_as_dict['unp_transcript'] = unp_transcript
        calculation_results_loader_as_dict['project'] =project

    # Provide the file name now.  We'll create the file's contents
    # (the domain graphics) below.  However, this filename ONLY makes sense
    # if the calculations have been launched....
    if 'unp' in parent_report_row and parent_report_row['unp']:
        pfamGraphicsIframe_fname = "%s_%s_PfamGraphicsIframe.html" % (
            parent_report_row['unp'],
            variant)

    template_render_dict = {
        # Convert calculation_results to a dictionary because it is infinitely more compatable with
        # jinja2 template
        'calculation_results': calculation_results_loader_as_dict,
        # 'structure_report_df_appended': structure_report_df_appended,
        "today_date": time.strftime("%Y-%m-%d"),
        "pfamGraphicsIframe_fname": pfamGraphicsIframe_fname,
        # html template references disease variant strings like clinvar/cosmic...
        'config_pathprox_dict': config_pathprox_dict,
        'rate4site_dict': rate4site_dict,
        'cosmis_dict': cosmis_dict
    }


    # Write out the individual structure.html files
    # website_filelist.append('color_chains_consecutively.js')
    if 'structure_graphics_dicts' in calculation_results_loader_as_dict:
        for structure_graphics_dict in calculation_results_loader_as_dict['structure_graphics_dicts']:
            LOGGER.info("Creating NGL Viewer for %s", structure_graphics_dict['html_div_id'])
            structureIframe_out = structure_iframe_template.render(
                {'structure': structure_graphics_dict,
                 'config_pathprox_dict': config_pathprox_dict}
            )
            structureIframe_base_filename = os.path.join(variant_directory_segment,
                                                         "%s_viewer.html" % structure_graphics_dict['html_div_id'])
            website_filelist.append(structureIframe_base_filename)
            with open(structureIframe_base_filename, "w") as html_f:
                html_f.write(structureIframe_out)

    if LOGGER.isEnabledFor(logging.DEBUG):
        pp = pprint.PrettyPrinter(indent=1)
        LOGGER.debug("Dictionary to render:\n%s" % pp.pformat(template_render_dict))

    template_render_dict['variant'] = variant
    template_render_dict['uniprot_id'] = uniprot_id
    template_render_dict['alphamissense_score'] = alphamissense_score

    variant_report_html_out = variant_isoform_template.render(template_render_dict)

    variant_report_base_filename = os.path.join(variant_directory_segment,
                                                "%s_report" % (variant_directory_segment))

    variant_isoform_summary['html_filename'] = variant_report_base_filename + ".html"
    website_filelist.append(variant_isoform_summary['html_filename'])
    with open(variant_isoform_summary['html_filename'], "w") as html_f:
        html_f.write(variant_report_html_out)

    pfg = PfamDomainGraphics(uniprot_id, variant)
    pfam_graphics_legend, pfam_domain_graphics_dict = pfg.create_graphics_legend_and_graphics_dict()
    pfam_domain_graphics_dict = pfg.add_our_variant_to_graphics_dict(pfam_domain_graphics_dict)
    pfam_domain_graphics_dict = pfg.update_graphics_dict_hrefs_to_point_to_pfam(pfam_domain_graphics_dict)
    template_render_dict["PfamResultTableHTML"] = pfg.create_pfam_html(pfam_graphics_legend, pfam_domain_graphics_dict)

    # Encapsulate in a list
    template_render_dict["DomainGraphicsJSON"] = '[' + json.dumps(pfam_domain_graphics_dict) + ']'

    template_pfam_graphics = jinja2_environment.get_template("html/pfamGraphicsIframeTemplate.html")
    pfam_graphics_html_full_iframe = template_pfam_graphics.render(template_render_dict)

    # The variant_report refers to the Pfam graphics without a directory indicator.
    # So we prefix the filename, here, with the report directory to write it.
    pfamGraphicsIframe_fname_with_directory = os.path.join(variant_directory_segment, pfamGraphicsIframe_fname)
    website_filelist.append(pfamGraphicsIframe_fname_with_directory)
    with open(pfamGraphicsIframe_fname_with_directory, "w") as html_f:
        html_f.write(pfam_graphics_html_full_iframe)

    # With the report fully written to disk, now populate summary details for return to caller.
    variant_isoform_summary['AA_seq'] = unp_transcript.aa_seq
    variant_isoform_summary['AA_len'] = len(variant_isoform_summary['AA_seq'])

    # gathered_info['gene_id']  PDBMapProtein.unp2gene_id(gethered_info['unp'])

    # template_vars = gathered_info.copy()

    # mutation = row0['mutation']  # Sorry for alias - but this is all over the code

    variant_isoform_summary['Error'] = ''  # < This is False for truth but looks fine on the report
    variant_isoform_summary['ddG Monomer Max'] = variant_isoform_summary['ddG Monomer Min'] = None
    variant_isoform_summary['ddG Cartesian Max'] = variant_isoform_summary['ddG Cartesian Min'] = None
    variant_isoform_summary['disease1_pp Min'] = variant_isoform_summary['disease1_pp Max'] = None
    variant_isoform_summary['disease2_pp Min'] = variant_isoform_summary['disease2_pp Max'] = None
    variant_isoform_summary['MusiteiDeepPTM'] = None
    variant_isoform_summary['ScanNetPPI'] = None

    variant_isoform_details['ddG Monomer'] = {}
    variant_isoform_details['ddG Cartesian'] = {}
   

    if vustruct_pipeline_launched and calculation_results_loader.ddG_monomer_results_dict_of_dfs:
        ddG_list = [calculation_results_loader.ddG_monomer_results_dict_of_dfs[ddg_results_key].ddG \
                    for ddg_results_key in calculation_results_loader.ddG_monomer_results_dict_of_dfs]
        variant_isoform_summary['ddG Monomer Max'] = max(ddG_list)
        variant_isoform_summary['ddG Monomer Min'] = min(ddG_list)

        # Added 2024 May to return individual structure results to include on the final spreadsheet
        _ddG_details = {ddg_results_key: 
            calculation_results_loader.ddG_monomer_results_dict_of_dfs[ddg_results_key].ddG \
                    for ddg_results_key in calculation_results_loader.ddG_monomer_results_dict_of_dfs}
        variant_isoform_details['ddG Monomer'] = _ddG_details
         

    if vustruct_pipeline_launched and calculation_results_loader.ddG_cartesian_results_dict_of_dfs:
        ddG_list = [calculation_results_loader.ddG_cartesian_results_dict_of_dfs[ddg_results_key].total \
                    for ddg_results_key in calculation_results_loader.ddG_cartesian_results_dict_of_dfs]
        variant_isoform_summary['ddG Cartesian Max'] = max(ddG_list)
        variant_isoform_summary['ddG Cartesian Min'] = min(ddG_list)
        # Added 2024 May to return individual structure results to include on the final spreadsheet
        _ddG_details = {ddg_results_key: 
            calculation_results_loader.ddG_cartesian_results_dict_of_dfs[ddg_results_key].total \
                    for ddg_results_key in calculation_results_loader.ddG_cartesian_results_dict_of_dfs}
        variant_isoform_details['ddG Cartesian'] = _ddG_details

    if vustruct_pipeline_launched and calculation_results_loader.scannet_prediction_dict:
        variant_isoform_summary['ScanNetPPI'] = calculation_results_loader.scannet_prediction_dict
    if vustruct_pipeline_launched and calculation_results_loader.musite_deep_neighborhood_dict:
        variant_isoform_summary['MusiteDeepPTM'] = calculation_results_loader.musite_deep_neighborhood_dict

    def nan_to_None(x: np.float64):
        return None if np.isnan(x) else x

    if vustruct_pipeline_launched and calculation_results_loader.pathprox_disease1_results_dict_of_dfs:
        pathprox_list = [
            calculation_results_loader.pathprox_disease1_results_dict_of_dfs[pathprox_results_key].iloc[0].pathprox \
            for pathprox_results_key in calculation_results_loader.pathprox_disease1_results_dict_of_dfs]
        # Get rid of Nans
        variant_isoform_summary['disease1_pp Max'] = nan_to_None(np.nanmax(np.array(pathprox_list)))
        variant_isoform_summary['disease1_pp Min'] = nan_to_None(np.nanmin(np.array(pathprox_list)))

    if vustruct_pipeline_launched and calculation_results_loader.pathprox_disease2_results_dict_of_dfs:
        pathprox_list = [
            calculation_results_loader.pathprox_disease2_results_dict_of_dfs[pathprox_results_key].iloc[0].pathprox \
            for pathprox_results_key in calculation_results_loader.pathprox_disease2_results_dict_of_dfs]

        variant_isoform_summary['disease2_pp Max'] = nan_to_None(np.nanmax(np.array(pathprox_list)))
        variant_isoform_summary['disease2_pp Min'] = nan_to_None(np.nanmin(np.array(pathprox_list)))

    if alphamissense_score: 
        variant_isoform_summary['alphamissense_score'] = alphamissense_score 

    return variant_isoform_summary, variant_isoform_details, website_filelist

    """jinja2_environment = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)))))
    isoform_variant_template = jinja2_environment.get_template("psb_report.html")

    if LOGGER.isEnabledFor(logging.DEBUG):
        pp = pprint.PrettyPrinter(indent=1)
        LOGGER.debug("Dictionary template_vars to render:\n%s", pp.pformat(template_vars))

    html_out = isoform_variant_template.render(template_vars)

    templatePfamGraphics = jinja2_environment.get_template("html/pfamGraphicsIframeTemplate.html")
    htmlPfamGraphics = templatePfamGraphics.render(template_vars)

    save_cwd = os.getcwd()
    os.chdir(variant_report_directory)

    # WE ARE NOW OPERATING FROM ..../UDN/CaseName target directory
    LOGGER.info("Now working in %s", variant_report_directory)

    html_fname = base_fname % "html"
    website_filelist.append(os.path.join(web_dir, html_fname))
    with open(html_fname, "w") as html_f:
        html_f.write(html_out)

    with open(pfamGraphicsIframe_fname, "w") as html_f:
        html_f.write(htmlPfamGraphics)
    """

    """write_pdfs = None
    if write_pdfs:
        pdf_fname = base_fname % "pdf"
        wkhtmltopdf_fname = base_fname % "wkhtml.pdf"
        print("\nWriting final reports to %s directory:\n  %s\n  %s\n  %s\n  %s\n" % (
            variant_report_directory, pdf_fname, wkhtmltopdf_fname, html_fname, pfamGraphicsIframe_fname))

        LOGGER.warning("Temporarily disabling all logging prior to calling weasyprint write_pdf()")

        logging.disable(sys.maxsize)
        HTML(string=html_out).write_pdf(pdf_fname, stylesheets=[
            "../html/css/typography.css"])  # ,CSS(string="@page { size: A3 landscape; }")])
        logging.disable(logging.NOTSET)
        LOGGER.warning("weasyprint write_pdf() has returned.  Restoring logging")

        # Write out another .pdf using the wkhtmltopdf tool

        wkhtml_command_list = ["wkhtmltopdf", "--print-media-type", "--no-outline", "--minimum-font-size", "12",
                               html_fname, wkhtmltopdf_fname]
        process = subprocess.Popen(wkhtml_command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (output, err) = process.communicate()

        # So unfortunately, the stderr from wkhtmltopdf is full of ==================
        # things, which are a mess for the LOGGER.
        # Get rid of stuff that ends in \r  Keep all other lines and ending in \r
        err_nl_array = err.split('\n')
        err_legit_list = []
        for ena in err_nl_array:
            crs = ena.split('\r')
            if len(crs) > 0:  ## Pick off the last string, which means all ending with Carriage Return are ignored
                cr_last = crs[-1].strip()
                if cr_last:
                    err_legit_list.append(cr_last)

        err_legit = '\n'.join(err_legit_list)

        # ContentNotFound really just info - not warning
        if (process.returncode != 0) and ("Exit with code 1 due to network error: ContentNotFoundError" not in err):
            print("Unable to complete %s exitstatus=%d due to error %s\n  output: %s\n" % (
                str(wkhtml_command_list), process.returncode, err_legit, output))
            LOGGER.warning("wkhtmltopdf stderr:\n%s", err_legit)
        else:
            LOGGER.info("wkhtmltopdf stderr:\n%s", err_legit)
    """

    """
    # WE HAVE NOW RETURNED TO the psb_pipeline/bin directory
    os.chdir(save_cwd)
    # Close out the local log file for this mutation
    _end_local_logging(local_logger_fh)
    gathered_info['variant_report_directory'] = os.path.basename(variant_report_directory)
    gathered_info['html_fname'] = html_fname
    return gathered_info  # End of function report_one_variant_one_isoform()
    """

if __name__ == '__main__':
    print("Only for use as library")
    sys.exit(1)
