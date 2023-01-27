#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# Filename       : psb_rep.py
# Authors        : Chris Moth and R. Michael Sivley
# project Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2018-01-22
# Description    : Output final html and pdf reports from completed jobs
#                : (output from psb_launch.py)
# =============================================================================#

from typing import Dict
import os
import pycurl
from io import  BytesIO
import json
from xml2json import Cxml2json
from lib import PDBMapGlobals

import logging
LOGGER = logging.getLogger(__name__)

class PfamDomainGraphics:
    """
    Manages loading of graphics in form of received json string from pfam.xfam.org.
    Important: In cases where pfam.xfam.org lacks JSON graphics, or fails to return
    a graphic in timely fashion, we build the JSON ourselves.

    I've not documented the format because you can see it online:
    http://pfam.xfam.org/protein/P45381/graphic

    with examples of how to display the graphic and the table format, not far from that.  Ex:
    http://pfam.xfam.org/protein/P45381/
    """

    def __init__(self, uniprot_id, variant):
        """
        Download (or build) a PFAM-graphics-library-compatiable dictionary of predicted protein domain
        annotations, so that pipeline case users can see a linear graphic layout of each protein

        @param uniprot_id: The canonical, or non-canonical uniprot identifier for the protein to diagram
        @param variant:    Our variant of interest, to add to the graphics and the table beneath.
        """
        self.unp = uniprot_id
        self.variant = variant

        assert self.unp is not None

    def _fetch_canonical_pfam_graphic_json_from_xfam(self, timeout_seconds=60) -> str:
        """
        Fetch a canonical PFAM graphic JSON string from the internet, using curl, for the
        unp passed to the class constructor.

        @return: PFAM graphic JSON string or None if pfam.xfam server down or bad data returned

        Requires config_dict['interpro_dir'] to be defined
        Inside that directory must be 'match_humanonly.xml' file
        """

        # xfam lacks graphics URLs for non-canonical isoforms.  We assemble those ourselves, elsewhere.
        canonical_pfam_xfam_url = 'http://pfam.xfam.org/protein/{0}/graphic'.format(
            self.unp.split(' ')[0].split('-')[0])

        LOGGER.debug("Attempting curl connection to %s to fetch canonical pfam graphics.  Timeout=%d\n",
                     canonical_pfam_xfam_url, timeout_seconds)
        buffer = BytesIO()

        try:
            c = pycurl.Curl()
            c.setopt(pycurl.CONNECTTIMEOUT, timeout_seconds)
            c.setopt(pycurl.TIMEOUT, timeout_seconds)
            # c.setopt(pycurl.WRITEFUNCTION, lambda x: None)
            c.setopt(pycurl.WRITEFUNCTION, buffer.write)
            c.setopt(pycurl.HEADERFUNCTION, lambda x: None)
            c.setopt(pycurl.NOSIGNAL, 1)
            c.setopt(pycurl.URL, canonical_pfam_xfam_url)
            c.setopt(pycurl.HTTPGET, 1)

            c.perform()
            c.close()

        except Exception as ex:
            LOGGER.exception('Failed to get JSON string from %s: %s' % (canonical_pfam_xfam_url, str(ex)))
            return None

        domain_graphics_json = buffer.getvalue().decode('latin')

        if len(domain_graphics_json) > 3 and domain_graphics_json[0] == '[':
            domain_graphics_json = domain_graphics_json[1:-1]
        elif domain_graphics_json[0] == '<':
            # Sometimes we get back html in the response, to tell us that the PFAM website is down
            LOGGER.warning("XML, not json, was returned from pfam:\n%s", domain_graphics_json)
            return None
        LOGGER.debug("pfam sent us:\n%s", domain_graphics_json)
        return domain_graphics_json

    def create_graphics_legend_and_graphics_dict(self) -> (str, Dict):
        """
        From self.unp, create a graphics Legend, initialize self._domain_graphics_dict.
        This involves opening a large xml file, and either getting the JSON from the internet, or constructing it manually
        
        @return: domain_graphics_legend,domain_graphics_dict
        """
        domain_graphics_legend = "Pfam Domain Graphic unavailable"
        # Get Uniprot Graphics string or put in a blank
        domain_graphics_json = None
        if self.unp:  # Fetch the JSON - wait up to 30 seconds
            interpro_xml_handle = Cxml2json(os.path.join(PDBMapGlobals.config['interpro_dir'], 'match_humanonly.xml'))
            if interpro_xml_handle.isCanonical(self.unp):  # Then we need to manually create the xml ourselves
                # Go for the web link because this is a canonical UNP - however that _can_ fail.
                domain_graphics_legend = "Downloaded Pfam Domain Graphic for Canonical Isoform %s" % self.unp
                LOGGER.info("Attempting download of Domain Graphics for %s from xfam.pfam" % self.unp)
                domain_graphics_json = self._fetch_canonical_pfam_graphic_json_from_xfam(30)
                if not domain_graphics_json:
                    LOGGER.info("Download returned nothing")
                    domain_graphics_legend = "Pfam Domain Graphic for Canonical Isoform %s" % self.unp
            else:
                domain_graphics_legend = "Pfam Domain Graphic for Non-canonical Isoform %s" % self.unp

            # _Either_ communications failure OR non-canonical isoform
            # So we create our own graphic from our local xml database
            if not domain_graphics_json:  # _Either_ communications failure OR non-canonical (no communication)
                LOGGER.info("Creating Domain Graphic for %s from xml", self.unp)
                nonCanonicalGraphicsJSON = interpro_xml_handle.PFAMgraphicsJSONfromUnp(self.unp)
                domain_graphics_json = json.dumps(nonCanonicalGraphicsJSON)

        domain_graphics_dict = json.loads(domain_graphics_json)

        return domain_graphics_legend, domain_graphics_dict

    def add_our_variant_to_graphics_dict(self, domain_graphics_dict: Dict) -> Dict:
        """
        Specifically add _our_ variant of interest to the pfam graphics, as a table entry, and as a diamond

        @return: reference to same domain_graphics_dict, with the ['markups'] altered.
        """
        # Add our variant point of interest to the PFAM-like domain graphics.
        # as a simple diamond in the graphics, and an entry in the table.
        variant_site_markup_dict = {'colour': '#e469fe',
                                    'display': True,
                                    'headStyle': 'diamond',
                                    'lineColour': '#333333',
                                    'metadata': {'database': 'UDN Case',
                                                 'description': '%s' % self.variant,
                                                 'start': self.variant[1:-1],
                                                 'type': 'Mutation'},
                                    'residue': 'X',
                                    'start': self.variant[1:-1],
                                    'type': 'UDN Mutation site',
                                    'v_align': 'top'}

        # Grab any _existing_ markups in our domain graphics.
        # Else default to an empty list.  Add _our_ variant markup to that list of markups.
        domain_graphics_dict['markups'] = domain_graphics_dict.get('markups', []) + [variant_site_markup_dict]

        # Now repackage in a list, so the "domain_graphics_json" is back to original form....
        return domain_graphics_dict

    def update_graphics_dict_hrefs_to_point_to_pfam(selfself, domain_graphics_dict: Dict) -> Dict:
        """
        Downloaded href's embedded in the graphics dictionary assume that the domains are described locally.
        We must convert them to point to pages on the UK webserver.

        @param domain_graphics_dict:
        @return:
        """

        for region in domain_graphics_dict.get('regions', []):
            if 'href' in region:
                region['href'] = 'http://pfam.xfam.org' + region['href']

        return domain_graphics_dict

    def create_pfam_html(self, doman_graphics_legend: str, domain_graphics_dict: Dict) -> str:
        """
        @param doman_graphics_legend: 
        @param domain_graphics_json: 
        @return html string:
        """

        # Build a table to explaine the pfam graphics
        PfamResultTableHTML = '''<table class="resultTable details"
                 id="imageKey"
                 summary="Key for the Pfam domain image">
            <thead>
              <tr>
                <th class="dh" rowspan="2">Source</th>
                <th class="dh" rowspan="2">Domain</th>
                <th class="dh" rowspan="2">Start</th>
                <th class="dh" rowspan="2">End</th>
                <th class="sh" style="display: none" colspan="2">Gathering threshold (bits)</th>
                <th class="sh" style="display: none" colspan="2">Score (bits)</th>
                <th class="sh" style="display: none" colspan="2">E-value</th>
              </tr>
              <tr>
                <th class="sh" style="display: none">Sequence</th>
                <th class="sh" style="display: none">Domain</th>
                <th class="sh" style="display: none">Sequence</th>
                <th class="sh" style="display: none">Domain</th>
                <th class="sh" style="display: none">Sequence</th>
                <th class="sh" style="display: none">Domain</th>
              </tr>
            </thead>'''

        PfamResultTableHTML += '<tbody>'

        table_rows = []
        for motif in domain_graphics_dict.get('motifs', []):
            table_rows.append({
                'type': motif['metadata']['type'],
                'domain': 'n/a',
                'start': motif['start'],
                'end': motif['end']})

        for region in domain_graphics_dict.get('regions', []):
            # It's nice to say something about a domain when we can
            # However not all subsequences of a protein is a domain of course
            domain_table_entry = ''
            if 'text' in region:
                region_text = region['text']
            else:
                region_text = 'Text Missing'

            if 'href' in region:
                domain_table_entry = "<a href=" + region['href'] + ">" + region_text + "</a>"
            else:
                domain_table_entry = region_text
            table_rows.append({
                'type': region['metadata']['type'],
                'domain': domain_table_entry,
                'start': region['start'],
                'end': region['end']})

        table_rows.append({
            'type': 'Variant',
            'domain': 'n/a',
            'start': self.variant[1:-1],
            'end': self.variant[1:-1]})

        sorted_table_rows = sorted(table_rows, key=lambda d: int(d['start']))

        for sr in sorted_table_rows:
            table_row = '<tr class="odd">'
            table_row += '<td class="domain">%s</td>' % sr['type']
            table_row += '<td><span class="inactive">%s</span></td>' % sr['domain']
            table_row += '<td>%s</td>' % sr['start']
            table_row += '<td>%s</td>' % sr['end']
            table_row += 6 * '<td class="sh" style="display: none"><span class="inactive">n/a</span></td>'
            table_row += '</tr>\n'
            PfamResultTableHTML += table_row

        PfamResultTableHTML += '</tbody></table>'

        return PfamResultTableHTML

