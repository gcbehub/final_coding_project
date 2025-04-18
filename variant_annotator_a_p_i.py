# Section 1. Import and Logging Setup
import os   # library required for operating sustem functionality
import logging  # library required for logging functionality
import requests # library required for making HTTP requests
from flask import Flask, jsonify, request   # web framework components
from flask_restful import Resource, Api # RESTful API support
import json # library required for JSON handling

# Set up logging configuaration
current_directory = os.path.dirname(os.path.abspath(__file__))  # gets current file's directory
log_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')   # defines log format
info_file_handler = logging.FileHandler(os.path.join(current_directory, 'info.log'))    # creates log file handler
info_file_handler.setFormatter(log_formatter)   # applies formatter to handler

logger = logging.getLogger(__name__)    # creates logger instance
logger.setLevel(logging.INFO)   # sets logging level
logger.addHandler(info_file_handler)    # adds file handler to logger

# Section 2. Main VariantAnnotatorSprint Class
class VariantAnnotatorSprint:
    """
    Handles variant annotation and validation processes using various external
    databases such as ClinVar and Ensembl.

    This class provides methods to fetch, validate, and process data from external
    sources for genetic variant information. The purpose of the class is to
    facilitate interaction with APIs related to genetic variants, offering tools
    to validate variants or access classification details.

    :ivar base_url: Base URL for NCBI variation API.
    :type base_url: str
    :ivar clinvar_base_url: Base URL for ClinVar database.
    :type clinvar_base_url: str
    :ivar ensembl_base_url: Base URL for Ensembl variation API (GRCh38).
    :type ensembl_base_url: str
    :ivar grch37_base_url: Base URL for Ensembl variation API (GRCh37).
    :type grch37_base_url: str
    """

    # 2-0 Main Class Initialization
    def __init__(self): # initializes base URLs for different services
        self.base_url = "https://api.ncbi.nlm.nih.gov/variation/v0/"    # NCBI API
        self.clinvar_base_url = "https://www.ncbi.nlm.nih.gov/clinvar/" # ClinVar database
        self.ensembl_base_url = "https://rest.ensembl.org"  # Ensembl API (GRCh38)
        self.grch37_base_url = "https://grch37.rest.ensembl.org"    # Ensembl API (GRCh37)

    # 2-1. Key method 1: HTTP Request Handler
    def _make_validation_request(self, url, params=None, headers=None): # private method to handle HTTP requests
        try:
            if headers is None:
                headers = {"Content-Type": "application/json"}
            response = requests.get(url, params=params, headers=headers)
            if response.status_code == 200:
                return response.json()
            logger.warning(f"Request failed with status code {response.status_code}")
            return None
        except Exception as e:
            logger.error(f"Error making validation request: {str(e)}")
            return None

    # 2-2. Key method 2: ClinVar Data Methods
    def get_clinvar_data(self, variant_id): # fetch data from ClinVar databases
        try:
            url = f"{self.clinvar_base_url}{variant_id}"
            response = requests.get(url)
            if response.status_code == 200:
                return response.text
            return None
        except Exception as e:
            logger.error(f"Error fetching ClinVar data: {str(e)}")
            return None

    def extract_classification(self, clinvar_data):
        """
        Extract clinical significance classification from ClinVar data.

        Args:
            clinvar_data (str): HTML or JSON response from ClinVar

        Returns:
            dict: Dictionary containing the classification, or None if extraction fails
        """
        try:
            if not clinvar_data:
                logger.warning("No ClinVar data provided")
                return None

            # If the response is in JSON format
            if isinstance(clinvar_data, dict):
                # Try to get classification from JSON structure
                if 'clinical_significance' in clinvar_data:
                    return {
                        "classification": clinvar_data['clinical_significance']
                    }

            # If the response is HTML (common for ClinVar web responses)
            elif isinstance(clinvar_data, str):
                # Look for common classification terms
                classification_terms = [
                    'Pathogenic',
                    'Likely pathogenic',
                    'Uncertain significance',
                    'Likely benign',
                    'Benign',
                    'Conflicting interpretations',
                    'Not provided'
                ]

                # Convert to lowercase for case-insensitive comparison
                clinvar_data_lower = clinvar_data.lower()

                # Find the first matching classification
                for term in classification_terms:
                    if term.lower() in clinvar_data_lower:
                        return {
                            "classification": term
                        }

            logger.warning("Could not find classification in ClinVar data")
            return {"classification": "Unknown"}

        except Exception as e:
            logger.error(f"Error extracting classification: {str(e)}")
            return None

    # 2-3. Key method 3: Variant Validation Methods (RefSeq and Ensembl databases)
    def validate_variant_refseq(self, variant_id, genome_build="GRCh38"):
        try:
            if not validate_genome_build(genome_build):
                return {"error": f"Invalid genome build: {genome_build}"}

            url = f"{self.base_url}beta/variation/{variant_id}"
            result = self._make_validation_request(url)

            if result:
                clinvar_data = self.get_clinvar_data(variant_id)
                classification = self.extract_classification(clinvar_data)

                response = {
                    "validated": True,
                    "variant_id": variant_id,
                    "genome_build": genome_build,
                    "details": result
                }

                if classification:
                    response["classification"] = classification

                return response
            return {"validated": False, "error": "Validation failed"}
        except Exception as e:
            logger.error(f"Error in RefSeq validation: {str(e)}")
            return {"error": str(e)}

    # 2-4. Key method 4: Ensembl Validation Methods
    def validate_variant_ensembl(self, variant_id, genome_build="GRCh38"):
        try:
            if not validate_genome_build(genome_build):
                return {"error": f"Invalid genome build: {genome_build}"}

            if genome_build == "GRCh38":
                base_url = self.ensembl_base_url
            else:
                base_url = self.grch37_base_url

            url = f"{base_url}/variation/human/{variant_id}"
            headers = {"Content-Type": "application/json"}
            result = self._make_validation_request(url, headers=headers)

            if result:
                return {
                    "validated": True,
                    "variant_id": variant_id,
                    "genome_build": genome_build,
                    "details": result
                }
            return {"validated": False, "error": "Validation failed"}
        except Exception as e:
            logger.error(f"Error in Ensembl validation: {str(e)}")
            return {"error": str(e)}

# Section 3. Utility Functions
def validate_genome_build(genome_build):
    """
    Validates if the provided genome build is supported.
    Args:
        genome_build (str): The genome build to validate
    Returns:
        bool: True if valid, False otherwise
    """
    valid_builds = ["GRCh37", "GRCh38"]
    return genome_build in valid_builds

# 3-1. Standalone function 1: Checks if a genome build is valid (GRCh37 or GRCh38)
def validate_variant_ensembl(variant_id, genome_build="GRCh38"):
    """
    Validate an Ensembl variant ID.
    Args:
        variant_id (str): The variant ID to validate
        genome_build (str): The genome build to use (default: GRCh38)
    Returns:
        dict: Validation results
    """
    try:
        if not validate_genome_build(genome_build):
            return {"error": f"Invalid genome build: {genome_build}"}

        if genome_build == "GRCh38":
            result = query_ensembl_variant_grch38(variant_id)
        else:
            result = query_ensembl_variant_grch37(variant_id)

        if result:
            return {
                "validated": True,
                "variant_id": variant_id,
                "genome_build": genome_build,
                "details": result
            }
        return {"validated": False, "error": "Validation failed"}

    except Exception as e:
        logger.error(f"Error in Ensembl validation: {str(e)}")
        return {"error": str(e)}

# 3-2. Standalone function 2: Queries Ensembl API for GRCh38 variants
def query_ensembl_variant_grch38(variant_id):
    """
    Query Ensembl API for GRCh38 variants.
    """
    try:
        base_url = "https://rest.ensembl.org/variation/human"
        headers = {"Content-Type": "application/json"}
        response = requests.get(f"{base_url}/{variant_id}", headers=headers)

        if response.status_code == 200:
            return response.json()
        logger.warning(f"Ensembl GRCh38 query failed with status code {response.status_code}")
        return None
    except Exception as e:
        logger.error(f"Error querying Ensembl GRCh38: {str(e)}")
        return None

# 3-3. Standalone function 3: Queries Ensembl API for GRCh37 variants
def query_ensembl_variant_grch37(variant_id):
    """
    Query Ensembl API for GRCh37 variants.
    """
    try:
        base_url = "https://grch37.rest.ensembl.org/variation/human"
        headers = {"Content-Type": "application/json"}
        response = requests.get(f"{base_url}/{variant_id}", headers=headers)

        if response.status_code == 200:
            return response.json()
        logger.warning(f"Ensembl GRCh37 query failed with status code {response.status_code}")
        return None
    except Exception as e:
        logger.error(f"Error querying Ensembl GRCh37: {str(e)}")
        return None

# Section 4. API Resource Classes
class ValidateRefSeqAPI(Resource):  # handles GET requests for RefSeq validation
    def get(self, variant_id):
        genome_build = request.args.get('genome_build', 'GRCh38')
        validator = VariantAnnotatorSprint()
        result = validator.validate_variant_refseq(variant_id, genome_build)
        return result

class ValidateEnsemblAPI(Resource): # handles GET requests for Ensemble validation
    """
    REST API resource for validating variants using Ensembl.
    Endpoints:
        GET /api/validate/ensembl/<variant_id>: Validate a variant using Ensembl database
    Parameters:
        variant_id (str): The identifier of the variant to validate
    """

    def __init__(self):
        self.validator = VariantAnnotatorSprint()

    def get(self, variant_id):
        """
        Validate a variant using Ensembl database.

        Args:
            variant_id (str): The identifier of the variant to validate

        Returns:
            dict: JSON response containing validation results

        Raises:
            HTTPException: If validation fails or variant is invalid
        """
        try:
            # Default to GRCh38 as per VariantAnnotatorSprint default
            result = self.validator.validate_variant_ensembl(variant_id)
            return result, 200
        except ValueError as e:
            return {"error": str(e)}, 400
        except Exception as e:
            return {"error": "Internal server error", "details": str(e)}, 500

# Section 5. Flask Application Setup
application = Flask(__name__)
"""Flask Application Configuration

This module is responsible for creating and configuring a Flask application
with REST API routes. It handles the initial setup of the Flask instance
and defines all API endpoints.

Features:
    - Flask application initialization
    - RESTful API route definitions
    - Endpoint configuration and mapping
"""
api = Api(application)
api.add_resource(ValidateRefSeqAPI, '/api/validate/refseq/<string:variant_id>')
api.add_resource(ValidateEnsemblAPI, '/api/validate/ensembl/<string:variant_id>')

@application.route('/health')
def health_check():
    return jsonify({"status": "healthy"})

def main():
    return application

if __name__ == '__main__':
    application = main()
    application.run(debug=True, host='0.0.0.0', port=5000)

# Section 6: Retrieving the sequence in FASTA format
from typing import Optional, Dict, Union

class GeneSequenceFetcher:
    def __init__(self):
        """Initialize the sequence fetcher with necessary API endpoints"""
        self.ensembl_base_url = "https://rest.ensembl.org"
        self.ncbi_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.logger = logging.getLogger(__name__)

    def get_gene_fasta(self, gene_id: str, species: str = "human") -> Optional[Dict[str, Union[str, dict]]]:
        """
        Fetch FASTA sequence for a given gene ID

        Args:
            gene_id (str): Gene identifier (can be Ensembl or HGNC symbol)
            species (str): Species name (default: "human")

        Returns:
            dict: Dictionary containing gene information and FASTA sequence
        """
        try:
            # First try to get Ensembl ID if gene symbol is provided
            ensembl_id = self._get_ensembl_id(gene_id)
            if not ensembl_id:
                self.logger.warning(f"Could not find Ensembl ID for gene: {gene_id}")
                return None

            # Construct URL for sequence request
            url = f"{self.ensembl_base_url}/sequence/id/{ensembl_id}"

            headers = {
                "Content-Type": "application/json",
                "Accept": "application/json"
            }

            params = {
                "type": "genomic",
                "expand_5prime": 0,
                "expand_3prime": 0,
                "format": "fasta"
            }

            response = requests.get(url, headers=headers, params=params)

            if response.status_code != 200:
                self.logger.error(f"Failed to fetch sequence. Status code: {response.status_code}")
                return None

            # Parse the response
            sequence_data = response.json()

            # Format the FASTA sequence
            fasta_sequence = self._format_fasta(
                header=f">{ensembl_id} {sequence_data.get('desc', '')}",
                sequence=sequence_data.get('seq', '')
            )

            return {
                "gene_id": gene_id,
                "ensembl_id": ensembl_id,
                "fasta": fasta_sequence,
                "metadata": {
                    "species": species,
                    "sequence_type": "genomic",
                    "length": len(sequence_data.get('seq', '')),
                    "molecule_type": sequence_data.get('molecule', 'dna')
                }
            }

        except Exception as e:
            self.logger.error(f"Error fetching gene sequence: {str(e)}")
            return None

    def _get_ensembl_id(self, gene_identifier: str) -> Optional[str]:
        """
        Convert gene symbol to Ensembl ID if needed

        Args:
            gene_identifier (str): Gene symbol or Ensembl ID

        Returns:
            str: Ensembl ID if found, None otherwise
        """
        try:
            # Check if already an Ensembl ID
            if gene_identifier.startswith('ENS'):
                return gene_identifier

            # Look up Ensembl ID using symbol
            url = f"{self.ensembl_base_url}/lookup/symbol/homo_sapiens/{gene_identifier}"
            headers = {"Content-Type": "application/json"}

            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                data = response.json()
                return data.get('id')
            return None

        except Exception as e:
            self.logger.error(f"Error converting gene symbol: {str(e)}")
            return None

    def _format_fasta(self, header: str, sequence: str, width: int = 60) -> str:
        """
        Format a sequence in FASTA format with specified line width

        Args:
            header (str): FASTA header line
            sequence (str): Sequence string
            width (int): Number of characters per line (default: 60)

        Returns:
            str: Formatted FASTA sequence
        """
        # Format sequence with specified width
        formatted_seq = '\n'.join(sequence[i:i + width]
                                  for i in range(0, len(sequence), width))
        return f"{header}\n{formatted_seq}"

    def save_fasta(self, gene_id: str, output_file: str) -> bool:
        """
        Fetch and save gene sequence to a FASTA file

        Args:
            gene_id (str): Gene identifier
            output_file (str): Path to output file

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            result = self.get_gene_fasta(gene_id)
            if not result:
                return False

            with open(output_file, 'w') as f:
                f.write(result['fasta'])

            self.logger.info(f"Successfully saved FASTA sequence to {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"Error saving FASTA file: {str(e)}")
            return False

# Section 7: Output
# 7-1. FASTA to JSON
def fasta_to_json(fasta_content):
    """
    Convert FASTA format to JSON with validation
    """
    try:
        if not isinstance(fasta_content, str):
            raise ValueError("Input must be a string")

        if not fasta_content.strip():
            raise ValueError("Empty FASTA content")

        sequences = []
        current_sequence = None
        line_number = 0

        lines = fasta_content.strip().split('\n')

        for line in lines:
            line_number += 1
            line = line.strip()

            if not line:
                continue

            if line.startswith('>'):
                if current_sequence:
                    if not current_sequence["sequence"]:
                        raise ValueError(f"No sequence found for header at line {line_number}")
                    sequences.append(current_sequence)

                header = line[1:].split(maxsplit=1)
                if not header[0]:
                    raise ValueError(f"Empty sequence identifier at line {line_number}")

                current_sequence = {
                    "id": header[0],
                    "description": header[1] if len(header) > 1 else "",
                    "sequence": ""
                }
            elif current_sequence is not None:
                # Validate sequence characters
                if not all(c.upper() in 'ATGCN-' for c in line):
                    raise ValueError(f"Invalid sequence characters at line {line_number}")
                current_sequence["sequence"] += line
            else:
                raise ValueError(f"Sequence data found before header at line {line_number}")

        if current_sequence:
            if not current_sequence["sequence"]:
                raise ValueError("Last sequence has no content")
            sequences.append(current_sequence)

        if not sequences:
            raise ValueError("No valid sequences found")

        return {
            "sequences": sequences,
            "count": len(sequences)
        }

    except Exception as e:
        logger.error(f"Error converting FASTA to JSON: {str(e)}")
        return {
            "error": str(e),
            "status": "failed"
        }
# 7-2. FASTA to VCF
from datetime import datetime
from typing import List, Dict, Optional, Tuple

class FastaToVcf:
    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def convert_fasta_to_vcf(self,
                             reference_sequence: str,
                             query_sequence: str,
                             chromosome: str = "1",
                             reference_name: str = "reference",
                             sample_name: str = "sample") -> Optional[str]:
        """
        Convert FASTA format to VCF by comparing against a reference sequence

        Args:
            reference_sequence (str): The reference sequence
            query_sequence (str): The query sequence to compare
            chromosome (str): Chromosome number/name
            reference_name (str): Name of the reference
            sample_name (str): Name of the sample

        Returns:
            str: VCF format string or None if error
        """
        try:
            # Validate inputs
            if len(reference_sequence) != len(query_sequence):
                raise ValueError("Reference and query sequences must be the same length")

            if not reference_sequence or not query_sequence:
                raise ValueError("Empty sequences provided")

            # Find variants
            variants = self._find_variants(reference_sequence, query_sequence)

            # Generate VCF content
            vcf_content = self._generate_vcf_header(reference_name, sample_name)

            # Add variant records
            for variant in variants:
                vcf_content += self._format_variant_record(
                    chromosome,
                    variant['position'],
                    variant['ref'],
                    variant['alt'],
                    sample_name
                )

            return vcf_content

        except Exception as e:
            self.logger.error(f"Error converting FASTA to VCF: {str(e)}")
            return None

    def _find_variants(self, reference: str, query: str) -> List[Dict]:
        """
        Find variants between reference and query sequences

        Args:
            reference (str): Reference sequence
            query (str): Query sequence

        Returns:
            List[Dict]: List of variant dictionaries
        """
        variants = []
        position = 1  # VCF positions are 1-based

        for ref_base, query_base in zip(reference, query):
            if ref_base != query_base:
                variants.append({
                    'position': position,
                    'ref': ref_base,
                    'alt': query_base
                })
            position += 1

        return variants

    def _generate_vcf_header(self, reference_name: str, sample_name: str) -> str:
        """
        Generate VCF header

        Args:
            reference_name (str): Name of the reference
            sample_name (str): Name of the sample

        Returns:
            str: VCF header string
        """
        current_date = datetime.now().strftime("%Y%m%d")
        header = (
            f"##fileformat=VCFv4.2\n"
            f"##fileDate={current_date}\n"
            f"##source=FastaToVcf\n"
            f"##reference={reference_name}\n"
            "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type of variant\">\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )
        return header

    def _format_variant_record(self,
                               chromosome: str,
                               position: int,
                               ref: str,
                               alt: str,
                               sample_name: str) -> str:
        """
        Format a single variant record in VCF format

        Args:
            chromosome (str): Chromosome number/name
            position (int): Position of variant
            ref (str): Reference allele
            alt (str): Alternative allele
            sample_name (str): Sample name

        Returns:
            str: Formatted VCF record
        """
        return (f"{chromosome}\t{position}\t.\t{ref}\t{alt}\t.\t"
                f"PASS\tTYPE=SNP\tGT\t1/1\n")

    def convert_fasta_file_to_vcf(self,
                                  reference_file: str,
                                  query_file: str,
                                  output_file: str,
                                  chromosome: str = "1") -> bool:
        """
        Convert FASTA files to VCF

        Args:
            reference_file (str): Path to reference FASTA file
            query_file (str): Path to query FASTA file
            output_file (str): Path to output VCF file
            chromosome (str): Chromosome number/name

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Read reference sequence
            ref_seq = self._read_fasta_file(reference_file)
            if not ref_seq:
                raise ValueError(f"Could not read reference file: {reference_file}")

            # Read query sequence
            query_seq = self._read_fasta_file(query_file)
            if not query_seq:
                raise ValueError(f"Could not read query file: {query_file}")

            # Convert to VCF
            vcf_content = self.convert_fasta_to_vcf(
                ref_seq,
                query_seq,
                chromosome=chromosome,
                reference_name=os.path.basename(reference_file),
                sample_name=os.path.basename(query_file)
            )

            if not vcf_content:
                raise ValueError("Failed to convert sequences to VCF")

            # Write VCF file

# Section 8: Export output
import pandas as pd
from typing import Dict, List, Any
import tkinter as tk
from tkinter import messagebox
from datetime import datetime
from pathlib import Path

# Class Definition and Initialization
class ExcelExporter:
    # 8-1. Logging configuration
    def __init__(self):
        self.logger = logging.getLogger(__name__)   # defines available parameters that can be exported
        self.available_parameters = {
            'Gene ID': False,
            'Chromosome Location': False,
            'FASTA Sequence': False,
            'Classification': False,
            'VCF Data': False,
            'Clinical Significance': False,
            'Variant Type': False,
            'Reference Allele': False,
            'Alternative Allele': False,
            'Protein Change': False
        }
        self.data = {}  # initialize an empty data dictionary of parameters (initially all false) to store actual data
    # 8-2. GUI creation method
    def create_parameter_selection_gui(self):
        """
        Create a GUI window to ask users which parameters they want to export
        """
        self.root = tk.Tk() # creates main window
        self.root.title("Select Parameters to Export")
        self.root.geometry("400x500")

        description = tk.Label( # adds description label
            self.root,
            text="Please select the parameters you want to export:",
            pady=10
        )
        description.pack()

        self.checkboxes = {}    # creates checkboxes for each parameter
        for param in self.available_parameters.keys():
            var = tk.BooleanVar()
            checkbox = tk.Checkbutton(
                self.root,
                text=param,
                variable=var,
                command=lambda p=param, v=var: self._update_selection(p, v)
            )
            checkbox.pack(anchor='w', padx=20)
            self.checkboxes[param] = var

        export_button = tk.Button(  # adds export button
            self.root,
            text="Export to Excel",
            command=self._handle_export,
            pady=10
        )
        export_button.pack(pady=20)

        self.root.mainloop()    # starts GUI event loop

    def _update_selection(self, parameter: str, var: tk.BooleanVar):
        """Update the selection status of parameters"""
        self.available_parameters[parameter] = var.get()

    def _handle_export(self):
        """Handle the export button click"""
        selected = [k for k, v in self.available_parameters.items() if v]
        if not selected:
            messagebox.showwarning(
                "No Selection",
                "Please select at least one parameter to export."
            )
            return

        try:
            success = self.export_to_excel(selected)
            if success:
                messagebox.showinfo(
                    "Success",
                    f"Data exported successfully to {self.output_file}"
                )
                self.root.destroy()
            else:
                messagebox.showerror(
                    "Error",
                    "Failed to export data. Please check the logs."
                )
        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {str(e)}")

    def collect_data(self, gene_id: str) -> Dict[str, Any]:
        """
        Collect all the required data based on selected parameters

        Args:
            gene_id (str): The gene identifier

        Returns:
            Dict[str, Any]: Dictionary containing the collected data
        """
        try:
            data = {"Gene ID": gene_id}

            if self.available_parameters["Chromosome Location"]:
                # Add code to fetch chromosome location
                data["Chromosome Location"] = self._get_chromosome_location(gene_id)

            if self.available_parameters["FASTA Sequence"]:
                # Add code to fetch FASTA sequence
                data["FASTA Sequence"] = self._get_fasta_sequence(gene_id)

            if self.available_parameters["Classification"]:
                # Add code to fetch classification
                data["Classification"] = self._get_classification(gene_id)

            if self.available_parameters["VCF Data"]:
                # Add code to fetch VCF data
                data["VCF Data"] = self._get_vcf_data(gene_id)

            # Add more parameter collections as needed

            return data

        except Exception as e:
            self.logger.error(f"Error collecting data: {str(e)}")
            return {}

    def export_to_excel(self, selected_parameters: List[str]) -> bool:
        """
        Export the collected data to Excel

        Args:
            selected_parameters (List[str]): List of parameters to export

        Returns:
            bool: True if export successful, False otherwise
        """
        try:
            # Create output directory if it doesn't exist
            output_dir = Path("output")
            output_dir.mkdir(exist_ok=True)

            # Generate output filename with timestamp
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.output_file = output_dir / f"gene_data_{timestamp}.xlsx"

            # Create DataFrame from collected data
            df = pd.DataFrame([self.data])

            # Select only the requested columns
            selected_columns = [col for col in df.columns if col in selected_parameters]
            df_selected = df[selected_columns]

            # Create Excel writer object
            writer = pd.ExcelWriter(self.output_file, engine='xlsxwriter')

            # Write data to Excel
            df_selected.to_excel(writer, sheet_name='Gene Data', index=False)

            # Get workbook and worksheet objects
            workbook = writer.book
            worksheet = writer.sheets['Gene Data']

            # Add formatting
            header_format = workbook.add_format({
                {
                    'bold': True,
                    'bg_color': '#4F81BD',
                    'font_color': 'white',
                    'border': 1,
                    'align': 'center'
                }
            })

            # Format headers
            for col_num, value in enumerate(df_selected.columns.values):
                worksheet.write(0, col_num, value, header_format)
                worksheet.set_column(col_num, col_num, 20)  # Set column width

            # Save the Excel file
            writer.close()

            self.logger.info(f"Data exported successfully to {self.output_file}")
            return True

        except Exception as e:
            self.logger.error(f"Error exporting to Excel: {str(e)}")
            return False

    def _get_chromosome_location(self, gene_id: str) -> str:
        """Fetch chromosome location"""
        # Implement the actual data fetching logic
        return "chr1:1234567-1234789"

    def _get_fasta_sequence(self, gene_id: str) -> str:
        """Fetch FASTA sequence"""
        # Implement the actual data fetching logic
        return "ATGC..."

    def _get_classification(self, gene_id: str) -> str:
        """Fetch classification"""
        # Implement the actual data fetching logic
        return "Pathogenic"

    def _get_vcf_data(self, gene_id: str) -> str:
        """Fetch VCF data"""
        # Implement the actual data fetching logic
        return "VCF data..."
