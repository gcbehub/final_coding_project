import os # imports the library required for file/directory operation
import logging # imports the library required for application logging
from logging.handlers import RotatingFileHandler # imports the library required for log file rotation
import json # imports the library required for JSON data handling
import requests # imports the library required for making HTTP requests
from flask import Flask, request, jsonify # imports the library required for creating web API
from flask_restful import Resource, Api # imports the library required for implementing RESTful API

# Logging setup
current_directory = os.path.dirname(os.path.abspath(__file__)) # gets current directory path
log_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s') # creates a log formatter
info_file_handler = RotatingFileHandler(
    os.path.join(current_directory, 'info.log'),
    mode='a',
    maxBytes=5 * 1024 * 1024,
    backupCount=2,
    encoding=None,
    delay=0
) # sets up rotating file handler

info_file_handler.setFormatter(log_formatter) # configures the log handler with formatter
info_file_handler.setLevel(logging.INFO) # sets logging level to INFO
logger = logging.getLogger('root') # creates root logger and adds handler
logger.setLevel(logging.INFO)
logger.addHandler(info_file_handler)

def validate_genome_build(genome_build):
    # validates the genome build parameter
    # takes the genome build to validate as argument with genome_build (str)
    # returns boolian (True if valid, False otherwise)

    valid_builds = ['GRCh37', 'GRCh38'] # only accepts GRCh37 or GRCh38
    try:
        if genome_build in valid_builds:
            return True
        return False
    except Exception as e: # logs any errors during validation
        logger.error(f"Error validating genome build: {str(e)}")
        return False

def main(): # main function to run the variant annotator
    try:
        genome_build = "GRCh38"  # Default genome build
        if not validate_genome_build(genome_build):
            logger.error(f"Invalid genome build: {genome_build}")
            return

        validation = VariantAnnotatorSprint()
        logger.info("Variant annotator initialized successfully")

    except Exception as e:
        logger.error(f"Error in main function: {str(e)}")
        return

class VariantAnnotatorSprint: # main class for validating and annotating genetic variants
    def __init__(self, base_url="https://rest.variantvalidator.org"):
        # initialize the variant validator
        # takes the API URL as argument

        self.base_url = base_url

    def validate_variant_refseq(self, variant_id):
        """
        Validate a variant using RefSeq
        Args:
            variant_id (str): The variant identifier
        Returns:
            dict: Validation results
        """
        endpoint = f"/VariantValidator/variantvalidator/{variant_id}/all"
        return self._make_validation_request(endpoint)

    def _make_validation_request(self, endpoint):
        """
        Make a request to the validation API
        Args:
            endpoint (str): API endpoint
        Returns:
            dict: API response
        """
        try:
            response = requests.get(f"{self.base_url}{endpoint}")
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logger.error(f"API request failed: {str(e)}")
            return None

    def validate_hgvs_variant(self, hgvs_id):
        """
        Validate an HGVS variant identifier
        Args:
            hgvs_id (str): The HGVS variant identifier
        Returns:
            dict: Validation results
        """
        try:
            validation = self.validate_variant_refseq(hgvs_id)
            if validation:
                validation_result = {
                    "is_valid": True,
                    "message": "Variant validated successfully"
                }
            else:
                validation_result = {
                    "is_valid": False,
                    "message": "Invalid variant format"
                }
            return validation_result
        except Exception as e:
            logger.error(f"Error validating HGVS variant: {str(e)}")
            return {"is_valid": False, "message": str(e)}

    def annotate_variant(self, variant_id, genome_build="GRCh38"):
        """
        Annotate a variant with additional information
        Args:
            variant_id (str): The variant identifier
            genome_build (str): Genome build version
        Returns:
            dict: Annotation results
        """
        try:
            validation = self.validate_variant_refseq(variant_id)
            if not validation:
                return {"error": "Invalid variant"}

            annotator = self.get_annotator(genome_build)
            validation_result = annotator(variant_id)

            if validation_result:
                annotation = {
                    "variant_id": variant_id,
                    "validation": validation_result,
                    "genome_build": genome_build
                }
            else:
                annotation = {
                    "error": "Failed to annotate variant",
                    "variant_id": variant_id
                }
            return annotation
        except Exception as e:
            logger.error(f"Error annotating variant: {str(e)}")
            return {"error": str(e)}

    def search_clinvar_by_hgvs(self, hgvs_id):
        """
        Search ClinVar database using HGVS identifier
        Args:
            hgvs_id (str): HGVS variant identifier
        Returns:
            dict: ClinVar search results
        """
        try:
            endpoint = f"/clinvar/{hgvs_id}"
            return self._make_validation_request(endpoint)
        except Exception as e:
            logger.error(f"Error searching ClinVar: {str(e)}")
            return None

    def extract_classifications(self, clinvar_data):
        """
        Extract clinical classifications from ClinVar data
        Args:
            clinvar_data (dict): ClinVar response data
        Returns:
            list: List of clinical classifications
        """
        classifications = []
        try:
            if clinvar_data and 'classifications' in clinvar_data:
                for classification in clinvar_data['classifications']:
                    classifications.append({
                        'source': classification.get('source'),
                        'classification': classification.get('term'),
                        'date': classification.get('date')
                    })
        except Exception as e:
            logger.error(f"Error extracting classifications: {str(e)}")
        return classifications

# Flask API setup
application = Flask(__name__)
api = Api(application)

def validate_variant_refseq(variant_id, genome_build="GRCh38"):
    """
    Validate a variant using RefSeq through API
    Args:
        variant_id (str): The variant identifier
        genome_build (str): Genome build version
    Returns:
        dict: Validation results
    """
    try:
        if not validate_genome_build(genome_build):
            return {"error": f"Invalid genome build: {genome_build}"}

        validator = VariantAnnotatorSprint()
        result = validator.validate_variant_refseq(variant_id)

        if result:
            return {
                "validated": True,
                "variant_id": variant_id,
                "genome_build": genome_build,
                "details": result
            }
        return {"validated": False, "error": "Validation failed"}

    except Exception as e:
        logger.error(f"Error in RefSeq validation: {str(e)}")
        return {"error": str(e)}

def validate_variant_ensembl(variant_id, genome_build="GRCh38"):
    """
    Validate a variant using Ensembl
    Args:
        variant_id (str): The variant identifier
        genome_build (str): Genome build version
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

def query_ensembl_variant_grch38(variant_id):
    """
    Query Ensembl API for GRCh38
    Args:
        variant_id (str): The variant identifier
    Returns:
        dict: Query results
    """
    try:
        base_url = "https://rest.ensembl.org"
        endpoint = f"/variation/human/{variant_id}?"

        headers = {"Content-Type": "application/json"}
        response = requests.get(f"{base_url}{endpoint}", headers=headers)

        if response.status_code == 200:
            return response.json()
        return None

    except Exception as e:
        logger.error(f"Error querying Ensembl GRCh38: {str(e)}")
        return None

def query_ensembl_variant_grch37(variant_id):
    """
    Query Ensembl API for GRCh37
    Args:
        variant_id (str): The variant identifier
    Returns:
        dict: Query results
    """
    try:
        base_url = "https://grch37.rest.ensembl.org"
        endpoint = f"/variation/human/{variant_id}?"

        headers = {"Content-Type": "application/json"}
        response = requests.get(f"{base_url}{endpoint}", headers=headers)

        if response.status_code == 200:
            return response.json()
        return None

    except Exception as e:
        logger.error(f"Error querying Ensembl GRCh37: {str(e)}")
        return None

# API Resource classes
class ValidateRefSeqAPI(Resource):
    def get(self, variant_id):
        """
        API endpoint for RefSeq validation
        """
        genome_build = request.args.get('genome_build', 'GRCh38')
        return validate_refseq_api(variant_id, genome_build)

def validate_refseq_api(variant_id, genome_build="GRCh38"):
    """
    API handler for RefSeq validation
    Args:
        variant_id (str): The variant identifier
        genome_build (str): Genome build version
    Returns:
        dict: API response
    """
    try:
        result = validate_variant_refseq(variant_id, genome_build)
        if "error" in result:
            return {"status": "error", "message": result["error"]}, 400
        return {"status": "success", "data": result}, 200

    except Exception as e:
        logger.error(f"API Error in RefSeq validation: {str(e)}")
        return {"status": "error", "message": str(e)}, 500

class ValidateEnsemblAPI(Resource):
    def get(self, variant_id):
        """
        API endpoint for Ensembl validation
        """
        genome_build = request.args.get('genome_build', 'GRCh38')
        return validate_ensembl_api(variant_id, genome_build)

def validate_ensembl_api(variant_id, genome_build="GRCh38"):
    """
    API handler for Ensembl validation
    Args:
        variant_id (str): The variant identifier
        genome_build (str): Genome build version
    Returns:
        dict: API response
    """
    try:
        result = validate_variant_ensembl(variant_id, genome_build)
        if "error" in result:
            return {"status": "error", "message": result["error"]}, 400
        return {"status": "success", "data": result}, 200

    except Exception as e:
        logger.error(f"API Error in Ensembl validation: {str(e)}")
        return {"status": "error", "message": str(e)}, 500

# API routes
api.add_resource(ValidateRefSeqAPI, '/api/validate/refseq/<string:variant_id>')
api.add_resource(ValidateEnsemblAPI, '/api/validate/ensembl/<string:variant_id>')

def validate_variant_refseq_api(variant_id):
    """
    Wrapper function for RefSeq API validation
    Args:
        variant_id (str): The variant identifier
    Returns:
        dict: Validation results
    """
    try:
        return validate_variant_refseq(variant_id)
    except Exception as e:
        logger.error(f"Error in RefSeq API validation: {str(e)}")
        return {"error": str(e)}

if __name__ == '__main__':
    main()
