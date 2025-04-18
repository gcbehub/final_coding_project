# Section 1. Import and Logging Setup
import os
import logging
import requests
from flask import Flask, jsonify, request
from flask_restful import Resource, Api
import json

# Set up logging
current_directory = os.path.dirname(os.path.abspath(__file__))
log_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
info_file_handler = logging.FileHandler(os.path.join(current_directory, 'info.log'))
info_file_handler.setFormatter(log_formatter)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(info_file_handler)

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
    def __init__(self):
        self.base_url = "https://api.ncbi.nlm.nih.gov/variation/v0/"
        self.clinvar_base_url = "https://www.ncbi.nlm.nih.gov/clinvar/"
        self.ensembl_base_url = "https://rest.ensembl.org"
        self.grch37_base_url = "https://grch37.rest.ensembl.org"

    # 2-1. Key method 1: HTTP Request Handler
    def _make_validation_request(self, url, params=None, headers=None):
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
    def get_clinvar_data(self, variant_id):
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
        try:
            if not clinvar_data:
                return None
            # Add your classification extraction logic here
            return {"classification": "classification_placeholder"}
        except Exception as e:
            logger.error(f"Error extracting classification: {str(e)}")
            return None

    # 2-3. Key method 3: Variant Validation Methods
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
    Validate if the provided genome build is supported.

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
class ValidateRefSeqAPI(Resource):
    def get(self, variant_id):
        genome_build = request.args.get('genome_build', 'GRCh38')
        validator = VariantAnnotatorSprint()
        result = validator.validate_variant_refseq(variant_id, genome_build)
        return result

class ValidateEnsemblAPI(Resource):
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
