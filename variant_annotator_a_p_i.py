import os
import logging
import requests
from flask import Flask, jsonify, request

# Set up logging
current_directory = os.path.dirname(os.path.abspath(__file__))
log_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
info_file_handler = logging.FileHandler(os.path.join(current_directory, 'info.log'))
info_file_handler.setFormatter(log_formatter)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(info_file_handler)

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

def query_ensembl_variant_grch38(variant_id):
    """
    Query Ensembl API for GRCh38 variants.

    Args:
        variant_id (str): The variant ID to query

    Returns:
        dict: Query results
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

def query_ensembl_variant_grch37(variant_id):
    """
    Query Ensembl API for GRCh37 variants.

    Args:
        variant_id (str): The variant ID to query

    Returns:
        dict: Query results
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

def validate_refseq_api(variant_id):
    """
    Validate a variant using RefSeq API.

    Args:
        variant_id (str): The variant ID to validate

    Returns:
        dict: Validation results
    """
    try:
        base_url = "https://api.ncbi.nlm.nih.gov/variation/v0"
        response = requests.get(f"{base_url}/variation/{variant_id}")

        if response.status_code == 200:
            return response.json()
        logger.warning(f"RefSeq API validation failed with status code {response.status_code}")
        return None
    except Exception as e:
        logger.error(f"Error in RefSeq API validation: {str(e)}")
        return None

def validate_ensembl_api(variant_id, genome_build="GRCh38"):
    """
    Validate a variant using Ensembl API.

    Args:
        variant_id (str): The variant ID to validate
        genome_build (str): The genome build to use (default: GRCh38)

    Returns:
        dict: Validation results
    """
    try:
        if genome_build == "GRCh38":
            result = query_ensembl_variant_grch38(variant_id)
        else:
            result = query_ensembl_variant_grch37(variant_id)

        return result if result else {"error": "Validation failed"}
    except Exception as e:
        logger.error(f"Error in Ensembl API validation: {str(e)}")
        return {"error": str(e)}

[... Previous VariantAnnotatorSprint

class and other classes remain the same...]

# Create Flask application
application = Flask(__name__)
api = None  # You'll need to set up your API routes here

@ application.route('/api/validate/refseq/<variant_id>')

def validate_variant_refseq_api(variant_id):
    """
    API endpoint for RefSeq variant validation.
    """
    genome_build = request.args.get('genome_build', 'GRCh38')
    validator = VariantAnnotatorSprint()
    result = validator.validate_variant_refseq(variant_id, genome_build)
    return jsonify(result)

@application.route('/api/validate/ensembl/<variant_id>')
def validate_variant_ensembl_api(variant_id):
    """
    API endpoint for Ensembl variant validation.
    """
    genome_build = request.args.get('genome_build', 'GRCh38')
    result = validate_variant_ensembl(variant_id, genome_build)
    return jsonify(result)

if __name__ == '__main__':
    application = main()
    application.run(debug=True)
