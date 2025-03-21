"""
PubMed Utilities Module

This module provides functions for interacting with the PubMed API using the Biopython Entrez module.
It includes functions for searching PubMed, fetching paper details, and converting the results to a
pandas DataFrame. It also includes utilities for handling SSL certificate verification issues.
"""

import pandas as pd
import ssl
import urllib.request
import re
from Bio import Entrez

# Create a global unverified SSL context
unverified_context = ssl._create_unverified_context()

# Patch the default opener to use our unverified context
def patch_ssl():
    """
    Patch urllib to use unverified SSL context.
    
    This function modifies the default urllib opener to use an unverified SSL context,
    which bypasses SSL certificate verification.
    
    Warning:
    --------
    This is a security risk and should only be used for testing purposes.
    """
    opener = urllib.request.build_opener(urllib.request.HTTPSHandler(context=unverified_context))
    urllib.request.install_opener(opener)

def validate_mesh_term(term, email, verify_ssl=True):
    """
    Validates if a term exists as a MeSH term in PubMed.
    
    Parameters:
    -----------
    term : str
        The term to validate (without [Mesh] qualifier)
    email : str
        Email address for PubMed API usage tracking
    verify_ssl : bool, default=True
        Whether to verify SSL certificates
        
    Returns:
    --------
    tuple
        (bool, str) - (is_valid, suggested_term)
        is_valid: True if term exists as MeSH term
        suggested_term: If found, returns the correct MeSH term, otherwise None
    """
    Entrez.email = email
    
    if not verify_ssl:
        patch_ssl()
    
    # Remove any [Mesh] qualifier and trim
    clean_term = term.replace('[Mesh]', '').strip().strip('"')
    
    try:
        # Search in MeSH database
        handle = Entrez.esearch(db="mesh", term=f"{clean_term}[MESH]", retmax=1)
        result = Entrez.read(handle)
        handle.close()
        
        if result['Count'] != '0':
            # Get the correct MeSH term
            handle = Entrez.efetch(db="mesh", id=result['IdList'][0], rettype="full")
            mesh_data = handle.read()
            handle.close()
            
            # Extract the descriptor name
            match = re.search(r'<DescriptorName[^>]*>(.*?)</DescriptorName>', mesh_data)
            if match:
                return True, match.group(1)
            return True, clean_term
        
        # If not found, try to find similar terms
        handle = Entrez.esearch(db="mesh", term=clean_term, retmax=5)
        suggestions = Entrez.read(handle)
        handle.close()
        
        if suggestions['Count'] != '0':
            # Get the first suggestion
            handle = Entrez.efetch(db="mesh", id=suggestions['IdList'][0], rettype="full")
            mesh_data = handle.read()
            handle.close()
            
            match = re.search(r'<DescriptorName[^>]*>(.*?)</DescriptorName>', mesh_data)
            if match:
                return False, match.group(1)
        
        return False, None
        
    except Exception as e:
        print(f"Error validating MeSH term '{term}': {str(e)}")
        return False, None

def clean_search_query(query, email, verify_ssl=True):
    """
    Cleans and validates a PubMed search query by checking MeSH terms.
    
    Parameters:
    -----------
    query : str
        The search query to clean
    email : str
        Email address for PubMed API usage tracking
    verify_ssl : bool, default=True
        Whether to verify SSL certificates
        
    Returns:
    --------
    tuple
        (str, list) - (cleaned_query, list of warnings/suggestions)
    """
    warnings = []
    
    # Find all MeSH terms in the query
    mesh_terms = re.findall(r'"([^"]+)"\[Mesh\]', query)
    
    cleaned_parts = query
    for term in mesh_terms:
        is_valid, suggestion = validate_mesh_term(term, email, verify_ssl)
        if not is_valid:
            if suggestion:
                warnings.append(f'MeSH term "{term}" not found. Suggested: "{suggestion}"')
                # Replace invalid MeSH term with suggested one
                cleaned_parts = cleaned_parts.replace(
                    f'"{term}"[Mesh]',
                    f'"{suggestion}"[Mesh]'
                )
            else:
                warnings.append(f'MeSH term "{term}" not found. Using as text search.')
                # Convert to All Fields search
                cleaned_parts = cleaned_parts.replace(
                    f'"{term}"[Mesh]',
                    f'"{term}"[All Fields]'
                )
    
    return cleaned_parts, warnings

def search_pubmed(query, email, retmax=1000, verify_ssl=True):
    """
    Search PubMed using the given query.
    
    Parameters:
    -----------
    query : str
        The search query to submit to PubMed.
    email : str
        Email address for PubMed API usage tracking.
    retmax : int, default=1000
        Maximum number of results to return.
    verify_ssl : bool, default=True
        Whether to verify SSL certificates. Set to False to bypass SSL verification.
        
    Returns:
    --------
    tuple
        A tuple containing (list of PubMed IDs, total count of results, list of warnings, cleaned_query).
    """
    if not query or not query.strip():
        raise ValueError("Empty search query")
        
    if not email or not email.strip():
        raise ValueError("Email address is required")
        
    if not isinstance(retmax, int) or retmax <= 0:
        raise ValueError("retmax must be a positive integer")
    
    # Clean and validate the query first
    cleaned_query, warnings = clean_search_query(query, email, verify_ssl)
    
    Entrez.email = email
    
    # If SSL verification is disabled, patch urllib
    if not verify_ssl:
        patch_ssl()
    
    try:
        # Clean up the query by removing any double spaces or trailing/leading spaces
        cleaned_query = ' '.join(cleaned_query.split())
        
        # Print debug information
        print(f"Executing PubMed search with query: {cleaned_query}")
        
        handle = Entrez.esearch(
            db="pubmed", 
            term=cleaned_query, 
            retmax=retmax,
            usehistory='y',
            retmode='xml'
        )
        
        record = Entrez.read(handle, validate=False)
        handle.close()
        
        if 'ErrorList' in record and record['ErrorList']:
            error_details = []
            for error_type, errors in record['ErrorList'].items():
                error_details.append(f"{error_type}: {', '.join(str(e) for e in errors)}")
            error_msg = '; '.join(error_details)
            raise Exception(f"PubMed API error: {error_msg}")
            
        pmid_list = record.get("IdList", [])
        count = int(record.get("Count", 0))
        
        return pmid_list, count, warnings, cleaned_query
        
    except Exception as e:
        print(f"Error in search_pubmed: {str(e)}")
        print(f"Query that caused error: {cleaned_query}")
        raise Exception(f"PubMed search error: {str(e)}\nQuery: {cleaned_query}")

def fetch_details(id_list, email, verify_ssl=True):
    """
    Fetches details for a list of PubMed IDs.
    
    Parameters:
    -----------
    id_list : list
        List of PubMed IDs to fetch details for.
    email : str
        Email address for PubMed API usage tracking.
    verify_ssl : bool, default=True
        Whether to verify SSL certificates. Set to False to bypass SSL verification.
        
    Returns:
    --------
    dict
        A dictionary containing the paper details in PubMed XML format.
    """
    if not id_list:
        return None
    Entrez.email = email
    
    # If SSL verification is disabled, patch urllib
    if not verify_ssl:
        patch_ssl()
    
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
    papers = Entrez.read(handle, validate=False)
    handle.close()
    return papers

def papers_to_df(papers):
    """
    Converts the PubMed XML structure into a Pandas DataFrame.
    
    Parameters:
    -----------
    papers : dict
        Dictionary containing paper details in PubMed XML format,
        as returned by the fetch_details function.
        
    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the paper details with columns:
        Title, Authors, Journal, Year, PMID, Abstract.
    """
    if not papers or 'PubmedArticle' not in papers:
        return pd.DataFrame()

    rows = []
    for paper in papers['PubmedArticle']:
        article = paper['MedlineCitation']['Article']
        title = article.get('ArticleTitle', 'Not available')
        pubmed_id = paper['MedlineCitation']['PMID']
        journal = article['Journal'].get('Title', 'Not available')
        pubdate = article['Journal']['JournalIssue']['PubDate']
        year = pubdate.get('Year', 'Not available')
        authors_list = article.get('AuthorList', [])
        authors = ', '.join([f"{a.get('LastName', '')}, {a.get('ForeName', '')}" 
                             for a in authors_list if 'LastName' in a and 'ForeName' in a])
        abstract_text = 'Not available'
        if 'Abstract' in article and 'AbstractText' in article['Abstract']:
            abstract_text = " ".join(article['Abstract']['AbstractText'])
        rows.append([
            title, 
            authors, 
            journal, 
            year, 
            pubmed_id, 
            abstract_text
        ])
    
    return pd.DataFrame(rows, columns=[
        'Title', 
        'Authors', 
        'Journal', 
        'Year', 
        'PMID', 
        'Abstract'
    ]) 