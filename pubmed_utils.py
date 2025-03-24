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

def clean_search_query(query: str, email: str = None) -> tuple:
    """
    Clean and validate a PubMed search query.
    
    Args:
        query (str): The original search query
        email (str, optional): Email for PubMed API
    
    Returns:
        tuple: (cleaned_query, warnings)
    """
    warnings = []
    original_query = query
    
    # If query uses OVID/MEDLINE syntax, convert it
    if any(pattern in query for pattern in ['exp ', '/', '.ti,ab.', '[pt]']):
        try:
            query = convert_ovid_to_pubmed(query)
            warnings.append("OVID/MEDLINE syntax detected and automatically converted to PubMed format")
        except Exception as e:
            warnings.append(f"Error converting OVID syntax: {str(e)}")
    
    # Fix common date format issues
    if "[Date - Publication]" in query:
        # Fix spacing before date ranges
        query = re.sub(r'(\d{4})/(\d{2})/(\d{2}):', r'\1/\2/\3 :', query)
        
        # Fix incorrect date format with year immediately after AND
        date_year_pattern = re.compile(r'AND\s+(\d{4})\s*\[(?:All Fields|Text Word)\](\d{2}/\d{2}:\d{4}/\d{2}/\d{2}\[Date - Publication\])')
        if date_year_pattern.search(query):
            query = date_year_pattern.sub(r'AND \1/\2', query)
            warnings.append("Fixed date filter being incorrectly treated as MeSH term")
    
    # Additional date format fix: look for year numbers followed by date pattern
    year_date_pattern = re.compile(r'(\d{4})\s*\[(?:All Fields|Text Word)\](\d{2}/\d{2}:\d{4}/\d{2}/\d{2}\[Date - Publication\])')
    if year_date_pattern.search(query):
        query = year_date_pattern.sub(r'\1/\2', query)
        warnings.append("Fixed date format where year was separated from date range")
    
    # Fix another common pattern: YYYY[All Fields]MM/DD:YYYY/MM/DD
    broken_date_pattern = re.compile(r'(\d{4})\s*\[(?:All Fields|Text Word)\](\d{2}/\d{2}:\d{4}/\d{2}/\d{2}\[Date - Publication\])')
    if broken_date_pattern.search(query):
        query = broken_date_pattern.sub(r'\1/\2', query)
        warnings.append("Corrected broken date format")
        
    # Check for AND YYYY[All Fields] pattern
    year_mesh_pattern = re.compile(r'AND\s+(\d{4})\s*\[(?:All Fields|Text Word|Mesh(?:\:NoExp)?)\]')
    if year_mesh_pattern.search(query):
        # Only fix if followed by date pattern
        date_after = re.search(r'AND\s+(\d{4})\s*\[(?:All Fields|Text Word|Mesh(?:\:NoExp)?)\]([^[]*?\[Date - Publication\])', query)
        if date_after:
            year = date_after.group(1)
            date_part = date_after.group(2)
            # Replace with proper format
            fixed_part = f"AND {year}/{date_part}"
            query = re.sub(r'AND\s+(\d{4})\s*\[(?:All Fields|Text Word|Mesh(?:\:NoExp)?)\]([^[]*?\[Date - Publication\])', fixed_part, query)
            warnings.append(f"MeSH term \"AND {year}\" not found. Using as text search.")
            warnings.append("Fixed date filter being incorrectly treated as MeSH term")
    
    # Validate MeSH terms if email is provided
    if email:
        # Find all MeSH terms in the query
        mesh_terms = re.findall(r'"([^"]+)"(?:\s*\[Mesh[^\]]*\])', query)
        mesh_terms += re.findall(r'([^"\s]+)(?:\s*\[Mesh[^\]]*\])', query)  # Also catch unquoted MeSH terms
        
        invalid_terms = []
        suggestions = {}
        
        for term in mesh_terms:
            # Skip checking YYYY as MeSH
            if re.match(r'^\d{4}$', term):
                invalid_terms.append(term)
                suggestions[term] = None
                continue
                
            try:
                is_valid, suggestion = validate_mesh_term(term, email, verify_ssl=False)
                if not is_valid:
                    invalid_terms.append(term)
                    suggestions[term] = suggestion
            except Exception as e:
                warnings.append(f"Error checking MeSH term '{term}': {str(e)}")
        
        # Add warnings for invalid terms
        for term in invalid_terms:
            if re.match(r'^\d{4}$', term):
                warnings.append(f"MeSH term \"{term}\" not found. Using as text search.")
                # Replace MeSH tag with All Fields
                query = re.sub(r'"?' + re.escape(term) + r'"?\s*\[Mesh[^\]]*\]', f'"{term}"[All Fields]', query)
            elif suggestions[term]:
                warnings.append(f"MeSH term \"{term}\" not found. Did you mean \"{suggestions[term]}\"?")
            else:
                warnings.append(f"MeSH term \"{term}\" not found. Using as text search.")
                # Replace MeSH tag with All Fields
                query = re.sub(r'"?' + re.escape(term) + r'"?\s*\[Mesh[^\]]*\]', f'"{term}"[All Fields]', query)
    
    # Clean up duplicate date filters
    date_filters = re.findall(r'\d{4}/\d{2}/\d{2}:\d{4}/\d{2}/\d{2}\[Date - Publication\]', query)
    if len(date_filters) > 1:
        # Keep just the first instance
        query_parts = query.split(date_filters[0])
        rest_of_query = query_parts[1]
        for date_filter in date_filters[1:]:
            rest_of_query = rest_of_query.replace(f" AND {date_filter}", "")
        query = query_parts[0] + date_filters[0] + rest_of_query
        warnings.append("Removed duplicate date filters")
    
    return query, warnings

def convert_ovid_to_pubmed(query):
    """
    Convert OVID/MEDLINE style queries to PubMed syntax.
    
    Parameters:
    -----------
    query : str
        The query in OVID/MEDLINE format
        
    Returns:
    --------
    str
        The query converted to PubMed syntax
    """
    print(f"Converting query from OVID to PubMed: {query}")
    
    # Make a copy of the original query
    converted_query = query
    
    # Handle exp
    converted_query = re.sub(r'exp\s+([^/]+)/', r'"\1"[Mesh]', converted_query)
    converted_query = re.sub(r'exp\s+([^/\s]+)', r'"\1"[Mesh]', converted_query)
    
    # Handle mesh terms with trailing slash
    # Pattern for quoted terms with trailing slash e.g. "Heart Disease"/
    converted_query = re.sub(r'"([^"]+)"/', r'"\1"[Mesh]', converted_query)
    
    # Pattern for unquoted terms with trailing slash e.g. Heart Disease/
    converted_query = re.sub(r'([a-zA-Z][a-zA-Z0-9 -]*)/', r'"\1"[Mesh]', converted_query)
    
    # Handle publication types
    converted_query = converted_query.replace('[pt]', '[Publication Type]')
    
    # Handle tiab formatting
    converted_query = converted_query.replace('[tiab]', '[Title/Abstract]')
    
    # Fix missing quotes around MeSH terms
    converted_query = re.sub(r'([A-Za-z][A-Za-z0-9 -]+)\[Mesh\]', r'"\1"[Mesh]', converted_query)
    converted_query = re.sub(r'([A-Za-z][A-Za-z0-9 -]+)\[Mesh:NoExp\]', r'"\1"[Mesh:NoExp]', converted_query)
    
    # Ensure proper spacing around operators
    converted_query = re.sub(r'\s+AND\s+', ' AND ', converted_query)
    converted_query = re.sub(r'\s+OR\s+', ' OR ', converted_query)
    converted_query = re.sub(r'\s+NOT\s+', ' NOT ', converted_query)
    
    # Fix spacing issues with asterisks
    converted_query = re.sub(r'(\w+)\*\s+\[', r'\1*[', converted_query)
    
    # Fix date format issues - make sure date ranges are properly formatted
    # Look for date ranges like YYYY/MM/DD:YYYY/MM/DD[Date - Publication]
    date_pattern = r'(\d{4}/\d{2}/\d{2}):(\d{4}/\d{2}/\d{2})\[Date - Publication\]'
    if re.search(date_pattern, converted_query):
        # Ensure proper spacing before the date filter
        converted_query = re.sub(r'AND(\d{4}/\d{2}/\d{2})', r'AND \1', converted_query)
        converted_query = re.sub(r'\)(\d{4}/\d{2}/\d{2})', r') \1', converted_query)
    
    # Print the result for debugging
    print(f"Converted query: {converted_query}")
    
    return converted_query

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
    cleaned_query, warnings = clean_search_query(query, email)
    
    Entrez.email = email
    
    # If SSL verification is disabled, patch urllib
    if not verify_ssl:
        patch_ssl()
    
    try:
        # Clean up the query by removing any double spaces or trailing/leading spaces
        cleaned_query = ' '.join(cleaned_query.split())
        
        # Print debug information
        print(f"Executing PubMed search with query: {cleaned_query}")
        print(f"Query length: {len(cleaned_query)} characters")
        
        # Check for potential issues with date formatting
        if ':' in cleaned_query and '[Date - Publication]' in cleaned_query:
            date_part = re.search(r'(\d{4}/\d{2}/\d{2}):(\d{4}/\d{2}/\d{2})\[Date - Publication\]', cleaned_query)
            if date_part:
                print(f"Date range found: {date_part.group(0)}")
            else:
                print("WARNING: Date format might be incorrect")
                # Try to fix date formatting if it seems incorrect
                if 'AND ' not in cleaned_query and re.search(r'AND\d{4}/\d{2}/\d{2}', cleaned_query):
                    print("Fixing missing space before date range")
                    cleaned_query = re.sub(r'AND(\d{4}/\d{2}/\d{2})', r'AND \1', cleaned_query)
        
        # Check for [Mesh:NoExp] format with year, which can cause issues
        if re.search(r'"AND \d{4}"', cleaned_query):
            print("WARNING: Detected potential issue with 'AND YYYY' being treated as MeSH term")
            # Fix the issue by ensuring space before date filter
            cleaned_query = re.sub(r'"AND (\d{4})"(\[Mesh:NoExp\])?', r'AND \1', cleaned_query)
        
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
            
            # Try to fix specific errors
            fixed_query, fix_warnings = try_fix_query_errors(cleaned_query, record['ErrorList'])
            warnings.extend(fix_warnings)
            
            if fixed_query != cleaned_query:
                # Retry with fixed query
                print(f"Retrying with fixed query: {fixed_query}")
                handle = Entrez.esearch(
                    db="pubmed", 
                    term=fixed_query, 
                    retmax=retmax,
                    usehistory='y',
                    retmode='xml'
                )
                
                record = Entrez.read(handle, validate=False)
                handle.close()
                
                if 'ErrorList' in record and record['ErrorList']:
                    # Still have errors after fixing
                    error_msg = '; '.join(error_details)
                    raise Exception(f"PubMed API error: {error_msg}")
                else:
                    # Success with fixed query
                    cleaned_query = fixed_query
            else:
                # Could not fix the query
                error_msg = '; '.join(error_details)
                raise Exception(f"PubMed API error: {error_msg}")
            
        pmid_list = record.get("IdList", [])
        count = int(record.get("Count", 0))
        
        return pmid_list, count, warnings, cleaned_query
        
    except Exception as e:
        print(f"Error in search_pubmed: {str(e)}")
        print(f"Query that caused error: {cleaned_query}")
        
        # If the error is related to PubMed API and it's a field or phrase error,
        # try one last recovery attempt by converting everything to text search
        if "PubMed API error" in str(e) and ("FieldNotFound" in str(e) or "PhraseNotFound" in str(e)):
            try:
                # Convert everything to [All Fields] as a last resort
                emergency_query = emergency_fix_query(cleaned_query)
                warnings.append("Emergency query conversion: Converting all terms to text search.")
                
                print(f"Attempting emergency query: {emergency_query}")
                
                handle = Entrez.esearch(
                    db="pubmed", 
                    term=emergency_query, 
                    retmax=retmax,
                    usehistory='y',
                    retmode='xml'
                )
                
                record = Entrez.read(handle, validate=False)
                handle.close()
                
                if 'ErrorList' not in record or not record['ErrorList']:
                    # Emergency query worked
                    pmid_list = record.get("IdList", [])
                    count = int(record.get("Count", 0))
                    return pmid_list, count, warnings, emergency_query
            except Exception as e2:
                # Even the emergency fix failed
                print(f"Emergency query also failed: {str(e2)}")
                pass  # Fall through to the original error
                
        # If we get here, either it wasn't a PubMed API error or the emergency fix failed
        raise Exception(f"PubMed search error: {str(e)}\nQuery: {cleaned_query}")

def try_fix_query_errors(query, error_list):
    """
    Attempts to fix specific query errors based on the error list from PubMed.
    
    Parameters:
    -----------
    query : str
        The query that caused the error
    error_list : dict
        ErrorList returned from PubMed API
        
    Returns:
    --------
    tuple
        (fixed_query, warnings) - (the fixed query, list of warnings about changes)
    """
    fixed_query = query
    warnings = []
    
    # Handle FieldNotFound
    if 'FieldNotFound' in error_list:
        for field in error_list['FieldNotFound']:
            warnings.append(f"Field not found: {field}")
            
            # Replace common field errors
            if field == '':
                # This happens with malformed field tags, look for patterns
                if '[pt]' in fixed_query.lower():
                    fixed_query = fixed_query.replace('[pt]', '[Publication Type]')
                    warnings.append("Replaced [pt] with [Publication Type]")
    
    # Handle PhraseNotFound
    if 'PhraseNotFound' in error_list:
        for phrase in error_list['PhraseNotFound']:
            warnings.append(f"Phrase not found: {phrase}")
            
            # Handle specific phrases that cause problems
            if phrase == 'animal':
                fixed_query = fixed_query.replace('animal [pt]', '"animals"[MeSH Terms]')
                fixed_query = fixed_query.replace('NOT animal', 'NOT "animals"[MeSH Terms]')
                fixed_query = fixed_query.replace('NOT (animal', 'NOT ("animals"[MeSH Terms]')
                warnings.append("Replaced 'animal' with 'animals'[MeSH Terms]")
    
    # Basic syntax fixes regardless of specific errors
    # Fix truncation syntax
    if '*' in fixed_query:
        fixed_query = re.sub(r'([a-zA-Z]+)\*\s+\[([^\]]+)\]', r'\1*[\2]', fixed_query)
        warnings.append("Fixed truncation syntax")
    
    # Fix multiple adjacent brackets
    fixed_query = re.sub(r'\]\s*\[', '] [', fixed_query)
    
    # Remove double spaces
    fixed_query = ' '.join(fixed_query.split())
    
    return fixed_query, warnings

def emergency_fix_query(query):
    """
    Emergency fix for a query by converting everything to [All Fields]
    
    Parameters:
    -----------
    query : str
        The query that caused errors
        
    Returns:
    --------
    str
        A simplified version of the query using only [All Fields]
    """
    # Replace any field qualifiers with [All Fields]
    # This will replace [Mesh], [Publication Type], etc.
    simplified_query = re.sub(r'\[[^\]]+\]', '[All Fields]', query)
    
    # Fix truncation if present
    simplified_query = re.sub(r'([a-zA-Z]+)\*\s+\[', r'\1*[', simplified_query)
    
    # Fix quoting issues - ensure terms with qualifiers are quoted
    simplified_query = re.sub(r'([A-Za-z][A-Za-z0-9 -]+)\[All Fields\]', r'"\1"[All Fields]', simplified_query)
    
    # Remove double spaces
    simplified_query = ' '.join(simplified_query.split())
    
    return simplified_query

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