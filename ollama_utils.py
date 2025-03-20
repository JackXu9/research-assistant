"""
Ollama Utilities Module

This module provides functions for interacting with the Ollama API.
It includes functions for generating search strategies and summarizing papers.
"""

import requests
import json
import re

def suggest_pubmed_search_strategy(research_question, search_type="broad", model_name="deepseek-r1:8b", ollama_url="http://localhost:11434/api"):
    """
    Generate a PubMed search strategy based on a research question.
    
    Parameters:
    -----------
    research_question : str
        The research question to generate a search strategy for.
    search_type : str, default="broad"
        Type of search. Can be "broad" (prioritizes recall) or "narrow" (prioritizes precision).
    model_name : str, default="deepseek-r1:8b"
        The Ollama model to use for generation.
    ollama_url : str, default="http://localhost:11434/api"
        The URL of the Ollama API.
        
    Returns:
    --------
    dict
        A dictionary containing the search string, explanation, and estimated results.
    """
    prompt = f"""
    You are an information specialist with expertise in literature searches for medical and scientific research.
    Generate a PubMed search strategy for the following research question: "{research_question}"
    
    The search type is: {search_type} ({"prioritize recall by being inclusive" if search_type == "broad" else "prioritize precision by being specific"})
    
    Format your response as a structured JSON with the following fields:
    - search_string: The PubMed search string using MeSH terms and Boolean operators
    - explanation: A brief explanation of your search strategy
    - estimated_results: Your estimate of how many results this will return (rough numerical estimate)
    
    Use proper PubMed syntax with field tags like [Title/Abstract], [MeSH Terms], etc. and Boolean operators AND, OR, NOT.
    """
    
    try:
        response = requests.post(
            f"{ollama_url}/generate",
            json={
                "model": model_name,
                "prompt": prompt,
                "stream": False
            }
        )
        response.raise_for_status()
        result = response.json()
        response_text = result.get("response", "")
        
        # Extract JSON from response
        json_match = re.search(r'```json\s*(.*?)\s*```', response_text, re.DOTALL)
        if json_match:
            json_str = json_match.group(1)
        else:
            json_str = response_text
        
        # Clean up any remaining markdown and try to parse JSON
        json_str = json_str.replace('```', '').strip()
        try:
            result_dict = json.loads(json_str)
            return result_dict
        except json.JSONDecodeError:
            # If JSON parsing fails, try to extract key fields manually
            search_match = re.search(r'"search_string"\s*:\s*"([^"]*)"', json_str)
            explanation_match = re.search(r'"explanation"\s*:\s*"([^"]*)"', json_str)
            estimated_match = re.search(r'"estimated_results"\s*:\s*"?([^",]*)"?', json_str)
            
            return {
                "search_string": search_match.group(1) if search_match else "Error parsing search string",
                "explanation": explanation_match.group(1) if explanation_match else "Error parsing explanation",
                "estimated_results": estimated_match.group(1) if estimated_match else "Unknown"
            }
            
    except Exception as e:
        return {
            "search_string": "",
            "explanation": f"Error generating search strategy: {str(e)}",
            "estimated_results": "Error"
        }

def summarize_papers(df, research_question, model_name="deepseek-r1:8b", ollama_url="http://localhost:11434/api"):
    """
    Summarize a collection of papers based on a research question.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the papers to summarize.
    research_question : str
        The research question to summarize the papers for.
    model_name : str, default="deepseek-r1:8b"
        The Ollama model to use for summarization.
    ollama_url : str, default="http://localhost:11434/api"
        The URL of the Ollama API.
        
    Returns:
    --------
    str
        A summary of the papers in relation to the research question.
    """
    prompt = f"""
    You are a scientific research assistant tasked with synthesizing findings from multiple papers.
    
    Research Question: "{research_question}"
    
    Please synthesize the findings from the following papers:
    
    """
    
    # Add paper details to the prompt
    for i, row in df.iterrows():
        prompt += f"""
        Paper {i+1}:
        Title: {row['Title']}
        Authors: {row['Authors']}
        Journal: {row['Journal']} ({row['Year']})
        Abstract: {row['Abstract']}
        
        """
    
    prompt += """
    Your synthesis should:
    1. Identify key findings relevant to the research question
    2. Note any consistencies or inconsistencies in the findings
    3. Comment on the strength of the evidence
    4. Highlight gaps in knowledge or areas for future research
    5. Provide a concise conclusion based on these papers
    
    Format your response in clear paragraphs with appropriate headings.
    """
    
    try:
        response = requests.post(
            f"{ollama_url}/generate",
            json={
                "model": model_name,
                "prompt": prompt,
                "stream": False
            }
        )
        response.raise_for_status()
        result = response.json()
        return result.get("response", "Error generating summary")
    except Exception as e:
        return f"Error generating summary: {str(e)}" 