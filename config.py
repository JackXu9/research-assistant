"""
Configuration Module

This module handles environment variables and configuration settings for the application.
It loads variables from a .env file and provides default values for optional settings.
"""

import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

class Config:
    """
    Configuration class that stores all application settings.
    
    This class loads settings from environment variables with fallbacks to default values.
    It also provides a validation method to ensure required settings are present.
    
    Attributes:
    -----------
    OLLAMA_API_URL : str
        The URL for the Ollama API, defaults to 'http://localhost:11434/api'.
    OLLAMA_MODEL : str
        The name of the Ollama model to use, defaults to 'deepseek-r1:8b'.
    PUBMED_EMAIL : str
        The email address to use for PubMed API requests (required).
    """
    
    # Ollama API URL - default to localhost if not specified
    OLLAMA_API_URL = os.getenv('OLLAMA_API_URL', 'http://localhost:11434/api')
    
    # Ollama model to use - default to deepseek-r1:8b if not specified
    OLLAMA_MODEL = os.getenv('OLLAMA_MODEL', 'deepseek-r1:8b')
    
    # Email for PubMed from environment variable
    PUBMED_EMAIL = os.getenv('PUBMED_EMAIL')
    
    @staticmethod
    def validate():
        """
        Validate required configuration variables.
        
        Raises:
        -------
        ValueError
            If any required environment variables are missing.
        """
        missing = []
        if not Config.PUBMED_EMAIL:
            missing.append('PUBMED_EMAIL')
        if missing:
            raise ValueError(f"Missing required environment variables: {', '.join(missing)}") 